"""Fixed-covariance calibrated-flux fits and conditional prediction."""

import json
import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares
from .flux_engine import BUNDLE


def jacobian(fun, x):
    h = np.array([2e-5, 2e-4, 2e-5, 2e-3] + ([2e-5] if len(x) == 5 else []))
    return np.column_stack(
        [
            (fun(x + np.eye(len(x))[j] * h[j]) - fun(x - np.eye(len(x))[j] * h[j]))
            / (2 * h[j])
            for j in range(len(x))
        ]
    )


def hessian(fun, x, scale=1.0):
    h = np.array([2e-4, 2e-3, 2e-4, 1e-2] + ([2e-4] if len(x) == 5 else [])) * scale
    f = fun(x)
    n = len(x)
    out = np.empty((n, n))
    eye = np.eye(n)
    for i in range(n):
        out[i, i] = (fun(x + h[i] * eye[i]) - 2 * f + fun(x - h[i] * eye[i])) / h[
            i
        ] ** 2
        for j in range(i):
            out[i, j] = out[j, i] = (
                fun(x + h[i] * eye[i] + h[j] * eye[j])
                - fun(x + h[i] * eye[i] - h[j] * eye[j])
                - fun(x - h[i] * eye[i] + h[j] * eye[j])
                + fun(x - h[i] * eye[i] - h[j] * eye[j])
            ) / (4 * h[i] * h[j])
    return out


def normal_logpdf(residual, cov):
    l = np.linalg.cholesky(cov)
    white = solve_triangular(l, residual, lower=True)
    return float(
        -0.5
        * (
            white @ white
            + 2 * np.log(np.diag(l)).sum()
            + len(residual) * np.log(2 * np.pi)
        )
    )


class Case:
    def __init__(self, engine, cid):
        self.engine, self.cid = engine, cid
        self.q = {k: v for k, v in np.load(BUNDLE / f"objective_{cid}.npz").items()}
        self.z = float(self.q["zHEL"][0])
        self.ebv = float(self.q["MWEBV"][0])
        self.p0 = self.q["parameters_x0_x1_c_t0"].copy()
        self.p0[0] = np.log(self.p0[0])
        self.tref = self.p0[3]
        self.x0 = self.p0.copy()
        self.x0[3] = 0.0
        self.t = self.q["MJD"]
        self.b = self.q["band"]
        self.y = self.q["data_flux"]
        self.grid = engine.prepare(self.z, self.ebv, 5.0)
        self.priors = json.loads((BUNDLE / "exact_prior_setup.json").read_text())[
            "objects"
        ][cid]
        self.priorcentre = self.priors["actual_t0_prior_centre"] - self.tref
        self.cov = self.q["frozen_flux_covariance"]

    def flux(self, x, family="salt", interpolation="sncosmo"):
        p = x.copy()
        p[3] += self.tref
        return self.engine.flux(
            p, self.b, self.t, self.z, self.grid, family, interpolation
        )

    def prior_residuals(self, x):
        out = [(x[3] - self.priorcentre) / 10.0]
        if len(x) == 5:
            out.append(x[4] / 0.1)
        return np.array(out)

    def fit(
        self,
        family="salt",
        train=None,
        cov=None,
        y=None,
        start=None,
        interpolation="sncosmo",
        two_starts=False,
    ):
        n = len(self.y)
        train = np.arange(n) if train is None else np.asarray(train)
        cov = self.cov if cov is None else cov
        y = self.y if y is None else y
        ctrain = cov[np.ix_(train, train)]
        l = np.linalg.cholesky(ctrain)
        pinit = np.r_[self.x0, 0.0] if family == "phase_colour" else self.x0.copy()
        if start is not None:
            pinit = np.array(start)
        low = np.array(
            [self.x0[0] - 3.0, -5.0, -0.5, -10.0]
            + ([-0.3] if family == "phase_colour" else [])
        )
        high = np.array(
            [self.x0[0] + 3.0, 5.0, 0.5, 10.0]
            + ([0.3] if family == "phase_colour" else [])
        )
        # Preserve model domain during all optimizer evaluations.
        low[3] = max(
            low[3],
            float(np.max(self.t - self.tref - (1 + self.z) * self.engine.tmax) + 1e-6),
        )
        high[3] = min(
            high[3],
            float(np.min(self.t - self.tref - (1 + self.z) * self.engine.tmin) - 1e-6),
        )

        def residual(x):
            r = solve_triangular(
                l, y[train] - self.flux(x, family, interpolation)[train], lower=True
            )
            return np.r_[r, self.prior_residuals(x)]

        starts = [np.clip(pinit, low + 1e-8, high - 1e-8)]
        if two_starts:
            starts.append(
                np.clip(
                    pinit
                    + np.array(
                        [0.05, 0.15, -0.01, 0.15]
                        + ([0.02] if family == "phase_colour" else [])
                    ),
                    low + 1e-8,
                    high - 1e-8,
                )
            )
        sols = [
            least_squares(
                residual,
                s,
                bounds=(low, high),
                diff_step=1e-4,
                xtol=1e-10,
                ftol=1e-10,
                gtol=1e-8,
                max_nfev=250,
            )
            for s in starts
        ]
        sol = min(sols, key=lambda s: s.fun @ s.fun)
        x = sol.x
        objective = lambda a: 0.5 * np.sum(residual(a) ** 2)
        h = hessian(objective, x)
        h2 = hessian(objective, x, 0.5)
        eig = np.linalg.eigvalsh(h)
        hc = None
        if eig[0] > 0:
            hc = np.linalg.inv(h)
        boundary = bool(np.any(np.minimum(x - low, high - x) < 1e-4))
        j = jacobian(lambda a: self.flux(a, family, interpolation), x)
        out = {
            "CID": self.cid,
            "family": family,
            "ntrain": len(train),
            "success": bool(sol.success),
            "message": sol.message,
            "nfev": int(sol.nfev),
            "x": x.tolist(),
            "chi2_total": float(sol.fun @ sol.fun),
            "chi2_data": float(
                sol.fun[: -len(self.prior_residuals(x))]
                @ sol.fun[: -len(self.prior_residuals(x))]
            ),
            "hessian_min_eigenvalue": float(eig[0]),
            "boundary": boundary,
            "hessian_stepsize_relative_difference": float(
                np.linalg.norm(h - h2) / np.linalg.norm(h)
            ),
            "hessian_covariance": None if hc is None else hc.tolist(),
            "gauss_newton_covariance": np.linalg.inv(sol.jac.T @ sol.jac).tolist(),
            "two_start_chi2_gap": float(
                max(s.fun @ s.fun for s in sols) - min(s.fun @ s.fun for s in sols)
            ),
        }
        return out, self.flux(x, family, interpolation), j, hc

    def predictive(self, fit, mu, j, hc, train, test, cov=None, y=None):
        cov = self.cov if cov is None else cov
        y = self.y if y is None else y
        train, test = np.asarray(train), np.asarray(test)
        if hc is None or fit["boundary"] or not fit["success"]:
            return {"valid": False}
        a = np.linalg.solve(cov[np.ix_(train, train)], cov[np.ix_(train, test)]).T
        condmean = mu[test] + a @ (y[train] - mu[train])
        noise = cov[np.ix_(test, test)] - a @ cov[np.ix_(train, test)]
        jp = j[test] - a @ j[train]
        predictive_cov = noise + jp @ hc @ jp.T
        r = y[test] - condmean
        return {
            "valid": True,
            "n_test": len(test),
            "logpdf_plugin": normal_logpdf(r, noise),
            "logpdf_laplace": normal_logpdf(r, predictive_cov),
            "test_conditional_chi2": float(r @ np.linalg.solve(predictive_cov, r)),
            "test_marginal_standardized_residuals": (
                r / np.sqrt(np.diag(predictive_cov))
            ).tolist(),
            "test_model": condmean.tolist(),
            "test_variance": np.diag(predictive_cov).tolist(),
        }
