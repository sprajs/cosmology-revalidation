"""Flat, nonaccelerating BAO shape cone with a free distance scale."""

import hashlib
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.linalg import solve_triangular
from scipy.optimize import nnls, minimize, LinearConstraint
from scipy.stats import beta
from .paths import DATA, ROOT

INPUT = DATA / "bao"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def load_data():
    names = [
        "desi_gaussian_bao_ALL_GCcomb_mean.txt",
        "desi_gaussian_bao_ALL_GCcomb_cov.txt",
    ]
    paths = [INPUT / name for name in names]
    hashes = {str(p.relative_to(ROOT)): digest(p.read_bytes()) for p in paths}
    df = pd.read_csv(paths[0], sep=r"\s+", comment="#", names=["z", "value", "kind"])
    cov = np.loadtxt(paths[1])
    zs = np.sort(df.loc[df.kind == "DH_over_rs", "z"].unique())
    order = []
    for kind in ["DM_over_rs", "DH_over_rs"]:
        for z in zs:
            idx = np.flatnonzero((df.z == z) & (df.kind == kind))
            assert len(idx) == 1
            order.append(int(idx[0]))
    scale = np.r_[np.ones(len(zs)), 1 + zs]
    y = df.value.to_numpy()[order] * scale
    C = cov[np.ix_(order, order)] * np.outer(scale, scale)
    return zs, y, C, hashes, order, df, cov


def cone_matrices(t):
    n = len(t)
    rows, labels, rays = [], [], []
    for j in range(n):
        r = np.zeros(2 * n)
        r[n + j] = 1
        rows.append(r)
        labels.append(f"positive_g_{j}")
        r = np.zeros(2 * n)
        r[j] = 1
        r[n + j] = -(t[j] - (t[j - 1] if j else 0))
        if j:
            r[j - 1] = -1
        rows.append(r)
        labels.append(f"lower_integral_{j}")
        if j:
            r = np.zeros(2 * n)
            r[n + j - 1] = 1
            r[n + j] = -1
            rows.append(r)
            labels.append(f"monotone_g_{j}")
            r = np.zeros(2 * n)
            r[j - 1] = 1
            r[j] = -1
            r[n + j - 1] = t[j] - t[j - 1]
            rows.append(r)
            labels.append(f"upper_integral_{j}")
        for side in ["before", "after"]:
            g = (np.arange(n) < j if side == "before" else np.arange(n) <= j).astype(
                float
            )
            rays.append(np.r_[np.minimum(t, t[j]), g])
    return np.array(rows), np.array(rays).T, labels


class Cone:
    def __init__(self, cov, rays, inequalities):
        self.L = np.linalg.cholesky(cov)
        self.rays = rays
        self.A = inequalities
        W = solve_triangular(self.L, rays, lower=True)
        self.norm = np.linalg.norm(W, axis=0)
        self.W = W / self.norm
        self.B = inequalities @ self.L
        self.B /= np.linalg.norm(self.B, axis=1)[:, None]

    def fit(self, y):
        v = solve_triangular(self.L, y, lower=True)
        coef, norm = nnls(self.W, v, maxiter=1000)
        fitted = self.rays @ (coef / self.norm)
        assert np.min(self.A @ fitted) > -1e-8
        resid = self.W @ coef - v
        grad = self.W.T @ resid
        kkt = max(
            float(np.max(np.abs(grad[coef > 1e-8]), initial=0)),
            float(np.max(-grad[coef <= 1e-8], initial=0)),
        )
        return norm * norm, fitted, kkt

    def independent(self, y):
        v = solve_triangular(self.L, y, lower=True)
        fit = minimize(
            lambda x: 0.5 * np.sum((x - v) ** 2),
            np.zeros(len(y)),
            jac=lambda x: x - v,
            method="SLSQP",
            constraints=[LinearConstraint(self.B, 0, np.inf)],
            options={"ftol": 1e-11, "maxiter": 2000},
        )
        feasibility = float(np.min(self.B @ fit.x))
        if not fit.success and feasibility < -1e-7:
            raise RuntimeError(str(fit.message))
        return 2 * float(fit.fun), feasibility, bool(fit.success)

    def statistics(self, means, standard_noise):
        return np.array(
            [
                self.fit(mean + self.L @ eps)[0]
                for mean, eps in zip(
                    np.broadcast_to(means, standard_noise.shape), standard_noise
                )
            ]
        )


def tail(values, threshold):
    n = len(values)
    k = int(np.sum(values >= threshold - 1e-10))
    lo = 0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return {
        "exceedances": k,
        "draws": n,
        "fraction": k / n,
        "binomial_95_interval": [lo, hi],
        "plus_one_tail": (k + 1) / (n + 1),
    }


def lcdm_vector(z, om, scale):
    e = lambda x: np.sqrt(om * (1 + x) ** 3 + 1 - om)
    d = np.array(
        [
            scale * quad(lambda x: 1 / e(x), 0, zz, epsabs=1e-12, epsrel=1e-12)[0]
            for zz in z
        ]
    )
    g = scale * (1 + z) / e(z)
    return np.r_[d, g]
