"""Released-distance fitting with explicit priors and sampler qualification."""

import numpy as np
import pandas as pd
import emcee
from scipy.optimize import minimize
from scipy.interpolate import PchipInterpolator
from astropy.cosmology import Flatw0waCDM
from lib.cosmology import Pantheon, BAO, mu, qvalue, integral
from lib.records import write_json, sha256
from lib.paths import ROOT

DEFAULTS = {
    "model": "lcdm",
    "data": "sn",
    "seed": 41001,
    "steps": 7000,
    "burn": 1500,
    "walkers": 40,
    "tau": 0.5,
    "correction_csv": None,
    "correction_column": "delta_mu",
    "correction_scale": 1.0,
    "amplitude": "none",
    "optimize_only": False,
}


def describe(values):
    return {
        "mean": float(np.mean(values)),
        "sd": float(np.std(values)),
        "q025_median_q975": np.quantile(values, [0.025, 0.5, 0.975]).tolist(),
    }


def run(out, cfg):
    model, dataset = cfg["model"], cfg["data"]
    choices = {
        "lcdm": ([0.33], ["Om"], [(0.01, 0.99)]),
        "wcdm": ([0.3, -0.9], ["Om", "w"], [(0.01, 0.99), (-3.0, 1.0)]),
        "cpl": (
            [0.3, -0.8, -0.5],
            ["Om", "w0", "wa"],
            [(0.01, 0.99), (-3.0, 1.0), (-3.0, 2.0)],
        ),
        "kinematic": ([-0.4, 0.8], ["q0", "q1"], [(-3.0, 2.0), (-8.0, 8.0)]),
        "qbins": (
            [-0.4, -0.4, -0.2, 0.1, 0.3],
            [f"q{i}" for i in range(5)],
            [(-3.0, 2.0)] * 5,
        ),
    }
    if model not in choices or dataset not in {"sn", "bao", "joint"}:
        raise ValueError("Unknown model or dataset")
    if dataset != "sn" and model != "cpl":
        raise ValueError("This preserved BAO likelihood is parameterized for CPL")
    if (
        cfg["amplitude"] not in {"none", "fixed", "normal", "uniform"}
        or cfg["tau"] <= 0
    ):
        raise ValueError("Invalid amplitude prior or q-bin smoothing scale")
    sn = Pantheon() if dataset in {"sn", "joint"} else None
    bao = BAO() if dataset in {"bao", "joint"} else None
    correction = cfg["correction_csv"]
    if bool(correction) != (cfg["amplitude"] != "none"):
        raise ValueError(
            "An explicit template and non-none amplitude setting must be supplied together"
        )
    if correction:
        if sn is None:
            raise ValueError("A magnitude correction requires supernova data")
        path = (ROOT / correction).resolve()
        if not path.is_relative_to(ROOT):
            raise ValueError("Correction must be inside this research folder")
        frame = pd.read_csv(path)
        values = (
            PchipInterpolator(
                frame.z, frame[cfg["correction_column"]], extrapolate=False
            )(sn.z)
            * cfg["correction_scale"]
        )
        if not np.isfinite(values).all():
            raise ValueError("Correction does not cover all input redshifts")
        sn.reset_template(values)
        write_json(
            out / "correction-input.json",
            {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)},
        )
    start, names, bounds = choices[model]
    ncos = len(start)
    if bao:
        start += [10000.0]
        names += ["H0_rd"]
        bounds += [(5000.0, 15000.0)]
    if cfg["amplitude"] in {"normal", "uniform"}:
        start += [1.0]
        names += ["evolution_amplitude"]
        bounds += [(-2.0, 4.0)]
    bounds = np.asarray(bounds)
    write_json(
        out / "priors.json",
        {
            "bounds": dict(zip(names, bounds.tolist())),
            "CPL_constraint": "w0+wa<0" if model == "cpl" else None,
            "normal_amplitude": [1.0, 4 / 30] if cfg["amplitude"] == "normal" else None,
            "qbin_difference_sd": cfg["tau"] if model == "qbins" else None,
        },
    )

    def chi(theta):
        theta = np.atleast_2d(theta)
        amp = (
            theta[:, -1]
            if cfg["amplitude"] in {"normal", "uniform"}
            else (1.0 if cfg["amplitude"] == "fixed" else 0.0)
        )
        value = np.zeros(len(theta))
        if sn is not None:
            value += sn.chisq(theta[:, :ncos], model, amp)
        if bao is not None:
            value += bao.chisq(theta[:, :ncos], theta[:, ncos])
        return value

    def logpost(theta):
        theta = np.atleast_2d(theta)
        lp = np.full(len(theta), -np.inf)
        keep = np.all((theta > bounds[:, 0]) & (theta < bounds[:, 1]), axis=1)
        if model == "cpl":
            keep &= theta[:, 1] + theta[:, 2] < 0
        values = theta[keep]
        if len(values):
            lp[keep] = -0.5 * chi(values)
            if cfg["amplitude"] == "normal":
                lp[keep] -= 0.5 * ((values[:, -1] - 1) / (4 / 30)) ** 2
            if model == "qbins":
                lp[keep] -= (
                    0.5
                    * np.sum(np.diff(values[:, :ncos], axis=1) ** 2, axis=1)
                    / cfg["tau"] ** 2
                )
        return lp

    rng = np.random.default_rng(cfg["seed"])
    opts = []
    for i in range(4):
        initial = np.array(start) + (
            rng.normal(size=len(start)) * 0.02 * np.diff(bounds).ravel() if i else 0
        )
        initial = np.clip(initial, bounds[:, 0] + 1e-5, bounds[:, 1] - 1e-5)
        fit = minimize(
            lambda x: -logpost(x)[0],
            initial,
            method="Nelder-Mead",
            options={"maxiter": 15000, "xatol": 1e-7, "fatol": 1e-7},
        )
        opts.append(
            {
                "x": fit.x.tolist(),
                "objective": float(fit.fun),
                "success": bool(fit.success),
            }
        )
    valid = [x for x in opts if x["success"] and np.isfinite(x["objective"])]
    if not valid:
        raise RuntimeError("No optimizer start converged")
    best = min(valid, key=lambda x: x["objective"])
    point = np.array(best["x"])
    result = {
        "scope": "Conditional inference from released corrected distances; no independent calibration or selection reconstruction.",
        "parameters": names,
        "optimizer_starts": opts,
        "mode": point.tolist(),
        "mode_chisq": float(chi(point)[0]),
        "mode_q0": float(qvalue([0.0], point[:ncos], model)[0, 0]),
        "q0_meaning": "first bin 0<z<0.1"
        if model == "qbins"
        else "model-extrapolated q(0)",
    }
    z = np.array([0.01, 0.1, 0.8, 2.3])
    ref = (
        Flatw0waCDM(H0=70, Om0=0.353, w0=-0.42, wa=-1.75, Tcmb0=0).distmod(z).value
        - 5 * np.log10(299792.458 / 70)
        - 25
    )
    numerical = {
        "astropy_distance_max_mag": float(
            np.max(abs(mu(z, [0.353, -0.42, -1.75])[0] - ref))
        ),
        "constant_q_zero_integral_error": float(
            np.max(abs(integral(z, [0.0] * 5, "qbins")[0] - np.log1p(z)))
        ),
    }
    if sn:
        amplitude = (
            point[-1]
            if cfg["amplitude"] in {"normal", "uniform"}
            else float(cfg["amplitude"] == "fixed")
        )
        numerical["mode_compressed_full_chisq_gap"] = float(
            abs(
                sn.chisq(point[:ncos], model, amplitude)[0]
                - sn.chisq_full(point[:ncos], model, amplitude)[0]
            )
        )
        if numerical["mode_compressed_full_chisq_gap"] > 0.01:
            raise RuntimeError(
                "Full-covariance likelihood does not match compressed likelihood"
            )
        result["data"] = sn.diagnostics()
    if (
        numerical["astropy_distance_max_mag"] > 1e-7
        or numerical["constant_q_zero_integral_error"] > 1e-12
    ):
        raise RuntimeError("Independent distance checks failed")
    result["numerical_checks"] = numerical
    if cfg["optimize_only"]:
        result["status"] = "optimizer_only_no_posterior"
        return result
    steps, burn, walkers = cfg["steps"], cfg["burn"], cfg["walkers"]
    if burn < 0 or steps - burn < 100 or walkers < 2 * len(point):
        raise ValueError(
            "Need at least 100 retained steps and twice the dimension in walkers"
        )
    initial = []
    for _ in range(100000):
        x = point + rng.normal(size=len(point)) * 0.005 * np.diff(bounds).ravel()
        if np.isfinite(logpost(x)[0]):
            initial.append(x)
        if len(initial) == walkers:
            break
    if len(initial) != walkers:
        raise RuntimeError("Could not initialize walkers inside the prior")
    np.random.seed(cfg["seed"])
    sampler = emcee.EnsembleSampler(walkers, len(point), logpost, vectorize=True)
    sampler.run_mcmc(np.array(initial), steps, progress=False)
    chain = sampler.get_chain(discard=burn)
    tau = sampler.get_autocorr_time(discard=burn, tol=0)
    good = bool(
        np.isfinite(tau).all()
        and (tau > 0).all()
        and np.min((steps - burn) / tau) >= 50
    )
    np.savez_compressed(
        out / "chains.npz",
        chain=chain,
        log_prob=sampler.get_log_prob(discard=burn),
        names=names,
    )
    result.update(
        status="posterior_diagnostics_pass"
        if good
        else "insufficient_autocorrelation_length",
        valid_for_posterior_summary=good,
        autocorrelation_time=tau.tolist(),
        acceptance_fraction=float(sampler.acceptance_fraction.mean()),
        convergence_scope="Within-ensemble autocorrelation length; independent seeds remain a separate robustness check.",
    )
    if good:
        flat = chain.reshape(-1, len(point))
        q = qvalue([0.0], flat[:, :ncos], model)[:, 0]
        result["posterior"] = {
            name: describe(flat[:, j]) for j, name in enumerate(names)
        }
        result["q0"] = describe(q)
        result["P_q0_negative"] = float(np.mean(q < 0))
        select = rng.choice(len(flat), min(5000, len(flat)), replace=False)
        grid = np.linspace(0, 2.3, 231)
        values = qvalue(grid, flat[select, :ncos], model)
        pd.DataFrame(
            {
                "z": grid,
                "q025": np.quantile(values, 0.025, axis=0),
                "median": np.median(values, axis=0),
                "q975": np.quantile(values, 0.975, axis=0),
            }
        ).to_csv(out / "q_history.csv", index=False)
    return result
