"""Validate bounded sign likelihoods and reproduce a seeded synthetic recovery."""

import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr
from scipy.stats import multivariate_normal
from scipy.optimize import minimize_scalar
from lib.selection import log_orthant, log_gaussian, log_censored
from lib.records import write_rows

DEFAULTS = {"replicates": 128, "seed": 26092703, "true_amplitude": 0.4, "epochs": 24}


def run(out, cfg):
    if (
        cfg["replicates"] < 1
        or cfg["epochs"] < 3
        or cfg["epochs"] > 100
        or not -1 < cfg["true_amplitude"] < 1
    ):
        raise ValueError("Unsupported recovery configuration")
    mean = np.array([-0.8, 0.2, 0.5])
    variance = np.array([1.0, 1.2, 0.9])
    loading = np.array([0.1, -0.1, 0.08])
    signs = np.array([1.0, -1.0, 1.0])
    numerical = quad(
        lambda z: (
            np.exp(-z * z / 2)
            / np.sqrt(2 * np.pi)
            * np.prod(ndtr(signs * (mean + loading * z) / np.sqrt(variance)))
        ),
        -12,
        12,
        epsabs=1e-12,
        epsrel=1e-12,
    )[0]
    gap = abs(log_orthant(mean, variance, loading, signs) - np.log(numerical))
    y = np.array([0.1, 0.3, -0.2])
    dense = multivariate_normal.logpdf(
        y, mean=mean, cov=np.diag(variance) + np.outer(loading, loading)
    )
    dense_gap = abs(log_gaussian(y, mean, variance, loading) - dense)
    if gap > 1e-9 or dense_gap > 1e-10:
        raise RuntimeError("Independent sign/density integration check failed")
    rng = np.random.default_rng(cfg["seed"])
    x = np.linspace(-2, 2, cfg["epochs"])
    template = np.exp(-x * x / 2)
    d = np.ones(len(x))
    v = np.full(len(x), 0.06)
    rows = []
    for replicate in range(cfg["replicates"]):
        y = (
            cfg["true_amplitude"] * template
            + v * rng.normal()
            + rng.normal(size=len(x))
        )
        keep = y > 0
        objectives = {
            "all_signed": lambda a: -log_gaussian(y, a * template, d, v),
            "positive_naive": lambda a: (
                -log_gaussian(y[keep], a * template[keep], d[keep], v[keep])
            ),
            "censored_with_known_schedule": lambda a: (
                -log_censored(
                    y[keep], a * template, d, v, keep, -np.ones((~keep).sum())
                )
            ),
        }
        for name, objective in objectives.items():
            fit = minimize_scalar(
                objective, bounds=(-3.0, 3.0), method="bounded", options={"xatol": 1e-8}
            )
            if not fit.success or not np.isfinite(fit.fun):
                raise RuntimeError("Synthetic optimizer failure")
            rows.append(
                {
                    "replicate": replicate,
                    "method": name,
                    "retained_positive": int(keep.sum()),
                    "amplitude": float(fit.x),
                    "boundary": abs(fit.x) > 2.999,
                }
            )
    write_rows(out / "recovery.csv", rows)
    summary = {}
    for name in objectives:
        values = np.array([r["amplitude"] for r in rows if r["method"] == name])
        summary[name] = {
            "mean": float(values.mean()),
            "bias": float(values.mean() - cfg["true_amplitude"]),
            "sd": float(values.std()),
            "mcse_of_mean": float(values.std(ddof=1) / np.sqrt(len(values)))
            if len(values) > 1
            else None,
            "boundary_fits": sum(r["boundary"] for r in rows if r["method"] == name),
        }
    return {
        "orthant_quad_log_gap": gap,
        "dense_gaussian_log_gap": dense_gap,
        "generator_factor_precision": float(np.sum(v * v / d)),
        "recovery": summary,
        "scope": "A clean, explicitly specified synthetic experiment, not a replay of every historical recovery design. All-signed, naive positive-only and full-schedule censored likelihoods use the same draws. Knowing omitted signs and the eligible schedule is an assumption; no real survey noise or selection law is inferred.",
    }
