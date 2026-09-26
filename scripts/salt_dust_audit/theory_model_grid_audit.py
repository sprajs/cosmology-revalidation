"""Bounded inspection of the pinned DES5YR SALT3 grids, without fitting data.

Run with phase2/env-official/bin/python. This is a support-domain diagnostic,
not an instrumented execution of SNANA's interpolation or covariance fallback.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
MODEL = ROOT / "sources/repos/des-science__DES-SN5YR@1.3/2_LCFIT_MODEL/SALT3.DES5YR"
OBS = ROOT / "phase2/official/results/recovered_refit_observables.csv"
OUT = ROOT / "docs/salt-dust-audit/theory-model-grid-audit.json"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name):
    path = MODEL / name
    return np.loadtxt(path)


def main():
    names = ["salt3_lc_variance_0.dat.gz", "salt3_lc_variance_1.dat.gz",
             "salt3_lc_covariance_01.dat.gz", "salt3_template_0.dat.gz",
             "salt3_template_1.dat.gz", "salt3_color_dispersion.dat.gz"]
    v0, v1, cov, m0, m1, disp = [load(name) for name in names]
    for arr in [v1, cov, m0, m1]:
        if not np.array_equal(v0[:, :2], arr[:, :2]):
            raise ValueError("Model grids differ; do not silently align them")
    if not all(np.isfinite(a).all() for a in [v0, v1, cov, m0, m1, disp]):
        raise ValueError("Nonfinite model value")
    obs = pd.read_csv(OBS)
    xlo, xhi = obs.x1_double.min(), obs.x1_double.max()
    clo, chi = obs.c_double.min(), obs.c_double.max()
    if not np.isfinite([xlo, xhi, clo, chi]).all():
        raise ValueError("Nonfinite fitted parameter range")
    results = {}
    for label, mask in {
        "whole_tabulated_domain": np.ones(len(v0), dtype=bool),
        "fit_phase_and_filter_mean_wavelength_domain": (
            (v0[:, 0] >= -15) & (v0[:, 0] <= 45)
            & (v0[:, 1] >= 3500) & (v0[:, 1] <= 8000)),
    }.items():
        a, b, d = v0[mask, 2], v1[mask, 2], cov[mask, 2]
        matrices = np.empty((len(a), 2, 2))
        matrices[:, 0, 0], matrices[:, 1, 1] = a, b
        matrices[:, 0, 1] = matrices[:, 1, 0] = d
        eig = np.linalg.eigvalsh(matrices)
        r = {"nodes": int(mask.sum()),
             "negative_min_eigenvalue_nodes": int((eig[:, 0] < 0).sum()),
             "minimum_eigenvalue": float(eig[:, 0].min()),
             "maximum_absolute_correlation": float(np.max(np.abs(d / np.sqrt(a*b)))),
             "negative_v0_nodes": int((a < 0).sum()),
             "negative_v1_nodes": int((b < 0).sum()),
             "parameter_domains": {}}
        for domain, (low, high) in {
            "published_basic_BBC_x1_interval": (-3., 3.),
            "envelope_of_all_recovered_refits": (xlo, xhi),
        }.items():
            # Endpoints plus an in-domain quadratic stationary point suffice.
            candidates = [np.full_like(b, low), np.full_like(b, high)]
            stationary = np.divide(-d, b, out=np.zeros_like(d), where=b != 0)
            candidates.append(np.clip(stationary, low, high))
            q = np.stack([a + 2*x*d + x*x*b for x in candidates])
            qmin = q.min(axis=0)
            # The mean flux is affine in x1: extrema occur at endpoints.
            fmin = np.minimum(m0[mask, 2] + low*m1[mask, 2],
                             m0[mask, 2] + high*m1[mask, 2])
            r["parameter_domains"][domain] = {
                "x1_min": float(low), "x1_max": float(high),
                "negative_quadratic_variance_nodes": int((qmin < 0).sum()),
                "minimum_quadratic_variance": float(qmin.min()),
                "nodes_with_negative_SED_for_some_x1_in_envelope": int((fmin < 0).sum()),
                "minimum_unscaled_SED_in_envelope": float(fmin.min()),
            }
        results[label] = r
    # Read the release law parameters and independently reproduce its C polynomial.
    params = np.loadtxt(MODEL / "salt3_color_correction.dat.gz", max_rows=6)[1:]
    wave = np.array([2800., 3500., 4302.57, 5428.55, 6500., 8000.])
    rl = (wave - 4302.57)/(5428.55 - 4302.57)
    polynomial = (1-params.sum())*rl
    for i, value in enumerate(params):
        polynomial += value*rl**(i+2)
    if not np.allclose(polynomial[2:4], [0, 1], atol=1e-12):
        raise ValueError("SALT B/V color-law anchors are inconsistent")
    grid_c = np.array([clo, chi])[:, None]
    multiplier = np.exp(np.log(10)/2.5 * grid_c * polynomial[None, :])
    # A finite color multiplier is positive and cannot change variance/SED signs.
    if not ((multiplier > 0).all() and np.isfinite(multiplier).all()):
        raise ValueError("Invalid finite-color multiplier")
    files = [MODEL/name for name in names]
    files += [MODEL/"SALT3.INFO", MODEL/"salt3_color_correction.dat.gz", OBS,
              Path(__file__).resolve(),
              ROOT/"sources/repos/RickKessler__SNANA/src/genmag_SALT2.c",
              ROOT/"sources/repos/RickKessler__SNANA/src/genmag_SALT2.h"]
    report = {
        "kind": "bounded_source_grid_diagnostic",
        "model": str(MODEL.relative_to(ROOT)),
        "refit_objects": len(obs), "refit_color_range": [clo, chi],
        "surface_tabulation": {"phase_days": [float(v0[:, 0].min()), float(v0[:, 0].max())],
                               "wavelength_angstrom": [float(v0[:, 1].min()), float(v0[:, 1].max())]},
        "grid_results": results,
        "color_dispersion": {"negative_nodes": int((disp[:, 1]<0).sum()),
                             "raw_maximum": float(disp[:, 1].max()),
                             "raw_nodes_above_INFO_cap_of_1": int((disp[:, 1]>1).sum())},
        "color_law_checks": {"wavelength_angstrom": wave.tolist(),
                             "C_polynomial": polynomial.tolist(),
                             "delta_monochromatic_magnitude_for_c_0p1": (-.1*polynomial).tolist(),
                             "refit_color_endpoint_flux_multipliers": multiplier.tolist()},
        "limitations": [
            "Uses original tabulation; does not reproduce SNANA spline resampling of mean surfaces.",
            "The rectangular phase/wavelength domain is not actual per-observation passband coverage.",
            "Filter wings can leave the allowed mean-wavelength interval.",
            "A negative monochromatic SED does not establish a negative integrated flux or an affected observation.",
            "Refit parameter envelope contains every successful recovered fit, including objects outside final BBC cuts.",
            "No execution counter for SNANA variance/flux sign fallback was instrumented.",
            "Positive covariance cannot establish an unbiased mean model or complete physical covariance.",
        ],
        "sha256": {str(p.relative_to(ROOT)): digest(p) for p in files},
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    print(json.dumps({"output": str(OUT.relative_to(ROOT)), "grid_results":results}, indent=2))


if __name__ == "__main__":
    main()
