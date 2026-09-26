"""Paired Gaussian control for nonlinear fitting; no selection or bias inference."""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.salt_dust_audit.flux_response import build_model, K

BASE = ROOT / "runs/research_2026_09_26/astra_design"
OUT = ROOT / "runs/research_2026_09_26/gaussian_fit_control"


def main():
    OUT.mkdir(exist_ok=False)
    source = Path(__file__).read_bytes()
    (OUT / "executed_source.py").write_bytes(source)
    cohort_path = ROOT / "runs/research_2026_09_26/calibration_pattern_refit/cohort.csv"
    coefficient_path = BASE / "validation1020/frozen-discovery-coefficients.npz"
    cohort = pd.read_csv(cohort_path, dtype={"CID": str})
    coefficients = np.load(coefficient_path)
    band_vector = coefficients["gauge_griz"] @ coefficients["basis_mean"]
    protocol = {"scope": "Conditional Gaussian nonlinear-fit mechanism control, no measured calibration bias or selection validation.",
                "cohort": "The same six previously fixed discovery redshift ranks; no new residual-based selection.",
                "draws_per_object": 128, "seed": 2026092619,
                "generator": "Native SALT mean at exact43 exported coordinates plus Gaussian noise from frozen native SNANA C, same accepted epochs.",
                "fit": "Four native SALT coordinates [delta_mB,delta_x1,delta_c,delta_t0], frozen C, no priors, bounds[-2,-8,-1,-30] to[2,8,1,30].",
                "score": "Unchanged exact43 observer vector; compute matched product and information at generating and fitted means/tangents with paired noise.",
                "verification": "Every16th draw also uses fixed shifted start; retain all failures/boundary hits and paired objective differences.",
                "limitations": "Not native SNANA optimizer/priors/clipping or model-C feedback; no object/epoch/classifier selection regeneration or physical-scatter truth."}
    (OUT / "protocol.json").write_text(json.dumps(protocol, indent=2) + "\n")
    model, bands, paths, zp = build_model()
    offsets = {str(row['Filter Name'])[-1]: float(row['Primary Mag']) for row in zp}
    inputs = [cohort_path, coefficient_path, ROOT / "scripts/salt_dust_audit/flux_response.py"] + paths
    rng = np.random.default_rng(protocol["seed"])
    steps = np.array([1e-4, 1e-3, 1e-4, .01])
    draws = []
    metadata = []
    bound_lo = np.array([-2., -8., -1., -30.])
    bound_hi = -bound_lo
    for row in cohort.itertuples():
        path = BASE / f"exact43/objectives/objective_{row.CID}.npz"
        inputs.append(path)
        data = np.load(path)
        times, labels = data["MJD"], data["band"]
        bandpasses = np.array([bands[b] for b in labels], dtype=object)
        conversion = np.array([10**(-.4 * (.27 + offsets[b])) for b in labels])
        x0, x1, color, t0 = data["parameters_x0_x1_c_t0"]
        redshift, ebv = float(data["zHEL"][0]), float(data["MWEBV"][0])
        chol = np.linalg.cholesky(data["frozen_flux_covariance"])
        dimming = np.array([band_vector['griz'.index(b)] for b in labels])

        def flux(theta):
            model.set(z=redshift, t0=t0 + theta[3], x0=x0 * np.exp(-K * theta[0]),
                      x1=x1 + theta[1], c=color + theta[2], mwebv=ebv, mwrv=3.1,
                      hostebv=0, hostrv=3.1)
            return conversion * model.bandflux(bandpasses, times, zp=27.5, zpsys='ab')

        def tangent(theta, mean):
            jac = np.column_stack([(flux(theta + np.eye(4)[j] * h) -
                                    flux(theta - np.eye(4)[j] * h)) / (2 * h)
                                   for j, h in enumerate(steps)])
            jac[:, 0] = -K * mean
            jw = solve_triangular(chol, jac, lower=True)
            u, singular, _ = np.linalg.svd(jw, full_matrices=False)
            assert np.sum(singular > singular[0] * 1e-10) == 4
            raw = solve_triangular(chol, -K * mean * dimming, lower=True)
            prediction = raw - u @ (u.T @ raw)
            assert np.linalg.norm(u.T @ prediction) < 1e-9
            return prediction, u

        truth = np.zeros(4)
        mean0 = flux(truth)
        vector0, u0 = tangent(truth, mean0)
        info0 = float(vector0 @ vector0)
        metadata.append({"CID": row.CID, "zHEL": redshift, "epochs": len(times),
                         "generating_information": info0,
                         "native_vs_exported_mean_max_fraction": float(max(abs(mean0 / data['model_flux'] - 1)))})
        for repeat in range(128):
            noise = rng.normal(size=len(times))
            y = mean0 + chol @ noise
            fun = lambda theta: solve_triangular(chol, y - flux(theta), lower=True)
            fit = least_squares(fun, truth, bounds=(bound_lo, bound_hi),
                                xtol=1e-10, ftol=1e-10, gtol=1e-9, max_nfev=300)
            mean = flux(fit.x)
            vector, fitted_u = tangent(fit.x, mean)
            residual = solve_triangular(chol, y - mean, lower=True)
            linear_residual = noise - u0 @ (u0.T @ noise)
            linear_a = float(vector0 @ linear_residual)
            assert abs(linear_a - vector0 @ noise) < 1e-9
            linear_chi2 = float(linear_residual @ linear_residual)
            nonlinear_a = float(vector @ residual)
            fixed_tangent_a = float(vector0 @ residual)
            info = float(vector @ vector)
            second_difference = None
            second_success = None
            if repeat % 16 == 0:
                second = least_squares(fun, np.array([.01, .1, .005, .1]), bounds=(bound_lo, bound_hi),
                                       xtol=1e-10, ftol=1e-10, gtol=1e-9, max_nfev=300)
                second_difference = float(second.fun @ second.fun - fit.fun @ fit.fun)
                second_success = bool(second.success)
            draws.append({"CID": row.CID, "replicate": repeat,
                          "success": bool(fit.success), "status": int(fit.status), "nfev": int(fit.nfev),
                          "boundary_hit": bool(np.any(fit.x - bound_lo < 1e-5) or np.any(bound_hi - fit.x < 1e-5)),
                          "optimality": float(fit.optimality), "dof": len(times) - 4,
                          "linear_matched_product": linear_a, "linear_information": info0,
                          "linear_gain": linear_a - .5 * info0, "linear_chi2": linear_chi2,
                          "nonlinear_matched_product": nonlinear_a, "nonlinear_information": info,
                          "nonlinear_gain": nonlinear_a - .5 * info,
                          "nonlinear_chi2": float(residual @ residual),
                          "nonlinear_residual_fixed_tangent_product": fixed_tangent_a,
                          "second_start_objective_difference": second_difference,
                          "second_start_success": second_success,
                          **{f"fit_delta_{k}": float(v) for k, v in zip(["mB", "x1", "c", "t0"], fit.x)}})
        print(row.CID, "completed128 Gaussian fits", flush=True)
    frame = pd.DataFrame(draws)
    frame.to_csv(OUT / "draws.csv", index=False)
    pd.DataFrame(metadata).to_csv(OUT / "cohort.csv", index=False)
    sums = frame.groupby("replicate").sum(numeric_only=True)
    summaries = {}
    metrics = ["linear_matched_product", "nonlinear_matched_product", "nonlinear_residual_fixed_tangent_product",
               "linear_information", "nonlinear_information", "linear_gain", "nonlinear_gain"]
    sums["paired_product_difference"] = sums.nonlinear_matched_product - sums.linear_matched_product
    sums["paired_information_difference"] = sums.nonlinear_information - sums.linear_information
    sums["paired_gain_difference"] = sums.nonlinear_gain - sums.linear_gain
    for metric in metrics + ["paired_product_difference", "paired_information_difference", "paired_gain_difference"]:
        value = sums[metric].to_numpy()
        summaries[metric] = {"mean": float(value.mean()), "Monte_Carlo_SE": float(value.std(ddof=1) / np.sqrt(len(value))),
                             "draw_sd": float(value.std(ddof=1))}
    report = {"scope": protocol["scope"], "objects": len(cohort), "draws_per_object": 128,
              "fits": len(frame), "all_success": bool(frame.success.all()),
              "boundary_hits": int(frame.boundary_hit.sum()),
              "second_start_checks": int(frame.second_start_success.notna().sum()),
              "second_start_all_success": bool(frame.second_start_success.dropna().all()),
              "second_start_max_absolute_objective_difference": float(frame.second_start_objective_difference.abs().max()),
              "linear_and_nonlinear_pooled_chi2_per_dof": [float(frame.linear_chi2.sum() / frame.dof.sum()),
                                                          float(frame.nonlinear_chi2.sum() / frame.dof.sum())],
              "six_object_sum_statistics": summaries,
              "normalization": "Each replicate sums six independently drawn objects. Monte Carlo SE describes this finite conditional experiment, not uncertainty in a real survey bias.",
              "limitations": protocol["limitations"]}
    (OUT / "result.json").write_text(json.dumps(report, indent=2) + "\n")
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    outputs = list(OUT.iterdir())
    (OUT / "manifest.json").write_text(json.dumps({"script_sha256": sha(__file__),
        "inputs_sha256": {str(p.relative_to(ROOT)): sha(p) for p in inputs},
        "outputs_sha256": {str(p.relative_to(ROOT)): sha(p) for p in outputs if p.is_file()}}, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    assert report["all_success"] and report["second_start_all_success"] and not report["boundary_hits"]
    assert report["second_start_max_absolute_objective_difference"] < 1e-5


if __name__ == "__main__":
    main()
