#!/usr/bin/env python3
"""CALSPEC-observer null and saved-contrast prior pushforwards; no object fits."""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = HERE.parent / "astra_design"
EXP = BASE / "expanded12"
DIST = HERE.parent / "shared_distance_response_12"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    # Analytically predicted null coordinate, without estimating it from data.
    units = json.loads((EXP / "calspec-unit-map.json").read_text())["bands"]
    magnitude = -np.array([1 / units[b]["flux_scale_float32"] - 1 for b in "griz"]) / (.4 * np.log(10))
    magnitude -= magnitude.mean()
    n = np.zeros(15)
    n[9] = 1
    n[12:] = -magnitude[[0, 2, 3]] / .02
    n /= np.linalg.norm(n)
    originalpost = np.load(EXP / "posterior-models.npz", allow_pickle=False)
    null = {"CALSPEC_equivalent_linear_zero_sum_griz_mag": magnitude.tolist(),
            "prior_whitened_null_vector": n.tolist(), "prior_null_variance": 1.}
    for name in ("discovery", "validation"):
        a = np.load(EXP / f"{name}-projected-modes.npz", allow_pickle=False)
        x = np.column_stack([a["calibration_design"], .02 * a["observer_design"]])
        null[name + "_null_prediction_norm"] = float(np.linalg.norm(x @ n))
        null[name + "_smallest_design_singular_value"] = float(np.linalg.svd(x, compute_uv=False)[-1])
        assert null[name + "_null_prediction_norm"] < 1e-10
    null["discovery_null_posterior_variance"] = float(n @ originalpost["systematics_plus_observer_discovery_posterior_covariance"] @ n)
    null["discovery_null_posterior_mean"] = float(n @ originalpost["systematics_plus_observer_discovery_posterior_mean"])
    assert abs(null["discovery_null_posterior_variance"] - 1) < 1e-10
    null["source_sha256"] = sha(__file__)
    (HERE / "calspec-observer-null-check.json").write_text(json.dumps(null, indent=2) + "\n")
    arrays = np.load(DIST / "responses.npz", allow_pickle=False)
    original = json.loads((DIST / "result.json").read_text())
    with (DIST / "contrast-membership.csv").open() as file:
        membership = list(csv.DictReader(file))
    assert [r["CID"] for r in membership] == list(arrays["validation_CID"])
    low = np.array([r["low_quartile"].lower() == "true" for r in membership])
    high = np.array([r["high_quartile"].lower() == "true" for r in membership])
    assert low.sum() == high.sum() == 255 and not (low & high).any()
    contrast = arrays["validation"][high].mean(0) - arrays["validation"][low].mean(0)
    contrast_error = float(np.max(np.abs(contrast - arrays["contrast"])))
    assert contrast_error < 1e-12
    paperpost = np.load(HERE / "colorlaw-weight-posterior.npz", allow_pickle=False)
    ciso = 2 * .02**2 * (np.eye(3) - np.ones((3, 3)) / 4)
    iso_map = np.linalg.cholesky(ciso) / .02
    results, errors = {}, []
    for weight_name, multiplier, post in (("config_amplitude_0.3", 1., originalpost),
                                          ("paper_amplitude_sqrt_one_third", np.sqrt(1 / 3) / .3, paperpost)):
        results[weight_name] = {}
        for model in ("systematics_only", "systematics_plus_observer", "systematics_plus_isotropic_griz_observer"):
            transform = np.eye(15)
            transform[11, 11] = multiplier
            if "isotropic" in model:
                transform[12:, 12:] = iso_map
            if model == "systematics_only":
                transform = transform[:, :12]
            h = arrays["contrast"] @ transform
            all_responses = arrays["validation"] @ transform
            independently_rebuilt = all_responses[high].mean(0) - all_responses[low].mean(0)
            errors.append(float(np.max(np.abs(h - independently_rebuilt))))
            suffix = "_discovery_posterior_" if weight_name == "config_amplitude_0.3" else "_"
            mean, covariance = post[model + suffix + "mean"], post[model + suffix + "covariance"]
            value, variance = float(h @ mean), float(h @ covariance @ h)
            assert variance >= 0
            chol = np.linalg.cholesky(covariance)
            independent_value = float(independently_rebuilt @ mean)
            independent_sd = float(np.linalg.norm(independently_rebuilt @ chol))
            assert abs(value - independent_value) < 1e-12
            assert abs(np.sqrt(variance) - independent_sd) < 1e-12
            results[weight_name][model] = {"discovery_conditioned_mean_mag": value,
                                          "discovery_conditioned_sd_mag": float(np.sqrt(variance))}
            if weight_name == "config_amplitude_0.3" and model in original["models"]:
                assert abs(value - original["models"][model]["mean_mag"]) < 1e-12
                assert abs(np.sqrt(variance) - original["models"][model]["discovery_conditioned_sd_mag"]) < 1e-12
    paths = [Path(__file__), HERE / "distance-prior-sensitivity-protocol.md", DIST / "responses.npz",
             DIST / "contrast-membership.csv", DIST / "result.json", EXP / "posterior-models.npz",
             HERE / "colorlaw-weight-posterior.npz", HERE / "colorlaw-weight-sensitivity.json"]
    output = {"status": "PASS; fixed-target arithmetic pushforward, no new fits or bin choices",
              "target": "Unweighted same saved high255 minus low255 validation zHEL ranks; pre-BBC local compensation, fixed alpha0.16087/beta3.1178",
              "contrast_reconstruction_max_error": contrast_error,
              "transformed_full_response_vs_contrast_max_error": max(errors),
              "conditional_contrast": results,
              "scope": "These means and SDs are reference-coordinate/prior-dependent finite-mode contributions, not measured biases, total statistical contrast errors or a cosmology posterior.",
              "sha256": {str(p.relative_to(ROOT)): sha(p) for p in paths}}
    (HERE / "distance-prior-sensitivity.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
