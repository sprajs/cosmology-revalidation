#!/usr/bin/env python3
"""Independent saved-score arithmetic check; no SN fitting or model changes."""
from pathlib import Path
import csv
import hashlib
import json

import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "astra_design"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def spd_solve(a, b):
    return cho_solve(cho_factor(a, lower=True), b)


def logdet(a):
    return float(2 * np.log(np.diag(np.linalg.cholesky(a))).sum())


def blocks_from_rows(a, t, y, offsets):
    fs, gs = [], []
    for left, right in zip(offsets[:-1], offsets[1:]):
        x = np.column_stack([a[left:right], t[left:right]])
        fs.append(x.T @ x)
        gs.append(x.T @ y[left:right])
    return np.array(fs), np.array(gs)


def main():
    paths = {
        "discovery": BASE / "shared43/projected-modes.npz",
        "validation": BASE / "shared1020/projected-modes.npz",
        "reported": BASE / "shared1020/result.json",
        "coefficients": BASE / "validation1020/frozen-discovery-coefficients.npz",
        "analysis": BASE / "shared1020_analyze.py",
        "protocol": BASE / "shared1020-protocol.md",
        "sensitivity_amendment": BASE / "shared1020/prior-amendment-manifest.json",
    }
    d = np.load(paths["discovery"], allow_pickle=False)
    v = np.load(paths["validation"], allow_pickle=False)
    reported = json.loads(paths["reported"].read_text())
    assert reported["source_sha256"] == sha(paths["analysis"])
    assert reported["protocol_sha256"] == sha(paths["protocol"])
    assert len(d["CID"]) == 43 and len(v["CID"]) == 1020
    assert not (set(d["CID"]) & set(v["CID"]))
    with (BASE.parents[1] / "salt_dust_audit/flux_response/selected_objects.csv").open() as f:
        original_ids = {r["CID"] for r in csv.DictReader(f)}
    assert len(original_ids) == 64 and not (original_ids & set(v["CID"]))
    ad = np.vstack([d[f"{cid}__calibration_design"] for cid in d["CID"]])
    td = np.vstack([d[f"{cid}__observer_design"] for cid in d["CID"]])
    yd = np.concatenate([d[f"{cid}__r"] for cid in d["CID"]])
    offsets_d = np.r_[0, np.cumsum([len(d[f"{cid}__r"]) for cid in d["CID"]])]
    av, tv, yv = v["calibration_design"], v["observer_design"], v["residual"]
    for arr in (ad, td, yd, av, tv, yv):
        assert np.isfinite(arr).all()
    fd, gd = blocks_from_rows(ad, td, yd, offsets_d)
    fv, gv = blocks_from_rows(av, tv, yv, v["row_offsets"])
    block_errors = {}
    for name, rebuilt, saved in (("discovery_F", fd, d["joint_F"]),
                                 ("discovery_g", gd, d["joint_u"]),
                                 ("validation_F", fv, v["joint_F"]),
                                 ("validation_g", gv, v["joint_u"])):
        block_errors[name] = float(np.max(np.abs(rebuilt - saved)))
        assert np.allclose(rebuilt, saved, rtol=3e-11, atol=1e-7), name
    # Build priors independently; equal total physical griz variance.
    band = np.array([[1, 0, 0], [-1, -1, -1], [0, 1, 0], [0, 0, 1.]])
    tau = 0.02
    c_iso = 2 * tau**2 * (np.eye(3) - np.ones((3, 3)) / 4)
    expected_physical = 2 * tau**2 * (np.eye(4) - np.ones((4, 4)) / 4)
    assert np.allclose(band @ c_iso @ band.T, expected_physical, rtol=0, atol=1e-18)
    assert abs(np.trace(expected_physical) - 6 * tau**2) < 1e-18
    transforms = {
        "calibration_only": np.eye(12)[:, :9],
        "calibration_plus_observer": np.diag(np.r_[np.ones(9), np.full(3, tau)]),
        "calibration_plus_isotropic_griz_observer": np.eye(12),
    }
    transforms["calibration_plus_isotropic_griz_observer"][9:, 9:] = np.linalg.cholesky(c_iso)
    results = {}
    posterior = {}
    xdraw, xvraw = np.column_stack([ad, td]), np.column_stack([av, tv])
    for name, transform in transforms.items():
        xd, xv = xdraw @ transform, xvraw @ transform
        sd = np.eye(xd.shape[1]) + xd.T @ xd
        pd = spd_solve(sd, np.eye(len(sd)))
        md = spd_solve(sd, xd.T @ yd)
        # Direct conditional density in a rank-aware SVD observation basis.
        # Independent from the producer's QR and R joint-evidence differences.
        u, s, vt = np.linalg.svd(xv, full_matrices=False)
        keep = s > s[0] * 1e-12
        u, s, vt = u[:, keep], s[keep], vt[keep]
        response = s[:, None] * vt
        cy = np.eye(len(s)) + response @ pd @ response.T
        projected_y = u.T @ yv
        e = projected_y - response @ md
        gain = float(-0.5 * (e @ spd_solve(cy, e) + logdet(cy)
                            - projected_y @ projected_y))
        expected = reported["models"][name]["joint_validation_predictive_gain_over_zero"]
        assert abs(gain - expected) < 1e-8, (name, gain, expected)
        expected_cov = v[name + "_discovery_posterior_covariance"]
        expected_mean = v[name + "_discovery_posterior_mean"]
        assert np.allclose(pd, expected_cov, rtol=1e-10, atol=1e-11)
        assert np.allclose(md, expected_mean, rtol=1e-10, atol=1e-10)
        results[name] = {
            "svd_direct_conditional_gain": gain,
            "saved_R_difference_gain": expected,
            "absolute_gain_difference": abs(gain - expected),
            "validation_design_rank": int(np.sum(keep)),
            "singular_values": s.tolist(),
            "posterior_mean_max_difference": float(np.max(np.abs(md - expected_mean))),
            "posterior_covariance_max_difference": float(np.max(np.abs(pd - expected_cov))),
        }
        posterior[name] = (md, pd)
    # Direct low-dimensional inverse in the calibration observation span,
    # independently of the producer's Woodbury inversion.
    md, pd = posterior["calibration_only"]
    u, s, vt = np.linalg.svd(av, full_matrices=False)
    response = s[:, None] * vt
    cy = np.eye(len(s)) + response @ pd @ response.T
    inv_correction = spd_solve(cy, np.eye(len(s))) - np.eye(len(s))
    pattern = tv @ np.load(paths["coefficients"], allow_pickle=False)["basis_mean"]
    centered = yv - av @ md
    pv, pe = u.T @ pattern, u.T @ centered
    information = float(pattern @ pattern + pv @ inv_correction @ pv)
    inner = float(pattern @ centered + pv @ inv_correction @ pe)
    z = inner / np.sqrt(information)
    pattern_result = {"information": information, "matched_filter_centered": inner,
                      "score_gain_for_conditional_mean_shift": inner - information / 2,
                      "conditional_Gaussian_Z": z,
                      "conditional_Gaussian_one_sided_tail": float(norm.sf(z))}
    expected = reported["original_frozen_observer_direction_after_calibration_null_conditioning"]
    for key, value in pattern_result.items():
        assert np.isclose(value, expected[key], rtol=1e-8, atol=1e-10), key
    ratios = {
        name: row["svd_direct_conditional_gain"] - results["calibration_only"]["svd_direct_conditional_gain"]
        for name, row in results.items() if name != "calibration_only"
    }
    output = {
        "status": "PASS; independent saved-score arithmetic, no observed-data refits",
        "discovery_objects": len(d["CID"]), "validation_objects": len(v["CID"]),
        "discovery_projected_rows": len(yd), "validation_projected_rows": len(yv),
        "all_64_discovery_ids_excluded_from_validation": True,
        "rebuilt_per_object_sufficient_statistic_max_errors": block_errors,
        "models": results, "conditional_log_ratios_vs_calibration": ratios,
        "fixed_observer_pattern_conditional_check": pattern_result,
        "isotropic_prior": {"coefficient_covariance": c_iso.tolist(),
                            "griz_covariance": expected_physical.tolist(),
                            "equal_total_griz_variance": float(np.trace(expected_physical))},
        "scope": "Same fixed Gaussian nine-mode approximation; does not independently re-export native predictions or calibrate full systematic/selection/adaptive-research significance.",
        "source_sha256": sha(__file__),
        "input_sha256": {str(p.relative_to(BASE.parents[2])): sha(p) for p in paths.values()},
    }
    (HERE / "observed-score-check.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
