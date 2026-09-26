#!/usr/bin/env python3
"""Frozen paper-versus-config COLORLAW amplitude sensitivity; saved arrays only."""
from pathlib import Path
import hashlib
import json

import numpy as np
from scipy.stats import norm

from check_observed_scores import spd_solve, logdet

HERE = Path(__file__).resolve().parent
BASE = HERE.parent / "astra_design"
ROOT = HERE.parents[2]
EXP = BASE / "expanded12"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def r_stat(x, y):
    s = np.eye(x.shape[1]) + x.T @ x
    g = x.T @ y
    return float(.5 * (g @ spd_solve(s, g) - logdet(s)))


def main():
    out = HERE / "colorlaw-weight-sensitivity.json"
    assert not out.exists(), "Preserve completed sensitivity"
    paths = [Path(__file__), HERE / "colorlaw-weight-protocol.md", HERE / "check_observed_scores.py",
             EXP / "discovery-projected-modes.npz", EXP / "validation-projected-modes.npz",
             EXP / "result.json", BASE / "validation1020/frozen-discovery-coefficients.npz",
             ROOT / "papers/text/2401.02945v2.txt",
             ROOT / "sources/repos/des-science__DES-SN5YR@1.3/7_PIPPIN_FILES/base_files/lcfit/fitopts.yml",
             ROOT / "phase2/official/build/SNANA-2fe0f56/util/create_covariance.py"]
    inputs = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    d = np.load(EXP / "discovery-projected-modes.npz", allow_pickle=False)
    v = np.load(EXP / "validation-projected-modes.npz", allow_pickle=False)
    original = json.loads((EXP / "result.json").read_text())
    assert original["mode_order"][11] == "COLORLAW" and original["mode_amplitude_scales"][11] == .3
    multiplier = np.sqrt(1 / 3) / .3
    ad, av = d["calibration_design"].copy(), v["calibration_design"].copy()
    ad[:, 11] *= multiplier
    av[:, 11] *= multiplier
    assert np.array_equal(ad[:, :11], d["calibration_design"][:, :11])
    assert np.array_equal(av[:, :11], v["calibration_design"][:, :11])
    xd0, xv0 = np.column_stack([ad, d["observer_design"]]), np.column_stack([av, v["observer_design"]])
    ciso = 2 * .02**2 * (np.eye(3) - np.ones((3, 3)) / 4)
    transforms = {"systematics_only": np.eye(15)[:, :12],
                  "systematics_plus_observer": np.diag(np.r_[np.ones(12), np.full(3, .02)]),
                  "systematics_plus_isotropic_griz_observer": np.eye(15)}
    transforms["systematics_plus_isotropic_griz_observer"][12:, 12:] = np.linalg.cholesky(ciso)
    scores, post, posterior_output = {}, {}, {}
    for name, transform in transforms.items():
        xd, xv = xd0 @ transform, xv0 @ transform
        sd = np.eye(xd.shape[1]) + xd.T @ xd
        pd = spd_solve(sd, np.eye(len(sd)))
        md = spd_solve(sd, xd.T @ d["residual"])
        u, s, vt = np.linalg.svd(xv, full_matrices=False)
        keep = s > s[0] * 1e-12
        u, s, vt = u[:, keep], s[keep], vt[keep]
        response = s[:, None] * vt
        covariance = np.eye(len(s)) + response @ pd @ response.T
        yp = u.T @ v["residual"]
        e = yp - response @ md
        gain = float(-.5 * (e @ spd_solve(covariance, e) + logdet(covariance) - yp @ yp))
        independent = r_stat(np.vstack([xd, xv]), np.r_[d["residual"], v["residual"]]) - r_stat(xd, d["residual"])
        assert abs(gain - independent) < 1e-8
        old = original["models"][name]["joint_validation_predictive_gain_over_zero"]
        scores[name] = {"SVD_conditional_gain": gain, "R_difference_gain": independent,
                        "independent_absolute_difference": abs(gain - independent),
                        "original_gain": old, "change_from_original": gain - old,
                        "validation_design_rank": int(keep.sum())}
        post[name] = md, pd
        posterior_output[name + "_mean"] = md
        posterior_output[name + "_covariance"] = pd
    md, pd = post["systematics_only"]
    u, s, vt = np.linalg.svd(av, full_matrices=False)
    response = s[:, None] * vt
    cov = np.eye(len(s)) + response @ pd @ response.T
    correction = spd_solve(cov, np.eye(len(s))) - np.eye(len(s))
    c = np.load(BASE / "validation1020/frozen-discovery-coefficients.npz", allow_pickle=False)["basis_mean"]
    pattern = v["observer_design"] @ c
    residual = v["residual"] - av @ md
    vp, ep = u.T @ pattern, u.T @ residual
    info = float(pattern @ pattern + vp @ correction @ vp)
    product = float(pattern @ residual + vp @ correction @ ep)
    z = product / np.sqrt(info)
    fixed = {"information": info, "centered_matched_product": product, "conditional_Gaussian_Z": z,
             "conditional_Gaussian_one_sided_tail": float(norm.sf(z)),
             "full_fixed_direction_shift_gain": product - .5 * info}
    fixed_changes = {key: value - original["original_fixed_observer_after_systematic_null_conditioning"][key]
                     for key, value in fixed.items()}
    ratios = {name: row["SVD_conditional_gain"] - scores["systematics_only"]["SVD_conditional_gain"]
              for name, row in scores.items() if name != "systematics_only"}
    result = {"status": "PASS; separately frozen post-result sensitivity, no native fits or validation coefficient tuning",
              "COLORLAW_amplitude": float(np.sqrt(1 / 3)), "COLORLAW_covariance_factor": 1 / 3,
              "original_amplitude": .3, "original_covariance_factor": .09,
              "column_multiplier": float(multiplier), "unchanged_first_eleven_mode_columns": True,
              "discovery_objects": len(d["CID"]), "validation_objects": len(v["CID"]),
              "models": scores, "observer_log_ratios": ratios,
              "observer_log_ratio_changes": {
                  "inherited": ratios["systematics_plus_observer"] - original["conditional_predictive_log_ratio_extra_observer"],
                  "isotropic": ratios["systematics_plus_isotropic_griz_observer"] - original["isotropic_prior_conditional_predictive_log_ratio_extra_observer"]},
              "fixed_original_observer": fixed, "fixed_original_observer_changes": fixed_changes,
              "scope": "This bounds one paper/config weighting sensitivity in a finite-mode model; it does not identify the original executed weight, calibrate the full systematic budget, or measure a correction.",
              "inputs_sha256": inputs}
    np.savez_compressed(HERE / "colorlaw-weight-posterior.npz", **posterior_output)
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
