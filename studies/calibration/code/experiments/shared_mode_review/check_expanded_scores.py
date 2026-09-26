#!/usr/bin/env python3
"""Native exported-quantity and independent 12/15-mode score checks; no fits."""
from pathlib import Path
import ast
import csv
import hashlib
import json

import numpy as np
from astropy.io import fits
from scipy.linalg import solve_triangular
from scipy.stats import norm

from check_observed_scores import spd_solve, logdet

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
BASE = HERE.parent / "astra_design"
EXP = BASE / "expanded12"


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def native_checks(d, v):
    kcor = ROOT / "phase2/official/inputs/SNDATA_ROOT/kcor/DES/DES-SN5YR/calib_DES-SN5YR_DES.fits.gz"
    mapping = json.loads((EXP / "calspec-unit-map.json").read_text())
    assert sha(kcor) == mapping["kcor_sha256"]
    means, scales = {}, {}
    with fits.open(kcor) as hdus:
        table = hdus["FilterTrans"].data
        wave = table["wavelength (A)"].astype(float)
        for band in "griz":
            trans = table["DES-" + band].astype(float)
            occupied = np.flatnonzero(trans > 1e-6)
            sl = slice(max(0, occupied[0] - 1), min(len(wave), occupied[-1] + 2))
            mean = float(np.dot(wave[sl], trans[sl]) / trans[sl].sum())
            f = np.float32
            delta = f(f(.00714) * f(f(mean) / f(1e4)))
            scale = f(f(10) ** f(-f(.4) * delta))
            assert mean == mapping["bands"][band]["transmission_weighted_mean_A"]
            assert float(scale) == mapping["bands"][band]["flux_scale_float32"]
            assert float(delta) == mapping["bands"][band]["delta_m_float32"]
            means[band] = mean
            scales[band] = float(scale)
    dcache = np.load(BASE / "exact43/comparison/matched-matrices.npz", allow_pickle=False)
    stats = {"native_object_variant_pairs": 0, "native_epoch_variant_pairs": 0,
             "calspec_flux_map_max_absolute_error": 0.,
             "calspec_error_map_max_absolute_error": 0.,
             "calspec_inverse_rounding_max_measurement_sigma": 0.,
             "MWEBV_float32_map_max_absolute_error": 0.,
             "new_projected_column_max_absolute_difference": 0.,
             "calspec_native_model_max_absolute_change": 0.,
             "nominal_projected_residual_max_absolute_difference": 0.}
    for cohort, sample in (("discovery", d), ("validation", v)):
        for index, cid in enumerate(sample["CID"]):
            nominal_dir = "exact43" if cohort == "discovery" else "validation1020"
            nominal = np.load(BASE / nominal_dir / "objectives" / f"objective_{cid}.npz", allow_pickle=False)
            order = np.array(sorted(range(len(nominal["MJD"])),
                                    key=lambda i: (nominal["MJD"][i], nominal["band"][i])))
            if cohort == "discovery":
                covariance = dcache[f"{cid}__exact_covariance"]
                tangent = dcache[f"{cid}__jacobian_flux"].copy()
            else:
                cache = np.load(BASE / "validation1020/analysis/objects" / f"{cid}.npz", allow_pickle=False)
                covariance = cache["exact_covariance"]
                tangent = cache["jacobian_flux"].copy()
            assert np.array_equal(covariance, nominal["frozen_flux_covariance"][np.ix_(order, order)])
            f0 = nominal["model_flux"][order]
            tangent[:, 0] = -.4 * np.log(10) * f0
            chol = np.linalg.cholesky(covariance)
            uw, sw, _ = np.linalg.svd(solve_triangular(chol, tangent, lower=True), full_matrices=True)
            assert np.sum(sw > sw[0] * 1e-10) == 4
            q = uw[:, 4:]
            lo, hi = sample["row_offsets"][index:index+2]
            r = q.T @ solve_triangular(chol, (nominal["data_flux"] - nominal["model_flux"])[order], lower=True)
            stats["nominal_projected_residual_max_absolute_difference"] = max(
                stats["nominal_projected_residual_max_absolute_difference"],
                float(np.max(np.abs(r - sample["residual"][lo:hi]))))
            for column, (name, weight) in enumerate((("CALSPEC", 1.), ("MWEBV", 1.), ("COLORLAW", .3)), 9):
                variant = np.load(EXP / cohort / name / "objectives" / f"objective_{cid}.npz", allow_pickle=False)
                for key in ("MJD", "band", "parameters_x0_x1_c_t0", "zHEL"):
                    assert np.array_equal(variant[key], nominal[key]), (cohort, cid, name, key)
                stats["native_object_variant_pairs"] += 1
                stats["native_epoch_variant_pairs"] += len(order)
                scale = np.ones(len(order))
                if name == "CALSPEC":
                    scale = np.array([scales[b] for b in nominal["band"]])
                    for key, metric in (("data_flux", "calspec_flux_map_max_absolute_error"),
                                        ("data_fluxerr", "calspec_error_map_max_absolute_error")):
                        predicted = (nominal[key].astype(np.float32) * scale.astype(np.float32)).astype(float)
                        error = float(np.max(np.abs(predicted - variant[key])))
                        stats[metric] = max(stats[metric], error)
                        assert np.array_equal(predicted, variant[key]), (cid, key)
                    rounding = np.abs(variant["data_flux"] / scale - nominal["data_flux"]) / nominal["data_fluxerr"]
                    stats["calspec_inverse_rounding_max_measurement_sigma"] = max(
                        stats["calspec_inverse_rounding_max_measurement_sigma"], float(rounding.max()))
                    change = float(np.max(np.abs(variant["model_flux"] - nominal["model_flux"])))
                    stats["calspec_native_model_max_absolute_change"] = max(
                        stats["calspec_native_model_max_absolute_change"], change)
                    assert np.array_equal(variant["MWEBV"], nominal["MWEBV"])
                else:
                    for key in ("data_flux", "data_fluxerr"):
                        assert np.array_equal(variant[key], nominal[key])
                    expected_mw = ((nominal["MWEBV"].astype(np.float32) * np.float32(.95)).astype(float)
                                   if name == "MWEBV" else nominal["MWEBV"])
                    assert np.array_equal(variant["MWEBV"], expected_mw)
                    if name == "MWEBV":
                        stats["MWEBV_float32_map_max_absolute_error"] = max(
                            stats["MWEBV_float32_map_max_absolute_error"],
                            float(np.max(np.abs(variant["MWEBV"] - expected_mw))))
                difference = (variant["model_flux"] / scale)[order] - f0
                rebuilt = weight * q.T @ solve_triangular(chol, difference, lower=True)
                saved = sample["calibration_design"][lo:hi, column]
                error = float(np.max(np.abs(rebuilt - saved)))
                assert np.allclose(rebuilt, saved, rtol=1e-10, atol=1e-10), (cid, name, error)
                stats["new_projected_column_max_absolute_difference"] = max(
                    stats["new_projected_column_max_absolute_difference"], error)
            if index and index % 250 == 0:
                print(f"Independent native map/projection check: {cohort} {index}/{len(sample['CID'])}", flush=True)
    assert stats["nominal_projected_residual_max_absolute_difference"] < 1e-10
    # Execute only the inspected pure label-filter function from archived source.
    grouping_file = ROOT / "phase2/official/build/SNANA-2fe0f56/util/create_covariance.py"
    definition = next(n for n in ast.parse(grouping_file.read_text()).body
                      if isinstance(n, ast.FunctionDef) and n.name == "apply_filter")
    namespace = {}
    exec(compile(ast.Module(body=[definition], type_ignores=[]), str(grouping_file), "exec"), namespace)
    labels = [f"cal_{j}" for j in range(1, 10)] + ["CALSPEC", "MWEBV", "COLORLAW"]
    grouped = [label for label in labels if namespace["apply_filter"](label, "+cal")]
    assert grouped == labels[:10]
    stats.update({"KCOR_transmission_weighted_mean_A": means, "source_float32_scales": scales,
                  "archived_calibration_distance_group_members": grouped,
                  "expanded_unique_native_modes": labels,
                  "no_distance_group_matrix_added": True})
    return stats


def main():
    d = np.load(EXP / "discovery-projected-modes.npz", allow_pickle=False)
    v = np.load(EXP / "validation-projected-modes.npz", allow_pickle=False)
    oldd = np.load(BASE / "shared43/projected-modes.npz", allow_pickle=False)
    oldv = np.load(BASE / "shared1020/projected-modes.npz", allow_pickle=False)
    reported = json.loads((EXP / "result.json").read_text())
    assert reported["source_sha256"] == sha(BASE / "expanded12_analyze.py")
    assert reported["protocol_sha256"] == sha(BASE / "expanded12-protocol.md")
    assert reported["mode_amplitude_scales"] == [.3] * 9 + [1., 1., .3]
    assert len(d["CID"]) == 43 and len(v["CID"]) == 1020
    with (ROOT / "runs/salt_dust_audit/flux_response/selected_objects.csv").open() as file:
        excluded = {row["CID"] for row in csv.DictReader(file)}
    assert len(excluded) == 64 and not (excluded & set(v["CID"]))
    assert np.array_equal(d["CID"], oldd["CID"]) and np.array_equal(v["CID"], oldv["CID"])
    assert np.array_equal(d["calibration_design"][:, :9], np.vstack([oldd[f"{cid}__calibration_design"] for cid in d["CID"]]))
    assert np.array_equal(v["calibration_design"][:, :9], oldv["calibration_design"])
    per_object_errors = {}
    for label, sample in (("discovery", d), ("validation", v)):
        x = np.column_stack([sample["calibration_design"], sample["observer_design"]])
        y = sample["residual"]
        assert np.isfinite(x).all() and np.isfinite(y).all()
        fs, gs = [], []
        for left, right in zip(sample["row_offsets"][:-1], sample["row_offsets"][1:]):
            fs.append(x[left:right].T @ x[left:right])
            gs.append(x[left:right].T @ y[left:right])
        for name, actual, target in (("F", np.array(fs), sample["joint_F"]), ("g", np.array(gs), sample["joint_u"])):
            error = float(np.max(np.abs(actual - target)))
            assert np.allclose(actual, target, rtol=1e-10, atol=1e-8)
            per_object_errors[label + "_" + name] = error
    native = native_checks(d, v)
    ciso = 2 * .02**2 * (np.eye(3) - np.ones((3, 3)) / 4)
    transforms = {"systematics_only": np.eye(15)[:, :12],
                  "systematics_plus_observer": np.diag(np.r_[np.ones(12), np.full(3, .02)]),
                  "systematics_plus_isotropic_griz_observer": np.eye(15)}
    transforms["systematics_plus_isotropic_griz_observer"][12:, 12:] = np.linalg.cholesky(ciso)
    xd0 = np.column_stack([d["calibration_design"], d["observer_design"]])
    xv0 = np.column_stack([v["calibration_design"], v["observer_design"]])
    savedpost = np.load(EXP / "posterior-models.npz", allow_pickle=False)
    scores, post = {}, {}
    for name, transform in transforms.items():
        xd, xv = xd0 @ transform, xv0 @ transform
        sd = np.eye(xd.shape[1]) + xd.T @ xd
        pd = spd_solve(sd, np.eye(len(sd)))
        md = spd_solve(sd, xd.T @ d["residual"])
        u, singular, vt = np.linalg.svd(xv, full_matrices=False)
        keep = singular > singular[0] * 1e-12
        u, singular, vt = u[:, keep], singular[keep], vt[keep]
        response = singular[:, None] * vt
        cy = np.eye(len(singular)) + response @ pd @ response.T
        yp = u.T @ v["residual"]
        e = yp - response @ md
        gain = float(-.5 * (e @ spd_solve(cy, e) + logdet(cy) - yp @ yp))
        expected = reported["models"][name]["joint_validation_predictive_gain_over_zero"]
        assert abs(gain - expected) < 1e-8, (name, gain, expected)
        mean_error = float(np.max(np.abs(md - savedpost[name + "_discovery_posterior_mean"])))
        cov_error = float(np.max(np.abs(pd - savedpost[name + "_discovery_posterior_covariance"])))
        assert max(mean_error, cov_error) < 1e-10
        scores[name] = {"direct_SVD_conditional_gain": gain, "saved_gain": expected,
                        "absolute_difference": abs(gain - expected), "validation_design_rank": int(keep.sum()),
                        "posterior_mean_max_difference": mean_error, "posterior_covariance_max_difference": cov_error}
        post[name] = md, pd
    md, pd = post["systematics_only"]
    a = v["calibration_design"]
    u, singular, vt = np.linalg.svd(a, full_matrices=False)
    response = singular[:, None] * vt
    cov = np.eye(len(singular)) + response @ pd @ response.T
    correction = spd_solve(cov, np.eye(len(singular))) - np.eye(len(singular))
    c = np.load(BASE / "validation1020/frozen-discovery-coefficients.npz", allow_pickle=False)["basis_mean"]
    pattern = v["observer_design"] @ c
    residual = v["residual"] - a @ md
    vp, rp = u.T @ pattern, u.T @ residual
    information = float(pattern @ pattern + vp @ correction @ vp)
    score = float(pattern @ residual + vp @ correction @ rp)
    z = score / np.sqrt(information)
    fixed = {"information": information, "centered_matched_product": score,
             "conditional_Gaussian_Z": z, "conditional_Gaussian_one_sided_tail": float(norm.sf(z)),
             "full_fixed_direction_shift_gain": score - information / 2}
    for key, value in fixed.items():
        assert np.isclose(value, reported["original_fixed_observer_after_systematic_null_conditioning"][key], rtol=1e-8, atol=1e-10), key
    paths = [Path(__file__), HERE / "check_observed_scores.py", EXP / "result.json", EXP / "input-manifest.json",
             EXP / "discovery-projected-modes.npz", EXP / "validation-projected-modes.npz", EXP / "posterior-models.npz",
             EXP / "calspec-unit-map.json", BASE / "expanded12_analyze.py", BASE / "expanded12-protocol.md"]
    result = {"status": "PASS; independent native export/unit/projection and 12/15-mode score check, no fits",
              "native_checks": native, "rebuilt_per_object_Gram_score_errors": per_object_errors,
              "models": scores,
              "observer_log_ratio": scores["systematics_plus_observer"]["direct_SVD_conditional_gain"] - scores["systematics_only"]["direct_SVD_conditional_gain"],
              "isotropic_observer_log_ratio": scores["systematics_plus_isotropic_griz_observer"]["direct_SVD_conditional_gain"] - scores["systematics_only"]["direct_SVD_conditional_gain"],
              "fixed_observer_conditional_check": fixed,
              "scope": "Expanded finite Gaussian modes with nominal local C/J, shared discovery conditioning; neither a full calibration prior nor full systematic/selection/cosmology validation.",
              "sha256": {str(p.relative_to(ROOT)): sha(p) for p in paths}}
    (HERE / "expanded-score-check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
