"""Independent arithmetic verification of the frozen validation1020 analysis.

Reads saved sufficient statistics and frozen coefficients; does not refit objects.
"""
from __future__ import annotations

import csv
import hashlib
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "runs/research_2026_09_26/astra_design/validation1020"
ANALYSIS = BASE / "analysis"
OUT = ROOT / "runs/research_2026_09_26/validation1020_independent_verification"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dense_shared_log_gain(u: np.ndarray, F: np.ndarray, c: np.ndarray, V: np.ndarray) -> float:
    """Log E_{beta~N(c,V)} exp(u.beta - beta.F.beta/2), dense 3x3 form."""
    P = np.linalg.inv(V)
    A = P + F
    b = P @ c + u
    sign_v, logdet_v = np.linalg.slogdet(V)
    sign_a, logdet_a = np.linalg.slogdet(A)
    assert sign_v == sign_a == 1
    return float((b @ np.linalg.solve(A, b) - c @ P @ c - logdet_v - logdet_a) / 2)


def main() -> None:
    assert not OUT.exists(), "preserve verification output"
    result = json.loads((ANALYSIS / "result.json").read_text())
    manifest = json.loads((BASE / "input-manifest.json").read_text())
    rest_manifest = json.loads((BASE / "rest-amendment-manifest.json").read_text())
    executed_source = ROOT / "runs/research_2026_09_26/astra_design/validation1020.py"
    assert sha(executed_source) == result["source_sha256"]
    observer = np.load(BASE / "frozen-discovery-coefficients.npz")
    rest = np.load(BASE / "frozen-rest-discovery-coefficients.npz")
    c, V = observer["basis_mean"], observer["basis_covariance"]
    cr, Vr = rest["basis_mean"], rest["basis_covariance"]
    assert sha(BASE / "frozen-discovery-coefficients.npz") == manifest["sha256"]["runs/research_2026_09_26/astra_design/validation1020/frozen-discovery-coefficients.npz"]
    assert sha(BASE / "frozen-discovery-coefficients.npz") == result["inputs_sha256"]["runs/research_2026_09_26/astra_design/validation1020/frozen-discovery-coefficients.npz"]
    assert sha(BASE / "frozen-rest-discovery-coefficients.npz") == rest_manifest["rest_coefficient_sha256"]
    assert np.array_equal(c, np.array(result["discovery_coefficient_mean"]))
    assert np.array_equal(V, np.array(result["discovery_coefficient_covariance"]))
    ridge = float(observer["ridge_sigma"])
    P = observer["discovery_F"] + np.eye(3) / ridge**2
    assert np.max(np.abs(np.linalg.inv(P) - V)) < 1e-15
    assert np.max(np.abs(np.linalg.solve(P, observer["discovery_u"]) - c)) < 1e-14
    assert np.array_equal(cr, np.array(result["frozen_rest_coefficient_mean"]))
    assert np.array_equal(Vr, np.array(result["frozen_rest_coefficient_covariance"]))
    with (ANALYSIS / "object-scores.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    with (ANALYSIS / "field-scores.csv").open(newline="") as handle:
        field_rows = list(csv.DictReader(handle))
    saved = np.load(ANALYSIS / "sufficient-arrays.npz")
    checks = {}
    global_max_score_error = 0.0
    for arm in ("published_mask", "first_exact_duplicate"):
        rr = [row for row in rows if row["arm"] == arm]
        ff = [row for row in field_rows if row["arm"] == arm]
        assert len(rr) == 1020 and len(ff) == 10
        assert len({row["CID"] for row in rr}) == 1020
        u6, F6 = saved[f"{arm}_u"], saved[f"{arm}_F"]
        assert u6.shape == (1020, 6) and F6.shape == (1020, 6, 6)
        a = u6[:, :3] @ c
        I = np.einsum("i,nij,j->n", c, F6[:, :3, :3], c)
        ar = u6[:, 3:] @ cr
        Ir = np.einsum("i,nij,j->n", cr, F6[:, 3:, 3:], cr)
        score_errors = [
            np.max(np.abs(a - np.array([float(row["matched_filter"]) for row in rr]))),
            np.max(np.abs(I - np.array([float(row["information"]) for row in rr]))),
            np.max(np.abs(a - I / 2 - np.array([float(row["fixed_prediction_gain"]) for row in rr]))),
            np.max(np.abs(ar - Ir / 2 - np.array([float(row["rest_fixed_prediction_gain"]) for row in rr]))),
        ]
        max_score_error = max(map(float, score_errors))
        global_max_score_error = max(global_max_score_error, max_score_error)
        fa = float(np.sum(a)); fi = float(np.sum(I)); gain = fa - fi / 2
        fra = float(np.sum(ar)); fri = float(np.sum(Ir)); rest_gain = fra - fri / 2
        # Group fields from object rows rather than reading their reported totals.
        groups = defaultdict(list)
        for idx, row in enumerate(rr):
            groups[row["field"]].append(idx)
        fields = sorted(groups)
        field_a = np.array([np.sum(a[groups[f]]) for f in fields])
        field_ar = np.array([np.sum(ar[groups[f]]) for f in fields])
        sign_matrix = np.array(list(itertools.product((-1.0, 1.0), repeat=10)))
        signed = sign_matrix @ field_a
        rest_signed = sign_matrix @ field_ar
        paired_signed = signed - rest_signed
        fixed_exceed = int(np.count_nonzero(signed >= fa - 1e-10))
        rest_exceed = int(np.count_nonzero(rest_signed >= fra - 1e-10))
        paired_exceed = int(np.count_nonzero(paired_signed >= fa - fra - 1e-10))
        u = u6[:, :3].sum(axis=0)
        F = F6[:, :3, :3].sum(axis=0)
        ur = u6[:, 3:].sum(axis=0)
        Fr = F6[:, 3:, 3:].sum(axis=0)
        joint = dense_shared_log_gain(u, F, c, V)
        rest_joint = dense_shared_log_gain(ur, Fr, cr, Vr)
        field_u = np.array([u6[groups[f], :3].sum(axis=0) for f in fields])
        signed_joint = np.array([dense_shared_log_gain(s @ field_u, F, c, V) for s in sign_matrix])
        joint_exceed = int(np.count_nonzero(signed_joint >= joint - 1e-10))
        reported = result["arms"][arm]
        assert abs(gain - reported["fixed_template_gain"]) < 1e-9
        assert abs(rest_gain - reported["frozen_rest_secondary"]["fixed_gain"]) < 1e-9
        assert abs(joint - reported["joint_discovery_posterior_predictive_gain"]) < 1e-8
        assert abs(rest_joint - reported["frozen_rest_secondary"]["joint_discovery_posterior_predictive_gain"]) < 1e-8
        assert fixed_exceed == reported["field_sign_exceed_count"]
        assert fixed_exceed / 1024 == reported["field_sign_exact_tail"]
        assert rest_exceed / 1024 == reported["frozen_rest_secondary"]["field_sign_exact_tail"]
        assert paired_exceed / 1024 == reported["frozen_rest_secondary"]["paired_field_sign_exact_tail_about_zero_residual"]
        assert joint_exceed / 1024 == reported["secondary_joint_predictive_field_sign_tail"]
        assert np.max(np.abs(signed - saved[f"{arm}_field_sign_matched_filter"])) < 1e-9
        assert np.max(np.abs(signed_joint - saved[f"{arm}_field_sign_joint_predictive"])) < 1e-8
        assert max(abs(float(row["matched_filter"]) - field_a[fields.index(row["field"])]) for row in ff) < 1e-9
        checks[arm] = {
            "objects": len(rr), "fields": len(fields), "epochs": sum(int(row["epoch_count"]) for row in rr),
            "fixed_observer_gain": gain, "fixed_rest_gain": rest_gain,
            "observer_matched_filter": fa, "observer_information": fi,
            "joint_observer_log_gain_dense": joint, "joint_rest_log_gain_dense": rest_joint,
            "field_sign_assignments": len(sign_matrix),
            "fixed_observer_sign_exceed": fixed_exceed, "fixed_rest_sign_exceed": rest_exceed,
            "paired_sign_exceed": paired_exceed, "joint_observer_sign_exceed": joint_exceed,
            "max_object_score_error": max_score_error,
            "max_saved_joint_sign_gain_error": float(np.max(np.abs(signed_joint - saved[f"{arm}_field_sign_joint_predictive"]))),
            "duplicate_objects": sum(int(row["exact_duplicate_groups"]) > 0 for row in rr),
            "duplicate_groups": sum(int(row["exact_duplicate_groups"]) for row in rr),
            "extra_duplicate_epochs": sum(int(row["exact_duplicate_extra_epochs"]) for row in rr),
        }
    # Independently recount exact duplicate tuples in all saved object exports.
    duplicate_objects = duplicate_groups = duplicate_extra = 0
    retained_epochs = published_epochs = 0
    for path in sorted((ANALYSIS / "objects").glob("*.npz")):
        obj = np.load(path)
        keys = list(zip(obj["MJD"], obj["band"], obj["observed_flux"], obj["quoted_error"]))
        counts = Counter(keys)
        groups = sum(n > 1 for n in counts.values())
        extra = sum(n - 1 for n in counts.values())
        duplicate_objects += groups > 0
        duplicate_groups += groups
        duplicate_extra += extra
        published_epochs += len(keys)
        retained_epochs += len(obj["first_duplicate_keep"])
        assert extra == len(keys) - len(obj["first_duplicate_keep"])
    assert (duplicate_objects, duplicate_groups, duplicate_extra) == (37, 59, 147)
    assert (published_epochs, retained_epochs) == (39606, 39459)
    assert global_max_score_error < 1e-9
    out = {
        "scope": "Independent arithmetic check of saved conditional validation1020 scores; no light-curve refit or shared calibration marginalization.",
        "frozen_prediction_verified": True,
        "frozen_observer_coefficient_sha256": sha(BASE / "frozen-discovery-coefficients.npz"),
        "frozen_rest_coefficient_sha256": sha(BASE / "frozen-rest-discovery-coefficients.npz"),
        "executed_analysis_source_sha256": sha(executed_source),
        "no_validation_coefficient_fit_in_primary_score": True,
        "dense_integral": "A=V^{-1}+F; b=V^{-1}c+u; log gain=(b^T A^{-1} b-c^T V^{-1}c-log det V-log det A)/2. A single 3x3 integration uses summed u,F over all objects.",
        "arms": checks,
        "duplicate_recount_from_saved_epochs": {"objects": duplicate_objects, "groups": duplicate_groups, "extra_epochs": duplicate_extra, "published_epochs": published_epochs, "first_exact_duplicate_epochs": retained_epochs},
        "qualification": "Conditional Gaussian and field-sign calibrations omit shared cross-object calibration/SALT3 covariance and do not express a full-systematics significance.",
        "source_sha256": sha(Path(__file__)),
        "input_sha256": {str(p.relative_to(ROOT)): sha(p) for p in [ANALYSIS / "result.json", ANALYSIS / "object-scores.csv", ANALYSIS / "field-scores.csv", ANALYSIS / "sufficient-arrays.npz", BASE / "input-manifest.json", BASE / "rest-amendment-manifest.json", BASE / "frozen-discovery-coefficients.npz", BASE / "frozen-rest-discovery-coefficients.npz", executed_source]},
    }
    OUT.mkdir(parents=True)
    (OUT / "verification.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({"arms": checks, "duplicate_recount": out["duplicate_recount_from_saved_epochs"]}, indent=2))


if __name__ == "__main__":
    main()
