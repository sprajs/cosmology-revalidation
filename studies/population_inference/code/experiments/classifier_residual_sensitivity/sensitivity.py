"""Declared post-validation intersection-cut sensitivity; no fit or inference."""
from __future__ import annotations

import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
ASTRA = ROOT / "runs/research_2026_09_26/astra_design"
VAL = ASTRA / "validation1020"
ANALYSIS = VAL / "analysis"
CLF = ROOT / "phase2/classification/reconstruction_20260926"
INPUTS = [
    OUT / "protocol.md",
    VAL / "cohort.csv",
    VAL / "cids.txt",
    VAL / "frozen-discovery-coefficients.npz",
    ANALYSIS / "object-scores.csv",
    ANALYSIS / "sufficient-arrays.npz",
    ANALYSIS / "result.json",
    ASTRA / "exact43/cids.txt",
    CLF / "des_clump_diagnostic.csv",
    CLF / "official_validate_comparison.csv",
]
FIELDS = ("C1", "C2", "C3", "E1", "E2", "S1", "S2", "X1", "X2", "X3")
THRESHOLD = 0.999


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def unique(rows: list[dict[str, str]], key: str) -> dict[str, dict[str, str]]:
    result = {row[key]: row for row in rows}
    assert len(result) == len(rows), f"duplicate {key}"
    return result


def totals(rows: list[dict]) -> dict:
    m = float(np.sum([r["matched_filter"] for r in rows], dtype=np.float64))
    i = float(np.sum([r["information"] for r in rows], dtype=np.float64))
    g = float(np.sum([r["fixed_prediction_gain"] for r in rows], dtype=np.float64))
    assert abs(g - (m - i / 2)) < 1e-9
    return {"objects": len(rows), "epochs": sum(r["epoch_count"] for r in rows),
            "matched_filter": m, "information": i, "fixed_prediction_gain": g}


def main() -> None:
    cohort = read_csv(VAL / "cohort.csv")
    scores = [r for r in read_csv(ANALYSIS / "object-scores.csv")
              if r["arm"] == "published_mask"]
    clump = unique(read_csv(CLF / "des_clump_diagnostic.csv"), "CID")
    official = unique(read_csv(CLF / "official_validate_comparison.csv"), "CID")
    ids = (VAL / "cids.txt").read_text().splitlines()
    discovery = (ASTRA / "exact43/cids.txt").read_text().splitlines()
    assert len(ids) == len(set(ids)) == len(cohort) == len(scores) == 1020
    assert len(discovery) == len(set(discovery)) == 43
    assert not set(ids) & set(discovery)
    assert ids == [r["CID"] for r in cohort] == [r["CID"] for r in scores]
    assert set(ids) | set(discovery) <= clump.keys() & official.keys()
    assert all(float(r["PROB_SNNV19"]) > THRESHOLD for r in cohort)
    assert all((float(clump[cid]["pIa"]) > THRESHOLD) ==
               (float(official[cid]["official_validate_pIa"]) > THRESHOLD)
               for cid in set(ids) | set(discovery))
    saved_result = json.loads((ANALYSIS / "result.json").read_text())
    with np.load(VAL / "frozen-discovery-coefficients.npz", allow_pickle=False) as z:
        coefficient = z["basis_mean"].copy()
    with np.load(ANALYSIS / "sufficient-arrays.npz", allow_pickle=False) as z:
        u = z["published_mask_u"][:, :3].copy()
        F = z["published_mask_F"][:, :3, :3].copy()
    assert u.shape == (1020, 3) and F.shape == (1020, 3, 3)
    assert np.array_equal(coefficient, np.array(saved_result["discovery_coefficient_mean"]))
    independently_matched = u @ coefficient
    independently_information = np.einsum("i,nij,j->n", coefficient, F, coefficient)

    ledger = []
    for idx, (co, score) in enumerate(zip(cohort, scores)):
        cid = score["CID"]
        p_release = float(co["PROB_SNNV19"])
        p_reconstructed = float(clump[cid]["pIa"])
        matched = float(score["matched_filter"])
        information = float(score["information"])
        gain = float(score["fixed_prediction_gain"])
        assert abs(matched - independently_matched[idx]) < 1e-9
        assert abs(information - independently_information[idx]) < 1e-9
        assert abs(gain - (matched - information / 2)) < 1e-9
        assert score["field"] == co["field"] in FIELDS
        ledger.append({"CID": cid, "field": score["field"],
                       "released_PROB_SNNV19": p_release,
                       "reconstructed_pIa": p_reconstructed,
                       "retained_both_gt_0p999": p_release > THRESHOLD and p_reconstructed > THRESHOLD,
                       "epoch_count": int(score["epoch_count"]),
                       "matched_filter": matched, "information": information,
                       "fixed_prediction_gain": gain})
    retained = [r for r in ledger if r["retained_both_gt_0p999"]]
    removed = [r for r in ledger if not r["retained_both_gt_0p999"]]
    full = totals(ledger)
    both = totals(retained)
    discarded = totals(removed)
    reported = saved_result["arms"]["published_mask"]
    assert abs(full["matched_filter"] - reported["matched_filter"]) < 1e-9
    assert abs(full["information"] - reported["fixed_prediction_information"]) < 1e-9
    assert abs(full["fixed_prediction_gain"] - reported["fixed_template_gain"]) < 1e-9
    assert abs(full["matched_filter"] - both["matched_filter"] - discarded["matched_filter"]) < 1e-9
    assert abs(full["information"] - both["information"] - discarded["information"]) < 1e-9
    assert abs(full["fixed_prediction_gain"] - both["fixed_prediction_gain"] - discarded["fixed_prediction_gain"]) < 1e-9
    grouped = defaultdict(list)
    for row in ledger:
        grouped[row["field"]].append(row)
    field_rows = []
    for field in FIELDS:
        all_rows = grouped[field]
        keep = [r for r in all_rows if r["retained_both_gt_0p999"]]
        drop = [r for r in all_rows if not r["retained_both_gt_0p999"]]
        field_rows.append({"field": field, "primary": totals(all_rows),
                           "intersection": totals(keep), "removed": totals(drop),
                           "removed_CIDs": [r["CID"] for r in drop]})
    assert sum(x["intersection"]["objects"] for x in field_rows) == len(retained)
    discovery_disagree = [cid for cid in discovery
                          if float(clump[cid]["pIa"]) <= THRESHOLD]
    assert all(float(official[cid]["PROB_SNNV19"]) > THRESHOLD for cid in discovery)
    output = {
        "scope": "Declared post-validation intersection sensitivity, published-mask scores only; no coefficient/template refit",
        "threshold": ">0.999 strict on both released and reconstructed probabilities",
        "reconstructed_source": rel(CLF / "des_clump_diagnostic.csv"),
        "official_validator_threshold_membership_agrees": True,
        "primary": full, "intersection": both, "removed_contribution": discarded,
        "delta_intersection_minus_primary": {
            key: both[key] - full[key]
            for key in ("objects", "epochs", "matched_filter", "information", "fixed_prediction_gain")},
        "retained_CIDs_in_frozen_order": [r["CID"] for r in retained],
        "removed_CIDs_in_frozen_order": [r["CID"] for r in removed],
        "fields": field_rows,
        "discovery43": {"objects": len(discovery), "disagree": len(discovery_disagree),
                        "disagree_CIDs": discovery_disagree,
                        "unchanged_frozen_coefficient_sha256": sha(VAL / "frozen-discovery-coefficients.npz")},
        "limitations": "Post-validation selection on same released high-probability cohort; neither classifier is ground truth. Fixed conditional residual score omits global calibration/SALT3 marginalization and selection regeneration. No new significance calculation.",
    }
    ledger_path = OUT / "membership_and_scores.csv"
    with ledger_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(ledger[0]))
        writer.writeheader()
        writer.writerows(ledger)
    result_path = OUT / "result.json"
    result_path.write_text(json.dumps(output, indent=2) + "\n")
    manifest = {"schema": "classifier_residual_sensitivity_manifest_v1",
                "generator": rel(Path(__file__).resolve()),
                "generator_sha256": sha(Path(__file__).resolve()),
                "inputs_sha256": {rel(p): sha(p) for p in INPUTS},
                "outputs_sha256": {p.name: sha(p) for p in (result_path, ledger_path)}}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({"primary": full, "intersection": both,
                      "removed_CIDs": output["removed_CIDs_in_frozen_order"],
                      "discovery43_disagree": discovery_disagree}, indent=2))


if __name__ == "__main__":
    main()
