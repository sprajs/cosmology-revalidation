"""Recover frozen Stage B primary summary from already saved CSV ledgers.

No FITS image is opened and no aperture extraction is repeated. The original
run.py, failure log and CSVs remain immutable. A separate independent checker
verifies the recovered statistics using its own arithmetic.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np

import run  # frozen numerical summary routine; imported without stage execution

HERE = Path(__file__).resolve().parent


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(name: str) -> list[dict]:
    with (HERE / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def plain(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [plain(v) for v in value]
    return value


def main():
    protocol = json.loads((HERE / "recovery-protocol.json").read_text())
    for record in protocol["inputs"]:
        p = Path(record["path"])
        if sha(p) != record["sha256"] or p.stat().st_size != record["bytes"]:
            raise ValueError("recovery frozen input mismatch: " + str(p))
    assert protocol["recovery_script_sha256"] == sha(Path(__file__))
    state = json.loads((HERE / "stage_a.json").read_text())
    design = json.loads((HERE.parent / "pixel_design" / "protocol.json").read_text())
    release = HERE / "root-stage-b-release.json"
    assert json.loads(release.read_text())["stage_a_sha256"] == sha(HERE / "stage_a.json")
    exposure, pairs, template = [rows(name) for name in
                                 ("signed_exposure_ledger.csv", "signed_pair_ledger.csv",
                                  "common_template_ledger.csv")]
    cohorts = {"strict": [r for r in state["rows"] if r.get("strict")],
               "secondary": [r for r in state["rows"] if r.get("secondary")]}
    expected = sum(len(cohort) * 8 for cohort in cohorts.values() if len(cohort) >= 30)
    assert len(exposure) == expected
    assert len(pairs) == expected // 4
    assert len(template) == expected // 8
    assert len({(r["branch"], r["candidate"], r["visit"], r["ordinal"]) for r in exposure}) == len(exposure)
    assert len({(r["branch"], r["candidate"], r["visit"]) for r in pairs}) == len(pairs)
    output = {"status": "stage_b_primary_recovered_from_saved_ledgers_after_json_type_failure",
              "no_fits_reopen_or_aperture_reextraction": True,
              "failed_stage_b_log_sha256": sha(HERE / "stage-b-execution.log"),
              "recovery_protocol_sha256": sha(HERE / "recovery-protocol.json"),
              "release_sha256": sha(release),
              "ledger_sha256": {name: sha(HERE / name) for name in
                               ("signed_exposure_ledger.csv", "signed_pair_ledger.csv",
                                "common_template_ledger.csv")},
              "branches": {}}
    for branch, cohort in cohorts.items():
        tile_counts = {str(t): sum(c["tile"] == t for c in cohort)
                       for t in sorted({c["tile"] for c in cohort})}
        entry = {"support_n": len(cohort), "support_tiles": sorted({c["tile"] for c in cohort}),
                 "tile_counts": tile_counts,
                 "primary_adequacy": ("full" if len(cohort) >= 100 and len(tile_counts) >= 12
                                      and all(n >= 3 for n in tile_counts.values()) else
                                      "limited" if len(cohort) >= 30 else "stop_below_30")}
        if len(cohort) < 30:
            output["branches"][branch] = entry
            continue
        visits = {}
        for visit in ("search", "template"):
            chosen = [r for r in pairs if r["branch"] == branch and r["visit"] == visit]
            assert len(chosen) == len(cohort)
            parsed = [{k: (int(r[k]) if k in ("candidate", "tile") else float(r[k]))
                       for k in ("candidate", "tile", "east", "north", "design_background_mean", "d", "V")}
                      for r in chosen]
            assert [r["candidate"] for r in parsed] == [c["id"] for c in cohort]
            visits[visit] = plain(run._score_pair(parsed, design["bootstrap"]["seed"],
                                                    design["bootstrap"]["spatial_tile_resamples"]))
        entry["visits"] = visits
        entry["count_weighted_pooled_uncentered_mean_z2"] = float(np.average(
            [visits[v]["uncentered_mean_z2"] for v in ("search", "template")],
            weights=[visits[v]["n"] for v in ("search", "template")]))
        tr = [r for r in template if r["branch"] == branch]
        assert len(tr) == len(cohort)
        entry["common_template"] = {
            "n": len(tr),
            "mean_crossproduct_difference": float(np.mean([float(r["crossproduct_difference"]) for r in tr])),
            "mean_predicted_offdiag": float(np.mean([float(r["predicted_common_cov_offdiag"]) for r in tr]))}
        output["branches"][branch] = entry
    path = HERE / "stage_b_recovered_primary.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps({"result": str(path), "result_sha256": sha(path),
                      "strict": output["branches"]["strict"]["primary_adequacy"],
                      "secondary": {v: {k: visits[k] for k in ("n", "signed_standardized_mean",
                           "uncentered_mean_z2", "centered_variance_ratio", "sum_d2_minus_sum_V",
                           "sum_d2_over_sum_V")} for v, visits in
                           output["branches"]["secondary"]["visits"].items()}}, indent=2))


if __name__ == "__main__":
    main()
