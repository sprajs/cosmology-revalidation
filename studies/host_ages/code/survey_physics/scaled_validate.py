#!/usr/bin/env python3
"""Integrity and scientific-gate audit for the separate enlarged native campaign."""

import argparse
import json
from pathlib import Path
import re
import numpy as np
from common import ROOT, HERE, WORK, RESULTS, DATA, sha
from scaled_bbc import SCALE, BBC_EXE, read_attempts, fitres
from scaled_response import analyze


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-historical-replay", action="store_true",
                        help="Validate newly regenerated final science without requiring the ignored historical k012 pilot; explicitly record that regression audit as not rerun.")
    parser.add_argument("--restored-fixtures", action="store_true",
                        help="Verify the independently reconstructed capacity-fixtures record instead of requiring original ignored pilot files.")
    args = parser.parse_args()
    repair_path = RESULTS / ("scaled-native-capacity-restoration.json" if args.restored_fixtures
                             else "scaled-native-capacity-replay.json")
    campaign = json.loads((RESULTS / "scaled-campaign.json").read_text())
    shards = campaign["paired_shards"]
    execution = json.loads((RESULTS / "scaled-native-execution.json").read_text())
    assert execution["completed_shards"] == 2 * shards
    assert len({v["name"] for v in execution["records"]}) == 2 * shards
    checks = {}
    times = {
        "summed_generation_wall_seconds": 0.0,
        "summed_fit_wall_seconds": 0.0,
        "classifier_inference_seconds": 0.0,
    }
    populated = 0
    generated = {}
    for job in campaign["jobs"]:
        assert sha(job["input"]) == job["input_sha256"]
        run = json.loads((SCALE / "fits" / job["name"] / "run.json").read_text())
        assert (
            run["fit_graceful"]
            and run["generation_returncode"] == 0
            and run["fit_returncode"] == 0
        )
        for path, digest in run["hashes"].items():
            assert sha(path) == digest
        times["summed_generation_wall_seconds"] += run["generation_seconds"]
        times["summed_fit_wall_seconds"] += run["fit_seconds"]
        classifier = json.loads(
            (SCALE / "fits" / job["name"] / "classifier.json").read_text()
        )
        times["classifier_inference_seconds"] += classifier["inference_seconds"]
        assert classifier["truth_inputs_used"] is False
        d = read_attempts(job["name"])
        for suffix in [".DUMP", "_HEAD.FITS", "_PHOT.FITS"]:
            path = SCALE / "simulations" / job["name"] / (job["name"] + suffix)
            generated[str(path.relative_to(ROOT))] = sha(path)
        assert len(d) == job["attempts"]
        valid = d.RV.ne(-9)
        assert np.max(abs(d.loc[valid, "RV"] - 3.1)) < 1e-5
        assert d.loc[valid, "AV"].ge(0).all()
        populated += int(valid.sum())
    checks["all_attempts_and_populated_dust_support"] = True
    checks["populated_dust_occurrences"] = populated
    evaluation_seed = int(
        re.search(
            r"(?m)^RANSEED:\s*(\d+)",
            (WORK / "inputs/SPP_EVAL_NOMINAL.input").read_text(),
        ).group(1)
    )
    assert evaluation_seed not in [v["seed"] for v in campaign["jobs"]]
    checks["training_evaluation_seed_separation"] = True
    final_inputs = json.loads(
        (RESULTS / f"scaled-inputs-k{shards:03d}.json").read_text()
    )
    for record in final_inputs["tables"].values():
        assert sha(record["path"]) == record["sha256"]
    literal = fitres(final_inputs["tables"]["age_literal"]["path"])
    retained = fitres(final_inputs["tables"]["age_retained_age"]["path"])
    assert literal.index.equals(retained.index)
    changed = [k for k in literal if not literal[k].equals(retained[k])]
    assert changed == ["SIM_gammaDM"] and retained.SIM_gammaDM.eq(0).all()
    assert np.array_equal(literal.SIM_gammaDM, retained.TRUE_AGE_MAGSHIFT)
    checks["training_target_changed_columns"] = changed
    for v in final_inputs["paired_attempt_validation"].values():
        assert all(v["latent_identity"].values())
    checks["all_paired_training_latents_identical"] = True
    variants = {}
    for source in sorted(RESULTS.glob(f"scaled-bbc-k{shards:03d}*.json")):
        bbc = json.loads(source.read_text())
        response = analyze(bbc["cases"])
        variants[source.name] = {
            "native_completed": {k: v["graceful"] for k, v in bbc["cases"].items()},
            "native_fit_identified": {
                k: v.get("native_fit_identified", False)
                for k, v in bbc["cases"].items()
            },
            "analyzable_slope_support": response["analyzable_slope_support"],
            "response_status": response["status"],
            "native_identity_audit": response["native_identity_audit"],
        }
        for table in bbc["tables"].values():
            assert sha(table["path"]) == table["sha256"]
        if "-common" in source.name:
            assert bbc.get("same_nominal_map_geometry_verified") is True
    boots = {
        p.name: json.loads(p.read_text()).get("status", "executed")
        for p in RESULTS.glob(f"scaled-bootstrap-*-k{shards:03d}*.json")
    }
    expected = {"scaled-bootstrap-" + mode + "-" + name.removeprefix("scaled-bbc-")
                for name in variants for mode in ["joint", "training_only"]}
    assert set(boots) == expected, "Missing final bootstrap execution/support-gate record."
    for name in expected:
        b = json.loads((RESULTS / name).read_text())
        planned = 200 if b["mode"] == "joint" else 100
        if b.get("status") == "not_run_support_gate":
            assert b["replicates"] == 0 and b["planned_replicates"] == planned
        else:
            assert b["replicates"] == planned
            assert sum(b["replicate_status_counts"].values()) == planned
            assert sha(b["executable"]["path"]) == b["executable"]["sha256"]
            for i, digest in b["replicate_result_hashes"].items():
                rp = SCALE / "bootstrap" / (f"k{shards:03d}" + b["variant"]) / b["mode"] / f"r{int(i):03d}" / "result.json"
                assert sha(rp) == digest
    if args.skip_historical_replay:
        checks["historical_capacity_replay"] = "not_rerun; archival regression evidence only"
    else:
        repair = json.loads(repair_path.read_text())
        failure_types = []
        for r in repair["records"]:
            for key in ["original_log", "original_input"]:
                assert sha(r[key]["path"]) == r[key]["sha256"]
            if r["original_success"]:
                assert r["science_rows_byte_identity"] and r["covariance_byte_identity"]
            else:
                failure_types.append(r["original_failure_type"])
                assert r["replay"]["graceful"] == (r["original_failure_type"] == "scan_capacity")
        assert failure_types.count("scan_capacity") == 1 and failure_types.count("zero_MAD_pull") == 5
        checks["historical_capacity_failure_and_five_native_zero_MAD_failures_preserved"] = True
    source_files = sorted(HERE.glob("scaled*.py")) + sorted(
        HERE.glob("scaled-*-design.json")
    ) + [HERE / "scaled_scatter_capacity.patch"]
    result_files = [
        p
        for p in sorted(RESULTS.glob("scaled-*.json"))
        if p.name != "scaled-validation.json"
    ]
    deps = [
        BBC_EXE,
        DATA / "sample_input_files/DES-SN5YR/base_files/bbc/BBC_des5yr.input",
        WORK / "SNANA/src/SALT2mu.c",
        HERE / "bbc_bounds.patch",
        WORK / "SNANA-legacy-eff/src/sntools.c",
        SCALE / "native-capacity/src/sntools.c",
        WORK / "SNANA-legacy-eff/bin/SALT2mu.exe",
        RESULTS / "positive-dust-campaign.json",
    ]
    input_tables = {v["path"]: v["sha256"] for v in final_inputs["tables"].values()}
    for source in RESULTS.glob(f"scaled-bbc-k{shards:03d}*.json"):
        input_tables.update({v["path"]: v["sha256"] for v in json.loads(source.read_text())["tables"].values()})
    out = {
        "status": "passed",
        "meaning": "Executed integrity checks pass; this does not mean all scientific support or identification gates pass.",
        "final_paired_shards": shards,
        "attempts_total": campaign["attempts_total"],
        "checks": checks,
        "native_runtime": times,
        "final_variants": variants,
        "bootstrap_status": boots,
        "code_sha256": {str(p.relative_to(ROOT)): sha(p) for p in source_files},
        "result_sha256": {str(p.relative_to(ROOT)): sha(p) for p in result_files},
        "dependencies_sha256": {str(p.relative_to(ROOT)): sha(p) for p in deps},
        "generated_simulation_sha256": generated,
        "final_input_table_sha256": input_tables,
        "capacity_replay_record_sha256": sha(repair_path),
    }
    (RESULTS / "scaled-validation.json").write_text(json.dumps(out, indent=2) + "\n")
    print(
        json.dumps(
            {"status": "passed", "variants": variants, "bootstrap_status": boots},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
