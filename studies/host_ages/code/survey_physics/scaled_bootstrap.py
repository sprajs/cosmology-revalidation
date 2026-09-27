#!/usr/bin/env python3
"""Paired attempt bootstrap of native correction maps and nuisance fits."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
from pathlib import Path
import time
import numpy as np
import pandas as pd
from common import RESULTS, sha
from scaled_bbc import SCALE, read_attempts, fitres, write_table, execute_bbc
from scaled_response import analyze
from scaled_text_tables import cache_table, write_clones
import scaled_bbc


def lookup(frame, attempts):
    f = frame.copy()
    f["CID"] = f.index.astype(int)
    f["LIBID"] = f.SIM_LIBID.astype(int)
    f = f.reset_index(drop=True)
    keys = attempts[["CID", "LIBID", "SPH_SHARD"]].copy()
    keys["ATTEMPT"] = np.arange(len(keys))
    f = f.merge(keys, on=["CID", "LIBID", "SPH_SHARD"], validate="one_to_one")
    assert len(f) == len(frame)
    return f.drop(columns=["CID", "LIBID"]).set_index("ATTEMPT")


def draw_table(frame, draws, cid_base):
    mapping = pd.DataFrame({"ATTEMPT": draws, "CID": np.arange(len(draws)) + cid_base})
    table = mapping.merge(
        frame, left_on="ATTEMPT", right_index=True, how="inner", validate="many_to_one"
    )
    table = table.drop(columns=["ATTEMPT"]).set_index("CID")
    table.index = table.index.astype(str)
    assert table.index.is_unique
    return table


def prepare(shards, common=False, gauge=False):
    rec = json.loads((RESULTS / f"scaled-inputs-k{shards:03d}.json").read_text())
    pairs = {}
    attempts = {}
    for pool, prefix in [("train", "SPB_"), ("eval", "SPP_EVAL_")]:
        if pool == "train":
            a = pd.concat(
                [read_attempts(f"SPB_N{i:03d}") for i in range(shards)],
                ignore_index=True,
            )
        else:
            a = read_attempts("SPP_EVAL_NOMINAL")
        attempts[pool] = len(a)
        for arm in ["nominal", "age"]:
            key = arm + "_literal" if pool == "train" else "eval_" + arm
            path = rec["tables"][key]["path"]
            assert sha(path) == rec["tables"][key]["sha256"]
            pairs[pool + "_" + arm] = lookup(fitres(path), a)
    if common or gauge:
        ids = pairs["eval_nominal"].index.intersection(pairs["eval_age"].index)
        if gauge:
            ids = ids[(pairs["eval_nominal"].loc[ids].HOST_LOGMASS >= 10.04)
                      & (pairs["eval_age"].loc[ids].HOST_LOGMASS >= 10.04)]
        for arm in ["nominal", "age"]:
            pairs["eval_" + arm] = pairs["eval_" + arm].loc[ids]
        assert np.array_equal(pairs["eval_nominal"].zHD, pairs["eval_age"].zHD)
    pairs["_caches"] = {k: cache_table(v) for k,v in rec["tables"].items()}
    return pairs, attempts


def one(index, mode, shards, tables, attempts, variant="", fixed=False, gauge=False):
    try:
        return one_native(index, mode, shards, tables, attempts, variant, fixed, gauge)
    except Exception as error:
        result = {
            "index": index,
            "seed": 972000 + index,
            "mode": mode,
            "status": "numerical_validation_failure",
            "contrasts": {},
            "error": repr(error),
        }
        path = (
            SCALE / "bootstrap" / (f"k{shards:03d}" + variant) / mode / f"r{index:03d}"
        )
        path.mkdir(parents=True, exist_ok=True)
        (path / "result.json").write_text(json.dumps(result, indent=2) + "\n")
        return result


def one_native(index, mode, shards, tables, attempts, variant="", fixed=False, gauge=False):
    seed = 972000 + index
    rng = np.random.default_rng(seed)
    path = SCALE / "bootstrap" / (f"k{shards:03d}" + variant) / mode / f"r{index:03d}"
    path.mkdir(parents=True, exist_ok=True)
    drawn = {
        pool: (
            rng.integers(n, size=n)
            if mode == "joint" or pool == "train"
            else np.arange(n)
        )
        for pool, n in attempts.items()
    }
    inputs = {}
    for pool in ["train", "eval"]:
        for arm in ["nominal", "age"]:
            f = draw_table(
                tables[pool + "_" + arm],
                drawn[pool],
                10000000 if pool == "train" else 20000000,
            )
            for target in (
                ["literal", "retained_age"]
                if (pool == "train" and arm == "age")
                else ["literal"]
            ):
                key = f"{pool}_{arm}_{target}"
                source_key = arm + "_" + target if pool == "train" else "eval_" + arm
                inputs[key] = write_clones(f, path / (key + ".FITRES"), tables["_caches"][source_key], target)
    cases = {}
    for name, ev, tr in [
        ("nominal", "nominal", "nominal_literal"),
        ("age_nominal", "age", "nominal_literal"),
        ("age_literal", "age", "age_literal"),
        ("age_retained", "age", "age_retained_age"),
    ]:
        cases[name] = execute_bbc(
            inputs[f"eval_{ev}_literal"]["path"],
            inputs["train_" + tr]["path"],
            path / name,
            fixed,
            gauge,
        )
    if "-common" in variant and all(
        cases[k]["graceful"] for k in ["nominal", "age_nominal"]
    ):
        assert (
            cases["nominal"]["native_map_geometry"]
            == cases["age_nominal"]["native_map_geometry"]
        )
    result = analyze(cases)
    for val in result["contrasts"].values():
        val.pop("common_occurrence_ids", None)
    result.update(
        index=index,
        seed=seed,
        mode=mode,
        inputs=inputs,
        native_cases={
            k: {
                q: v[q]
                for q in [
                    "graceful",
                    "seconds",
                    "output_rows",
                    "parameters",
                    "log_sha256",
                    "native_support_lines",
                    "native_fit_identified",
                    "nuisance_design",
                    "executable_sha256",
                    "scatter_capacity_patch_sha256",
                    "scatter_capacity_scan_diagnostics",
                    "nuisance_native_bounds",
                    "gamma_gauge_max_host_factor_deviation",
                    "failure_kind",
                    "native_warning_lines",
                ]
                if q in v
            }
            for k, v in cases.items()
        },
    )
    (path / "result.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--shards", type=int, required=True)
    p.add_argument("--replicates", type=int, default=200)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--mode", choices=["joint", "training_only"], default="joint")
    p.add_argument("--common", action="store_true")
    p.add_argument("--fixed", action="store_true")
    p.add_argument("--highmass-gauge", action="store_true")
    a = p.parse_args()
    assert not (a.fixed and a.highmass_gauge)
    variant = ("-common-highmass-gauge" if a.highmass_gauge else
               ("-common" if a.common else "") + ("-fixed" if a.fixed else ""))
    tables, attempts = prepare(a.shards, a.common, a.highmass_gauge)
    start = time.monotonic()
    results = []
    baseline = json.loads(
        (RESULTS / (f"scaled-response-k{a.shards:03d}" + variant + ".json")).read_text()
    )
    destination = RESULTS / (
        f"scaled-bootstrap-{a.mode}-k{a.shards:03d}" + variant + ".json"
    )
    if not baseline["analyzable_slope_support"]:
        destination.write_text(
            json.dumps(
                {
                    "status": "not_run_support_gate",
                    "planned_replicates": a.replicates,
                    "replicates": 0,
                    "mode": a.mode,
                    "shards": a.shards,
                    "variant": variant,
                    "reason": "No baseline paired slope passes sample, nuisance-identifiability and native output gates.",
                    "code_sha256": sha(__file__),
                },
                indent=2,
            )
            + "\n"
        )
        return
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        jobs = [
            pool.submit(one, i, a.mode, a.shards, tables, attempts, variant, a.fixed, a.highmass_gauge)
            for i in range(a.replicates)
        ]
        for job in as_completed(jobs):
            r = job.result()
            results.append(r)
            print(r["index"], r["status"], len(results), flush=True)
    summaries = {}
    for label in sorted(set(k for r in results for k in r["contrasts"])):
        valid = [
            r["contrasts"][label]
            for r in results
            if label in r["contrasts"]
            and r["contrasts"][label]["slope_mag_per_Gyr"] is not None
        ]
        slope = np.array([r["slope_mag_per_Gyr"] for r in valid])
        bins = []
        for i in range(5):
            means = [
                r["centered_redshift_means"][i]["mean_mag"]
                for r in valid
                if r["centered_redshift_means"][i]["mean_mag"] is not None
            ]
            bins.append(
                {
                    "bin": i,
                    "supported_replicates": len(means),
                    "sd_mag": float(np.std(means, ddof=1)) if len(means) > 1 else None,
                    "percentile_95_mag": (
                        np.quantile(means, [0.025, 0.975]).tolist()
                        if len(means) > 1
                        else None
                    ),
                }
            )
        summaries[label] = {
            "supported_replicates": len(valid),
            "slope_sd_mag_per_Gyr": (
                float(np.std(slope, ddof=1)) if len(slope) > 1 else None
            ),
            "percentile_95_mag_per_Gyr": (
                np.quantile(slope, [0.025, 0.975]).tolist() if len(slope) > 1 else None
            ),
            "median_common_n": (
                float(np.median([r["n"] for r in valid])) if valid else None
            ),
            "redshift_means": bins,
        }
        active = [
            i
            for i, v in enumerate(
                baseline["contrasts"].get(label, {}).get("centered_redshift_means", [])
            )
            if v["mean_mag"] is not None
        ]
        complete = np.array(
            [
                [r["centered_redshift_means"][i]["mean_mag"] for i in active]
                for r in valid
                if active
                and all(
                    r["centered_redshift_means"][i]["mean_mag"] is not None
                    for i in active
                )
            ]
        )
        summaries[label]["centered_bin_covariance"] = {
            "baseline_supported_bin_indices": active,
            "complete_replicates": len(complete),
            "covariance_mag2": (
                np.atleast_2d(np.cov(complete, rowvar=False, ddof=1)).tolist()
                if len(complete) > 1
                else None
            ),
        }
        summaries[label]["fraction_of_injected_slope_95"] = (
            np.quantile(slope / -0.03, [0.025, 0.975]).tolist()
            if len(slope) > 1
            else None
        )
    differences = {}
    for before, after in [
        ("frozen_nominal", "age_nominal"),
        ("age_nominal", "age_literal"),
        ("age_literal", "age_retained"),
        ("frozen_nominal", "age_retained"),
    ]:
        vals = []
        for r in results:
            b = r["contrasts"].get(before + "_all_cases", {}).get("slope_mag_per_Gyr")
            c = r["contrasts"].get(after + "_all_cases", {}).get("slope_mag_per_Gyr")
            if b is not None and c is not None:
                vals.append(c - b)
        differences[after + "_minus_" + before] = {
            "complete_paired_replicates": len(vals),
            "sd_mag_per_Gyr": float(np.std(vals, ddof=1)) if len(vals) > 1 else None,
            "percentile_95_mag_per_Gyr": (
                np.quantile(vals, [0.025, 0.975]).tolist() if len(vals) > 1 else None
            ),
        }
    out = {
        "mode": a.mode,
        "shards": a.shards,
        "variant": variant,
        "replicates": a.replicates,
        "workers": a.workers,
        "seconds": time.monotonic() - start,
        "code_sha256": sha(__file__),
        "attempt_pool_sizes": attempts,
        "executable": {"path": str(scaled_bbc.BBC_EXE), "sha256": sha(scaled_bbc.BBC_EXE)},
        "row_cloning_code_sha256": sha(Path(__file__).with_name("scaled_text_tables.py")),
        "bbc_driver_code_sha256": sha(Path(__file__).with_name("scaled_bbc.py")),
        "seed_base": 972000,
        "summary": summaries,
        "paired_estimand_differences": differences,
        "uncertainty_scope": "Attempt resampling conditional on released cadence, W22 mock host population, frozen SALT/classifier and input physics. Not clustered survey-systematic uncertainty.",
        "interval_scope": "Percentile intervals summarize supported replicates only; failures and unsupported replicates are reported, not imputed. They are not established frequentist coverage guarantees, especially near native map-support boundaries.",
        "native_case_failure_counts": pd.Series([v.get("failure_kind", "none")
                 for r in results for v in r.get("native_cases", {}).values()
                 if not v.get("graceful", False)]).value_counts().to_dict(),
        "expanded_scan_case_count": sum(bool(v.get("scatter_capacity_scan_diagnostics", []))
                 for r in results for v in r.get("native_cases", {}).values()),
        "native_boundary_case_count": sum(any(b.get("at_bound",False) for b in v.get("nuisance_native_bounds",{}).values())
                 for r in results for v in r.get("native_cases", {}).values()),
        "replicate_status_counts": pd.Series([r["status"] for r in results])
        .value_counts()
        .to_dict(),
        "replicate_result_hashes": {
            str(r["index"]): sha(
                SCALE
                / "bootstrap"
                / (f"k{a.shards:03d}" + variant)
                / a.mode
                / f"r{r['index']:03d}"
                / "result.json"
            )
            for r in results
        },
    }
    destination.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
