#!/usr/bin/env python3
"""Released-dimensional BBC support on enlarged, independently generated mocks."""

import argparse
import datetime
import io
import json
import re
import sys
from pathlib import Path
import subprocess
import time
import numpy as np
import pandas as pd
from scipy.special import expit
from common import ROOT, HERE, WORK, RESULTS, DATA, sha, native_env

sys.path.insert(0, str(ROOT))
from lib.records import fitres

SCALE = WORK / "scaled-bbc"
BBC_EXE = SCALE / "native-capacity/bin/SALT2mu.exe"
EDGES = np.array([0.05, 0.3, 0.5, 0.7, 0.9, 1.2])


def quality(f):
    errors = f[["x0ERR", "x1ERR", "cERR"]].to_numpy()
    covariance = np.zeros((len(f), 3, 3))
    covariance[:, range(3), range(3)] = errors**2
    for i, j, k in [(0, 1, "COV_x1_x0"), (0, 2, "COV_c_x0"), (1, 2, "COV_x1_c")]:
        covariance[:, i, j] = covariance[:, j, i] = f[k]
    corr = covariance / errors[:, :, None] / errors[:, None, :]
    finite = np.isfinite(corr).all(axis=(1, 2))
    mineig = np.full(len(f), -np.inf)
    mineig[finite] = np.linalg.eigvalsh(corr[finite])[:, 0]
    keep = (
        f.ERRFLAG_FIT.eq(0)
        & f.x1.between(-3, 3)
        & f.c.between(-0.3, 0.3)
        & f.x1ERR.between(0, 1, inclusive="neither")
        & f.cERR.between(0, 1.5, inclusive="neither")
        & f.PKMJDERR.between(0, 2, inclusive="neither")
        & f.FITPROB.gt(0.001)
        & f.zHD.between(0.05, 1.2)
        & (mineig > 0)
        & f.x0.gt(0)
    )
    return keep


def read_selected(name):
    base = SCALE if name.startswith("SPB_") else WORK
    folder = base / "fits" / name
    source = folder / "fit.FITRES.TEXT"
    run = json.loads((folder / "run.json").read_text())
    assert run["fit_graceful"]
    assert sha(source) == run["fit.FITRES.TEXT"]["sha256"]
    f = fitres(source)
    assert f.index.is_unique
    cl = json.loads((folder / "classifier.json").read_text())
    assert cl["native_fit_sha256"] == sha(source)
    probability = (
        pd.read_csv(folder / "classifier.csv", dtype={"CID": str}).set_index("CID").pIa
    )
    probability = probability.reindex(f.index)
    assert probability.notna().all()
    q = quality(f)
    selected = f[q & probability.gt(0.5)].copy()
    assert np.max(abs(selected.SIM_RV - 3.1)) < 1e-6 and selected.SIM_AV.ge(0).all()
    # Preserve original observed/truth columns; extra labels never used as predictors.
    selected["TRUE_AGE_MAGSHIFT"] = selected.SIM_gammaDM
    selected["ORIGINAL_CID"] = selected.index.astype(int)
    if name.startswith("SPB_"):
        selected["SPH_SHARD"] = int(name[-3:])
    else:
        selected["SPH_SHARD"] = -1
    counts = {
        "fitres": len(f),
        "quality": int(q.sum()),
        "pIa_quality": len(selected),
        "input_sha256": sha(source),
        "classifier_sha256": sha(folder / "classifier.csv"),
    }
    return selected, counts


def read_attempts(name):
    base = SCALE if name.startswith("SPB_") else WORK
    path = base / "simulations" / name / (name + ".DUMP")
    lines = path.read_text().splitlines()
    header = next(line.split()[1:] for line in lines if line.startswith("VARNAMES:"))
    d = pd.read_csv(
        io.StringIO("\n".join(x[3:] for x in lines if x.startswith("SN:"))),
        sep=r"\s+",
        names=header,
    )
    assert not d.duplicated(["CID", "LIBID"]).any()
    d["SPH_SHARD"] = int(name[-3:]) if name.startswith("SPB_") else -1
    return d


def write_table(f, path, target="literal"):
    output = f.copy()
    if target == "retained_age":
        output["SIM_gammaDM"] = 0.0
    assert np.array_equal(output.SIM_DLMAG, f.SIM_DLMAG)
    changed = (
        [
            c
            for c in f
            if not np.array_equal(output[c].to_numpy(), f[c].to_numpy(), equal_nan=True)
        ]
        if all(f[c].dtype.kind not in "OU" for c in f)
        else [c for c in f if not output[c].equals(f[c])]
    )
    assert set(changed) <= ({"SIM_gammaDM"} if target == "retained_age" else set())
    output = output.reset_index(names="CID")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as stream:
        stream.write(
            "# Explicit survey-physics derived table; native source and target mapping in scaled-campaign/scaled-analysis-design.\n"
        )
        stream.write("VARNAMES: " + " ".join(output.columns) + "\n")
        output.insert(0, "ROW", "SN:")
        output.to_csv(stream, sep=" ", index=False, header=False, float_format="%.17g")
    return {
        "path": str(path),
        "rows": len(f),
        "sha256": sha(path),
        "target": target,
        "changed_columns": changed,
    }


def merge(shards):
    manifest = json.loads((RESULTS / "scaled-campaign.json").read_text())
    assert shards <= manifest["paired_shards"]
    records = {}
    counts = {}
    pairs = {}
    for arm, label in [("nominal", "N"), ("age", "A")]:
        tables = []
        for i in range(shards):
            name = f"SPB_{label}{i:03d}"
            f, c = read_selected(name)
            tables.append(f)
            counts[name] = c
        joined = pd.concat(tables)
        assert joined.index.is_unique
        if arm == "age":
            assert (
                np.max(abs(joined.SIM_gammaDM + 0.03 * (joined.SIM_HOSTLIB_SN_age - 3)))
                < 2e-7
            )
        else:
            assert joined.SIM_gammaDM.eq(0).all()
        for target in ["literal"] if arm == "nominal" else ["literal", "retained_age"]:
            path = SCALE / "tables" / f"train-{arm}-{target}-k{shards:03d}.FITRES"
            records[arm + "_" + target] = write_table(joined, path, target)
    for label in ["NOMINAL", "AGE"]:
        f, c = read_selected("SPP_EVAL_" + label)
        counts["SPP_EVAL_" + label] = c
        records["eval_" + label.lower()] = write_table(
            f, SCALE / "tables" / ("eval-" + label.lower() + ".FITRES")
        )
    for i in range(shards):
        a = read_attempts(f"SPB_N{i:03d}")
        b = read_attempts(f"SPB_A{i:03d}")
        common = a.merge(
            b,
            on=["CID", "LIBID", "SPH_SHARD"],
            validate="one_to_one",
            suffixes=("_nom", "_age"),
        )
        assert len(a) == len(b) == len(common) == manifest["attempts_per_shard"]
        cols = [
            "GENZ",
            "GALID",
            "SN_age",
            "PEAKMJD",
            "SALT2x1",
            "SALT2c",
            "AV",
            "RV",
            "MU",
        ]
        identity = {
            x: bool((common[x + "_nom"] == common[x + "_age"]).all()) for x in cols
        }
        assert all(identity.values())
        pairs[str(i)] = {
            "attempts": len(a),
            "latent_identity": identity,
            "selected_nominal": int(a.FLAG_ACCEPT.sum()),
            "selected_age": int(b.FLAG_ACCEPT.sum()),
            "selection_disagreements": int(
                (common.FLAG_ACCEPT_nom != common.FLAG_ACCEPT_age).sum()
            ),
        }
    result = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "paired_shards": shards,
        "code_sha256": sha(__file__),
        "analysis_design_sha256": sha(HERE / "scaled-analysis-design.json"),
        "counts": counts,
        "tables": records,
        "paired_attempt_validation": pairs,
    }
    (RESULTS / f"scaled-inputs-k{shards:03d}.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(
        json.dumps(
            {
                "shards": shards,
                "table_rows": {k: v["rows"] for k, v in records.items()},
            },
            indent=2,
        )
    )
    return result


def source_config(data, train, fixed=False, gauge=False):
    source = DATA / "sample_input_files/DES-SN5YR/base_files/bbc/BBC_des5yr.input"
    text = source.read_text().split("#END_YAML", 1)[1]
    for key in [
        "cid_reject_file",
        "varname_pIa",
        "simfile_ccprior",
        "idsurvey_list_probcc0",
    ]:
        text = re.sub(r"(?m)^" + key + r"=.*\n", "", text)
    changes = {
        "datafile": data,
        "simfile_biascor": train,
        "prefix": "bbc",
        "surveygroup_biascor": "'DES(zbin=0.075)'",
        "u13": "0",
        "p13": "0",
        "zmin": ".05",
        "ndump_nobiascor": "0",
    }
    for key, value in changes.items():
        text = re.sub(r"(?m)^" + key + r"=.*$", key + "=" + str(value), text)
    assert "opt_biascor=4336" in text and "u2=3" in text
    assert "snrmin_sigint_biascor" not in text
    if fixed:
        for key, value in {"u1": "0", "u2": "0", "u5": "0", "p5": "0"}.items():
            text = re.sub(r"(?m)^" + key + r"=.*$", key + "=" + value, text)
            changes[key] = value
    elif gauge:
        for key, value in {"u5": "0", "p5": "0"}.items():
            text = re.sub(r"(?m)^" + key + r"=.*$", key + "=" + value, text)
            changes[key] = value
    return text, source, changes


def execute_bbc(data, train, folder, fixed=False, gauge=False):
    folder.mkdir(parents=True, exist_ok=True)
    text, source, changes = source_config(data, train, fixed, gauge)
    inp = folder / "bbc.input"
    inp.write_text(text)
    start = time.monotonic()
    with (folder / "bbc.log").open("w") as out:
        p = subprocess.run(
            [str(BBC_EXE), str(inp)],
            cwd=folder,
            env=native_env(),
            stdout=out,
            stderr=subprocess.STDOUT,
        )
    log = (folder / "bbc.log").read_text()
    graceful = (
        p.returncode == 0 and "Done." in log[-1000:] and "FATAL ERROR ABORT" not in log
    )
    r = {
        "returncode": p.returncode,
        "seconds": time.monotonic() - start,
        "graceful": graceful,
        "input_sha256": sha(inp),
        "data_table_sha256": sha(data),
        "training_table_sha256": sha(train),
        "reference_input_sha256": sha(source),
        "reference_config_overrides": {k: str(v) for k, v in changes.items()},
        "executable_sha256": sha(BBC_EXE),
        "capacity_patch_sha256": sha(HERE / "bbc_bounds.patch"),
        "scatter_capacity_patch_sha256": sha(HERE / "scaled_scatter_capacity.patch"),
        "scatter_capacity_scan_diagnostics": [x.strip() for x in log.splitlines() if "CAPACITY_SCAN:" in x],
        "log_sha256": sha(folder / "bbc.log"),
        "folder": str(folder),
        "nuisance_mode": (
            "fixed_gamma_highmass_width_colour_refit"
            if gauge
            else ("fixed_released_construction_values" if fixed else "native_refit")
        ),
        "native_support_lines": [
            x.strip()
            for x in log.splitlines()
            if any(
                s in x
                for s in [
                    "Stored list of",
                    "NBIASCOR_CUTS",
                    "Rejected ",
                    "No BIASCOR",
                    "NSNFIT",
                    "MINPERCELL",
                    "min_per_cell",
                    "BCOR_LOSS",
                    "sigInt[",
                    "NCELL_INTERP_USE",
                ]
            )
        ],
        "native_warning_lines": [
            x.strip()
            for x in log.splitlines()
            if any(
                v in x
                for v in [
                    "WARNING(exec_mnpout",
                    "CrazyERR",
                    "BAD_OUTPUT",
                    "HAS_COV_PROBLEM",
                    "HAS_CRAZY_ERRORS",
                ]
            )
        ],
    }
    if not graceful:
        r["failure_tail"] = log[-2500:]
        r["failure_kind"] = ("scan_storage_capacity" if "exceeds bound MXSTORE_PULL" in log
                             else "zero_MAD_pull" if "Invalid stdPull = 0.000000" in log
                             else "no_native_bias_support" if "No BIASCOR" in log
                             else "other_native_failure")
    else:
        f = fitres(folder / "bbc.FITRES")
        r["output_rows"] = len(f)
        text = (folder / "bbc.FITRES").read_text()
        r["parameters"] = {
            m.group(1): {
                "value": float(m.group(2)),
                "conditional_error": float(m.group(3)),
            }
            for m in re.finditer(
                r"#\s+(alpha0|beta0|gamma0)\s+=\s+([\deE+.-]+)\s+\+-\s+([\deE+.-]+)",
                text,
            )
        }
        if fixed:
            r["parameters"] = {
                k: {"value": v, "conditional_error": 0.0, "fixed": True}
                for k, v in [("alpha0", 0.15), ("beta0", 2.87), ("gamma0", 0.0)]
            }
        if gauge:
            r["parameters"]["gamma0"] = {
                "value": 0.0,
                "conditional_error": 0.0,
                "fixed": True,
            }
            deviation = float(
                np.max(abs(0.5 - expit((f.HOST_LOGMASS - 10) / 0.001) + 0.5))
            )
            assert f.HOST_LOGMASS.ge(10.04).all() and deviation < 1e-15
            r["gamma_gauge_max_host_factor_deviation"] = deviation
        for key in ["M0avg", "sigint"]:
            hit = re.search(r"#\s+" + key + r"\s+=\s+([\deE+.-]+)", text)
            if hit:
                r[key] = float(hit.group(1))
        r["global_gamma_offset_lines"] = [
            x.strip()
            for x in log.splitlines()
            if "gammadm" in x.lower() and "OFFSET" in x
        ]
        r["outputs"] = {
            name: sha(folder / name) for name in ["bbc.FITRES", "bbc.M0DIF", "bbc.COV"]
        }
        r["MUCOV"] = {
            key: {
                "min": float(f[key].min()),
                "median": float(f[key].median()),
                "max": float(f[key].max()),
            }
            for key in ["biasCor_muCOVSCALE", "biasCor_muCOVADD"]
            if key in f
        }
        values = [
            line
            for line in (folder / "bbc.COV").read_text().splitlines()
            if line and not line.startswith("#")
        ]
        n = int(values[0])
        c = np.array([float(v) for v in values[1:]]).reshape(n, n)
        assert np.isfinite(c).all() and np.max(abs(c - c.T)) < 1e-8
        r["covariance_dimension"] = n
        r["covariance_min_eigenvalue"] = float(np.linalg.eigvalsh(c).min())
        r["populated_IZBIN_counts"] = {
            str(k): int(v) for k, v in f.IZBIN.value_counts().sort_index().items()
        }
        informative = sorted(
            int(v) for v in f.loc[f.M0DIFERR.lt(998), "IZBIN"].unique()
        )
        r["informative_covariance_bins"] = informative
        r["informative_covariance_min_eigenvalue"] = (
            float(np.linalg.eigvalsh(c[np.ix_(informative, informative)]).min())
            if informative
            else None
        )
        r["covariance_scope"] = (
            "Native conditional raw fitted-bin Hesse covariance, with free global normalization; "
            "not the exactly centered M0DIF covariance, correction-training Monte Carlo or survey systematics. Empty bins "
            "with MUDIFERR=999 are not informative, regardless of matrix positivity."
        )
        bins = sorted(f.IZBIN.unique())
        design = np.column_stack(
            [f.x1, f.c, 0.5 - expit((f.HOST_LOGMASS - 10) / 0.001)]
            + [(f.IZBIN == b).astype(float) for b in bins]
        )
        if fixed:
            design = design[:, 3:]
        elif gauge:
            design = np.column_stack([design[:, :2], design[:, 3:]])
        weighted = design / f.MUERR.to_numpy()[:, None]
        norms = np.linalg.norm(weighted, axis=0)
        scaled = weighted / np.where(norms > 0, norms, 1)
        _, singular, right = np.linalg.svd(scaled, full_matrices=True)
        rank = int(np.sum(singular > singular[0] * 1e-10))
        null = right[rank:] / np.where(norms > 0, norms, 1)[None, :]
        nfree = 2 if gauge else 3
        centered_shape = (
            design[:, :nfree] @ null[:, :nfree].T
            if not fixed
            else np.zeros((len(f), len(null)))
        )
        centered_shape -= centered_shape.mean(axis=0, keepdims=True)
        null_shape_max = (
            float(np.max(abs(centered_shape))) if centered_shape.size else 0.0
        )
        r["nuisance_design"] = {
            "rank": rank,
            "columns": design.shape[1],
            "singular_values": singular.tolist(),
            "rank_relative_threshold": 1e-10,
            "null_direction_centered_mu_shape_max": null_shape_max,
            "centered_mu_shape_identified_by_design": bool(null_shape_max < 1e-8),
            "low_mass_count": int(f.HOST_LOGMASS.lt(10).sum()),
            "high_mass_count": int(f.HOST_LOGMASS.ge(10).sum()),
            "mixed_mass_native_bins": [
                int(b)
                for b in bins
                if f.loc[f.IZBIN.eq(b), "HOST_LOGMASS"].lt(10).any()
                and f.loc[f.IZBIN.eq(b), "HOST_LOGMASS"].ge(10).any()
            ],
        }
        r["native_fit_identified"] = bool(
            rank == design.shape[1]
            and "BAD_OUTPUT ERROR detected" not in log
            and r["informative_covariance_min_eigenvalue"] is not None
            and r["informative_covariance_min_eigenvalue"] > 0
        )
        r["host_mass_overlap_supported"] = bool(
            fixed
            or gauge
            or (
                r["nuisance_design"]["low_mass_count"] > 0
                and r["nuisance_design"]["high_mass_count"] > 0
                and len(r["nuisance_design"]["mixed_mass_native_bins"]) > 0
            )
        )
        r["gamma_at_native_bound"] = bool(
            not fixed
            and not gauge
            and abs(r["parameters"]["gamma0"]["value"]) >= 0.5 - 1e-5
        )
        r["nuisance_native_bounds"] = {
            k: {
                "min": lo,
                "max": hi,
                "at_bound": bool(
                    not fixed
                    and not (gauge and k == "gamma0")
                    and (
                        r["parameters"][k]["value"] <= lo + 1e-5
                        or r["parameters"][k]["value"] >= hi - 1e-5
                    )
                ),
            }
            for k, lo, hi in [
                ("alpha0", 0.02, 0.30),
                ("beta0", 1.0, 6.0),
                ("gamma0", -0.5, 0.5),
            ]
        }
        r["native_fit_identified"] = bool(
            r["native_fit_identified"]
            and r["host_mass_overlap_supported"]
            and not any(v["at_bound"] for v in r["nuisance_native_bounds"].values())
        )
        r["native_map_geometry"] = [
            x.strip()
            for x in log.splitlines()
            if any(
                v in x
                for v in ["NBINz =", "NBIASCOR_CUTS =", "set_MAPCELL_biasCor : malloc"]
            )
        ]
    (folder / "run.json").write_text(json.dumps(r, indent=2) + "\n")
    return r


def run_stage(shards, common=False, fixed=False, gauge=False):
    record = json.loads((RESULTS / f"scaled-inputs-k{shards:03d}.json").read_text())
    t = record["tables"]
    if gauge:
        assert not fixed
        common = True
    suffix = (
        ("-common" if common else "")
        + ("-fixed" if fixed else "")
        + ("-highmass-gauge" if gauge else "")
    )
    if common:
        n, a = [fitres(t["eval_" + arm]["path"]) for arm in ["nominal", "age"]]
        ids = n.index.intersection(a.index)
        original_common = len(ids)
        if gauge:
            ids = ids[
                n.loc[ids, "HOST_LOGMASS"].ge(10.04)
                & a.loc[ids, "HOST_LOGMASS"].ge(10.04)
            ]
        assert np.array_equal(n.loc[ids, "zHD"], a.loc[ids, "zHD"])
        assert np.array_equal(n.loc[ids, "FIELD"], a.loc[ids, "FIELD"])
        for arm, f in [("nominal", n), ("age", a)]:
            t["eval_" + arm] = write_table(
                f.loc[ids], SCALE / "tables" / f"eval-{arm}{suffix}.FITRES"
            )
    cases = [
        ("nominal", "eval_nominal", "nominal_literal"),
        ("age_nominal", "eval_age", "nominal_literal"),
        ("age_literal", "eval_age", "age_literal"),
        ("age_retained", "eval_age", "age_retained_age"),
    ]
    results = {}
    for name, e, tr in cases:
        for which in [e, tr]:
            assert sha(t[which]["path"]) == t[which]["sha256"]
        folder = SCALE / "bbc" / (f"k{shards:03d}" + suffix) / name
        results[name] = execute_bbc(t[e]["path"], t[tr]["path"], folder, fixed, gauge)
        print(
            name,
            results[name]["graceful"],
            results[name].get("output_rows"),
            round(results[name]["seconds"], 2),
            flush=True,
        )
    out = {
        "paired_shards": shards,
        "attempts_per_arm": json.loads((RESULTS / "scaled-campaign.json").read_text())[
            "attempts_per_shard"
        ]
        * shards,
        "code_sha256": sha(__file__),
        "analysis_design_sha256": sha(HERE / "scaled-analysis-design.json"),
        "cases": results,
        "tables": t,
        "evaluation_cohort": (
            "common_preBBC_quality_classifier_occurrences"
            if common
            else "arm_specific_selected_occurrences"
        ),
    }
    if common and all(results[k]["graceful"] for k in ["nominal", "age_nominal"]):
        assert (
            results["nominal"]["native_map_geometry"]
            == results["age_nominal"]["native_map_geometry"]
        )
        out["same_nominal_map_geometry_verified"] = True
    if gauge:
        out["highmass_stratum"] = {
            "min_logmass": 10.04,
            "pre_stratum_common_candidates": original_common,
            "retained_common_candidates": len(ids),
            "candidate_loss": original_common - len(ids),
        }
    (RESULTS / (f"scaled-bbc-k{shards:03d}" + suffix + ".json")).write_text(
        json.dumps(out, indent=2) + "\n"
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("stage", choices=["merge", "bbc"])
    p.add_argument("--shards", type=int, required=True)
    p.add_argument("--common", action="store_true")
    p.add_argument("--fixed", action="store_true")
    p.add_argument("--highmass-gauge", action="store_true")
    a = p.parse_args()
    if a.stage == "merge":
        merge(a.shards)
    else:
        run_stage(a.shards, a.common, a.fixed, a.highmass_gauge)


if __name__ == "__main__":
    main()
