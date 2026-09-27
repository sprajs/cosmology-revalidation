#!/usr/bin/env python3
"""Audit nebular dust, infrared information/catalogues and enlarged native BBC.

This preserves the first physical campaign as a fixed historical record. An
executed scientific support failure can be recorded honestly; missing or still
running calculations cannot pass this audit.
"""
from __future__ import annotations
import ast
import difflib
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "studies/host_ages/code"
RESULTS = ROOT / "studies/host_ages/results"
NOTES = ROOT / "studies/host_ages/notes"
BASE_COMMIT = "ca8cfaabba5db7d438aa8ced24b1e947cf7b9eec"
REPORT = ROOT / "validation/reports/physical-extension.json"
sys.path.insert(0, str(ROOT))
from lib.records import fitres as parse_fitres
CACHE = {}
TABLE_CACHE = {}
FAILURES = []
SECTIONS = {}


def rel(path):
    p = Path(path).resolve()
    return str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p)


def sha(path):
    p = Path(path).resolve()
    s = p.stat()
    if p not in CACHE:
        with p.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        CACHE[p] = (digest, s.st_size, s.st_mtime_ns)
    else:
        assert CACHE[p][1:] == (s.st_size, s.st_mtime_ns), "Changed during audit: " + str(p)
    return CACHE[p][0]


def read(path):
    sha(path)
    return json.loads(Path(path).read_text())


def fitres(path):
    path = Path(path).resolve()
    sha(path)
    if path not in TABLE_CACHE:
        TABLE_CACHE[path] = parse_fitres(path)
    return TABLE_CACHE[path]


def check(condition, label, detail=None):
    if not bool(condition):
        FAILURES.append({"check": label, "detail": detail})


def identity(path, digest):
    check(sha(path) == digest, "sha256", rel(path))


def identities(mapping, base=ROOT):
    for path, digest in mapping.items():
        identity(base / path, digest)


def close(actual, expected, label, atol=1e-9, rtol=1e-9):
    check(np.allclose(actual, expected, atol=atol, rtol=rtol, equal_nan=False), label)


def quantiles(series):
    return dict(zip(["p05", "median", "p95"], map(float, np.quantile(series, [.05, .5, .95]))))


def first_extension():
    path = ROOT / "provenance/physical-program-execution.json"
    previous = subprocess.check_output(
        ["git", "show", BASE_COMMIT + ":provenance/physical-program-execution.json"], cwd=ROOT)
    identity(path, hashlib.sha256(previous).hexdigest())
    record = json.loads(previous)
    counts, excluded = {}, []
    for name in ["code_and_design_sha256", "compact_result_sha256", "figure_sha256",
                 "environment_sha256", "validation_sha256"]:
        count = 0
        for filename, digest in record[name].items():
            if filename.endswith(".md"):
                excluded.append(filename)
                continue
            identity(ROOT / filename, digest)
            count += 1
        counts[name] = count
    SECTIONS["first_extension_preservation"] = {
        "base_commit": BASE_COMMIT, "execution_manifest_sha256": sha(path),
        "verified_historical_files": sum(counts.values()), "counts": counts,
        "excluded_branch_documentation": excluded,
        "scope": "All frozen scientific code, designs, results, figure, environment and validation identities preserved. Shared narrative documentation and edition metadata may advance."}


def legacy_catalogue_label():
    record = read(RESULTS/"host_transport/des-deep-interface.json")
    for source in record["source_catalogues"]:identity(ROOT/source["path"],source["sha256"])
    program = '''
import json,sys
from pathlib import Path
import pyarrow.parquet as pq
root=Path(sys.argv[1]);rows=[]
for field in ['C3','E2','X3']:
 p=root/f'.work/host-transport/des/Y3_DEEP_FIELDS_PHOTOM-SN-{field}-0000.parquet'
 q=pq.read_table(p,columns=['KNN_CLASS']).column('KNN_CLASS').to_numpy()
 rows.append(dict(field=field,all_catalogue_objects=len(q),galaxy_class1=int((q==1).sum())))
print(json.dumps(rows))
'''
    rows = json.loads(subprocess.check_output([str(ROOT/".work/host-transport/.venv/bin/python"),"-c",program,str(ROOT)]))
    for row in rows:
        source = next(item for item in record["source_catalogues"] if item["field"] == row["field"])
        check(row["all_catalogue_objects"] == source["galaxies"],"historical galaxies key actually counts all catalogue rows")
    cross = pd.read_csv(ROOT/".work/host-transport/des-deep-crosswalk.csv",dtype={"SNID":str})
    selected = cross[cross.primary_match & cross.selected_Dovekie & cross.deep_quality & cross.host_dlr_lt4]
    check(len(selected) == 331,"host-quality selected sample unaffected by label correction")
    SECTIONS["legacy_catalogue_label_clarification"] = {
        "fields":rows,"total_objects":sum(row["all_catalogue_objects"] for row in rows),
        "total_galaxy_class1":sum(row["galaxy_class1"] for row in rows),"selected_hosts":len(selected),
        "historical_field":"des-deep-interface.json source_catalogues[].galaxies",
        "meaning":"The preserved historical key contains total raw catalogue rows, not galaxy-classified counts. Original science code/results remain unchanged; current narrative clarifies the denominator. KNN_CLASS=1 is a catalogue classification, not ground-truth galaxy identity."}


def nebular():
    folder, out = CODE / "nebular_dust", RESULTS / "nebular_dust"
    result, validation = read(out / "summary.json"), read(out / "validation.json")
    check(validation["passed"], "nebular component validation passed")
    identity(out / "summary.json", validation["summary_sha256"])
    identities(validation["code_sha256"])
    for key in ["input_sha256", "dependency_sha256", "output_sha256"]:
        identities(result[key])
    identity(folder / "analyze.py", result["code_sha256"])
    identity(folder / "design.json", result["design_sha256"])
    identity(out / "sources.json", result["source_context_sha256"])
    for item in read(out / "sources.json")["files"]:
        identity(ROOT / item["path"], item["sha256"])
    design = read(folder / "design.json")
    ledger = pd.read_csv(ROOT / ".work/nebular-dust/cohort-ledger.csv", dtype={"specObjID": str})
    parent = pd.read_csv(ROOT / ".work/galaxy-validation/brightness-cohort.csv", dtype={"specObjID": str})
    check(ledger.source_id.tolist() == parent.source_id.tolist(), "nebular parent sample identity")
    check(ledger.specObjID.tolist() == parent.specObjID.tolist(), "nebular exact spectrum identifiers")
    check(ledger.physical_host.tolist() == parent.physical_host.tolist(), "nebular exact physical hosts")
    check(ledger.groupby("physical_host").fold.nunique().max() == 1, "nebular no host leakage")
    a, b = ledger.h_alpha_flux.to_numpy(), ledger.h_beta_flux.to_numpy()
    sa, sb = ledger.h_alpha_flux_err.to_numpy(), ledger.h_beta_flux_err.to_numpy()
    valid = np.isfinite(np.c_[a,b,sa,sb]).all(axis=1) & (sa > 0) & (sb > 0)
    valid &= ~np.isin(np.c_[a,b], [-9999., 9999.]).any(axis=1)
    masks = {"all_valid_lines": valid, "balmer_snr3": valid & (a > 3*sa) & (b > 3*sb),
             "balmer_snr5": valid & (a > 5*sa) & (b > 5*sb)}
    for key, mask in masks.items():
        check(np.array_equal(mask, ledger[key]), "nebular independent selection " + key)
    for key, value in result["selection_counts"].items():
        check(int(ledger[key].sum()) == value, "nebular saved cohort size " + key)
    h = pd.read_csv(ROOT / ".work/nebular-dust/heldout.csv")
    recomputed = -.5*(np.log(2*np.pi*h.variance)+(h.y-h.prediction)**2/h.variance)
    score_error = float(abs(recomputed-h.score).max())
    check(score_error < 1e-12 and (h.variance > 0).all(), "nebular direct predictive likelihood")
    primary = []
    model_checks = 0
    for run in result["runs"]:
        check(not run["gates"], "nebular no unresolved fit/prediction gates", run["gates"])
        scores = {}
        for name, fit in run["fits"].items():
            if "host_age" in fit["parameters"]:
                value, se = fit["parameters"]["host_age"], fit["standard_errors"]["host_age"]
                close([value-1.96*se, value+1.96*se], fit["conditional_age_slope_95"], "nebular age interval arithmetic")
            if "heldout_log_score" not in fit:
                continue
            q = h[(h.cohort == run["cohort"]) & (h.SED_controls == run["SED_controls"])
                  & (h.setting == run["setting"]) & (h.model == name)].sort_values("source_id")
            check(len(q) == run["n"] and q.source_id.is_unique, "nebular heldout sample completeness")
            close(q.score.sum(), fit["heldout_log_score"], "nebular log-score total")
            close(np.sqrt(np.mean((q.y-q.prediction)**2)), fit["heldout_RMSE_mag"], "nebular heldout RMSE")
            scores[name] = q.set_index("source_id")
            model_checks += 1
        rng = np.random.default_rng(design["seed"])
        for added, base in [("age","base"), ("nebular","base"), ("nebular_age","nebular"), ("nebular_age","age")]:
            if added not in scores or base not in scores:
                continue
            q, p = scores[added], scores[base]
            check(q.index.equals(p.index), "nebular common-cohort prediction comparison")
            difference = q.score - p.score
            comparison = run["comparisons"][added+"_minus_"+base]
            close(difference.sum(), comparison["delta_heldout_log_score"], "nebular score difference")
            units = difference.groupby(q.physical_host).sum().sort_index().to_numpy()
            draws = rng.choice(units, (2000, len(units)), replace=True).sum(axis=1)
            close(np.quantile(draws,[.025,.975]), comparison["host_bootstrap95_conditional_on_folds"], "nebular conditional host bootstrap")
        if run["setting"] == "formal":
            f = run["fits"]["nebular_age"]
            primary.append({"cohort": run["cohort"], "SED_controls": run["SED_controls"],
                            "hosts": run["n"], "age_slope_mag_per_Gyr": f["parameters"]["host_age"],
                            "conditional_SE_mag_per_Gyr": f["standard_errors"]["host_age"],
                            "age_heldout_gain_after_nebular": run["comparisons"]["nebular_age_minus_nebular"]})
    independent = result["independent_fit_checks"]
    check(all(r["passed"] and r["independent_optimizer_success"] for r in independent), "nebular independent optimization")
    SECTIONS["nebular_dust"] = {"selected_counts": result["selection_counts"], "recomputed_heldout_models": model_checks,
        "predictive_score_identity_max_error": score_error, "independently_checked_full_fits": len(independent),
        "primary_results": primary,
        "scope": "Conditional, selected-host association with first-order line-error propagation. Weak-line nonlinearity, missing joint SED errors, aperture and unknown ZTF amplitude blinding remain explicit."}


def infrared():
    folder, out = CODE / "infrared_resolution", RESULTS / "infrared_resolution"
    acquisition = read(out / "acquisition.json")
    identity(folder / "acquire.py", acquisition["code_sha256"])
    for item in acquisition["files"]:
        identity(ROOT / item["path"], item["sha256"])
        check((ROOT / item["path"]).stat().st_size == item["bytes"], "infrared source byte count", item["path"])
    result = read(out / "information.json")
    validation = read(out / "validation.json")
    identity(folder / "analyze.py", result["code_sha256"])
    identities(result["dependencies_sha256"])
    identities(validation["dependencies_sha256"])
    identity(ROOT / result["output"], result["output_sha256"])
    check(validation["status"] == "passed", "infrared component validator passed")
    frame = pd.read_csv(ROOT / result["output"], dtype={"sn_id": str})
    check(len(frame) == result["rows"] and frame.sn_id.nunique() == result["objects"], "infrared result dimensions")
    check(not frame.duplicated(["sn_id","band_um"]).any(), "infrared object/band uniqueness")
    maps = pd.read_csv(ROOT / ".work/host-transport/des-spire-photometry.csv", dtype={"sn_id":str})
    maps = maps[maps.eightband_331 & maps.valid_map_measurement & (maps.common_sky_offsets >= 64)]
    counts = maps.groupby("sn_id").band_um.nunique()
    check(set(frame.sn_id) == set(counts[counts == 3].index), "infrared exact source-derived selection")
    check(set(frame.band_um) == {250,350,500}, "infrared native bands")
    check(frame.all_candidates_noise_inflation.notna().all(), "infrared unresolved templates remain explicit")
    close(frame.all_candidates_noise_inflation * frame.host_residual_template_norm_fraction, 1,
          "infrared residual-template information identity", rtol=1e-10)
    pair = 1/np.sqrt(1-frame.nearest_pair_template_correlation**2)
    close(pair, frame.nearest_pair_radial_noise_inflation, "infrared direct pair correlation identity", rtol=1e-8)
    resolution = frame.nearest_separation_arcsec * np.sqrt(4*np.log(2)/(-np.log(.75)))
    close(resolution, frame.Gaussian_FWHM_for_pair_inflation2_arcsec, "infrared independent Gaussian pair criterion")
    check((frame.all_candidates_noise_inflation >= frame.nearest_pair_radial_noise_inflation-1e-7).all(),
          "adding unconstrained neighbours cannot improve information")
    recomputed = []
    for row in result["summaries"]:
        group = frame[frame.band_um == row["band_um"]]
        check(len(group) == row["hosts"], "infrared per-band object count")
        check(int(group.all_candidates_noise_inflation.isna().sum()) == row["numerically_unresolved"], "infrared failed-template accounting")
        for key, value in row.items():
            if isinstance(value, dict) and set(value) == {"p05","median","p95"}:
                calculated = quantiles(group[key])
                close(list(calculated.values()), list(value.values()), "infrared summary " + key, rtol=1e-8)
        check(int((group.all_candidates_noise_inflation > 10).sum()) == row["all_candidates_inflation_gt10"], "infrared confusion count")
        check(int((group.nearest_pair_radial_noise_inflation > 2).sum()) == row["nearest_pair_inflation_gt2"], "infrared nearest-pair count")
        recomputed.append({"band_um": row["band_um"], "hosts": len(group),
            "median_pair_inflation": float(group.nearest_pair_radial_noise_inflation.median()),
            "median_all_candidate_inflation": float(group.all_candidates_noise_inflation.median())})
    checks = validation["checks"]
    check(len(checks) == validation["checked_templates"], "infrared independent template count")
    check(len(checks) == 3*len(validation["selected_hosts"]), "infrared all selected hosts/bands tested")
    for item in checks:
        check(item["primary_relative_error"] < 2e-7 and item["QR_relative_error"] < 2e-7, "infrared independent projection agreement")
        check(item["original_unit_orthogonality"] < 2e-6, "infrared original-unit orthogonality")
        check(abs(item["noise_mean_standard_errors"]) < 5 and .91 < item["noise_sd_to_expected"] < 1.09, "infrared white-noise estimator check")
    for band, scenarios in validation["sensitivity_ratios"].items():
        rows = [item for item in checks if item["band_um"] == int(band)]
        for scenario, summary in scenarios.items():
            values = [item["sensitivity_ratio"][scenario] for item in rows]
            finite = [value for value in values if value is not None]
            check(all(np.isfinite(value) and value > 0 for value in finite),
                  "infrared identified sensitivity is positive finite", [band, scenario])
            check(summary["identified"] == len(finite) and
                  summary["numerically_unresolved"] == len(values)-len(finite),
                  "infrared sensitivity unresolved counts", [band, scenario])
            if finite:
                close(list(quantiles(finite).values()),
                      [summary[key] for key in ["p05", "median", "p95"]],
                      "infrared sensitivity conditional quantiles " + band + ":" + scenario)
            else:
                check(all(summary[key] is None for key in ["p05", "median", "p95"]),
                      "unresolved sensitivity has no finite precision", [band, scenario])
    method = read(out / "method-review.json")
    identities(method["source_sha256"])
    check(not method["findings"]["material_arithmetic_or_identity_defect_found"], "infrared independent method review")
    figure = read(out/"figure.json")
    identities(figure["dependencies_sha256"])
    check((ROOT/figure["figure"]).exists(), "infrared validated figure exists")
    SECTIONS["infrared_resolution"] = {"objects": result["objects"], "rows": len(frame),
        "independent_templates": len(checks), "recomputed_summaries": recomputed,
        "sensitivity_ratios": validation["sensitivity_ratios"],
        "scope": "Conditional unregularized signed-amplitude information under assumed beam/noise. No deblended infrared flux, host luminosity or universal angular-resolution guarantee."}


def scaled_survey():
    folder, out = CODE/"survey_physics", RESULTS/"survey_physics"
    work = ROOT/".work/survey-physics/scaled-bbc"
    validation = read(out/"scaled-validation.json")
    check(validation["status"] == "passed", "expanded survey component validation passed")
    identities(validation["code_sha256"])
    identities(validation["result_sha256"])
    for key in ["input_sha256", "dependencies_sha256"]:
        if key in validation:
            identities(validation[key])
    identities(validation["generated_simulation_sha256"])
    identities(validation["final_input_table_sha256"])
    identity(out/"scaled-native-capacity-replay.json",validation["capacity_replay_record_sha256"])
    campaign = read(out/"scaled-campaign.json")
    k = validation["final_paired_shards"]
    check(campaign["paired_shards"] == k, "final survey stage equals completed campaign")
    jobs = campaign["jobs"]
    check(len(jobs) == 2*k and len({j["name"] for j in jobs}) == 2*k, "all paired survey jobs unique")
    check(len(validation["generated_simulation_sha256"]) == 3*len(jobs),"all generated truth/head/photon identities present")
    check(sum(j["attempts"] for j in jobs) == campaign["attempts_total"], "generated attempt total")
    for arm in ["nominal","age"]:
        selected = [j for j in jobs if j["arm"] == arm]
        check({j["shard"] for j in selected} == set(range(k)), "complete training shards", arm)
        check(sum(j["attempts"] for j in selected) == campaign["attempts_per_arm"], "per-arm attempt total", arm)
    native = read(out/"scaled-native-execution.json")
    # This legacy field counts arm-specific shard jobs, not paired shards.
    check(native["completed_shards"] == len(jobs) and len(native["records"]) == len(jobs), "native execution completeness")
    completed = {r["name"]:r for r in native["records"]}
    check(set(completed) == {j["name"] for j in jobs}, "all declared native runs recorded")
    for job in jobs:
        path = work/"fits"/job["name"]
        entry = completed[job["name"]]
        identity(path/"run.json", entry["native_run_sha256"])
        run = read(path/"run.json")
        check(run["fit_graceful"] and run["generation_returncode"] == run["fit_returncode"] == 0,
              "native generation and light-curve fit completed", job["name"])
        identities(run["hashes"])
        identity(job["input"], job["input_sha256"])
        check(run["generated_attempts_requested"] == job["attempts"], "native declared attempt count", job["name"])
        for filename in ["fit.FITRES.TEXT","fit.SNANA.TEXT"]:
            identity(path/filename, run[filename]["sha256"])
        classifier = entry["classifier"]
        identity(path/"fit.FITRES.TEXT", classifier["native_fit_sha256"])
        identity(path/"classifier.csv", classifier["output_sha256"])
        check(classifier["truth_inputs_used"] is False and classifier["single_batch_max_probability_difference"] < 1e-6,
              "scaled classifier batch identity and no truth predictors", job["name"])
    stages = {"selected_input": audit_survey_stage(folder,out,k,"")}
    check((out/f"scaled-bbc-k{k:03d}-common.json").exists(),
          "declared common-input map-geometry sensitivity executed")
    prefix = f"scaled-bbc-k{k:03d}"
    for path in sorted(out.glob(prefix+"-*.json")):
        suffix = path.stem.removeprefix(prefix)
        stages[suffix.removeprefix("-")] = audit_survey_stage(folder,out,k,suffix)
    expected_variants = {prefix+("" if name == "selected_input" else "-"+name)+".json" for name in stages}
    check(set(validation["final_variants"]) == expected_variants,"component validator covers every final variant")
    # Completed uncertainty or an explicitly demonstrated baseline-support gate
    # is required; an unfinished campaign cannot pass by merely recording status.
    check("bootstrap_status" in validation, "expanded native bootstrap execution or support gate recorded")
    bootstrap = audit_bootstrap(folder,out,work,k,stages)
    capacity = audit_capacity_repair(folder,out)
    restored_capacity = audit_capacity_repair(folder,out,"scaled-native-capacity-restoration.json")
    redshift = audit_redshift_summary(folder,out,work,k)
    review = audit_survey_review(bootstrap,redshift)
    SECTIONS["scaled_survey"] = {"paired_shards":k,"attempts":campaign["attempts_total"],
        "native_execution_jobs":len(jobs),"stages":stages,
        "bootstrap_status":validation.get("bootstrap_status"),
        "bootstrap_independent_audit":bootstrap,
        "capacity_repair_audit":capacity,
        "restored_capacity_fixture_audit":restored_capacity,
        "independent_method_review":review,
        "exploratory_redshift_summary":redshift,
        "scope":"Hypothetical paired mock intervention and conditional native maps. Unsupported slopes/bins remain null; native covariance is not training-uncertainty or full-survey coverage."}


def infrared_photometry():
    folder, out = CODE/"infrared_photometry", RESULTS/"infrared_photometry"
    work = ROOT/".work/infrared-photometry"
    acquisition = read(out/"acquisition.json")
    identity(folder/"acquire.py",acquisition["code_sha256"])
    identity(folder/"design.json",acquisition["design_sha256"])
    identities(acquisition["inputs"])
    identity(work/"hosts.csv",acquisition["hosts_sha256"])
    hosts = pd.read_csv(work/"hosts.csv",dtype={"SNID":str,"deep_ID":str})
    source = pd.read_csv(ROOT/".work/host-transport/des-deep-crosswalk.csv",dtype={"SNID":str,"deep_ID":"Int64"})
    source = source[source.primary_match & source.selected_Dovekie & source.deep_quality & source.host_dlr_lt4]
    maps = pd.read_csv(ROOT/".work/host-transport/des-spire-photometry.csv",dtype={"sn_id":str})
    maps = maps[maps.eightband_331 & maps.valid_map_measurement & (maps.common_sky_offsets>=64)]
    nband = maps.groupby("sn_id").band_um.nunique()
    source = source[source.SNID.isin(nband[nband == 3].index)].set_index("SNID")
    check(len(hosts) == acquisition["hosts"] == len(source) and hosts.SNID.is_unique and hosts.deep_ID.is_unique,
          "infrared catalogue host sample dimensions")
    check(set(hosts.SNID) == set(source.index),"infrared catalogue exact frozen source cohort")
    source = source.loc[hosts.SNID]
    check(np.array_equal(hosts.deep_ID,source.deep_ID.astype(str)),"infrared catalogue physical host identity")
    close(hosts[["deep_RA","deep_DEC","redshift"]],source[["deep_RA","deep_DEC","redshift"]],
          "infrared catalogue host coordinates and redshift")
    catalogues = {}
    for item in acquisition["catalogues"]:
        identity(ROOT/item["path"],item["sha256"])
        parts = []
        for request in item["requests"]:
            identity(ROOT/request["path"],request["sha256"])
            check((ROOT/request["path"]).stat().st_size == request["bytes"],"infrared catalogue download bytes")
            if "SELECT *" in request["query"]:
                parts.append(pd.read_csv(ROOT/request["path"],dtype=str))
        raw = pd.concat(parts,ignore_index=True)
        check(raw.groupby("cntr",dropna=False).nunique(dropna=False).max().max()<=1,
              "infrared duplicate query rows identical",item["catalogue"])
        unique = raw.drop_duplicates("cntr").sort_values("cntr").reset_index(drop=True)
        table = pd.read_csv(ROOT/item["path"],dtype=str)
        check(table.equals(unique) and len(table) == item["rows"],"infrared native-table union exact",item["catalogue"])
        catalogues[item["catalogue"]] = table
    for item in acquisition["documentation"]:
        identity(ROOT/item["path"],item["sha256"])
    summary = read(out/"summary.json")
    check(summary["status"] == "catalogue_associations_only","infrared catalogue scope remains association-only")
    identity(folder/"analyze.py",summary["code_sha256"])
    identity(folder/"design.json",summary["design_sha256"])
    identity(out/"acquisition.json",summary["acquisition_sha256"])
    for key in ["dependencies_sha256","optical_input_sha256","output_sha256"]:
        identities(summary[key])
    candidates = pd.read_csv(work/"candidates.csv",dtype={"source_key":str,"cntr":str})
    check(len(candidates) == summary["candidate_sources"] == sum(map(len,catalogues.values())),
          "infrared candidate union dimensions")
    check(candidates.source_key.is_unique,"infrared source identifiers globally unique")
    for name,table in catalogues.items():
        rows = candidates[candidates.catalogue == name].set_index("cntr")
        check(set(rows.index) == set(table.cntr),"infrared candidate/native row identity",name)
        table = table.set_index("cntr").loc[rows.index]
        check(np.array_equal(rows.source_key,name+":"+rows.index),"infrared native source-key construction")
        close(rows[["ra","dec"]],table[["ra","dec"]].astype(float),"infrared native candidate coordinates")
    associations = pd.read_csv(work/"associations.csv",dtype={"SNID":str,"deep_ID":str,"source_key":str}).fillna({"source_key":""})
    check(len(associations) == summary["associations"] == len(hosts)*3 and
          not associations.duplicated(["SNID","catalogue"]).any(),"infrared complete unique host-family accounting")
    byhost = hosts.set_index("SNID")
    nearby = pd.read_csv(work/"host-candidates.csv",dtype={"SNID":str,"deep_ID":str,"source_key":str})
    check(not nearby.duplicated(["SNID","source_key"]).any(),"infrared full candidate-pair ledger unique")
    candidate_groups = {name:rows for name,rows in candidates.groupby("catalogue")}
    shared = associations[associations.primary_unique_ir].groupby("source_key").SNID.nunique()
    optical = read_optical_neighbours()
    for row in associations.itertuples(index=False):
        host = byhost.loc[row.SNID]
        check(row.deep_ID == host.deep_ID and row.field == host.field,"infrared association host identity")
        data = candidate_groups[row.catalogue]
        distances = vector_angles(sky_vectors([host.deep_RA],[host.deep_DEC])[0],sky_vectors(data.ra,data.dec))
        local = distances<=75+1e-7
        ledger = nearby[(nearby.SNID == row.SNID)&(nearby.catalogue == row.catalogue)].set_index("source_key")
        check(set(ledger.index) == set(data.loc[local,"source_key"]),"infrared full nearby candidate ledger")
        if len(ledger):
            close(ledger.loc[data.loc[local,"source_key"],"separation_arcsec"],distances[local],
                  "infrared independent candidate separations",atol=1e-7,rtol=0)
        check(row.candidates_within75 == int(local.sum()),"infrared 75-arcsec count")
        for radius in [1,2,3,6]:
            check(getattr(row,f"candidates_within{radius}") == int((distances<=radius).sum()),"infrared radius sensitivity count")
        primary = 2. if "_24_" in row.catalogue else 1.
        unique = int((distances<=primary).sum()) == 1
        check(row.primary_radius_arcsec == primary and row.primary_unique_ir == unique,"infrared primary-radius gate")
        if local.any():
            best = int(np.argmin(distances))
            check(row.source_key == data.iloc[best].source_key,"infrared nearest candidate identifier")
            close(row.separation_arcsec,distances[best],"infrared primary separation",atol=1e-7,rtol=0)
            others = optical[row.source_key]
            point = sky_vectors([data.iloc[best].ra],[data.iloc[best].dec])[0]
            angles = vector_angles(point,sky_vectors([v[1] for v in others],[v[2] for v in others]))
            is_target = np.array([v[0] == host.deep_ID for v in others])
            target_angle = float(angles[is_target][0]) if is_target.any() else float("inf")
            rival = angles[~is_target]
            nearest = not len(rival) or target_angle<rival.min()
            check(row.target_nearest_optical == nearest,"infrared nearest optical host")
            for radius,key in [(primary,"other_optical_within_primary"),(3,"other_optical_within3"),(6,"other_optical_within6")]:
                check(getattr(row,key) == int((rival<=radius).sum()),"infrared optical-competitor count")
        else:
            check(row.source_key == "" and pd.isna(row.separation_arcsec),"infrared unmatched source remains missing")
        is_shared = unique and shared.get(row.source_key,0)>1
        geometry = unique and row.target_nearest_optical and row.other_optical_within_primary == 0 and not is_shared
        check(row.source_shared_by_selected_hosts == is_shared and row.geometry_pass == geometry,"infrared conservative geometry gate")
    photometry = audit_infrared_fluxes(work,catalogues,associations,summary)
    controls = audit_matching_controls(out,work,hosts,candidates)
    SECTIONS["infrared_photometry"] = {"hosts":len(hosts),"catalogue_rows":{name:len(rows) for name,rows in catalogues.items()},
        "photometry":photometry,"independent_positional_controls":controls,
        "scope":"Catalogue detections and geometric association candidates; missing records are not forced-photometry upper limits or empirical dust corrections."}


def read_optical_neighbours():
    # The optional parquet reader stays in its existing isolated environment.
    # Return raw qualified optical positions; all angles/counts are recomputed
    # independently in the main numpy-only audit above.
    program = '''
import json,sys
from pathlib import Path
import numpy as np,pandas as pd
import pyarrow.parquet as pq
from scipy.spatial import cKDTree
root=Path(sys.argv[1]);work=root/'.work/infrared-photometry'
a=pd.read_csv(work/'associations.csv',dtype={'source_key':str}).dropna(subset=['source_key'])
c=pd.read_csv(work/'candidates.csv',dtype={'source_key':str}).set_index('source_key')
def xyz(ra,dec):
 r,d=np.deg2rad(ra),np.deg2rad(dec)
 return np.column_stack([np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)])
out={}
for field in ['C3','X3']:
 p=root/f'.work/host-transport/des/Y3_DEEP_FIELDS_PHOTOM-SN-{field}-0000.parquet'
 o=pq.read_table(p,columns=['ID','RA','DEC','FLAGS','MASK_FLAGS','KNN_CLASS']).to_pandas()
 o=o[(o.FLAGS==0)&(o.MASK_FLAGS==0)&(o.KNN_CLASS==1)]
 tree=cKDTree(xyz(o.RA,o.DEC))
 for key in a.loc[a.field==field,'source_key'].unique():
  source=c.loc[key];q=xyz([source.ra],[source.dec])[0]
  ids=tree.query_ball_point(q,2*np.sin(np.deg2rad(75.000001/3600)/2))
  out[key]=[[str(int(row.ID)),float(row.RA),float(row.DEC)] for row in o.iloc[ids].itertuples()]
print(json.dumps(out,allow_nan=False))
'''
    return json.loads(subprocess.check_output([str(ROOT/".work/host-transport/.venv/bin/python"),"-c",program,str(ROOT)],cwd=ROOT))


def audit_infrared_fluxes(work,catalogues,associations,summary):
    phot = pd.read_csv(work/"photometry.csv",dtype={"SNID":str,"deep_ID":str,"source_key":str,"detid_24":str,"native_flag_text":str}).fillna({"source_key":"","native_flag_text":""})
    check(len(phot) == summary["photometry_rows"] == len(associations)//3*8 and
          not phot.duplicated(["SNID","catalogue","band_um"]).any(),"infrared complete host/band ledger")
    tables = {name:table.set_index("cntr") for name,table in catalogues.items()}
    metadata = {}
    for name in catalogues:
        lines = (work/f"{name}-columns.csv").read_text().splitlines()
        check(lines[0] == "column_name,datatype,unit,description","infrared native metadata header")
        parsed = [re.fullmatch(r"([^,]+),([^,]*),([^,]*),(.*)",line).groups() for line in lines[1:] if line]
        table = pd.DataFrame(parsed,columns=["column_name","datatype","unit","description"])
        check(table.column_name.is_unique,"infrared metadata column identity")
        metadata[name] = table.set_index("column_name")
    amap = associations.set_index(["SNID","catalogue"])
    def native(key):
        if not key:return None
        name,index = key.split(":",1)
        return tables[name].loc[index]
    def numeric(value):
        try:return float(value)
        except (ValueError,TypeError):return float("nan")
    def flag(value):
        match = re.search(r"(-?\d+)$",str(value))
        return int(match[1]) if match else -1
    def same(a,b,label):
        check(np.allclose(a,b,atol=1e-12,rtol=1e-12,equal_nan=True),label)
    for row in phot.itertuples(index=False):
        assoc = amap.loc[(row.SNID,row.catalogue)]
        check(row.source_key == assoc.source_key and row.deep_ID == assoc.deep_ID and row.field == assoc.field,
              "infrared photometry association identity")
        servs = row.catalogue.startswith("servs")
        key = ({3.6:"1",4.5:"2"} if servs else {3.6:"36",4.5:"45",5.8:"58",8.:"80",24.:"24"})[row.band_um]
        fc,ec = (f"flux_aper_2_{key}",f"fluxerr_aper_2_{key}") if servs else (f"flux_ap2_{key}",f"uncf_ap2_{key}")
        check(row.native_flux_column == fc and row.native_error_column == ec,"infrared native aperture/column identity")
        check(metadata[row.catalogue].loc[fc,"unit"] == metadata[row.catalogue].loc[ec,"unit"] == "uJy",
              "infrared flux and uncertainty native units")
        source = native(row.source_key)
        value = numeric(source[fc]) if source is not None else float("nan")
        error = numeric(source[ec]) if source is not None else float("nan")
        same([value,error],[row.native_flux_uJy,row.native_error_uJy],"infrared native signed flux/error identity")
        valid = np.isfinite(value) and np.isfinite(error) and error>0 and value != -99 and error != -99
        check(row.flux_valid == valid,"infrared sentinel and signed-flux gate")
        coverage_column = f"cov_avg_{key}" if "_24_" in row.catalogue else f"cov_{key}"
        coverage_kind = "hundred_second_frames" if servs else ("mean_number_frames" if "_24_" in row.catalogue else "coverage_flag")
        check(row.coverage_column == coverage_column and row.coverage_kind == coverage_kind,"infrared catalogue-specific coverage units")
        coverage = numeric(source.get(coverage_column,float("nan"))) if source is not None else float("nan")
        flagtext = source[("flags_" if servs else "flgs_")+key] if source is not None else ""
        native_flag = flag(flagtext)
        mask = (flag(source[f"mask_{key}"]) if source is not None else -1) if servs else float("nan")
        extension = float("nan") if servs or source is None else numeric(source.get(f"ext_fl_{key}",float("nan")))
        same([coverage,extension],[row.coverage_native,row.native_extended_flag],"infrared coverage and extension metadata")
        check(row.native_flag == native_flag and row.native_flag_text == ("" if pd.isna(flagtext) else flagtext),"infrared native quality flags")
        same(mask,row.native_mask,"infrared native masks retain missingness")
        conflict = row.catalogue == "chandra_cat_f05" and row.band_um != 24
        check(row.aperture_metadata_conflict == conflict,"infrared release metadata conflict explicit")
        extension_conflict = not servs and extension == -1
        check(row.extension_semantics_conflict == extension_conflict,"infrared extension-flag contradiction explicit")
        clean = valid and native_flag == 0 and (not servs or mask == 0) and not(extension>0) and not conflict and not extension_conflict
        check(row.clean_photometry == clean,"infrared optional quality subset")
        within = np.isfinite(assoc.separation_arcsec) and assoc.separation_arcsec<=assoc.primary_radius_arcsec
        geometry = bool(assoc.geometry_pass)
        link = "not_applicable"
        if row.band_um == 24 and "_24_" not in row.catalogue:
            other_name = next(name for name in tables if "_24_" in name and
                              (name.startswith("chandra") == row.catalogue.startswith("chandra")))
            other = amap.loc[(row.SNID,other_name)]
            other_source = native(other.source_key)
            detid = str(source.get("detid_24","")) if source is not None else ""
            matched_id = other_source is not None and detid == str(other_source.detid_24) and detid not in ["","-99","nan"]
            geometry = geometry and other.geometry_pass and matched_id
            link = "same_detection_verified" if geometry else "separate_centroid_or_identity_not_verified"
        check(row.anchor_geometry_pass == assoc.geometry_pass and row.geometry_pass == geometry and
              row.bandmerged24_centroid_link == link and row.within_primary == within,"infrared band-specific centroid gate")
        if not within:state="no_catalogue_association_coverage_unknown"
        elif not valid:state="covered_flux_missing" if coverage>0 else ("no_coverage" if coverage == 0 else "flux_missing_coverage_unknown")
        elif not geometry:state="ambiguous_association"
        else:state="catalogue_measurement_geometry_pass"
        check(row.status == state,"infrared missingness and association semantics")
    check(phot.status.value_counts().to_dict() == summary["photometry_statuses"],"infrared complete status counts")
    for row in summary["bands"]:
        group = phot[(phot.catalogue == row["catalogue"])&(phot.band_um == row["band_um"])]
        good = group.geometry_pass & group.flux_valid
        check(int(good.sum()) == row["geometry_valid_measurements"] and
              int((good&group.clean_photometry).sum()) == row["clean_geometry_measurements"] and
              int((good&group.extension_semantics_conflict).sum()) == row["extension_semantics_conflict_valid_associations"] and
              group.status.value_counts().to_dict() == row["statuses"],"infrared per-band counts")
    for row in summary["unique_host_union_by_band"]:
        group = phot[(phot.band_um == row["band_um"])&phot.geometry_pass&phot.flux_valid]
        check(group.SNID.nunique() == row["unique_hosts_with_valid_geometric_association"],"infrared unique-host union avoids duplicate releases")
    for row in summary["matching"]:
        group = associations[(associations.family == row["family"])&(associations.field == row["field"])]
        check(len(group) == row["hosts"],"infrared family/field denominator")
        for radius,stats in row["radii_arcsec"].items():
            values = group[f"candidates_within{radius}"]
            check(int((values>0).sum()) == stats["any_match"] and int((values == 1).sum()) == stats["unique_ir"],
                  "infrared family radius counts")
        for key in ["primary_unique_ir","geometry_pass","source_shared_by_selected_hosts"]:
            check(int(group[key].sum()) == row[key],"infrared family geometry counts")
        check(int((group.primary_unique_ir & ~group.target_nearest_optical).sum()) == row["target_not_nearest_within_primary"] and
              int((group.primary_unique_ir & group.other_optical_within_primary.gt(0)).sum()) == row["optical_competition"],
              "infrared ambiguous optical-host counts")
        for radius in [3,6]:
            check(int((group.geometry_pass & group[f"other_optical_within{radius}"].gt(0)).sum()) == row[f"passing_with_other_optical_within{radius}"],
                  "infrared crowding beyond primary radius remains explicit")
    consistency = pd.read_csv(work/"release-consistency.csv",dtype={"SNID":str,"swire_source_key":str,"servs_source_key":str})
    check(not consistency.duplicated(["SNID","band_um"]).any(),"infrared release comparison uniqueness")
    expected_pairs = set()
    for field in ["C3","X3"]:
        for band in [3.6,4.5]:
            selected = phot[(phot.field == field)&(phot.band_um == band)&phot.geometry_pass]
            left = selected[~selected.catalogue.str.startswith("servs")]
            right = selected[selected.catalogue.str.startswith("servs")]
            expected_pairs |= {(snid,band) for snid in set(left.SNID)&set(right.SNID)}
    check(set(zip(consistency.SNID,consistency.band_um)) == expected_pairs,"infrared release-comparison cohort complete")
    for row in consistency.itertuples(index=False):
        left,right = native(row.swire_source_key),native(row.servs_source_key)
        angle = vector_angles(sky_vectors([float(left.ra)],[float(left.dec)])[0],sky_vectors([float(right.ra)],[float(right.dec)]))[0]
        close(angle,row.infrared_centroid_separation_arcsec,"independent inter-release centroid",atol=1e-7,rtol=0)
        check(row.centroids_agree_within1 == (angle<=1),"infrared release centroid gate")
        pair = phot[(phot.SNID == row.SNID)&(phot.band_um == row.band_um)].set_index("source_key")
        l,r = pair.loc[row.swire_source_key],pair.loc[row.servs_source_key]
        check(l.geometry_pass and r.geometry_pass,"infrared release comparison source geometry")
        valid,positive = bool(l.flux_valid and r.flux_valid),bool(l.flux_valid and r.flux_valid and l.native_flux_uJy>0 and r.native_flux_uJy>0)
        check(row.valid_native_pair == valid and row.positive_flux_pair == positive,"infrared ratio eligibility")
        check(row.both_clean_photometry == bool(l.clean_photometry and r.clean_photometry),"infrared release flags-clear comparison")
        same([row.swire_flux_uJy,row.servs_flux_uJy],[l.native_flux_uJy,r.native_flux_uJy],"infrared cross-release signed fluxes")
        same(row.signed_servs_minus_swire_uJy,r.native_flux_uJy-l.native_flux_uJy if valid else float("nan"),"infrared signed release difference")
        same(row.servs_over_swire,r.native_flux_uJy/l.native_flux_uJy if positive else float("nan"),"infrared positive-pair release ratio")
    for row in summary["release_consistency"]:
        group = consistency[(consistency.field == row["field"])&(consistency.band_um == row["band_um"])]
        good = group.centroids_agree_within1 & group.positive_flux_pair
        values = group.loc[good,"servs_over_swire"]
        check(len(group) == row["geometry_pairs"] and len(values) == row["ratio_pairs"] and
              int((~group.centroids_agree_within1).sum()) == row["centroid_disagreements_gt1"] and
              int((~group.valid_native_pair).sum()) == row["missing_or_invalid_pairs"] and
              int((group.valid_native_pair & ~group.positive_flux_pair).sum()) == row["nonpositive_valid_pairs"],
              "infrared release consistency denominator")
        check(int((good&group.both_clean_photometry).sum()) == row["clean_pairs"],"infrared strict release-comparison subset")
        if len(values):close(np.quantile(values,[.05,.5,.95]),row["ratio_p05_median_p95"],"infrared conditional release-ratio quantiles")
        else:check(row["ratio_p05_median_p95"] is None,"no unsupported infrared ratio summary")
    return {"rows":len(phot),"statuses":summary["photometry_statuses"],"bands":summary["bands"],
            "unique_host_union_by_band":summary["unique_host_union_by_band"],
            "release_consistency":summary["release_consistency"],"negative_valid_native_fluxes":int((phot.flux_valid & phot.native_flux_uJy.lt(0)).sum())}


def audit_survey_stage(folder,out,k,suffix):
    inputs = read(out/f"scaled-inputs-k{k:03d}.json")
    bbc = read(out/f"scaled-bbc-k{k:03d}{suffix}.json")
    response = read(out/f"scaled-response-k{k:03d}{suffix}.json")
    for record, script in [(inputs,"scaled_bbc.py"),(bbc,"scaled_bbc.py"),(response,"scaled_response.py")]:
        identity(folder/script,record["code_sha256"])
        identity(folder/"scaled-analysis-design.json",record["analysis_design_sha256"])
    check(inputs["paired_shards"] == bbc["paired_shards"] == k, "final native BBC stage identity")
    for item in bbc["tables"].values():
        identity(item["path"],item["sha256"])
    tables = {}
    for name,item in inputs["tables"].items():
        identity(item["path"],item["sha256"])
        table = fitres(item["path"])
        check(len(table) == item["rows"] and table.index.is_unique, "scaled training/evaluation table completeness", name)
        tables[name] = table
    literal, retained = tables["age_literal"], tables["age_retained_age"]
    check(literal.index.equals(retained.index), "literal and retained-age target identical occurrences")
    check(literal.drop(columns="SIM_gammaDM").equals(retained.drop(columns="SIM_gammaDM")),
          "retained-age target only changes SIM_gammaDM")
    check(retained.SIM_gammaDM.eq(0).all(), "retained-age declared target transformation")
    check(np.max(abs(literal.SIM_gammaDM+.03*(literal.SIM_HOSTLIB_SN_age-3))) < 2e-7,
          "training intervention sign and magnitude")
    check(tables["nominal_literal"].SIM_gammaDM.eq(0).all(), "nominal training has no age injection")
    frames, native_checks = {}, {}
    for name,case in bbc["cases"].items():
        path = Path(case["folder"])
        identity(path/"bbc.input",case["input_sha256"])
        identity(path/"bbc.log",case["log_sha256"])
        audit_bbc_configuration(case,path/"bbc.input")
        identity(case["reference_config_overrides"]["datafile"],case["data_table_sha256"])
        identity(case["reference_config_overrides"]["simfile_biascor"],case["training_table_sha256"])
        log = (path/"bbc.log").read_text()
        graceful = case["returncode"] == 0 and "Done." in log[-1000:] and "FATAL ERROR ABORT" not in log
        check(case["graceful"] == graceful, "native BBC completion status", name)
        if not case["graceful"]:
            check(bool(case.get("failure_tail")), "native BBC failure preserved", name)
            continue
        for filename,digest in case["outputs"].items():
            identity(path/filename,digest)
        f = fitres(path/"bbc.FITRES")
        check(f.index.is_unique and len(f) == case["output_rows"] and (f.MUERR>0).all(),
              "native BBC accepted rows", name)
        raw = native_distance_shape(f,case["parameters"])
        constant = np.median(raw-f.MU)
        shape_error = float(abs(raw-constant-f.MU).max())
        offset_error = float(abs(f.MU-f.MUMODEL-f.MURES-f.M0DIF).max())
        check(shape_error < 3e-4 and offset_error < 1.5e-4, "independent native distance identities", name)
        close(shape_error,response["native_identity_audit"][name]["native_Tripp_shape_max_error_mag"], "saved native shape audit")
        close(offset_error,response["native_identity_audit"][name]["MU_MUMODEL_MURES_M0DIF_identity_max_error_mag"], "saved native redshift-bin offset audit")
        populated = {str(key):int(value) for key,value in f.IZBIN.value_counts().sort_index().items()}
        check(populated == case["populated_IZBIN_counts"], "native covariance informative-bin membership", name)
        lines = [line for line in (path/"bbc.COV").read_text().splitlines() if line and not line.startswith("#")]
        size = int(lines[0])
        cov = np.array([float(v) for v in lines[1:]]).reshape(size,size)
        check(np.isfinite(cov).all() and np.max(abs(cov-cov.T)) < 1e-8 and np.linalg.eigvalsh(cov).min() > -1e-8,
              "native conditional covariance numerical validity", name)
        identified = independent_native_gate(f,cov,case,log,name)
        frames[name] = f
        native_checks[name] = {"accepted_rows":len(f),"distance_shape_max_error_mag":shape_error,
                              "bin_offset_identity_max_error_mag":offset_error,"populated_bins":populated,
                              "independently_identified":identified,"nuisance_design_rank":case["nuisance_design"]["rank"],
                              "nuisance_design_columns":case["nuisance_design"]["columns"],
                              "host_mass_overlap_supported":case["host_mass_overlap_supported"],
                              "nuisance_native_bounds":case["nuisance_native_bounds"]}
    if "nominal" not in frames:
        check(response["status"] == "unsupported_nominal_map" and not response["contrasts"] and
              not response["analyzable_slope_support"], "failed nominal gate remains explicit")
        return {"native_BBC":native_checks,"response_status":response["status"],
                "independently_recomputed_contrasts":{}}
    nominal = frames["nominal"]
    edges = np.array(read(folder/"scaled-analysis-design.json")["redshift_bins"])
    common = set(nominal.index)
    for f in frames.values():
        common &= set(f.index)
    if len(frames) != 4:
        common = set()
    recomputed = {}
    for label,saved in response["contrasts"].items():
        case_name,support = label.rsplit("_",1)
        if support == "cases":
            case_name = label.removesuffix("_all_cases")
            ids = common
        else:
            ids = set(nominal.index) & set(frames["age_nominal" if case_name == "frozen_nominal" else case_name].index)
        target_name = "age_nominal" if case_name == "frozen_nominal" else case_name
        ids = sorted(ids)
        check(set(saved["common_occurrence_ids"]) == set(ids), "native response exact common support", label)
        n, a = nominal.loc[ids], frames[target_name].loc[ids]
        for column in ["SIM_HOSTLIB_SN_age","SIM_DLMAG","SIM_LIBID","SIM_ZCMB"]:
            close(n[column],a[column],"paired latent identity "+column,atol=2e-6,rtol=0)
        mu = native_distance_shape(a,bbc["cases"]["nominal"]["parameters"]) if case_name == "frozen_nominal" else a.MU.to_numpy()
        delta = (mu-a.SIM_DLMAG.to_numpy())-(n.MU-n.SIM_DLMAG).to_numpy()
        identified = native_checks["nominal"]["independently_identified"] and (
            case_name == "frozen_nominal" or native_checks[target_name]["independently_identified"])
        check(saved["native_nuisance_identified"] == identified,"native response identification gate",label)
        got = compare_response(n,delta,saved,edges,label,identified=identified)
        shift = 5*np.log10(relative_luminosity_distance(a.zHD)/relative_luminosity_distance(a.SIM_ZCMB))
        nshift = 5*np.log10(relative_luminosity_distance(n.zHD)/relative_luminosity_distance(n.SIM_ZCMB))
        compare_response(n,delta-shift+nshift,saved["observed_z_truth_sensitivity"],edges,label+": redshift truth",identified=identified)
        if got["slope_mag_per_Gyr"] is None:
            check(saved["fraction_of_injected_slope"] is None, "unsupported attenuation fraction remains null",label)
        else:
            close(got["slope_mag_per_Gyr"]/-.03,saved["fraction_of_injected_slope"],"injection fraction arithmetic")
        recomputed[label] = got
    compare_response(nominal,(nominal.MU-nominal.SIM_DLMAG).to_numpy(),response["closure"]["nominal"],edges,"nominal closure",
                     identified=native_checks["nominal"]["independently_identified"])
    check(response["analyzable_slope_support"] == any(row["slope_mag_per_Gyr"] is not None for row in recomputed.values()),
          "native analyzable support status")
    if suffix.startswith("-common"):
        n = fitres(bbc["cases"]["nominal"]["reference_config_overrides"]["datafile"])
        a = fitres(bbc["cases"]["age_nominal"]["reference_config_overrides"]["datafile"])
        check(n.index.equals(a.index), "common input occurrence identity")
        expected = set(tables["eval_nominal"].index)&set(tables["eval_age"].index)
        if suffix.endswith("-highmass-gauge"):
            before = len(expected)
            expected &= set(tables["eval_nominal"].index[tables["eval_nominal"].HOST_LOGMASS.ge(10.04)])
            expected &= set(tables["eval_age"].index[tables["eval_age"].HOST_LOGMASS.ge(10.04)])
            check(n.HOST_LOGMASS.ge(10.04).all() and a.HOST_LOGMASS.ge(10.04).all(),"high-mass input stratum exact threshold")
            stratum = bbc["highmass_stratum"]
            check(stratum["min_logmass"] == 10.04 and stratum["pre_stratum_common_candidates"] == before and
                  stratum["retained_common_candidates"] == len(expected) and stratum["candidate_loss"] == before-len(expected),
                  "high-mass stratum input denominator")
        check(set(n.index) == expected,"common input cohort derived only from declared selection")
        close(n.zHD,a.zHD,"common input observed-redshift geometry")
        check(np.array_equal(n.FIELD,a.FIELD),"common input field identity")
        check(bbc["cases"]["nominal"]["native_map_geometry"] == bbc["cases"]["age_nominal"]["native_map_geometry"],
              "common input native map geometry identity")
    return {"native_BBC":native_checks,"response_status":response["status"],
            "independently_recomputed_contrasts":recomputed}


def independent_native_gate(frame,cov,case,log,label):
    """Separate native completion from nuisance and informative-covariance rank."""
    informative = sorted(frame.loc[frame.M0DIFERR.lt(998),"IZBIN"].astype(int).unique())
    check(informative == case["informative_covariance_bins"], "informative native-bin identity",label)
    mineig = float(np.linalg.eigvalsh(cov[np.ix_(informative,informative)]).min()) if informative else None
    if mineig is None:
        check(case["informative_covariance_min_eigenvalue"] is None,"uninformative covariance explicit",label)
    else:
        close(mineig,case["informative_covariance_min_eigenvalue"],"native informative covariance eigenvalue")
    bins = sorted(frame.IZBIN.unique())
    logistic = 1/(1+np.exp(np.clip(-(frame.HOST_LOGMASS.to_numpy()-10)/.001,-700,700)))
    mode = case["nuisance_mode"]
    fixed = mode == "fixed_released_construction_values"
    gauge = mode == "fixed_gamma_highmass_width_colour_refit"
    check(mode in ["native_refit","fixed_released_construction_values","fixed_gamma_highmass_width_colour_refit"],
          "declared native nuisance mode",label)
    free = [] if fixed else (["alpha0","beta0"] if gauge else ["alpha0","beta0","gamma0"])
    columns = [] if fixed else [frame.x1,frame.c] + ([] if gauge else [.5-logistic])
    design = np.column_stack(columns + [(frame.IZBIN == value).astype(float) for value in bins])
    weighted = design/frame.MUERR.to_numpy()[:,None]
    norms = np.linalg.norm(weighted,axis=0)
    normalized = weighted/np.where(norms>0,norms,1)
    _, singular, right = np.linalg.svd(normalized,full_matrices=True)
    rank = int(np.sum(singular>singular[0]*1e-10))
    recorded = case["nuisance_design"]
    check(rank == recorded["rank"] and design.shape[1] == recorded["columns"],"independent native nuisance rank",label)
    close(singular,recorded["singular_values"],"native scaled design singular values")
    null = right[rank:]/np.where(norms>0,norms,1)[None,:]
    shapes = np.zeros((len(frame),len(null))) if fixed else design[:,:len(free)]@null[:,:len(free)].T
    shapes -= shapes.mean(axis=0,keepdims=True)
    shape_max = float(abs(shapes).max()) if shapes.size else 0.0
    close(shape_max,recorded["null_direction_centered_mu_shape_max"],"native null-direction shape")
    low = int(frame.HOST_LOGMASS.lt(10).sum())
    high = int(frame.HOST_LOGMASS.ge(10).sum())
    mixed = [int(value) for value in bins if frame.loc[frame.IZBIN.eq(value),"HOST_LOGMASS"].lt(10).any()
             and frame.loc[frame.IZBIN.eq(value),"HOST_LOGMASS"].ge(10).any()]
    check(low == recorded["low_mass_count"] and high == recorded["high_mass_count"] and
          mixed == recorded["mixed_mass_native_bins"],"native within-bin host-mass overlap",label)
    overlap = fixed or gauge or (low>0 and high>0 and len(mixed)>0)
    bound = "gamma0" in free and abs(case["parameters"]["gamma0"]["value"])>=.5-1e-5
    boundary = {key:bool(key in free and (case["parameters"][key]["value"]<=lo+1e-5 or
                    case["parameters"][key]["value"]>=hi-1e-5))
                for key,(lo,hi) in {"alpha0":(.02,.30),"beta0":(1.,6.),"gamma0":(-.5,.5)}.items()}
    for key,(lo,hi) in {"alpha0":(.02,.30),"beta0":(1.,6.),"gamma0":(-.5,.5)}.items():
        value = case["nuisance_native_bounds"][key]
        check(value["min"] == lo and value["max"] == hi and value["at_bound"] == boundary[key],
              "native nuisance boundary audit",[label,key])
    check(overlap == case["host_mass_overlap_supported"] and bound == case["gamma_at_native_bound"],
          "native host-step support and boundary gates",label)
    identified = (rank == design.shape[1] and "BAD_OUTPUT ERROR detected" not in log and mineig is not None and mineig>0
                  and overlap and not any(boundary.values()))
    check(identified == case["native_fit_identified"],"native identification status",label)
    if fixed:
        close([case["parameters"][key]["value"] for key in ["alpha0","beta0","gamma0"]],[.15,2.87,0],"fixed diagnostic values")
    if gauge:
        check(frame.HOST_LOGMASS.ge(10.04).all() and np.all(logistic == 1.),"high-mass numerical constant-step gauge")
        close(np.max(abs(1-logistic)),case["gamma_gauge_max_host_factor_deviation"],"recorded high-mass gauge numerical variation")
        check(case["parameters"]["gamma0"]["value"] == 0,"high-mass gauge gamma fixed to zero")
        for key in free:
            error = case["parameters"][key]["conditional_error"]
            check(np.isfinite(error) and error>0,"high-mass free width/colour conditional errors",[label,key])
    return bool(identified)


def audit_bootstrap(folder,out,work,k,stages):
    reports = {}
    for stage in stages:
        suffix = "" if stage == "selected_input" else "-"+stage
        baseline = read(out/f"scaled-response-k{k:03d}{suffix}.json")
        for mode in ["joint","training_only"]:
            path = out/f"scaled-bootstrap-{mode}-k{k:03d}{suffix}.json"
            record = read(path)
            identity(folder/"scaled_bootstrap.py",record["code_sha256"])
            check(record["mode"] == mode and record["shards"] == k and record["variant"] == suffix,
                  "bootstrap target identity",rel(path))
            if record.get("status") == "not_run_support_gate":
                check(not baseline["analyzable_slope_support"] and record["replicates"] == 0 and
                      record["planned_replicates"] in ([200] if mode == "joint" else [100,200]) and bool(record["reason"]),
                      "bootstrap nonexecution is an explicit unsupported-baseline gate",rel(path))
                reports[stage+":"+mode] = {"status":"not_run_support_gate","planned":record["planned_replicates"],"executed":0}
                continue
            check(baseline["analyzable_slope_support"],"bootstrap baseline supported",rel(path))
            count = record["replicates"]
            planned = 200
            if stage in ["common-fixed","common-highmass-gauge"]:
                design_path = folder/("scaled-grid-audit-design.json" if stage == "common-fixed" else "scaled-highmass-design.json")
                plan = read(design_path)
                statement = plan["fixed_nuisance_conditional_control"]["bootstrap"] if stage == "common-fixed" else plan["uncertainty"]
                check("200 joint" in statement and "100 training-only" in statement,"frozen conditional bootstrap allocation")
                planned = 200 if mode == "joint" else 100
            check(count == planned and record.get("planned_replicates",planned) == planned,"bootstrap planned replicate count completed",rel(path))
            hashes = record["replicate_result_hashes"]
            check(set(hashes) == {str(i) for i in range(count)},"bootstrap no missing replicate records",rel(path))
            replicates = []
            for index,digest in hashes.items():
                directory = work/"bootstrap"/f"k{k:03d}{suffix}"/mode/f"r{int(index):03d}"
                identity(directory/"result.json",digest)
                item = read(directory/"result.json")
                check(item["index"] == int(index) and item["seed"] == 972000+int(index) and item["mode"] == mode,
                      "bootstrap replicate identity",rel(directory))
                for source in item.get("inputs",{}).values():
                    identity(source["path"],source["sha256"])
                for case,details in item.get("native_cases",{}).items():
                    identity(directory/case/"bbc.log",details["log_sha256"])
                replicates.append(item)
            actual = pd.Series([item["status"] for item in replicates]).value_counts().to_dict()
            check(actual == record["replicate_status_counts"],"bootstrap all failure/support statuses counted",rel(path))
            labels = {name for item in replicates for name in item["contrasts"]}
            check(labels == set(record["summary"]),"bootstrap summary contrast completeness",rel(path))
            for label, summary in record["summary"].items():
                valid = [item["contrasts"][label] for item in replicates if label in item["contrasts"] and
                         item["contrasts"][label]["slope_mag_per_Gyr"] is not None]
                slopes = np.array([item["slope_mag_per_Gyr"] for item in valid])
                check(len(valid) == summary["supported_replicates"],"bootstrap supported slope count",label)
                check_distribution(slopes,summary["slope_sd_mag_per_Gyr"],summary["percentile_95_mag_per_Gyr"],label)
                if len(slopes)>1:
                    close(np.quantile(slopes/-.03,[.025,.975]),summary["fraction_of_injected_slope_95"],"bootstrap signed response interval")
                else:
                    check(summary["fraction_of_injected_slope_95"] is None,"unsupported response fraction uncertainty null")
                if valid:
                    close(np.median([item["n"] for item in valid]),summary["median_common_n"],"bootstrap median common support")
                else:
                    check(summary["median_common_n"] is None,"bootstrap no supported common summary")
                for row in summary["redshift_means"]:
                    means = [item["centered_redshift_means"][row["bin"]]["mean_mag"] for item in valid
                             if item["centered_redshift_means"][row["bin"]]["mean_mag"] is not None]
                    check(len(means) == row["supported_replicates"],"bootstrap supported redshift count")
                    check_distribution(means,row["sd_mag"],row["percentile_95_mag"],label+":bin")
                covariance = summary["centered_bin_covariance"]
                active = [i for i,row in enumerate(baseline["contrasts"].get(label,{}).get("centered_redshift_means",[]))
                          if row["mean_mag"] is not None]
                check(active == covariance["baseline_supported_bin_indices"],"bootstrap centered covariance baseline bins")
                complete = [[item["centered_redshift_means"][i]["mean_mag"] for i in active] for item in valid if active and
                            all(item["centered_redshift_means"][i]["mean_mag"] is not None for i in active)]
                check(len(complete) == covariance["complete_replicates"],"bootstrap centered covariance completeness")
                if len(complete)>1:
                    close(np.atleast_2d(np.cov(complete,rowvar=False,ddof=1)),covariance["covariance_mag2"],"bootstrap centered-bin covariance")
                else:
                    check(covariance["covariance_mag2"] is None,"unsupported centered-bin covariance null")
            for label,summary in record["paired_estimand_differences"].items():
                after,before = label.split("_minus_")
                pairs = [(item["contrasts"].get(after+"_all_cases",{}).get("slope_mag_per_Gyr"),
                          item["contrasts"].get(before+"_all_cases",{}).get("slope_mag_per_Gyr")) for item in replicates]
                differences = [a-b for a,b in pairs if a is not None and b is not None]
                check(len(differences) == summary["complete_paired_replicates"],"paired estimand bootstrap completeness")
                check_distribution(differences,summary["sd_mag_per_Gyr"],summary["percentile_95_mag_per_Gyr"],label)
            reports[stage+":"+mode] = {"status":"executed","executed":count,"replicate_statuses":actual,
                                        "supported":{name:item["supported_replicates"] for name,item in record["summary"].items()}}
    return reports


def check_distribution(values,sd,interval,label):
    if len(values)>1:
        close(np.std(values,ddof=1),sd,"bootstrap standard deviation "+label)
        close(np.quantile(values,[.025,.975]),interval,"bootstrap percentile interval "+label)
    else:
        check(sd is None and interval is None,"unsupported bootstrap uncertainty null "+label)


def audit_redshift_summary(folder,out,work,k):
    record = read(out/f"scaled-redshift-summary-k{k:03d}.json")
    identity(folder/"scaled_redshift_summary.py",record["code_sha256"])
    identity(folder/"scaled-redshift-summary-design.json",record["design_sha256"])
    identities(record["dependencies_sha256"])
    check(record["status"] == "computed" and record["exploratory"] is True,"post-pattern redshift contrast scope explicit")
    design = read(folder/"scaled-redshift-summary-design.json")
    check("after inspection" in design["status"],"redshift-summary timing remains exploratory")
    labels = ["frozen_nominal","age_nominal","age_literal","age_retained"]
    def bins(item,label):
        entry = item.get("contrasts",{}).get(label+"_all_cases",{})
        if entry.get("slope_mag_per_Gyr") is None:return None
        return [entry["centered_redshift_means"][i]["mean_mag"] for i in range(4)]
    checked = {}
    for variant,saved in record["variants"].items():
        point = read(out/f"scaled-response-k{k:03d}-{variant}.json")
        hilo = {label:bins(point,label)[3]-bins(point,label)[0] for label in labels}
        for label,value in hilo.items():close(value,saved["point_high_minus_low_mag"][label],"exploratory point high-low difference")
        for key,before in [("frozen","frozen_nominal"),("nominal_refit","age_nominal")]:
            close(hilo["age_retained"]-hilo[before],saved[f"point_retained_minus_{key}_mag"],"paired point redshift contrast")
        counts = {}
        for mode,summary in saved["bootstrap"].items():
            boot = read(out/f"scaled-bootstrap-{mode}-k{k:03d}-{variant}.json")
            reps = []
            for index,digest in boot["replicate_result_hashes"].items():
                path = work/"bootstrap"/f"k{k:03d}-{variant}"/mode/f"r{int(index):03d}"/"result.json"
                identity(path,digest)
                reps.append(read(path))
            check(summary["replicates"] == len(reps) and summary["replicate_status_counts"] == boot["replicate_status_counts"],
                  "redshift bootstrap full denominator")
            highlow = {label:[] for label in labels}
            changes = {"retained_minus_frozen":[],"retained_minus_nominal_refit":[]}
            complete = []
            for item in reps:
                vectors = {label:bins(item,label) for label in labels}
                differences = {label:vector[3]-vector[0] for label,vector in vectors.items()
                               if vector is not None and vector[0] is not None and vector[3] is not None}
                for label,value in differences.items():highlow[label].append(value)
                for name,before in [("retained_minus_frozen","frozen_nominal"),("retained_minus_nominal_refit","age_nominal")]:
                    if "age_retained" in differences and before in differences:
                        changes[name].append(differences["age_retained"]-differences[before])
                if all(vector is not None and all(v is not None for v in vector) for vector in vectors.values()):
                    complete.append([value for label in labels for value in vectors[label]])
            for key,values in [("high_minus_low",highlow),("paired_changes",changes)]:
                for label,sample in values.items():
                    row = summary[key][label]
                    check(len(sample) == row["supported_replicates"],"redshift paired bootstrap support")
                    check_distribution(sample,row["sd_mag"],row["percentile_95_mag"],"redshift "+label)
            covariance = summary["joint_bin_covariance"]
            check(covariance["ordering"] == [f"{label}:bin{i}" for label in labels for i in range(4)] and
                  covariance["complete_replicates"] == len(complete),"joint redshift-estimand covariance ordering and support")
            if len(complete)>1:close(np.cov(complete,rowvar=False,ddof=1),covariance["covariance_mag2"],"joint paired redshift covariance")
            else:check(covariance["covariance_mag2"] is None,"unsupported joint redshift covariance remains null")
            counts[mode] = {"attempted":len(reps),"complete_joint_vectors":len(complete)}
        checked[variant] = {"point_high_minus_low_mag":hilo,"bootstrap":counts}
    check(set(checked) == {"common-fixed","common-highmass-gauge"},"declared exploratory redshift variants complete")
    return {"exploratory":True,"variants":checked,"scope":"Post-pattern contrast, conditional supported-replicate distributions; no cosmological correction or verified interval coverage."}


def sky_vectors(ra,dec):
    r,d = np.deg2rad(ra),np.deg2rad(dec)
    return np.column_stack([np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)])


def vector_angles(vector,others):
    return np.rad2deg(np.arctan2(np.linalg.norm(np.cross(vector,others),axis=1),others@vector))*3600


def audit_matching_controls(out,work,hosts,candidates):
    result = read(out/"matching-controls.json")
    identities(result["dependencies_sha256"])
    check(result["status"] == "passed","infrared positional controls complete")
    design = read(CODE/"infrared_photometry/matching-control-design.json")
    saved = pd.read_csv(work/"matching-controls.csv",dtype={"SNID":str,"deep_ID":str,"nearest_source_key":str})
    radii = np.array(design["radii_arcsec"])
    offsets = [(r,b) for r in design["position_controls"]["radii_arcsec"]
               for b in design["position_controls"]["bearings_degrees"]]
    check(len(saved) == result["rows"] == len(hosts)*3*len(radii),"infrared positional control dimensions")
    check(result["host_catalogue_configurations"] == len(hosts)*3 and
          result["shifted_positions_per_configuration"] == len(offsets),"infrared control count accounting")
    check(max(r for r,b in offsets)+max(radii)<design["position_controls"]["candidate_query_radius_arcsec"],
          "shifted control searches inside acquired cones")
    rng = np.random.default_rng(927061)
    max_angle_error = 0.0
    for catalogue in dict.fromkeys(row["catalogue"] for row in result["summaries"]):
        data = candidates[candidates.catalogue == catalogue]
        fields = data.field.unique()
        check(len(fields) == 1,"infrared catalogue single declared field",catalogue)
        parent = hosts[hosts.field == fields[0]]
        points = sky_vectors(data.ra.to_numpy(),data.dec.to_numpy())
        count = np.zeros((len(parent),len(radii)),dtype=int)
        shifted = np.zeros((len(parent),len(offsets),len(radii)),dtype=int)
        rows = saved[saved.catalogue == catalogue].set_index(["SNID","radius_arcsec"])
        for j,host in enumerate(parent.itertuples(index=False)):
            center = sky_vectors([host.deep_RA],[host.deep_DEC])[0]
            angles = vector_angles(center,points)
            keep = angles<=75+1e-6
            nearby = points[keep]
            angles = angles[keep]
            count[j] = (angles[:,None]<=radii).sum(axis=0)
            nearest = int(np.argmin(angles)) if len(angles) else None
            ra,dec = np.deg2rad(host.deep_RA),np.deg2rad(host.deep_DEC)
            east = np.array([-np.sin(ra),np.cos(ra),0.])
            north = np.array([-np.sin(dec)*np.cos(ra),-np.sin(dec)*np.sin(ra),np.cos(dec)])
            for q,(offset,bearing) in enumerate(offsets):
                distance,theta = np.deg2rad(offset/3600),np.deg2rad(bearing)
                displaced = np.cos(distance)*center+np.sin(distance)*(np.cos(theta)*north+np.sin(theta)*east)
                shifted[j,q] = (vector_angles(displaced,nearby)[:,None]<=radii).sum(axis=0)
            for q,radius in enumerate(radii):
                row = rows.loc[(host.SNID,radius)]
                check(row.candidates_in_radius == count[j,q] and row.shifted_with_match == np.sum(shifted[j,:,q]>0)
                      and row.shifted_with_multiple == np.sum(shifted[j,:,q]>1),"infrared independent observed/shifted counts")
                check(row.candidates_within75 == len(angles),"infrared local acquisition counts")
                if nearest is None:
                    check(pd.isna(row.nearest_source_key) and pd.isna(row.nearest_separation_arcsec),"empty infrared neighbourhood explicit")
                else:
                    check(row.nearest_source_key == data.loc[keep].iloc[nearest].source_key,"independent infrared nearest identifier")
                    error = abs(float(row.nearest_separation_arcsec)-angles[nearest])
                    max_angle_error = max(max_angle_error,error)
                    check(error<1e-7,"independent infrared spherical angle")
        resamples = rng.integers(0,len(parent),size=(2000,len(parent)))
        for q,radius in enumerate(radii):
            row = next(item for item in result["summaries"] if item["catalogue"] == catalogue and item["radius_arcsec"] == radius)
            detected = (count[:,q]>0).astype(float)
            control = (shifted[:,:,q]>0).mean(axis=1)
            check(row["hosts"] == len(parent) and row["actual_with_match"] == detected.sum() and
                  row["actual_unique"] == np.sum(count[:,q] == 1) and row["actual_multiple"] == np.sum(count[:,q]>1),
                  "infrared positional summary counts")
            close(control.mean(),row["shifted_match_fraction"],"infrared shifted match fraction")
            close((detected-control).mean(),row["actual_minus_shifted_fraction"],"infrared positional excess")
            close(np.quantile((detected-control)[resamples].mean(axis=1),[.025,.975]),row["conditional_host_bootstrap_95"],
                  "infrared positional host-bootstrap interval")
    return {"independent_geometry_max_error_arcsec":max_angle_error,"configurations":len(hosts)*3,
            "positions_per_configuration":len(offsets),"scope":"Conditional local catalogue coincidence; no calibrated association probability or footprint completeness."}


def native_distance_shape(frame, parameters):
    alpha,beta,gamma = [parameters[key]["value"] for key in ["alpha0","beta0","gamma0"]]
    logistic = 1/(1+np.exp(np.clip(-(frame.HOST_LOGMASS.to_numpy()-10)/.001,-700,700)))
    return -2.5*np.log10(frame.x0.to_numpy())+alpha*frame.x1.to_numpy()-beta*frame.c.to_numpy()-gamma*(.5-logistic)-frame.biasCor_mu.to_numpy()


def audit_bbc_configuration(case,path):
    reference = ROOT/".work/survey-physics/SNDATA_ROOT/sample_input_files/DES-SN5YR/base_files/bbc/BBC_des5yr.input"
    identity(reference,case["reference_input_sha256"])
    original = reference.read_text().split("#END_YAML",1)[1]
    actual = path.read_text()
    removed = {"cid_reject_file","varname_pIa","simfile_ccprior","idsurvey_list_probcc0"}
    allowed = {"datafile","simfile_biascor","prefix","surveygroup_biascor","u13","p13","zmin","ndump_nobiascor"}
    mode = case["nuisance_mode"]
    if mode == "fixed_released_construction_values":allowed |= {"u1","u2","u5","p5"}
    if mode == "fixed_gamma_highmass_width_colour_refit":allowed |= {"u5","p5"}
    declared = case["reference_config_overrides"]
    check(set(declared) == allowed,"only declared native scope/nuisance overrides")
    def parsed(text):
        output = {}
        retained = []
        for line in text.splitlines():
            match = re.match(r"^([A-Za-z0-9_]+)=(.*)$",line)
            key = match[1] if match else None
            if key is not None:output[key] = match[2].strip()
            if key not in allowed|removed:retained.append(line)
        return output,retained
    old,old_lines = parsed(original)
    new,new_lines = parsed(actual)
    check(old_lines == new_lines,"all unmodified native scientific configuration lines retained")
    check(not (removed&set(new)),"declared pure-Ia configuration scope")
    for key,value in declared.items():check(new.get(key) == str(value),"native override exact value",key)
    check(new.get("opt_biascor") == old.get("opt_biascor") == "4336","native multidimensional correction option retained")
    check("snrmin_sigint_biascor" not in new,"native intrinsic-scatter SNR minimum not overridden")
    for key,value in {"prefix":"bbc","surveygroup_biascor":"'DES(zbin=0.075)'","u13":"0","p13":"0","zmin":".05","ndump_nobiascor":"0"}.items():
        check(declared[key] == value,"declared native mock-scope fixed value",key)


def audit_survey_review(bootstrap,redshift):
    record = read(RESULTS/"galaxy_validation/scaled-bbc-method-review.json")
    identities(record["source_sha256"])
    check(record["read_only_source_review"] is True and record["large_run_review_pending"] == [],
          "independent native review completed with no pending calculation")
    check(bool(record["remaining_scientific_limitations"]),"independent review retains scientific limitations")
    final = record["k016_final_bootstrap_and_target_review"]
    identities(final["dependencies"])
    check(final["status"] == "passed","independent final native numerical review passed")
    executed = {key:item for key,item in bootstrap.items() if item["status"] == "executed"}
    check(set(final["bootstrap"]) == set(executed),"independent bootstrap review scope complete")
    for key,item in final["bootstrap"].items():
        actual = executed[key]
        variant,mode = key.split(":")
        check(item["replicates"] == actual["executed"] and item["status_counts"] == actual["replicate_statuses"],
              "independent review matches reconstructed replicate denominator")
        check(item["all_case_slope_supported"] == actual["supported"]["age_retained_all_cases"],
              "independent review matches reconstructed slope support")
        check(item["full_16bin_covariance_supported"] ==
              redshift["variants"][variant]["bootstrap"][mode]["complete_joint_vectors"],
              "independent review matches reconstructed covariance support")
        check(item["max_summary_difference"]<1e-10,"independent summary arithmetic agrees")
    raw = final["raw_output_audits"]
    check(len(raw) == 32 and all(item["max_slope_or_bin_difference"]<1e-10 for item in raw),
          "independent raw native-output audits complete")
    return {"completed":True,"source_identities":len(record["source_sha256"]),
            "final_dependencies":len(final["dependencies"]),"raw_output_audits":len(raw),
            "replicates_reconstructed":sum(item["replicates"] for item in final["bootstrap"].values()),
            "remaining_scientific_limitations":record["remaining_scientific_limitations"]}


def audit_capacity_repair(folder,out,record_name="scaled-native-capacity-replay.json"):
    build = read(out/"scaled-native-capacity-build.json")
    replay = read(out/record_name)
    identity(folder/"scaled_build.py",build["code_sha256"])
    identity(folder/"scaled_capacity_replay.py",replay["code_sha256"])
    identity(out/"scaled-native-capacity-build.json",replay["build_record_sha256"])
    identities(build["source"])
    for key in ["patch","previous_high_snr_capacity_patch","build_log","original_executable","executable"]:
        identity(build[key]["path"],build[key]["sha256"])
    for key in ["patch","executable"]:identity(replay[key]["path"],replay[key]["sha256"])
    check(build["status"] == "built" and replay["status"] == "passed","native storage repair built and replayed")
    if record_name == "scaled-native-capacity-restoration.json":
        check(replay["restored_fixtures"] is True,"capacity replay uses genuinely restored fixtures")
        identity(folder/"scaled_restore_capacity.py",replay["restoration_code_sha256"])
        manifest = replay["restoration_manifest"]
        identity(manifest["path"],manifest["sha256"])
        restoration = read(manifest["path"])
        identity(folder/"scaled_restore_capacity.py",restoration["code_sha256"])
        identity(restoration["original_executable"],restoration["original_executable_sha256"])
        check(restoration["original_native_runs"] == 24 and restoration["pilot_replicates"] == 20 and
              restoration["bootstrap_seed_base"] == 972000,"bounded restored original-run design")
        original_root = Path(manifest["path"]).parent/"original"
        originals = sorted(original_root.glob("*/original-run.json"))
        check(len(originals) == 24,"all original-binary restoration runs present")
        failed = []
        for path in originals:
            run = read(path)
            identity(path.parent/"bbc.input",run["input_sha256"])
            identity(path.parent/"bbc.log",run["log_sha256"])
            identity(restoration["original_executable"],run["executable_sha256"])
            log = (path.parent/"bbc.log").read_text()
            complete = "Done." in log[-1000:] and "FATAL ERROR ABORT" not in log
            check(complete == run["graceful"],"restored original native completion from log")
            if run["name"].startswith("pilot-"):
                check(run["seed"] == 972000+int(run["name"][-3:]),"unchanged restored bootstrap seed")
            else:check(run["seed"] is None and complete,"all restored original baseline cases completed")
            if not complete:failed.append(run["name"])
        check(failed == restoration["failures"] == [f"pilot-r{i:03d}" for i in [5,7,8,16,17,19]],
              "restoration reproduces all six original pilot failures")
    old = ROOT/".work/survey-physics/SNANA-legacy-eff/src"
    new = ROOT/".work/survey-physics/scaled-bbc/native-capacity/src"
    check((old/"SALT2mu.c").read_bytes() == (new/"SALT2mu.c").read_bytes(),"native BBC scientific source unchanged by storage repair")
    delta = "".join(difflib.unified_diff((old/"sntools.c").read_text().splitlines(keepends=True),
                 (new/"sntools.c").read_text().splitlines(keepends=True),fromfile="a/src/sntools.c",tofile="b/src/sntools.c"))
    check(delta == Path(build["patch"]["path"]).read_text(),"native repaired source exactly matches saved patch")
    removed = [line[1:].strip() for line in delta.splitlines() if line.startswith("-") and not line.startswith("---")]
    added = [line[1:].strip() for line in delta.splitlines() if line.startswith("+") and not line.startswith("+++")]
    check(removed == ["#define MXSTORE_PULL 100","else {"],"capacity patch removes no scientific calculation")
    permitted = {"", "// Storage covers the unchanged descending scan through sigint_min.",
        "// Two spare entries preserve the existing increment-before-bound check.",
        "const int MXSTORE_PULL = (int)fmax(100.0,", "ceil((sigTmp_hi - sigint_min) / sigint_bin) + 2.0);",
        "else {", "if ( NBIN_SIGINT >= 100 ) {", "}",
        'printf(" CAPACITY_SCAN: %s N=%d start=%.17g mean=%.17g STD=%.17g steps=%d result=%.17g lower_floor=1\\n",',
        'callFun, N_LIST, sigTmp_hi, AVG_MURES, STD_MURES_ORIG, NBIN_SIGINT, sigint_min);',
        'printf(" CAPACITY_SCAN: %s N=%d start=%.17g mean=%.17g STD=%.17g steps=%d result=%.17g lower_floor=0\\n",',
        'callFun, N_LIST, sigTmp_hi, AVG_MURES, STD_MURES_ORIG, NBIN_SIGINT, sigint);'}
    check(set(added)<=permitted,"capacity patch adds only storage sizing and diagnostic output")
    counts = {"unchanged_success":0,"repaired_storage_abort":0,"persistent_zero_MAD":0}
    for item in replay["records"]:
        for key in ["original_log","original_input"]:identity(item[key]["path"],item[key]["sha256"])
        original_log = Path(item["original_log"]["path"]).read_text()
        original = Path(item["original_input"]["path"]).parent
        case = item["replay"]
        path = Path(case["folder"])
        identity(path/"bbc.log",case["log_sha256"])
        identity(path/"bbc.input",case["input_sha256"])
        check((path/"bbc.input").read_bytes() == Path(item["original_input"]["path"]).read_bytes(),"storage repair identical input")
        audit_bbc_configuration(case,path/"bbc.input")
        log = (path/"bbc.log").read_text()
        if item["original_success"]:
            check(case["graceful"] and "Done." in original_log[-1000:],"previous successful native case remains successful")
            left = [line for line in (original/"bbc.FITRES").read_text().splitlines() if line.startswith(("SN:","VARNAMES:"))]
            right = [line for line in (path/"bbc.FITRES").read_text().splitlines() if line.startswith(("SN:","VARNAMES:"))]
            check(left == right and fitres(original/"bbc.FITRES").equals(fitres(path/"bbc.FITRES")),"storage repair exact science table identity")
            check((original/"bbc.COV").read_bytes() == (path/"bbc.COV").read_bytes(),"storage repair exact covariance identity")
            counts["unchanged_success"] += 1
        else:
            check("FATAL ERROR ABORT" in original_log,"original native pilot failure retained")
            if item["original_failure_type"] == "scan_capacity":
                check("NBIN_SIGINT=100 exceeds bound MXSTORE_PULL=100" in original_log and case["graceful"],
                      "native array-capacity failure repaired separately")
                check(bool(item["expanded_scan_cells"]),"extended scan diagnostics present")
                for scan in item["expanded_scan_cells"]:
                    capacity = max(100,int(np.ceil((scan["start"]+.3)/.01)+2))
                    check(capacity == scan["requested_capacity"] and 100<=scan["steps"]<capacity,"unchanged scan fits allocated capacity")
                    check(scan["lower_floor"] == 1 and scan["result"] == -.3,"repaired cell retains original lower-floor result")
                counts["repaired_storage_abort"] += 1
            elif item["original_failure_type"] == "zero_MAD_pull":
                check("Invalid stdPull = 0.000000" in original_log and "Invalid stdPull = 0.000000" in log
                      and not case["graceful"],"scientific zero-MAD gate remains failed")
                counts["persistent_zero_MAD"] += 1
            else:check(False,"unclassified original pilot failure",item["label"])
        if case["graceful"]:
            for filename,digest in case["outputs"].items():identity(path/filename,digest)
    check(counts == {"unchanged_success":4,"repaired_storage_abort":1,"persistent_zero_MAD":5},"complete pilot replay accounting")
    return counts


def relative_luminosity_distance(redshift):
    # Independent fixed Gauss quadrature, H0 cancels in the distance ratio.
    z = np.asarray(redshift)
    nodes,weights = np.polynomial.legendre.leggauss(32)
    x = z[:,None]*(nodes+1)/2
    return (1+z)*z/2*np.sum(weights/np.sqrt(.315*(1+x)**3+.685),axis=1)


def compare_response(frame, residual, record, edges, label, identified=True):
    weights = 1/frame.MUERR.to_numpy()**2
    bins = np.digitize(frame.zHD.to_numpy(),edges)-1
    center = np.average(residual,weights=weights)
    close(center,record["raw_weighted_global_offset_mag"],"native response global normalization "+label)
    allowed = []
    for i, row in enumerate(record["centered_redshift_means"]):
        keep = bins == i
        check(int(keep.sum()) == row["n"],"native response bin count "+label)
        if keep.sum() >= 10:
            allowed.append(i)
            close(np.average(residual[keep]-center,weights=weights[keep]),row["mean_mag"],"native response redshift mean "+label)
        else:
            check(row["mean_mag"] is None,"unsupported native mean remains null "+label)
    keep = np.isin(bins,allowed)
    check(int(keep.sum()) == record["slope_n"] and len(frame) == record["n"],"native response slope support "+label)
    slope = None
    if identified and keep.sum() >= 30:
        x = np.column_stack([frame.SIM_HOSTLIB_SN_age.to_numpy()[keep]]+[(bins[keep] == i).astype(float) for i in allowed])
        if np.linalg.matrix_rank(x) == x.shape[1]:
            fit = np.linalg.lstsq(x*np.sqrt(weights[keep,None]),residual[keep]*np.sqrt(weights[keep]),rcond=None)[0]
            slope = float(fit[0])
    if slope is None:
        check(record["slope_mag_per_Gyr"] is None,"unsupported native slope remains null "+label)
    else:
        close(slope,record["slope_mag_per_Gyr"],"native response slope "+label,atol=2e-9)
    return {"n":len(frame),"slope_n":int(keep.sum()),"slope_mag_per_Gyr":slope}


def navigation():
    scripts = list((CODE/"nebular_dust").glob("*.py")) + list((CODE/"infrared_resolution").glob("*.py"))
    scripts += list((CODE/"infrared_photometry").glob("*.py"))
    scripts += list((CODE/"survey_physics").glob("scaled*.py"))
    scripts += [Path(__file__), ROOT/"validation/record_physical_extension.py"]
    for path in scripts:
        ast.parse(path.read_text(), filename=str(path))
        sha(path)
    documents = list((CODE/"nebular_dust").glob("*.md")) + list((CODE/"infrared_resolution").glob("*.md"))
    documents += list((CODE/"infrared_photometry").glob("*.md"))
    documents += list((CODE/"survey_physics").glob("scaled*.md"))
    documents += [NOTES/"nebular-dust-results.md", NOTES/"infrared-resolution-results.md",
                  NOTES/"survey-bbc-support-results.md",NOTES/"infrared-photometry-results.md"]
    documents += list(NOTES.glob("scaled*bbc*.md")) + list(NOTES.glob("survey-bbc*.md"))
    documents = sorted(set(documents))
    links = 0
    for path in documents:
        check(path.exists(), "new branch documentation exists", rel(path))
        if not path.exists():
            continue
        sha(path)
        for target in re.findall(r"\]\(([^\s)]+)(?:\s+\"[^\"]*\")?\)", path.read_text()):
            if target.startswith(("http:","https:","mailto:","#","data:")):
                continue
            target = unquote(target.strip("<>").split("#",1)[0])
            check((path.parent/target).exists(), "new local Markdown link", {"document":rel(path),"target":target})
            links += 1
    SECTIONS["navigation"] = {"parsed_new_python_files": len(scripts), "new_markdown_files":len(documents), "local_links":links}


def main():
    for function in [first_extension, legacy_catalogue_label, nebular, infrared, infrared_photometry, scaled_survey, navigation]:
        try:
            function()
        except Exception as exc:
            FAILURES.append({"check": function.__name__, "exception":type(exc).__name__, "detail":str(exc)})
    code_files = list((CODE/"nebular_dust").glob("*")) + list((CODE/"infrared_resolution").glob("*")) + list((CODE/"survey_physics").glob("scaled*"))
    code_files += list((CODE/"infrared_photometry").glob("*"))
    result_files = list((RESULTS/"nebular_dust").glob("*.json")) + list((RESULTS/"infrared_resolution").glob("*.json")) + list((RESULTS/"survey_physics").glob("scaled*.json"))
    result_files += list((RESULTS/"galaxy_validation").glob("scaled*review.json"))
    result_files += list((RESULTS/"infrared_photometry").glob("*.json"))
    code_hashes = {rel(path):sha(path) for path in code_files if path.is_file()}
    result_hashes = {rel(path):sha(path) for path in result_files if path.is_file()}
    checker_hash = sha(Path(__file__))
    for path, (_,size,mtime) in CACHE.items():
        stat = path.stat()
        check((size,mtime) == (stat.st_size,stat.st_mtime_ns), "unchanged audit snapshot", rel(path))
    report = {"completed_utc":datetime.now(timezone.utc).isoformat(),
        "status":"passed" if not FAILURES else "failed", "passed":not FAILURES,
        "checker_sha256":checker_hash, "first_extension_preservation":SECTIONS.get("first_extension_preservation"),
        "sections":SECTIONS, "code_sha256":code_hashes, "result_sha256":result_hashes,
        "source_sha256":{rel(path):values[0] for path,values in CACHE.items()
                         if path.is_relative_to(ROOT/".work")},
        "files_hashed":len(CACHE), "bytes_hashed":sum(x[1] for x in CACHE.values()), "failures":FAILURES,
        "scope":"Identity, execution and independent arithmetic checks; numerical success is not empirical proof of an age correction or a new cosmological measurement."}
    REPORT.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":report["status"],"files_hashed":len(CACHE),"failures":FAILURES},indent=2))
    raise SystemExit(0 if report["passed"] else 1)


if __name__=="__main__":main()
