#!/usr/bin/env python3
"""Audit four physical-research branches using only NumPy and pandas.

Run after the branch validators. Source identities, catalogue membership,
interval summaries and matrix products are checked independently here. This
does not certify stellar physics or imply an observational cosmology correction.
"""
from __future__ import annotations

import ast
import hashlib
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "studies/host_ages"
CODE, RESULTS = BASE / "code", BASE / "results"
BRANCHES = ("galaxy_validation", "host_transport", "physical_ages", "survey_physics")
REPORT = ROOT / "validation/reports/physical-program.json"
CACHE: dict[Path, tuple[str, int, int]] = {}
FAILURES: list[dict] = []
SECTIONS: dict = {}


def relative(path):
    path = Path(path).resolve()
    return str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)


def sha(path):
    path = Path(path).resolve()
    stat = path.stat()
    if path not in CACHE:
        with path.open("rb") as stream:
            digest = hashlib.file_digest(stream, "sha256").hexdigest()
        CACHE[path] = (digest, stat.st_size, stat.st_mtime_ns)
    else:
        assert CACHE[path][1:] == (stat.st_size, stat.st_mtime_ns), f"Changed during audit: {path}"
    return CACHE[path][0]


def read(path):
    sha(path)
    return json.loads(Path(path).read_text())


def check(condition, label, detail=None):
    if not bool(condition):
        FAILURES.append({"check": label, "detail": detail})


def identity(path, expected):
    check(sha(path) == expected, "sha256", relative(path))


def identities(mapping, base=ROOT):
    for path, digest in mapping.items():
        identity(base / path, digest)


def quantiles(series):
    x = pd.to_numeric(series, errors="coerce").dropna()
    return dict(zip(("p05", "median", "p95"), map(float, np.quantile(x, [.05, .5, .95])))) if len(x) else None


def same_numbers(actual, expected, label, tolerance=1e-9):
    if isinstance(actual, dict):
        check(isinstance(expected, dict) and actual.keys() == expected.keys(), label + ":keys")
        if isinstance(expected, dict):
            for key in actual.keys() & expected.keys():
                same_numbers(actual[key], expected[key], label + ":" + key, tolerance)
    elif actual is None or expected is None:
        check(actual is expected, label, {"actual": actual, "expected": expected})
    else:
        check(np.isclose(actual, expected, atol=tolerance, rtol=1e-10), label,
              {"actual": actual, "expected": expected})


def galaxy():
    code, out, work = CODE / "galaxy_validation", RESULTS / "galaxy_validation", ROOT / ".work/galaxy-validation"
    validation = read(out / "validation.json")
    check(validation["passed"], "galaxy validator passed")
    identity(code / "validate.py", validation["code_sha256"])
    identities(validation["result_hashes"], out)
    lock = read(code / "source-lock.json")["files"]
    identities(lock, work)
    check(len(lock) == validation["source_bytes_verified"], "galaxy source count")
    provenance = read(out / "diagnostics-summary.json")["provenance"]
    for key in ("code_sha256", "design_sha256"):
        identities(provenance[key], code)
    identities(provenance["dependencies_sha256"], CODE / "environment_validation")
    for key in ("input_sha256", "output_sha256"):
        identities(provenance[key])
    d = pd.read_csv(work / "host-spectrum-age-ledger.csv", dtype={"specObjID": str, "source_id": str})
    check(d.source_id.is_unique, "unique galaxy host associations")
    check(len(d) == validation["matched_host_rows"], "galaxy row count")
    check(d.specObjID.nunique() == validation["distinct_spectra"], "distinct spectra")
    check(d.position_z_accepted.all() and d.unique_photometric_counterpart.all(), "host association gates")
    check((d.separation_arcsec <= 3).all() and (abs(d.delta_z) < .001).all(), "host sky and redshift separation")
    folds = pd.read_csv(work / "brightness-heldout.csv")
    check(folds.groupby("physical_host").fold.nunique().max() == 1, "no physical-host heldout leakage")
    means = pd.read_csv(work / "remeasured-dn4000.csv")
    means = means[means.valid]
    ratio_error = float(abs(means.raw_red_Fnu_uJy / means.raw_blue_Fnu_uJy - means.raw_Dn4000).max())
    check(ratio_error < 1e-12, "independent observed spectral ratio identity")
    check(abs(validation["independent_SLSQP_minus_Minuit_NLL"]) < 1e-5, "galaxy independent likelihood closure")
    SECTIONS["galaxy_observations"] = {"source_files": len(lock), "host_associations": len(d),
        "distinct_spectra": int(d.specObjID.nunique()), "recomputed_spectral_ratio_max_error": ratio_error,
        "physical_host_fold_leakage": False}


def transport():
    code, out = CODE / "host_transport", RESULTS / "host_transport"
    validation = read(out / "validation.json")
    check(validation["passed"] and validation["physical_operator_review_included"], "transport full validator passed")
    identities(validation["code_sha256"])
    identities(validation["record_sha256"], out)
    identity(out / "inputs.json", validation["source_registry_sha256"])
    pinned = read(out / "inputs.json")["files"]
    for item in pinned:
        identity(ROOT / item["path"], item["sha256"])
        check((ROOT / item["path"]).stat().st_size == item["bytes"], "source byte count", item["path"])
    # Follow identities rather than merely trusting the branch's passed field.
    for name in validation["record_sha256"]:
        record = read(out / name)
        for key in ("input_sha256", "code_sha256", "output_sha256"):
            if isinstance(record.get(key), dict):
                identities(record[key])
    roman = read(out / "roman-interface.json")
    identity(ROOT / roman["output_path"], roman["output_sha256"])
    identity(ROOT / roman["properties_path"], roman["properties_sha256"])
    filters = read(out / "filter-interface.json")
    identity(code / "filters.py", filters["code_sha256"])
    identity(ROOT / filters["npz_path"], filters["npz_sha256"])
    for band in filters["filters"]:
        for kind in ("input", "output"):
            identity(ROOT / band[kind + "_path"], band[kind + "_sha256"])
    review = read(out / "physical-operator-review.json")
    identity(code / "review_physical_operator.py", review["review_code_sha256"])
    identity(CODE / "physical_ages/model.py", review["model_sha256"])
    worst = max(row["max_abs_colour_error_mag"] for row in review["cases"])
    same_numbers(worst, review["max_abs_colour_error_mag"], "dense review maximum")
    check(review["passed"] and worst < 1e-6, "independent accurate photometry quadrature")
    infrared = read(out / "infrared-validation.json")
    check(infrared["passed"] and infrared["minimum_covariance_eigenvalue_Jy2"] >= -1e-15,
          "infrared covariance proof")
    SECTIONS["host_transport"] = {"pinned_sources": len(pinned),
        "signed_negative_local_fluxes_retained": validation["roman_negative_measurements_retained"],
        "dense_quadrature_cases": len(review["cases"]), "dense_quadrature_max_colour_error_mag": worst,
        "infrared_covariances_recomputed_by_branch_validator": infrared["covariances_recomputed"],
        "infrared_physical_status": "Confused beam measurements with empirical sky covariance; no deblended host dust luminosity or energy-balance constraint."}


def selected_inputs(sample):
    """Independent membership reconstruction, keeping signed measurements."""
    if sample == "sdss":
        table = pd.read_csv(ROOT / ".work/galaxy-validation/host-spectrum-age-ledger.csv", dtype={"specObjID": str})
        table = table[table.specObjID.notna()].drop_duplicates("specObjID")
        flux = table[[f"fiberFlux_{band}" for band in "ugriz"]].to_numpy(float)
        ivar = table[[f"fiberFluxIvar_{band}" for band in "ugriz"]].to_numpy(float)
        valid = np.isfinite(flux).all(axis=1) & np.isfinite(ivar).all(axis=1) & (ivar > 0).all(axis=1)
        selected = table.loc[valid]
        return dict(zip(selected.specObjID, selected.spectral_z)), int((~valid).sum())
    path = ROOT / (".work/host-transport/roman-photometry.csv" if sample == "roman" else ".work/host-transport/des-deep-photometry.csv")
    table = pd.read_csv(path, dtype={"sn_id": str, "host_id": str})
    if sample == "roman":
        table = table[(table.aperture == "local") & (table.filter_family == "SDSS")]
        identifier, bands, flux_col, error_col = "sn_id", list("ugriz"), "flux_ab_nanomaggies", "flux_error_ab_nanomaggies"
    else:
        table = table[table.selected_Dovekie & table.deep_quality & table.host_dlr_lt4]
        identifier, bands, flux_col, error_col = "host_id", ["u", "g", "r", "i", "z", "j", "h", "ks"], "flux_dered_ab_nanomaggies", "flux_error_dered_ab_nanomaggies"
    selected, excluded = {}, 0
    for name, group in table.groupby(identifier):
        for _, repeated in group.groupby("band"):
            check(repeated[[flux_col, error_col, "redshift"]].drop_duplicates().shape[0] == 1,
                  "duplicate host measurements agree", {"sample": sample, "id": name})
        group = group.drop_duplicates("band").set_index("band").reindex(bands)
        values = group[[flux_col, error_col]].to_numpy(float)
        valid = np.isfinite(values).all() and (values[:, 1] > 0).all()
        if valid:
            selected[str(name)] = float(group.redshift.iloc[0])
        else:
            excluded += 1
    return selected, excluded


def interval_checks(table, label):
    bounds = ["age_low_Gyr", "age_high_Gyr", "joint_age_low_Gyr", "joint_age_high_Gyr", "photometry_Dn_low", "photometry_Dn_high"]
    failed = table.numerical_failures > 0
    check(table.loc[failed, bounds].isna().all().all(), label + ": numerical failures fail closed")
    check(not failed.any(), label + ": no unresolved numerical failures", table.loc[failed, "id"].astype(str).tolist())
    for lower, upper in zip(bounds[::2], bounds[1::2]):
        check((table[lower].isna() == table[upper].isna()).all(), label + ": paired bounds " + lower)
        q = table[table[lower].notna()]
        check((q[lower] <= q[upper] + 2e-5).all(), label + ": ordered interval " + lower)
        if "Gyr" in lower:
            check((q[lower] >= .005 - 2e-5).all() and (q[upper] <= q.cosmic_age_ceiling_Gyr + 2e-5).all(),
                  label + ": physical age support " + lower)
    check(np.isfinite(table.minimum_chi2).all() and (table.minimum_chi2 >= 0).all(), label + ": finite nonnegative residuals")
    check(((table.nuisance_cells_feasible > 0) == table.age_low_Gyr.notna())[~failed].all(), label + ": photometric feasible-cell closure")
    check(((table.joint_nuisance_cells_feasible > 0) == table.joint_age_low_Gyr.notna())[~failed].all(), label + ": joint feasible-cell closure")


def physical():
    code, out = CODE / "physical_ages", RESULTS / "physical_ages"
    acquisition = read(out / "acquisition.json")
    identity(code / "acquire.py", acquisition["code_sha256"])
    identity(ROOT / acquisition["stellar_inventory_path"], acquisition["stellar_inventory_sha256"])
    stellar = read(ROOT / acquisition["stellar_inventory_path"])
    identities(stellar, ROOT / ".work/physical-ages/fsps")
    check(len(stellar) == acquisition["stellar_files"], "stellar inventory count")
    for item in acquisition["files"]:
        identity(ROOT / item["path"], item["sha256"])
    for name in ("stellar-library.json", "stellar-library-refined.json"):
        record = read(out / name)
        identity(ROOT / record["output"], record["output_sha256"])
        identity(code / "build_library.py", record["code_sha256"])
        identity(code / "design.json", record["design_sha256"])
        identity(out / "acquisition.json", record["acquisition_sha256"])
    summaries = {}
    for sample in ("sdss", "roman", "des"):
        record = read(out / (sample + "-age-bounds.json"))
        identities(record["source_sha256"])
        for required in (".work/physical-ages/ssp-library.npz", ".work/physical-ages/sfd/SFD_dust_4096_ngp.fits", ".work/physical-ages/sfd/SFD_dust_4096_sgp.fits"):
            check(required in record["source_sha256"], "physical source dependency recorded", {"sample": sample, "source": required})
        identities(record["code_sha256"], code)
        identity(code / "design.json", record["design_sha256"])
        identity(ROOT / record["output"], record["output_sha256"])
        d = pd.read_csv(ROOT / record["output"], dtype={"id": str})
        selected, excluded = selected_inputs(sample)
        check(set(d.id) == set(selected), sample + ": exact source-derived selected IDs")
        check(len(d) == 2 * len(selected) == 2 * record["objects"], sample + ": two scenarios per object")
        check(not d.duplicated(["id", "model_discrepancy_mag"]).any(), sample + ": unique object/scenario")
        check(set(d.model_discrepancy_mag) == {0., .03}, sample + ": declared discrepancy scenarios")
        check(sum(record["excluded"].values()) == excluded, sample + ": exclusion count")
        check(np.allclose(d.z, d.id.map(selected), atol=1e-14, rtol=0), sample + ": exact selected redshifts")
        interval_checks(d, sample)
        # Elementary flat LCDM cosmic time, independent of Astropy's routine.
        hubble_gyr = 3.0856775814913673e19 / 70 / (365.25 * 86400 * 1e9)
        ceiling = 2 * hubble_gyr / (3 * np.sqrt(.7)) * np.arcsinh(np.sqrt(.7 / .3) / (1 + d.z) ** 1.5)
        check(np.allclose(d.cosmic_age_ceiling_Gyr, ceiling, atol=2e-9, rtol=0), sample + ": independent cosmic clock")
        summaries[sample] = []
        for summary in record["summaries"]:
            q = d[d.model_discrepancy_mag == summary["model_discrepancy_mag"]]
            good, joint = q.age_low_Gyr.notna(), q.joint_age_low_Gyr.notna()
            tested = q.heldout_Dn_intervals_disjoint.notna()
            actual = {"objects": len(q), "photometry_feasible": int(good.sum()), "joint_feasible": int(joint.sum()),
                "heldout_Dn_tested": int(tested.sum()), "heldout_Dn_disjoint": int(q.loc[tested, "heldout_Dn_intervals_disjoint"].astype(bool).sum()),
                "numerical_failures": int(q.numerical_failures.sum()),
                "age_interval_width_Gyr": quantiles(q.age_high_Gyr - q.age_low_Gyr), "lower_age_Gyr": quantiles(q.age_low_Gyr),
                "upper_age_Gyr": quantiles(q.age_high_Gyr), "joint_interval_width_Gyr": quantiles(q.joint_age_high_Gyr - q.joint_age_low_Gyr),
                "minimum_chi2": quantiles(q.minimum_chi2)}
            for key, value in actual.items():
                same_numbers(value, summary[key], sample + ": summary " + key)
            summaries[sample].append({"model_discrepancy_mag": summary["model_discrepancy_mag"], **actual})
    SECTIONS["physical_age_bounds"] = {"source_stellar_files": len(stellar), "independent_membership_and_summaries": summaries,
        "inference_status": "Conditional compatibility intervals on finite SSP/dust/metallicity support and an assumed cosmological clock. No empirical progenitor-age posterior or distance correction.",
        "confidence_scope": "Only the fixed-error, included-model Gaussian construction has stated coverage. The observed-flux-dependent 0.03 mag floor is a sensitivity, not calibrated coverage."}


def physical_numerics():
    code, out, work = CODE / "physical_ages", RESULTS / "physical_ages", ROOT / ".work/physical-ages"
    record = read(out / "numerical-validation.json")
    identities(record["code_sha256"], code)
    check(record["passed"], "physical numerical validator passed")
    dense_error = max(row["max_colour_difference_mag"] for row in record["final_gauss_vs_dense_checks"])
    check(dense_error < 1e-6, "final quadrature versus dense physical integration")
    check(record["vectorized_operator_max_relative_difference"] < 1e-12, "operator linearity")
    check(record["conic_coverage"]["solver_failures_on_feasible_inputs"] == 0, "injected feasible cone failures")
    primal = record["independent_primal"]
    check(primal["verified_endpoints"] == primal["endpoints"] == len(primal["records"]), "all independent primal endpoints verified")
    differences = []
    for endpoint in primal["records"]:
        best = endpoint["independent_best"]
        check(endpoint["passed"] and best["success"] and best["feasible"], "independent feasible endpoint")
        check(best["chi2"] <= 11.070497693516351 + primal["criteria"]["maximum_chi2_excess"]
              and best["simplex_error"] <= primal["criteria"]["maximum_simplex_error"]
              and best["minimum_coordinate"] >= primal["criteria"]["minimum_coordinate"],
              "independent endpoint original-unit acceptance")
        check(not any(attempt["uses_conic_witness"] for attempt in endpoint["attempts"]), "independent starts exclude conic witness")
        differences.append(abs(best["age_Gyr"] - endpoint["conic_outward_bound_Gyr"]))
    same_numbers(max(differences), primal["maximum_best_difference_Gyr"], "independent primal maximum")
    check(max(differences) < 1e-5, "independent primal endpoint accuracy")
    for name, script in (("sensitivity.json", "sensitivity.py"), ("physical-injections.json", "injections.py")):
        r = read(out / name)
        identity(code / script, r["code_sha256"])
        identities(r["dependencies_sha256"])
        identity(ROOT / r["output"], r["output_sha256"])
        d = pd.read_csv(ROOT / r["output"], dtype={"id": str})
        interval_checks(d, name)
        grouping = ["sample", "variant" if name == "sensitivity.json" else "generator", "model_discrepancy_mag"]
        for summary in r["summaries"]:
            q = d[np.logical_and.reduce([d[key] == summary[key] for key in grouping])]
            count_key = "n" if name == "sensitivity.json" else "draws"
            feasible_key = "feasible" if name == "sensitivity.json" else "photometry_feasible"
            check(len(q) == summary[count_key], name + ": scenario count")
            check(int(q.age_low_Gyr.notna().sum()) == summary[feasible_key], name + ": feasible count")
            check(int(q.numerical_failures.sum()) == summary["numerical_failures"] == 0, name + ": solver failure accounting")
    benchmark = read(out / "compute-benchmark.json")
    identity(code / "gpu_benchmark.py", benchmark["code_sha256"])
    identities(benchmark["dependencies_sha256"])
    identity(work / "gpu-forward-fixture.npz", benchmark["fixture_sha256"])
    identity(work / "des-gpu-forward-flux.npy", benchmark["GPU_result_sha256"])
    fixture = np.load(work / "gpu-forward-fixture.npz")
    cpu = fixture["spectra"] @ fixture["kernels"]
    gpu = np.load(work / "des-gpu-forward-flux.npy")
    gpu_error = float(np.max(abs(cpu - gpu) / np.maximum(abs(cpu), 1e-300)))
    check(gpu_error < 1e-12, "independent saved GPU/CPU matrix product")
    selected, _ = selected_inputs("des")
    check(set(fixture["ids"].astype(str)) == set(selected), "GPU actual DES source IDs")
    retry = read(out / "endpoint-retry-regression.json")
    identity(code / "model.py", retry["current_model_sha256"])
    identity(ROOT / retry["fixture_path"], retry["fixture_sha256"])
    before, after = retry["before"]["rows"], retry["after"]["rows"]
    check(sum(r["numerical_failures"] for r in before) > 0 and sum(r["numerical_failures"] for r in after) == 0,
          "captured endpoint regression resolved")
    direct = retry["direct_validation"]
    check(direct["chi2"] <= direct["threshold"] + 2e-5 and direct["minimum_weight"] >= -1e-7
          and abs(direct["sum_weights"] - 1) < 1e-6 and direct["inequality_slack"] < 1e-7,
          "endpoint retry original-unit feasibility")
    check(retry["successful_cone"]["duality_gap"] < 1e-5, "endpoint retry primal-dual bound")
    audits = []
    for path in sorted(out.glob("solver-audit-*.json")):
        r = read(path)
        identity(code / "solver_audit.py", r["code_sha256"])
        identity(work / "nnls-regression-fixture.npz", r["fixture_sha256"])
        check(r["bounded_variable_success"] and r["bounded_variable_KKT"] < 1e-6, "BVLS regression certificate", path.name)
        audits.append({k: r[k] for k in ("scipy", "numpy", "nnls_residual_identity_passed", "bounded_variable_squared_residual")})
    SECTIONS["physical_numerics"] = {"independent_primal_attempts": primal["attempts"], "independently_verified_endpoints": primal["verified_endpoints"],
        "independent_primal_failed_attempts_retained": primal["attempts"] - primal["converged"],
        "final_gauss_vs_dense_max_colour_error_mag": dense_error,
        "saved_GPU_recomputed_max_relative_error": gpu_error, "endpoint_regression_resolved": True,
        "historical_NNLS_fixture": audits,
        "scope": "Historical failing NNLS fixture is retained; final physical fits use a different bounded solver and corrected quadrature. Failed historical checks are not relabelled as successful."}


def survey():
    code, out, work = CODE / "survey_physics", RESULTS / "survey_physics", ROOT / ".work/survey-physics"
    validation = read(out / "validation.json")
    identity(code / "validate.py", validation["code_sha256"])
    for key in ("source_code_current_sha256", "binary_sha256"):
        identities(validation[key])
    identities(validation["model_and_cadence_sha256"], work)
    review = validation["independent_physics_review"]
    identity(ROOT / review["path"], review["sha256"])
    r = read(ROOT / review["path"])
    identity(CODE / "galaxy_validation/review_survey_physics.py", r["code_sha256"])
    for name in ("main_W22_campaign", "fixed_RV31_campaign"):
        identities(r[name]["source_sha256"])
        check(r[name]["maximum_noiseless_fluxratio_relative_error"] < 1e-5, "native signed photon injection", name)
    jobs = [job for name in ("campaign.json", "positive-dust-campaign.json") for job in read(out / name)["jobs"]]
    native, classifiers, history = read(out / "native-execution.json"), read(out / "classifier-execution.json"), []
    check(len(jobs) == len(native) == len(classifiers) == validation["native_run_count"], "all native runs accounted for")
    for job, saved_run, saved_classifier in zip(jobs, native, classifiers):
        folder = work / "fits" / job["name"]
        run, classifier = read(folder / "run.json"), read(folder / "classifier.json")
        check(run == saved_run and classifier == saved_classifier, "native execution ledger equality", job["name"])
        check(run["fit_graceful"] and run["generation_returncode"] == run["fit_returncode"] == 0, "native successful execution", job["name"])
        identity(job["input"], job["input_sha256"])
        check(run["hashes"][job["input"]] == job["input_sha256"], "native input identity", job["name"])
        for filename in ("fit.FITRES.TEXT", "fit.SNANA.TEXT"):
            identity(folder / filename, run[filename]["sha256"])
        identity(folder / "fit.FITRES.TEXT", classifier["native_fit_sha256"])
        check(classifier["truth_inputs_used"] is False and classifier["single_batch_max_probability_difference"] < 1e-6,
              "classifier no truth leakage and batch identity", job["name"])
        # Preserve generation-time records; later driver review edits are not
        # grounds to relabel old generation hashes as current source hashes.
        history.append({"job": job["name"], "run_record_sha256": sha(folder / "run.json")})
    min_eigenvalues = {}
    for filename in ("correction-response.json", "positive-dust-response.json"):
        r = read(out / filename)
        identity(code / "analyze.py", r["code_sha256"])
        identity(code / "analysis-design.json", r["analysis_design_sha256"])
        for count in r["counts"].values():
            identities(count["files_sha256"])
        for contrast, row in r["variants"]["primary"]["differences"].items():
            c = np.array(row["zbin_covariance"], float)
            active = np.isfinite(np.diag(c))
            c = c[np.ix_(active, active)]
            if len(c):
                check(np.isfinite(c).all() and np.allclose(c, c.T, atol=1e-12), "supported covariance finite and symmetric", contrast)
                minimum = float(np.linalg.eigvalsh(c).min())
                check(minimum >= -1e-12, "supported covariance PSD", contrast)
                min_eigenvalues[filename + ":" + contrast] = minimum
    failed_bbc = read(out / "native-bbc.json")
    bbc_gates = {key: {"returncode": value["returncode"], "graceful": value.get("graceful")} for key, value in failed_bbc.items()}
    for key, value in bbc_gates.items():
        check(value["graceful"] == (value["returncode"] == 0), "native BBC status truthfulness", key)
    SECTIONS["survey_physics"] = {"native_runs": len(jobs), "generated_attempts": sum(job["attempts"] for job in jobs),
        "generation_history": history, "supported_covariance_min_eigenvalues": min_eigenvalues, "BBC_execution_gates": bbc_gates,
        "claim_boundary": "Mock population, pure Ia, fixed light-curve surface and classifier. Production 4D BBC failures remain explicit; diagnostic 1D/nearest-neighbor closure does not measure an observational age correction."}


def code_and_documentation():
    parsed, documents, links = [], [], []
    for branch in BRANCHES:
        for path in sorted((CODE / branch).glob("*.py")):
            ast.parse(path.read_text(), filename=str(path))
            parsed.append(relative(path))
        documents += list((CODE / branch).glob("*.md"))
    documents += [BASE / "notes" / (branch.replace("_", "-") + "-results.md") for branch in BRANCHES]
    for path in documents:
        check(path.exists(), "branch documentation exists", relative(path))
        if not path.exists():
            continue
        sha(path)
        for target in re.findall(r"\]\(([^\s)]+)(?:\s+\"[^\"]*\")?\)", path.read_text()):
            if target.startswith(("https:", "http:", "mailto:", "#", "data:")):
                continue
            target = unquote(target.strip("<>").split("#", 1)[0])
            dest = path.parent / target
            check(dest.exists(), "local Markdown target", {"document": relative(path), "target": target})
            links.append({"document": relative(path), "target": target})
    SECTIONS["authored_code_and_navigation"] = {"python_AST_files": len(parsed), "markdown_files": len(documents), "local_links_checked": len(links)}


def main():
    for function in (galaxy, transport, physical, physical_numerics, survey, code_and_documentation):
        try:
            function()
        except Exception as exc:
            FAILURES.append({"check": function.__name__, "exception": type(exc).__name__, "detail": str(exc)})
    result_hashes = {relative(p): sha(p) for branch in BRANCHES for p in sorted((RESULTS / branch).glob("*.json"))}
    code_hashes = {relative(p): sha(p) for branch in BRANCHES for p in sorted((CODE / branch).iterdir()) if p.is_file()}
    checker_hash = sha(Path(__file__))
    for path, (_, size, mtime) in list(CACHE.items()):
        stat = path.stat()
        check((size, mtime) == (stat.st_size, stat.st_mtime_ns), "snapshot unchanged during audit", relative(path))
    report = {"completed_utc": datetime.now(timezone.utc).isoformat(), "passed": not FAILURES,
        "status": "passed" if not FAILURES else "failed",
        "checker_sha256": checker_hash, "scope": "Cross-branch byte identities and independent arithmetic/selection checks; not validation of the physical model or an empirical cosmology correction.",
        "unique_files_hashed": len(CACHE), "bytes_hashed": sum(value[1] for value in CACHE.values()),
        "branch_result_sha256": result_hashes,
        "authored_code_and_design_sha256": code_hashes,
        "sections": SECTIONS, "failures": FAILURES,
        "legacy_baseline": "391-input/29-run verification remains in validation/record.py; not replaced by this extension."}
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"passed": report["passed"], "unique_files_hashed": report["unique_files_hashed"], "failures": FAILURES}, indent=2))
    raise SystemExit(0 if report["passed"] else 1)


if __name__ == "__main__":
    main()
