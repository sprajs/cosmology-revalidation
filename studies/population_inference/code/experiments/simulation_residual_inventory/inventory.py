"""Read-only metadata inventory of the nine existing forward simulation arms."""
from __future__ import annotations

import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
SIM = ROOT / "phase2/literature/simulations"
FWD = ROOT / "phase2/hierarchy/forward"
FIT = ROOT / "phase2/official/results"
NML = ROOT / "phase2/checkpoint/official-inputs"
ARMS = ("P21", "BS21", "G10", "P21_dmplus010", "P21_dmminus010",
        "P21_dmz020", "P21_rho000", "P21_rho090", "P21_noisetrue120")


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for part in iter(lambda: f.read(1024 * 1024), b""):
            h.update(part)
    return h.hexdigest()


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def counts(values) -> dict[str, int]:
    return {str(k).strip(): int(v) for k, v in zip(*np.unique(values, return_counts=True))}


def main() -> None:
    records = []
    inputs = [Path(__file__), FWD / "summary.json", FWD / "manifest.json"]
    for arm in ARMS:
        version = "PH2_pilot02_" + arm
        base = SIM / "outputs" / version
        head_path = base / f"{version}_HEAD.FITS"
        phot_path = base / f"{version}_PHOT.FITS"
        dump_path = base / f"{version}.DUMP"
        input_path = SIM / "inputs" / f"{version}.input"
        run_manifest = SIM / "manifests" / f"{version}.json"
        nml_path = NML / f"snana_forward_{arm.lower()}.nml"
        actual_nml_path = ROOT / "phase2/official/inputs" / f"snana_forward_{arm.lower()}.nml"
        fit_path = FIT / f"snana_forward_{arm.lower()}.FITRES.TEXT"
        fit_log_path = ROOT / "phase2/official/diagnostics" / f"snana-forward-{arm.lower().replace('_', '-')}.log"
        attempts_path = FWD / f"{arm}-attempts.csv.gz"
        fitted_path = FWD / f"{arm}-fitted.csv.gz"
        cov_path = FWD / f"{arm}-covariance.npz"
        paths = [head_path, phot_path, dump_path, input_path, run_manifest,
                 nml_path, actual_nml_path, fit_path, fit_log_path,
                 attempts_path, fitted_path, cov_path]
        assert all(p.exists() for p in paths), [rel(p) for p in paths if not p.exists()]
        inputs.extend(paths)
        assert sha(nml_path) == sha(actual_nml_path)
        run = json.loads(run_manifest.read_text())
        assert run["exit_code"] == 0 and run["success_marker"]
        for item in run["outputs"]:
            assert sha(ROOT / item["path"]) == item["sha256"]
        nml = nml_path.read_text()
        inp = input_path.read_text()
        fit_log_head = "\n".join(fit_log_path.read_text().splitlines()[:7])
        assert "SNANA-current/bin/snlc_fit.exe" in fit_log_head
        assert "SNANA_VERSION: 886408a4-output17g" in fit_log_head
        assert "SNTABLE_LIST      = 'FITRES(text:host)'" in nml
        assert "GENFILTERS:  griz" in inp
        with fits.open(head_path, memmap=True) as h, fits.open(phot_path, memmap=True) as p:
            heads = h[1].data
            phot = p[1].data
            hcols = list(h[1].columns.names)
            pcols = list(p[1].columns.names)
            cids = np.asarray(heads["SNID"]).astype(int)
            assert len(cids) == len(set(cids))
            first = np.asarray(heads["PTROBS_MIN"], dtype=int)
            last = np.asarray(heads["PTROBS_MAX"], dtype=int)
            nobs = np.asarray(heads["NOBS"], dtype=int)
            assert first[0] == 1 and np.array_equal(last - first + 1, nobs)
            assert np.array_equal(first[1:], last[:-1] + 2)
            assert last[-1] == len(phot) - 1
            sep_idx = last
            assert np.all(np.char.strip(np.asarray(phot["BAND"][sep_idx]).astype(str)) == "-")
            is_observation = np.ones(len(phot), dtype=bool)
            is_observation[sep_idx] = False
            assert int(is_observation.sum()) == int(nobs.sum())
            raw_fields = np.char.strip(np.asarray(phot["FIELD"][first - 1]).astype(str))
            head_field = dict(zip(cids, raw_fields))
            field_mult = sum(len(set(np.char.strip(np.asarray(phot["FIELD"][a-1:b]).astype(str)))) > 1
                             for a, b in zip(first, last))
            raw_flags = counts(phot["PHOTFLAG"][is_observation])
            raw_bands = counts(phot["BAND"][is_observation])
            raw_fields_count = counts(phot["FIELD"][is_observation])
            raw_bad_error = int(np.count_nonzero(phot["FLUXCALERR"][is_observation] <= 0))
            head_fields = list(head_field.values())
            head_types = counts(heads["SNTYPE"])
            head_sim_models = counts(heads["SIM_MODEL_NAME"])
            nhead = len(heads)
            nphot = len(phot)
            nlinked = int(nobs.sum())
        attempts = pd.read_csv(attempts_path)
        fitted = pd.read_csv(fitted_path)
        assert len(attempts) == 26518
        assert attempts.generated_attempt_index.is_unique
        assert np.array_equal(attempts.generated_attempt_index.to_numpy(), np.arange(1, 26519))
        written = attempts.loc[attempts.written]
        assert len(written) == nhead and set(written.CID.astype(int)) == set(cids)
        assert not written.CID.duplicated().any()
        assert len(fitted) == attempts.fit_pass.sum()
        assert set(fitted.CID.astype(int)) <= set(cids)
        assert int(fitted.basic_quality_pass.sum()) == int(attempts.basic_quality_pass.sum())
        fit_first_field_disagree = int(sum(str(row.FIELD).strip() != head_field[int(row.CID)]
                                           for row in fitted.itertuples(index=False)))
        with np.load(cov_path, allow_pickle=False) as cv:
            cov_keys = {key: {"shape": list(cv[key].shape), "dtype": str(cv[key].dtype)}
                        for key in cv.files}
            assert cv["mag_parameters"].shape == (len(fitted), 4)
            assert cv["mag_covariance"].shape == (len(fitted), 4, 4)
            assert np.array_equal(cv["CID"], fitted.CID.to_numpy(dtype=int))
        fit_header = next(line for line in fit_path.open() if line.startswith("VARNAMES:"))
        assert "FIELD" in fit_header and "SIM_LIBID" in fit_header
        fit_row_count = sum(line.startswith("SN:") for line in fit_path.open())
        assert fit_row_count == len(fitted)
        records.append({
            "arm": arm, "version": version,
            "paths": {"simulation_input": rel(input_path), "run_manifest": rel(run_manifest),
                      "head": rel(head_path), "phot": rel(phot_path), "dump": rel(dump_path),
                      "fit_input": rel(nml_path), "fitres": rel(fit_path),
                      "executed_fit_input": rel(actual_nml_path),
                      "fit_log": rel(fit_log_path),
                      "attempts": rel(attempts_path), "fitted": rel(fitted_path),
                      "parameter_covariance": rel(cov_path)},
            "counts": {"generated_attempts": len(attempts),
                       "generated_duplicate_cids": int(attempts.CID.duplicated().sum()),
                       "written_head_objects": nhead, "head_observations": nlinked,
                       "phot_rows_including_separators": nphot,
                       "separator_rows": nphot - nlinked,
                       "fitted_objects": len(fitted),
                       "basic_quality_objects": int(fitted.basic_quality_pass.sum()),
                       "raw_fluxerr_nonpositive": raw_bad_error,
                       "written_objects_with_multiple_phot_fields": int(field_mult),
                       "fitted_field_vs_first_phot_field_disagreements": fit_first_field_disagree},
            "schema": {"head_columns": hcols, "phot_columns": pcols,
                       "fitres_columns": fit_header.split()[1:],
                       "parameter_covariance": cov_keys},
            "raw_photflag_observation_counts": raw_flags,
            "raw_band_observation_counts": raw_bands,
            "raw_field_observation_counts": raw_fields_count,
            "first_phot_field_object_counts": dict(Counter(head_fields)),
            "head_sntype_counts": head_types,
            "head_sim_model_counts": head_sim_models,
            "simulation_recipe": {key: [line.strip() for line in inp.splitlines()
                                       if line.lstrip().startswith(key)]
                                  for key in ("GENMODEL:", "KCOR_FILE:", "FLUXERRMODEL_FILE:",
                                              "SMEARFLAG_FLUX:", "SMEARFLAG_ZEROPT:",
                                              "APPLY_SEARCHEFF_OPT:", "GENMAG_SMEAR_ADDPHASECOR:")},
            "input_hashes": {"head_sha256": sha(head_path), "phot_sha256": sha(phot_path),
                             "fitres_sha256": sha(fit_path), "nml_sha256": sha(nml_path),
                             "fit_log_sha256": sha(fit_log_path)},
        })
    result = {
        "schema": "simulation_residual_feasibility_inventory_v1",
        "scope": "Metadata/schema/counts only; no new simulation or residual score",
        "arms": records,
        "epoch_covariance_export_present": False,
        "per_epoch_fitter_accepted_mask_export_present": False,
        "simulated_classifier_probabilities_present": False,
        "provenance": "Existing run manifests and forward hierarchy summary; all five original sim outputs per arm checked against manifest hashes",
    }
    result_path = OUT / "inventory.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    # An unexecuted, single-CID native-export candidate. It deliberately keeps
    # the original P21 epoch cuts and uses a measured-fit seed, not truth.
    with np.load(FWD / "P21-covariance.npz", allow_pickle=False) as cv:
        idx = np.flatnonzero(cv["CID"] == 16)
        assert len(idx) == 1
        m_b, x_1, color, t_0 = cv["mag_parameters"][idx[0]]
    x_0 = 10 ** ((10.635 - m_b) / 2.5)
    seed_path = OUT / "candidate_p21_cid16.seed.FITRES"
    seed_path.write_text("VARNAMES: CID PKMJD x0 x1 c\n"
                         f"SN: 16 {t_0:.17g} {x_0:.17g} {x_1:.17g} {color:.17g}\n")
    original = (NML / "snana_forward_p21.nml").read_text()
    candidate = original.replace("&SNLCINP", "&SNLCINP\n"
                                 "    OPT_SNCID_LIST = 3\n"
                                 f"    SNCID_LIST_FILE = '{seed_path.resolve()}'", 1)
    candidate = candidate.replace("SNTABLE_LIST      = 'FITRES(text:host)'",
                                  "SNTABLE_LIST      = 'FITRES(text:host) LCPLOT(text:col)'", 1)
    candidate = re.sub(r"(?m)^    TEXTFILE_PREFIX\s*=.*$",
                       f"    TEXTFILE_PREFIX = '{(OUT / 'candidate_p21_cid16').resolve()}'",
                       candidate, count=1)
    candidate = candidate.replace("&FITINP", "&FITINP\n    LFIXPAR_ALL = T", 1)
    assert candidate.count("LFIXPAR_ALL = T") == 1
    candidate_path = OUT / "candidate_p21_cid16.nml"
    candidate_path.write_text(candidate)
    manifest = {"script_sha256": sha(Path(__file__)),
                "inputs_sha256": {rel(p): sha(p) for p in sorted(set(inputs))},
                "inventory_sha256": sha(result_path),
                "candidate_sha256": {rel(p): sha(p) for p in (seed_path, candidate_path)},
                "candidate_executed": False}
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({r["arm"]: r["counts"] for r in records}, indent=2))


if __name__ == "__main__":
    main()
