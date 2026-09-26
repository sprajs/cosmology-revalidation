"""Execute frozen P21/G10 native-fit provenance gate without scoring residuals.

The simulation design and law choices are frozen in Astra's hashed handoff.
This runner only prepares measured-initial-value native fits, runs an isolated
audit executable, and checks output/objective identity. It does not simulate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import linear_sum_assignment
from scipy.linalg import solve_triangular

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "runs/research_2026_09_26/astra_design/simulation_design"
OUT = ROOT / "runs/research_2026_09_26/simulation_residual_control"
PROTO = ROOT / "runs/research_2026_09_26/astra_design/simulation-residual-protocol.md"
LAW = ROOT / "runs/research_2026_09_26/astra_design/simulation-law-clarification.md"
EXE = ROOT / "phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe"
EXPORT = ROOT / "scripts/research_2026_09_26/export_flux_objectives.py"
PYTHON = ROOT / "phase2/env-official/bin/python"
EXPECTED_PROTOCOL = "764e9b1394b2a6c5df708f745bcc6c713640b9d553e57aed98c716c0dd121efe"
EXPECTED_HANDOFF = "8a19e0660ecdd9f8ccfc0504c4c9adf5887c01434df6400d695ac805ad023415"
EXPECTED_LAW = "fba384256d687eb3043b28450ec7a31dc89de27ddab52918859886c3655f1288"
SCORE_PLAN = OUT / "score-plan.md"
EXPECTED_SCORE_PLAN = "0be5fb8a8b34482bb516e07dcf0ca2877799c00edb72e689d185b2ec6de6ded6"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def preflight() -> dict:
    assert sha(PROTO) == EXPECTED_PROTOCOL
    assert sha(DESIGN / "handoff-manifest.json") == EXPECTED_HANDOFF
    assert sha(LAW) == EXPECTED_LAW
    handoff = json.loads((DESIGN / "handoff-manifest.json").read_text())
    for name, expected in handoff["sha256"].items():
        assert sha(ROOT / name) == expected, name
    return handoff


def make_input(arm: str, law: str, cohort: str) -> Path:
    assert arm in ("P21", "G10") and law in ("approx_minus99", "exact_99")
    assert cohort in ("engineering8", "full256")
    preflight()
    path = OUT / cohort / arm / law
    if path.exists():
        raise RuntimeError(f"Preserve existing run directory: {path}")
    path.mkdir(parents=True)
    ids_file = DESIGN / f"{arm}-{'engineering8-' if cohort == 'engineering8' else ''}cids.txt"
    ids = [int(value) for value in ids_file.read_text().split()]
    assert len(ids) == (8 if cohort == "engineering8" else 256)
    assert len(ids) == len(set(ids))
    source = DESIGN / f"{arm}-cohort.csv"
    table = pd.read_csv(source)
    table = table.set_index("CID").loc[ids]
    assert np.all(np.isfinite(table[["t0_double", "mB_double", "x1_double", "c_double"]]))
    assert np.all(np.isfinite(table.x0) & (table.x0 > 0))
    x0 = 10 ** ((10.635 - table.mB_double.to_numpy(float)) / 2.5)
    seed = path / "measured_initial_values.FITRES"
    seed.write_text("VARNAMES: CID PKMJD x0 x1 c\n" + "".join(
        f"SN: {cid} {t0:.17g} {amp:.17g} {x1:.17g} {color:.17g}\n"
        for cid, t0, amp, x1, color in zip(ids, table.t0_double, x0,
                                            table.x1_double, table.c_double)))
    base = ROOT / f"phase2/checkpoint/official-inputs/snana_forward_{arm.lower()}.nml"
    actual = ROOT / f"phase2/official/inputs/snana_forward_{arm.lower()}.nml"
    assert sha(base) == sha(actual)
    original = base.read_text()
    assert "SNTABLE_LIST      = 'FITRES(text:host)'" in original
    assert "OPT_MWCOLORLAW = 99" in original
    assert "OPT_SNCID_LIST" not in original and "SNCID_LIST_FILE" not in original
    assert "LFIXPAR_ALL" not in original
    choice = "-99" if law == "approx_minus99" else "99"
    text = original.replace("&SNLCINP", "&SNLCINP\n"
                            "    MXLC_PLOT = 1000\n"
                            "    OPT_SNCID_LIST = 2\n"
                            f"    SNCID_LIST_FILE = '{seed.resolve()}'", 1)
    text = text.replace("SNTABLE_LIST      = 'FITRES(text:host)'",
                        "SNTABLE_LIST      = 'FITRES(text:host) LCPLOT(text:col)'", 1)
    text = re.sub(r"(?m)^    TEXTFILE_PREFIX\s*=.*$",
                  f"    TEXTFILE_PREFIX = '{(path / 'fit').resolve()}'", text, count=1)
    text = text.replace("OPT_MWCOLORLAW = 99", f"OPT_MWCOLORLAW = {choice}", 1)
    nml = path / "fit.nml"
    nml.write_text(text)
    inputs = {str(p.relative_to(ROOT)): sha(p) for p in
              (PROTO, LAW, DESIGN / "handoff-manifest.json", ids_file,
               source, base, actual, EXE, EXPORT, Path(__file__))}
    info = {"arm": arm, "law": law, "law_option": int(choice), "cohort": cohort,
            "selected_ids": ids, "selected_count": len(ids),
            "opt_sncid_list": 2,
            "opt_sncid_source_reason": "SNANA bit1 reads measured initial values while bit0 remains off, preserving ordinary cuts",
            "fitter_parameters_free": True, "truth_coordinates_used": False,
            "inputs_sha256": inputs, "seed_sha256": sha(seed), "nml_sha256": sha(nml),
            "execution_started": False}
    (path / "prepared.json").write_text(json.dumps(info, indent=2) + "\n")
    return path


def run_fit(path: Path) -> dict:
    info = json.loads((path / "prepared.json").read_text())
    nml = path / "fit.nml"
    seed = path / "measured_initial_values.FITRES"
    assert sha(nml) == info["nml_sha256"] and sha(seed) == info["seed_sha256"]
    assert not (path / "fit.log").exists()
    env = os.environ.copy()
    env.update(SNANA_DIR=str(ROOT / "phase2/official/build/SNANA-audit-v3"),
               SNDATA_ROOT=str(ROOT / "phase2/official/inputs/SNDATA_ROOT"),
               LD_LIBRARY_PATH=str(ROOT / "phase2/official/build/sysroot/usr/lib"),
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (path / "fit.log").open("x") as log:
        result = subprocess.run([str(EXE), str(nml)], cwd=path, env=env,
                                stdout=log, stderr=subprocess.STDOUT, check=False)
    text = (path / "fit.log").read_text()
    fitres = path / "fit.FITRES.TEXT"
    ids = set(map(int, info["selected_ids"]))
    fitted = set()
    if fitres.exists():
        for line in fitres.read_text().splitlines():
            if line.startswith("SN:"):
                fitted.add(int(line.split()[1]))
    context = set()
    for line in text.splitlines():
        if line.lstrip().startswith("PHASE2_OBJECTIVE:"):
            context.add(int(line.split()[1]))
    assert fitted <= ids and context <= ids
    check = {"returncode": result.returncode,
             "graceful": "ENDING PROGRAM GRACEFULLY." in text,
             "fitres_count": len(fitted), "objective_count": len(context),
             "fitres_cids": sorted(fitted), "objective_cids": sorted(context),
             "missing_fitres_cids": sorted(ids - fitted),
             "missing_objective_cids": sorted(ids - context),
             "log_sha256": sha(path / "fit.log"),
             "fitres_sha256": sha(fitres) if fitres.exists() else None}
    (path / "fit-gate.json").write_text(json.dumps(check, indent=2) + "\n")
    return check


def export(path: Path) -> dict:
    fit = json.loads((path / "fit-gate.json").read_text())
    assert fit["returncode"] == 0 and fit["graceful"]
    assert not (path / "objectives").exists()
    result = subprocess.run([str(PYTHON), str(EXPORT), "--log", str(path / "fit.log"),
                             "--output", str(path / "objectives")], cwd=ROOT,
                            capture_output=True, text=True, check=False)
    (path / "export.stdout").write_text(result.stdout)
    (path / "export.stderr").write_text(result.stderr)
    report = {"returncode": result.returncode,
              "manifest_exists": (path / "objectives/manifest.json").exists(),
              "stdout": result.stdout[-1000:], "stderr": result.stderr[-2000:]}
    if report["manifest_exists"]:
        obj = json.loads((path / "objectives/manifest.json").read_text())
        report["objects"] = obj["count"]
        report["epochs"] = obj["epochs"]
        report["max_objective_error"] = max(abs(v["objective_residual"])
                                            for v in obj["checks"].values())
        assert set(obj["checks"]) == set(map(str, fit["objective_cids"]))
    (path / "export-gate.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def match_rows(left: np.ndarray, right: np.ndarray, both_flux: bool) -> tuple[float, float, float]:
    """One-to-one native/recorded epoch identity, preserving duplicates."""
    assert len(left) <= len(right)
    # Columns are MJD, band, flux, quoted error. Raw PHOT and output lists
    # can contain the same MJD in multiple bands and repeated exposures.
    cost = (np.abs(left[:, 0].astype(float)[:, None] - right[:, 0].astype(float)[None, :])
            + 1e6 * (left[:, 1][:, None] != right[:, 1][None, :]))
    if both_flux:
        cost += 1e-4 * np.abs(left[:, 2].astype(float)[:, None]
                                 - right[:, 2].astype(float)[None, :])
    rows, cols = linear_sum_assignment(cost)
    assert np.array_equal(rows, np.arange(len(left)))
    assert np.array_equal(left[:, 1], right[cols, 1])
    dt = float(np.max(np.abs(left[:, 0].astype(float) - right[cols, 0].astype(float))))
    df = float(np.max(np.abs(left[:, 2].astype(float) - right[cols, 2].astype(float))))
    de = float(np.max(np.abs(left[:, 3].astype(float) - right[cols, 3].astype(float))))
    return dt, df, de


def gate(path: Path) -> dict:
    """Before scores, establish CID, flux-unit, and accepted-epoch closure."""
    preflight()
    info = json.loads((path / "prepared.json").read_text())
    fit = json.loads((path / "fit-gate.json").read_text())
    exp = json.loads((path / "export-gate.json").read_text())
    assert fit["returncode"] == exp["returncode"] == 0
    assert fit["graceful"] and exp["manifest_exists"]
    assert exp["max_objective_error"] < 1e-7
    assert exp["objects"] == len(fit["objective_cids"])
    arm = info["arm"]
    cohort = pd.read_csv(DESIGN / f"{arm}-cohort.csv").set_index("CID")
    sys_path = ROOT / "scripts/phase2/official"
    import sys
    sys.path.insert(0, str(sys_path))
    from audit_fits import read_fit
    native = read_fit(path / "fit.FITRES.TEXT")
    native.CID = native.CID.astype(int)
    assert native.CID.is_unique
    native = native.set_index("CID")
    assert set(native.index) == set(map(int, fit["fitres_cids"]))
    version = f"PH2_pilot02_{arm}"
    base = ROOT / "phase2/literature/simulations/outputs" / version
    with fits.open(base / f"{version}_HEAD.FITS", memmap=True) as head_hdu, \
         fits.open(base / f"{version}_PHOT.FITS", memmap=True) as phot_hdu:
        heads = head_hdu[1].data
        phot = phot_hdu[1].data
        head_index = {int(row["SNID"]): row for row in heads}
        assert len(head_index) == len(heads)
        raw = {}
        for cid in fit["objective_cids"]:
            h = head_index[cid]
            start, stop = int(h["PTROBS_MIN"]), int(h["PTROBS_MAX"])
            q = phot[start - 1:stop]
            assert len(q) == int(h["NOBS"])
            raw[cid] = q.copy()
    lc = pd.read_csv(path / "fit.LCPLOT.TEXT", sep=r"\s+", comment="#", header=None)
    assert lc.shape[1] == 10
    scored_ids = sorted(set(fit["objective_cids"]) & set(fit["fitres_cids"]))
    rows = []
    for cid in scored_ids:
        obj = np.load(path / "objectives" / f"objective_{cid}.npz", allow_pickle=False)
        mjd, band, y, err = (obj[k] for k in ("MJD", "band", "data_flux", "data_fluxerr"))
        assert len(mjd) >= 5
        out = np.column_stack([mjd, band, y, err])
        q = raw[cid]
        phot_rows = np.column_stack([q["MJD"].astype(float),
                                     np.char.strip(q["BAND"].astype(str)),
                                     q["FLUXCAL"].astype(float),
                                     q["FLUXCALERR"].astype(float)])
        dt_raw, df_raw, de_raw = match_rows(out, phot_rows, True)
        assert dt_raw < 0.005 and df_raw < 0.001 and de_raw < 0.001, (cid, dt_raw, df_raw, de_raw)
        accepted = lc.loc[(lc[0] == cid) & (lc[6] == 1)]
        plot_rows = np.column_stack([accepted[1].to_numpy(float), accepted[7].to_numpy(str),
                                     accepted[4].to_numpy(float), accepted[5].to_numpy(float)])
        assert len(plot_rows) == len(out), (cid, len(plot_rows), len(out))
        dt_plot, _, _ = match_rows(out, plot_rows, False)
        assert dt_plot < 0.005, (cid, dt_plot)
        arch = cohort.loc[cid]
        fitrow = native.loc[cid]
        assert str(arch.field_first_phot) == str(q["FIELD"][0]).strip()
        assert abs(float(fitrow.PKMJDINI) - float(arch.PKMJDINI)) < 0.005
        assert abs(float(fitrow.zHEL) - float(arch.zHEL)) < 1e-6
        pars = obj["parameters_x0_x1_c_t0"]
        dp = np.array([float(fitrow[k]) for k in ("x0", "x1", "c", "PKMJD")]) - pars
        assert abs(dp[0]) < 2e-8 and abs(dp[1]) < 1e-5 and abs(dp[2]) < 1e-5 and abs(dp[3]) < 0.005, (cid, dp)
        C = obj["frozen_flux_covariance"]
        assert C.shape == (len(out), len(out))
        mineig = float(np.linalg.eigvalsh(C).min())
        assert mineig > 0
        rows.append(dict(CID=cid, arm=arm, law=info["law"], epochs=len(out),
                         raw_mjd_max_abs=dt_raw, raw_flux_max_abs=df_raw,
                         raw_error_max_abs=de_raw, lcplot_mjd_max_abs=dt_plot,
                         first_phot_field=str(arch.field_first_phot),
                         archived_basic_quality=bool(arch.basic_quality_pass),
                         original_pkmjdini=float(arch.PKMJDINI),
                         new_pkmjdini=float(fitrow.PKMJDINI),
                         archived_pkmjd=float(arch.PKMJD),
                         new_pkmjd=float(fitrow.PKMJD),
                         min_epoch_cov_eigenvalue=mineig))
    frame = pd.DataFrame(rows)
    frame.to_csv(path / "provenance-ledger.csv", index=False, float_format="%.17g")
    result = {"cohort": info["cohort"], "arm": arm, "law": info["law"],
              "selected": info["selected_count"], "exported": len(frame),
              "success_fraction": len(frame) / info["selected_count"],
              "full_pilot_adequate": len(frame) / info["selected_count"] >= 0.95,
              "epochs": int(frame.epochs.sum()),
              "max_raw_mjd_difference": float(frame.raw_mjd_max_abs.max()),
              "max_raw_flux_difference": float(frame.raw_flux_max_abs.max()),
              "max_raw_error_difference": float(frame.raw_error_max_abs.max()),
              "max_lcplot_mjd_difference": float(frame.lcplot_mjd_max_abs.max()),
              "max_pkmjdini_difference": float(np.max(np.abs(frame.new_pkmjdini-frame.original_pkmjdini))),
              "max_pkmjd_difference": float(np.max(np.abs(frame.new_pkmjd-frame.archived_pkmjd))),
              "min_epoch_cov_eigenvalue": float(frame.min_epoch_cov_eigenvalue.min()),
              "source_sha256": sha(Path(__file__)),
              "provenance_ledger_sha256": sha(path / "provenance-ledger.csv")}
    (path / "provenance-gate.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def tangent_gate(path: Path) -> dict:
    """Check native/independent flux units and local nuisance rank; no scores."""
    preflight()
    info = json.loads((path / "prepared.json").read_text())
    provenance = json.loads((path / "provenance-gate.json").read_text())
    assert provenance["full_pilot_adequate"]
    source = ROOT / "scripts/salt_dust_audit/flux_response.py"
    spec = importlib.util.spec_from_file_location("flux_response", source)
    fr = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(fr)
    model, bands, _, zp = fr.build_model()
    magoff = {str(row["Filter Name"])[-1]: float(row["Primary Mag"]) for row in zp}
    ids = pd.read_csv(path / "provenance-ledger.csv").CID.astype(int).tolist()
    rows = []
    for cid in ids:
        obj = np.load(path / "objectives" / f"objective_{cid}.npz", allow_pickle=False)
        pars = obj["parameters_x0_x1_c_t0"]
        row = SimpleNamespace(CID=str(cid), x0=float(pars[0]), x1=float(pars[1]),
                              c=float(pars[2]), PKMJD=float(pars[3]),
                              zHEL=float(obj["zHEL"][0]))
        points = pd.DataFrame({"MJD": obj["MJD"], "BAND": obj["band"],
                               "FLUXCAL": obj["data_flux"],
                               "FLUXCALERR": obj["data_fluxerr"]})
        order = points.sort_values(["MJD", "BAND"], kind="stable").index.to_numpy()
        _, audit, _, diag = fr.analyze(row, points, {"MWEBV": obj["MWEBV"][0]},
                                      model, bands, magoff, False)
        f = obj["model_flux"][order]
        assert np.allclose(obj["MJD"][order], audit["mjd"], atol=0, rtol=0)
        assert np.array_equal(obj["band"][order], audit["band"])
        C = obj["frozen_flux_covariance"][np.ix_(order, order)]
        J = audit["jacobian_flux"].copy()
        J[:, 0] = -fr.K * f
        singular = np.linalg.svd(np.linalg.solve(np.linalg.cholesky(C), J),
                                 compute_uv=False)
        rank = int(np.sum(singular > singular[0] * 1e-10))
        assert rank == 4
        sigma = obj["data_fluxerr"][order]
        delta = audit["flux_model"] - f
        rows.append(dict(CID=cid, epochs=len(f), nuisance_rank=rank,
                         nuisance_condition=float(singular[0] / singular[-1]),
                         max_native_minus_independent_flux_sigma=float(np.max(np.abs(delta / sigma))),
                         rms_native_minus_independent_flux_sigma=float(np.sqrt(np.mean((delta / sigma) ** 2))),
                         derivative_halfstep_relative_error=float(diag["derivative_relative_error"])))
    frame = pd.DataFrame(rows)
    frame.to_csv(path / "tangent-gate-ledger.csv", index=False, float_format="%.17g")
    result = {"arm": info["arm"], "law": info["law"], "objects": len(frame),
              "all_rank4": bool((frame.nuisance_rank == 4).all()),
              "max_native_minus_independent_flux_sigma": float(frame.max_native_minus_independent_flux_sigma.max()),
              "max_derivative_halfstep_relative_error": float(frame.derivative_halfstep_relative_error.max()),
              "flux_response_source_sha256": sha(source),
              "ledger_sha256": sha(path / "tangent-gate-ledger.csv"),
              "scope": "Independent sncosmo tangent diagnostic, no residual scores. The independent model uses modern exact F99 even for the native -99 arm."}
    (path / "tangent-gate.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def mean_edges(path: Path) -> dict:
    """Metadata-only model/flux discrepancy ledger, without residual scoring."""
    assert json.loads((path / "tangent-gate.json").read_text())["all_rank4"]
    source = ROOT / "scripts/salt_dust_audit/flux_response.py"
    spec = importlib.util.spec_from_file_location("flux_response", source)
    fr = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(fr)
    model, bands, _, zp = fr.build_model()
    magoff = {str(row["Filter Name"])[-1]: float(row["Primary Mag"]) for row in zp}
    ids = pd.read_csv(path / "provenance-ledger.csv").CID.astype(int)
    records = []
    for cid in ids:
        obj = np.load(path / "objectives" / f"objective_{cid}.npz", allow_pickle=False)
        pars = obj["parameters_x0_x1_c_t0"]
        row = SimpleNamespace(CID=str(cid), x0=float(pars[0]), x1=float(pars[1]),
                              c=float(pars[2]), PKMJD=float(pars[3]),
                              zHEL=float(obj["zHEL"][0]))
        points = pd.DataFrame({"MJD": obj["MJD"], "BAND": obj["band"],
                               "FLUXCAL": obj["data_flux"],
                               "FLUXCALERR": obj["data_fluxerr"]})
        order = points.sort_values(["MJD", "BAND"], kind="stable").index.to_numpy()
        _, audit, _, _ = fr.analyze(row, points, {"MWEBV": obj["MWEBV"][0]},
                                   model, bands, magoff, False)
        native = obj["model_flux"][order]
        independent = audit["flux_model"]
        err = obj["data_fluxerr"][order]
        delta_sigma = (independent-native)/err
        flag = (native <= 0) | (independent <= 0) | (np.abs(delta_sigma) > 0.1)
        for rank in np.flatnonzero(flag):
            records.append(dict(CID=cid, rank_after_MJD_band_sort=int(rank),
                                MJD=float(obj["MJD"][order][rank]),
                                band=str(obj["band"][order][rank]),
                                native_model_flux=float(native[rank]),
                                independent_model_flux=float(independent[rank]),
                                observed_flux=float(obj["data_flux"][order][rank]),
                                quoted_error=float(err[rank]),
                                independent_minus_native_in_quoted_sigma=float(delta_sigma[rank]),
                                native_nonpositive=bool(native[rank] <= 0),
                                independent_nonpositive=bool(independent[rank] <= 0),
                                model_discrepancy_gt_point1_sigma=bool(abs(delta_sigma[rank]) > 0.1)))
    frame = pd.DataFrame(records)
    frame.to_csv(path / "mean-edge-ledger.csv", index=False, float_format="%.17g")
    result = {"arm": json.loads((path / "prepared.json").read_text())["arm"],
              "law": json.loads((path / "prepared.json").read_text())["law"],
              "eligible_objects": len(ids), "flagged_epochs": len(frame),
              "flagged_objects": int(frame.CID.nunique()) if len(frame) else 0,
              "native_nonpositive_epochs": int(frame.native_nonpositive.sum()) if len(frame) else 0,
              "independent_nonpositive_epochs": int(frame.independent_nonpositive.sum()) if len(frame) else 0,
              "model_discrepancy_gt_point1_sigma_epochs": int(frame.model_discrepancy_gt_point1_sigma.sum()) if len(frame) else 0,
              "ledger_sha256": sha(path / "mean-edge-ledger.csv"),
              "note": "Flag-only; no epochs excluded from native objective or future primary score."}
    (path / "mean-edge-gate.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def score() -> dict:
    """Execute the frozen score and cell transport after all provenance gates."""
    preflight()
    assert sha(SCORE_PLAN) == EXPECTED_SCORE_PLAN
    dest = OUT / "score"
    if dest.exists():
        raise RuntimeError("Preserve existing score output")
    dest.mkdir()
    # Verify every objective gate and mean-edge ledger before reading residuals.
    paths = {}
    for arm in ("P21", "G10"):
        for law in ("approx_minus99", "exact_99"):
            path = OUT / "full256" / arm / law
            prepared = json.loads((path / "prepared.json").read_text())
            fit = json.loads((path / "fit-gate.json").read_text())
            exported = json.loads((path / "export-gate.json").read_text())
            provenance = json.loads((path / "provenance-gate.json").read_text())
            tangent = json.loads((path / "tangent-gate.json").read_text())
            edges = json.loads((path / "mean-edge-gate.json").read_text())
            assert prepared["selected_count"] == 256
            assert fit["returncode"] == exported["returncode"] == 0
            assert fit["graceful"] and provenance["full_pilot_adequate"]
            assert tangent["all_rank4"] and edges["eligible_objects"] == provenance["exported"]
            assert sha(path / "mean-edge-ledger.csv") == edges["ledger_sha256"]
            assert sha(path / "provenance-ledger.csv") == provenance["provenance_ledger_sha256"]
            paths[(arm, law)] = path
    source = ROOT / "scripts/salt_dust_audit/flux_response.py"
    spec = importlib.util.spec_from_file_location("flux_response", source)
    fr = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(fr)
    model, bands, _, zp = fr.build_model()
    magoff = {str(row["Filter Name"])[-1]: float(row["Primary Mag"]) for row in zp}
    coeff_file = ROOT / "runs/research_2026_09_26/astra_design/validation1020/frozen-discovery-coefficients.npz"
    with np.load(coeff_file, allow_pickle=False) as frozen:
        c = frozen["basis_mean"].astype(float)
    assert c.shape == (3,) and np.all(np.isfinite(c))
    support = pd.read_csv(DESIGN / "truth-support-ledger.csv")
    assert support.groupby(["arm", "support_class"]).size().to_dict() == {
        ("G10", "no_explicit_host_screen_sentinel"): 256,
        ("P21", "negative_passive_attenuation"): 14,
        ("P21", "nonnegative_declared_grid"): 242}
    records = []
    sufficient = {}
    for (arm, law), path in paths.items():
        cohort = pd.read_csv(DESIGN / f"{arm}-cohort.csv").set_index("CID")
        accepted = set(pd.read_csv(path / "provenance-ledger.csv").CID.astype(int))
        other = set(pd.read_csv(paths[(arm, "exact_99" if law == "approx_minus99" else "approx_minus99")]
                                / "provenance-ledger.csv").CID.astype(int))
        ids = sorted(accepted & other)
        assert len(ids) >= 244  # registered >=95% of 256
        sys_path = ROOT / "scripts/phase2/official"
        import sys
        sys.path.insert(0, str(sys_path))
        from audit_fits import read_fit
        newfit = read_fit(path / "fit.FITRES.TEXT")
        newfit.CID = newfit.CID.astype(int)
        newfit = newfit.set_index("CID")
        for cid in ids:
            obj = np.load(path / "objectives" / f"objective_{cid}.npz", allow_pickle=False)
            pars = obj["parameters_x0_x1_c_t0"]
            row = SimpleNamespace(CID=str(cid), x0=float(pars[0]), x1=float(pars[1]),
                                  c=float(pars[2]), PKMJD=float(pars[3]),
                                  zHEL=float(obj["zHEL"][0]))
            points = pd.DataFrame({"MJD": obj["MJD"], "BAND": obj["band"],
                                   "FLUXCAL": obj["data_flux"],
                                   "FLUXCALERR": obj["data_fluxerr"]})
            order = points.sort_values(["MJD", "BAND"], kind="stable").index.to_numpy()
            _, audit, _, diagnostic = fr.analyze(row, points,
                {"MWEBV": obj["MWEBV"][0]}, model, bands, magoff, False)
            f = obj["model_flux"][order]
            y = obj["data_flux"][order]
            b = obj["band"][order]
            C = obj["frozen_flux_covariance"][np.ix_(order, order)]
            J = audit["jacobian_flux"].copy()
            J[:, 0] = -fr.K * f
            L = np.linalg.cholesky(C)
            wj = solve_triangular(L, J, lower=True)
            U, s, _ = np.linalg.svd(wj, full_matrices=True)
            rank = int(np.sum(s > s[0] * 1e-10))
            assert rank == 4
            Q = U[:, 4:]
            residual = Q.T @ solve_triangular(L, y-f, lower=True)
            griz = np.column_stack([-fr.K*f*(b == band) for band in "griz"])
            observer = np.column_stack([griz[:, 0]-griz[:, 1],
                                        griz[:, 2]-griz[:, 1],
                                        griz[:, 3]-griz[:, 1]])
            T = Q.T @ solve_triangular(L, observer, lower=True)
            u = T.T @ residual
            F = T.T @ T
            a = float(c @ u)
            information = float(c @ F @ c)
            assert information > 0 and np.isfinite(a)
            gain = a - information / 2
            archived = cohort.loc[cid]
            fitrow = newfit.loc[cid]
            quality = (abs(float(fitrow.x1)) < 3 and abs(float(fitrow.c)) < .3
                       and float(fitrow.x1ERR) < 1 and float(fitrow.PKMJDERR) < 2
                       and float(fitrow.cERR) < 1.5 and float(fitrow.FITPROB) > .001
                       and .025 < float(fitrow.zHD) < 1.2)
            tag = support.loc[(support.arm == arm) & (support.CID == cid), "support_class"]
            assert len(tag) == 1
            records.append(dict(arm=arm, law=law, CID=cid, field=str(archived.field),
                                zHEL=float(archived.zHEL), zbin=int(archived.zbin),
                                cell=str(archived.cell), LIBID=int(archived.SIM_LIBID),
                                SNRMAX1_archived=float(archived.SNRMAX1),
                                support_class=str(tag.iloc[0]),
                                archived_basic_quality=bool(archived.basic_quality_pass),
                                new_basic_quality=bool(quality),
                                epochs=len(y), projected_dimension=len(residual),
                                projected_chi2=float(residual @ residual),
                                matched_filter=a, information=information,
                                fixed_gain=gain, nuisance_condition=float(s[0]/s[-1]),
                                max_native_independent_mean_sigma=float(np.max(np.abs((audit["flux_model"]-f)/obj["data_fluxerr"][order]))),
                                independent_derivative_halfstep_error=float(diagnostic["derivative_relative_error"])))
            sufficient[(arm, law, cid)] = (u, F)
    scored = pd.DataFrame(records).sort_values(["arm", "law", "CID"])
    scored.to_csv(dest / "object-scores.csv", index=False, float_format="%.17g")
    np.savez_compressed(dest / "sufficient-arrays.npz",
                        arm=scored.arm.to_numpy(), law=scored.law.to_numpy(), CID=scored.CID.to_numpy(int),
                        u=np.array([sufficient[(r.arm,r.law,r.CID)][0] for r in scored.itertuples()]),
                        F=np.array([sufficient[(r.arm,r.law,r.CID)][1] for r in scored.itertuples()]))
    assert len(scored) >= 4*244
    # Transport source is the frozen real validation1020, not simulated truth.
    real_path = ROOT / "runs/research_2026_09_26/astra_design/validation1020/analysis/object-scores.csv"
    real_cohort_path = ROOT / "runs/research_2026_09_26/astra_design/validation1020/cohort.csv"
    real_snr_path = ROOT / "runs/research_2026_09_26/astra_design/validation1020/fit.FITRES.TEXT"
    real = pd.read_csv(real_path).query("arm == 'published_mask'").copy()
    assert len(real) == 1020 and real.CID.is_unique
    real_cohort = pd.read_csv(real_cohort_path)[["CID", "field", "zHEL"]]
    assert set(real.CID) == set(real_cohort.CID)
    assert real.set_index("CID").field.sort_index().equals(real_cohort.set_index("CID").field.sort_index())
    from audit_fits import read_fit
    realfit = read_fit(real_snr_path)
    realfit.CID = realfit.CID.astype(int)
    real = real.merge(realfit[["CID", "SNRMAX1"]], on="CID", validate="one_to_one")
    edges = [.05, .2, .35, .5, .65, .8, 1.2]
    real["zbin"] = pd.cut(real.zHEL, edges, include_lowest=True, labels=False)
    assert real.zbin.notna().all()
    real["zbin"] = real.zbin.astype(int)
    real["cell"] = real.field + "_" + real.zbin.astype(str)
    real["snrbin"] = np.where(real.SNRMAX1 < 15, "lt15", "ge15")
    real["fixed_gain"] = real.fixed_prediction_gain
    scored["snrbin"] = np.where(scored.SNRMAX1_archived < 15, "lt15", "ge15")
    stages = ("common_success", "archived_quality", "new_quality", "p21_grid_nonnegative")
    rng = np.random.default_rng(26092691)
    libs = np.array(sorted(set(scored.LIBID.astype(int))))
    draws = rng.integers(0, len(libs), size=(1000, len(libs)))
    multiplicities = np.zeros((1000, len(libs)), dtype=np.int16)
    for j in range(1000):
        multiplicities[j] = np.bincount(draws[j], minlength=len(libs))
    libindex = {lib: i for i, lib in enumerate(libs)}
    summary_rows = []
    field_rows = []
    bootstrap_rows = []
    for law in ("approx_minus99", "exact_99"):
        for stage in stages:
            selected = {}
            for arm in ("P21", "G10"):
                q = scored[(scored.arm == arm) & (scored.law == law)].copy()
                if stage == "archived_quality":
                    q = q[q.archived_basic_quality]
                elif stage == "new_quality":
                    q = q[q.new_basic_quality]
                elif stage == "p21_grid_nonnegative" and arm == "P21":
                    q = q[q.support_class == "nonnegative_declared_grid"]
                selected[arm] = q
            for split in ("field_z", "field_z_snr15"):
                real_df = real.copy()
                for arm in selected:
                    selected[arm] = selected[arm].copy()
                    selected[arm]["transport_cell"] = (selected[arm].cell if split == "field_z"
                                                        else selected[arm].cell + "_" + selected[arm].snrbin)
                real_df["transport_cell"] = (real_df.cell if split == "field_z"
                                             else real_df.cell + "_" + real_df.snrbin)
                counts_p = selected["P21"].transport_cell.value_counts()
                counts_g = selected["G10"].transport_cell.value_counts()
                supported = sorted(set(counts_p[counts_p >= 2].index) &
                                   set(counts_g[counts_g >= 2].index))
                supported = [cell for cell in supported if (real_df.transport_cell == cell).any()]
                target = real_df[real_df.transport_cell.isin(supported)]
                target_a = float(target.matched_filter.sum())
                target_i = float(target.information.sum())
                target_w = len(target)
                assert target_i > 0 and target_w > 0
                for arm in ("P21", "G10"):
                    full = selected[arm]
                    q = full[full.transport_cell.isin(supported)].copy()
                    n_per = q.transport_cell.value_counts()
                    real_per = target.transport_cell.value_counts()
                    q["base_weight"] = q.transport_cell.map(real_per).to_numpy(float) / q.transport_cell.map(n_per).to_numpy(float)
                    w = q.base_weight.to_numpy(float)
                    assert abs(w.sum()-target_w) < 1e-8
                    a = q.matched_filter.to_numpy(float)
                    information = q.information.to_numpy(float)
                    gain = q.fixed_gain.to_numpy(float)
                    weighted_a = float(w @ a)
                    weighted_i = float(w @ information)
                    assert weighted_i > 0
                    suffix = dict(arm=arm, law=law, stage=stage, transport=split)
                    summary_rows.append(dict(**suffix, selected=len(full), supported_objects=len(q),
                        supported_cells=len(supported), supported_sim_fraction=len(q)/len(full),
                        supported_real_objects=target_w, supported_real_fraction=target_w/len(real),
                        raw_a=float(full.matched_filter.sum()), raw_I=float(full.information.sum()),
                        raw_G=float(full.fixed_gain.sum()),
                        supported_unweighted_a=float(a.sum()), supported_unweighted_I=float(information.sum()),
                        supported_unweighted_G=float(gain.sum()),
                        transported_amplitude=weighted_a/weighted_i,
                        transported_mean_gain=float(w @ gain / w.sum()),
                        real_supported_a=target_a, real_supported_I=target_i,
                        real_supported_G=float(target.fixed_gain.sum()),
                        real_supported_amplitude=target_a/target_i,
                        real_supported_mean_gain=float(target.fixed_gain.mean()),
                        supported_cell_ids="|".join(supported)))
                    for field, frows in q.groupby("field"):
                        ww = frows.base_weight.to_numpy(float)
                        field_rows.append(dict(**suffix, field=field, n=len(frows),
                            raw_a=float(frows.matched_filter.sum()),
                            raw_I=float(frows.information.sum()), raw_G=float(frows.fixed_gain.sum()),
                            weighted_a=float(ww @ frows.matched_filter.to_numpy(float)),
                            weighted_I=float(ww @ frows.information.to_numpy(float)),
                            weighted_G=float(ww @ frows.fixed_gain.to_numpy(float))))
                    mult = multiplicities[:, [libindex[int(x)] for x in q.LIBID]]
                    bw = mult * w[None, :]
                    den_i = bw @ information
                    den_w = bw.sum(axis=1)
                    good = (den_i > 0) & (den_w > 0)
                    amp = np.full(1000, np.nan)
                    mean_gain = np.full(1000, np.nan)
                    amp[good] = (bw[good] @ a)/den_i[good]
                    mean_gain[good] = (bw[good] @ gain)/den_w[good]
                    for rep in range(1000):
                        bootstrap_rows.append(dict(**suffix, replicate=rep,
                            valid=bool(good[rep]), amplitude=float(amp[rep]),
                            mean_gain=float(mean_gain[rep])))
    summary = pd.DataFrame(summary_rows)
    boot = pd.DataFrame(bootstrap_rows)
    for key, group in boot.groupby(["arm", "law", "stage", "transport"]):
        mask = np.logical_and.reduce([summary[k] == v for k, v in zip(("arm", "law", "stage", "transport"), key)])
        valid = group[group.valid]
        summary.loc[mask, "bootstrap_failures"] = len(group)-len(valid)
        for col in ("amplitude", "mean_gain"):
            if len(valid):
                lo, hi = np.quantile(valid[col], [.025,.975])
                summary.loc[mask, "bootstrap_"+col+"_ci_low"] = lo
                summary.loc[mask, "bootstrap_"+col+"_ci_high"] = hi
    summary.to_csv(dest / "transport-summary.csv", index=False, float_format="%.17g")
    pd.DataFrame(field_rows).to_csv(dest / "field-contributions.csv", index=False, float_format="%.17g")
    boot.to_csv(dest / "bootstrap.csv.gz", index=False, compression={"method": "gzip", "mtime": 0}, float_format="%.17g")
    inputs = [SCORE_PLAN, PROTO, LAW, DESIGN/"handoff-manifest.json",
              DESIGN/"truth-support-ledger.csv", coeff_file, source,
              real_path, real_cohort_path, real_snr_path, Path(__file__)]
    inputs += [p/f for p in paths.values() for f in ("prepared.json", "fit-gate.json", "export-gate.json",
                                               "provenance-gate.json", "tangent-gate.json", "mean-edge-gate.json",
                                               "provenance-ledger.csv", "mean-edge-ledger.csv", "objectives/manifest.json")]
    manifest = {"inputs_sha256": {str(p.relative_to(ROOT)): sha(p) for p in inputs},
                "outputs_sha256": {str(p.relative_to(ROOT)): sha(p) for p in dest.iterdir() if p.name != "manifest.json"},
                "frozen_coefficients": c.tolist(), "bootstrap_seed": 26092691,
                "bootstrap_replicates": 1000, "score_plan_sha256": EXPECTED_SCORE_PLAN,
                "scope": "Frozen diagnostic on P21/G10 native refits; no classifier or cosmology claim."}
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return {"object_rows": len(scored), "transport_rows": len(summary),
            "bootstrap_rows": len(boot), "max_bootstrap_failures": int(summary.bootstrap_failures.max()),
            "output": str(dest)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("prepare", "fit", "export", "gate", "tangent", "mean_edges", "score"))
    parser.add_argument("cohort", choices=("engineering8", "full256"), nargs="?")
    parser.add_argument("arm", choices=("P21", "G10"), nargs="?")
    parser.add_argument("law", choices=("approx_minus99", "exact_99"), nargs="?")
    args = parser.parse_args()
    if args.action == "score":
        print(json.dumps(score(), indent=2))
        return
    assert args.cohort and args.arm and args.law
    path = OUT / args.cohort / args.arm / args.law
    if args.action == "prepare":
        path = make_input(args.arm, args.law, args.cohort)
        result = {"prepared": str(path), "nml_sha256": sha(path / "fit.nml")}
    elif args.action == "fit":
        result = run_fit(path)
    elif args.action == "export":
        result = export(path)
    elif args.action == "gate":
        result = gate(path)
    elif args.action == "tangent":
        result = tangent_gate(path)
    else:
        result = mean_edges(path)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
