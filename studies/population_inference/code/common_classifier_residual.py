"""Frozen common measured-input SuperNNova sensitivity for P21/G10 pilots."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/research_2026_09_26/common_classifier_residual"
DESIGN = ROOT / "runs/research_2026_09_26/astra_design/simulation_design"
CONTROL = ROOT / "runs/research_2026_09_26/simulation_residual_control"
PROTO = ROOT / "runs/research_2026_09_26/astra_design/common-classifier-protocol.md"
MANIFEST = ROOT / "runs/research_2026_09_26/astra_design/common-classifier-manifest.json"
CLUMP = ROOT / "phase2/classification/reconstruction_20260926/clump_run"
CLASS = ROOT / "phase2/classification/reconstruction_20260926"
VAL = ROOT / "runs/research_2026_09_26/astra_design/validation1020"
MODEL = ROOT / "phase2/official/inputs/SNDATA_ROOT/models/classifiers/DES-SN5YR/SNNTRAINV19_z_TRAINDES_V19/model.pt"
PY = ROOT / "phase2/classification/env/bin/python"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def preflight() -> None:
    frozen = json.loads(MANIFEST.read_text())["sha256"]
    for rel, expected in frozen.items():
        assert digest(ROOT / rel) == expected, rel
    assert digest(PROTO) == "e060491bc41b217f0c78065e08ffc60009b7e12b914a617c197bbdc663a358b8"
    assert (CONTROL / "score/manifest.json").exists()
    assert (CONTROL / "score-independent-verification.json").exists()


def wrapper():
    path = ROOT / "scripts/phase2/classification/reconstruct.py"
    spec = importlib.util.spec_from_file_location("frozen_snn_wrapper", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(mod)
    return mod


def cohort_ids(arm: str) -> list[str]:
    ids = [s.strip() for s in (DESIGN / f"{arm}-cids.txt").read_text().splitlines() if s.strip()]
    assert len(ids) == 256 and len(set(ids)) == 256
    return ids


def prepare_clump(arm: str) -> Path:
    preflight()
    assert arm in ("P21", "G10")
    dest = OUT / "clump" / arm
    dest.mkdir(parents=True, exist_ok=False)
    version = f"PH2_pilot02_{arm}"
    base = (CLUMP / "clump.nml").read_text()
    old_path = "PRIVATE_DATA_PATH = '/home/szymon/Documents/ChatGPT/supernova/phase2/official/inputs/SNDATA_ROOT/lcmerge/DES-SN5YR'"
    new_path = f"PRIVATE_DATA_PATH = '{(ROOT / 'phase2/literature/simulations/outputs').resolve()}'"
    assert base.count(old_path) == 1
    base = base.replace(old_path, new_path)
    assert base.count("VERSION_PHOTOMETRY = 'DES-SN5YR_DES'") == 1
    base = base.replace("VERSION_PHOTOMETRY = 'DES-SN5YR_DES'", f"VERSION_PHOTOMETRY = '{version}'")
    assert base.count("TEXTFILE_PREFIX = 'DES-SN5YR_DES'") == 1
    base = base.replace("TEXTFILE_PREFIX = 'DES-SN5YR_DES'", f"TEXTFILE_PREFIX = '{version}'")
    assert "OPT_SETPKMJD = 16" in base and "PHOTFLAG_MSKREJ = 1016" in base
    (dest / "clump.nml").write_text(base)
    (dest / "prepared.json").write_text(json.dumps({
        "arm": arm, "version": version,
        "base_nml_sha256": digest(CLUMP / "clump.nml"),
        "nml_sha256": digest(dest / "clump.nml"),
        "binary_sha256": digest(CLUMP / "snana.exe"),
        "changed_fields": ["PRIVATE_DATA_PATH", "VERSION_PHOTOMETRY", "TEXTFILE_PREFIX"],
        "archived_fit_setpkmjd": 20, "common_clump_setpkmjd": 16,
        "truth_peak_used": False}, indent=2) + "\n")
    return dest


def run_clump(arm: str) -> dict:
    preflight()
    path = OUT / "clump" / arm
    prep = json.loads((path / "prepared.json").read_text())
    assert digest(path / "clump.nml") == prep["nml_sha256"]
    assert not (path / "clump.log").exists()
    env = os.environ.copy()
    env.update(SNANA_DIR=str(ROOT / "phase2/official/build/SNANA-current"),
               SNDATA_ROOT=str(ROOT / "phase2/official/inputs/SNDATA_ROOT"),
               LD_LIBRARY_PATH=str(ROOT / "phase2/official/build/sysroot/usr/lib"),
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (path / "clump.log").open("x") as stream:
        proc = subprocess.run([str(CLUMP / "snana.exe"), str(path / "clump.nml")],
                              cwd=path, env=env, stdout=stream, stderr=subprocess.STDOUT,
                              check=False)
    output = path / f"PH2_pilot02_{arm}.SNANA.TEXT"
    status = {"arm": arm, "returncode": proc.returncode,
              "graceful": "ENDING PROGRAM GRACEFULLY." in (path / "clump.log").read_text(),
              "output_exists": output.exists(), "log_sha256": digest(path / "clump.log")}
    if output.exists():
        status["output_sha256"] = digest(output)
    (path / "run-gate.json").write_text(json.dumps(status, indent=2) + "\n")
    return status


def read_peak(path: Path) -> pd.DataFrame:
    names = None
    records = []
    for line in path.read_text().splitlines():
        if line.startswith("VARNAMES:"):
            names = line.split()[1:]
        elif line.startswith("SN:"):
            assert names is not None
            values = line.split()[1:]
            assert len(values) == len(names)
            records.append(dict(zip(names, values)))
    assert names and "PKMJDINI" in names
    q = pd.DataFrame(records)
    assert q.CID.is_unique
    q["PKMJDINI"] = pd.to_numeric(q.PKMJDINI, errors="raise")
    assert np.isfinite(q.PKMJDINI).all()
    return q


def classify() -> dict:
    preflight()
    OUT.mkdir(parents=True, exist_ok=True)
    dest = OUT / "inference"
    dest.mkdir(exist_ok=False)
    wrapped = wrapper()
    records = []
    single_checks = []
    peak_rows = []
    for arm in ("P21", "G10"):
        path = OUT / "clump" / arm
        gate = json.loads((path / "run-gate.json").read_text())
        assert gate["returncode"] == 0 and gate["graceful"] and gate["output_exists"]
        ids = cohort_ids(arm)
        peak = read_peak(path / f"PH2_pilot02_{arm}.SNANA.TEXT")
        assert set(ids).issubset(set(peak.CID))
        peaks = peak.set_index("CID").PKMJDINI.loc[ids]
        archived = pd.read_csv(DESIGN / f"{arm}-cohort.csv", dtype={"CID": str}).set_index("CID")
        for cid in ids:
            peak_rows.append(dict(arm=arm, CID=cid, common_peak=float(peaks[cid]),
                                  archived_fit_peak=float(archived.loc[cid, "PKMJDINI"]),
                                  common_minus_archived=float(peaks[cid]-archived.loc[cid, "PKMJDINI"])))
        stem = ROOT / f"phase2/literature/simulations/outputs/PH2_pilot02_{arm}/PH2_pilot02_{arm}"
        df, meta = wrapped.load_curves(stem, set(ids), peaks.to_dict())
        assert len(meta) == 256 and set(meta.CID) == set(ids)
        # Declared engineering ranks, fixed before inference.
        check_ids = [ids[int(i)] for i in np.linspace(0, 255, 5, dtype=int)]
        subset = df.loc[df.SNID.isin(check_ids)]
        alone = wrapped.infer(subset, chunk_size=1).set_index("CID").pIa
        batched = wrapped.infer(df, chunk_size=128).set_index("CID").pIa
        delta = abs(alone - batched.loc[alone.index])
        assert delta.max() <= 1e-6
        assert ((alone > .999) == (batched.loc[alone.index] > .999)).all()
        single_checks.extend(dict(arm=arm, CID=cid, single=float(alone[cid]),
                                  batched=float(batched[cid]), abs_delta=float(delta[cid]))
                             for cid in alone.index)
        out = meta.merge(batched.rename("pIa"), left_on="CID", right_index=True,
                         validate="one_to_one")
        out.insert(0, "arm", arm)
        out["gt999"] = out.pIa > .999
        records.append(out)
    allsim = pd.concat(records, ignore_index=True).sort_values(["arm", "CID"])
    assert len(allsim) == 512
    allsim.to_csv(dest / "simulation-probabilities.csv", index=False, float_format="%.17g")
    pd.DataFrame(peak_rows).sort_values(["arm", "CID"]).to_csv(dest / "peak-ledger.csv", index=False, float_format="%.17g")
    pd.DataFrame(single_checks).sort_values(["arm", "CID"]).to_csv(dest / "single-batch-check.csv", index=False, float_format="%.17g")
    # Recompute real probabilities using the same preserved current clump peaks.
    real_ids = set(pd.read_csv(VAL / "cohort.csv", dtype={"CID": str}).CID)
    assert len(real_ids) == 1020
    realpeak = read_peak(CLUMP / "DES-SN5YR_DES.SNANA.TEXT").set_index("CID").PKMJDINI
    assert real_ids <= set(realpeak.index)
    stem = ROOT / "sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES/DES-SN5YR_DES"
    real_df, real_meta = wrapped.load_curves(stem, real_ids, realpeak.loc[sorted(real_ids)].to_dict())
    real_pred = wrapped.infer(real_df, chunk_size=128)
    real_out = real_meta.merge(real_pred, on="CID", validate="one_to_one")
    reference = pd.read_csv(CLASS / "des_clump_diagnostic.csv", dtype={"CID": str}).set_index("CID")
    chk = real_out.set_index("CID").join(reference[["pIa"]], rsuffix="_preserved")
    assert len(chk) == 1020 and chk.pIa_preserved.notna().all()
    delta = abs(chk.pIa - chk.pIa_preserved)
    assert delta.max() <= 1e-6
    assert ((chk.pIa > .999) == (chk.pIa_preserved > .999)).all()
    real_out["gt999"] = real_out.pIa > .999
    real_out.sort_values("CID").to_csv(dest / "real-probabilities.csv", index=False, float_format="%.17g")
    result = {"simulated": len(allsim), "real": len(real_out),
              "max_single_batch_abs": float(pd.DataFrame(single_checks).abs_delta.max()),
              "max_real_preserved_abs": float(delta.max()),
              "real_gt999": int(real_out.gt999.sum()),
              "threshold_near_1e6": allsim.loc[abs(allsim.pIa-.999) <= 1e-6,
                  ["arm", "CID", "pIa"]].to_dict("records"),
              "sim_gt999": allsim.groupby("arm").gt999.sum().astype(int).to_dict()}
    (dest / "inference-gate.json").write_text(json.dumps(result, indent=2) + "\n")
    sources = [MANIFEST, PROTO, MODEL, Path(__file__),
               CLASS / "des_clump_diagnostic.csv", CLUMP / "DES-SN5YR_DES.SNANA.TEXT"]
    sources += [OUT / "clump" / arm / f"PH2_pilot02_{arm}.SNANA.TEXT" for arm in ("P21", "G10")]
    (dest / "manifest.json").write_text(json.dumps({
        "inputs_sha256": {str(p.relative_to(ROOT)): digest(p) for p in sources},
        "outputs_sha256": {str(p.relative_to(ROOT)): digest(p) for p in dest.iterdir() if p.name != "manifest.json"}}, indent=2) + "\n")
    return result


def score() -> dict:
    """Apply frozen membership to existing, immutable native per-object scores."""
    preflight()
    inf = OUT / "inference"
    gate = json.loads((inf / "inference-gate.json").read_text())
    assert gate["simulated"] == 512 and gate["real"] == 1020
    assert gate["max_real_preserved_abs"] <= 1e-6 and gate["max_single_batch_abs"] <= 1e-6
    dest = OUT / "score"
    dest.mkdir(exist_ok=False)
    sim = pd.read_csv(inf / "simulation-probabilities.csv", dtype={"CID": int})
    realpred = pd.read_csv(inf / "real-probabilities.csv", dtype={"CID": int})
    old = pd.read_csv(ROOT / "runs/research_2026_09_26/classifier_residual_sensitivity/membership_and_scores.csv")
    good_real = set(realpred.loc[realpred.gt999, "CID"])
    assert len(good_real) == 1006
    assert good_real == set(old.loc[old.retained_both_gt_0p999, "CID"])
    scores = pd.read_csv(CONTROL / "score/object-scores.csv")
    assert len(scores) == 1020
    # The same one probability per CID/arm selects both native-law score outputs.
    scores = scores.merge(sim[["arm", "CID", "pIa", "gt999"]], on=["arm", "CID"],
                          validate="many_to_one")
    assert len(scores) == 1020 and scores.groupby(["arm", "CID"]).law.nunique().eq(2).all()
    real = pd.read_csv(VAL / "analysis/object-scores.csv")
    real = real.loc[(real.arm == "published_mask") & (real.CID.isin(good_real))].copy()
    assert len(real) == 1006 and real.CID.is_unique
    real["fixed_gain"] = real.fixed_prediction_gain
    edges = [.05, .2, .35, .5, .65, .8, 1.2]
    real["zbin"] = pd.cut(real.zHEL, edges, include_lowest=True, labels=False)
    assert real.zbin.notna().all()
    real["zbin"] = real.zbin.astype(int)
    real["cell"] = real.field + "_" + real.zbin.astype(str)
    # The real SNR split is the same measured FITRES source as the primary control.
    from sys import path as syspath
    syspath.insert(0, str(ROOT / "scripts/phase2/official"))
    from audit_fits import read_fit
    rf = read_fit(VAL / "fit.FITRES.TEXT")[["CID", "SNRMAX1"]]
    rf.CID = rf.CID.astype(int)
    real = real.merge(rf, on="CID", validate="one_to_one")
    real["snrbin"] = np.where(real.SNRMAX1 < 15, "lt15", "ge15")
    scores["snrbin"] = np.where(scores.SNRMAX1_archived < 15, "lt15", "ge15")
    # Additional unweighted stage ledger retains all pilot/classifier counts.
    stage_ledger = []
    for arm in ("P21", "G10"):
        pilot = sim.loc[sim.arm == arm]
        assert len(pilot) == 256
        for law in ("approx_minus99", "exact_99"):
            q = scores.loc[(scores.arm == arm) & (scores.law == law)]
            assert len(q) == 255
            for stage, sub in (
                ("native_success", q),
                ("archived_quality", q.loc[q.archived_basic_quality]),
                ("classifier_only", q.loc[q.gt999]),
                ("quality_and_classifier", q.loc[q.archived_basic_quality & q.gt999]),
                ("quality_classifier_p21_nonnegative", q.loc[q.archived_basic_quality & q.gt999 & (q.support_class == "nonnegative_declared_grid")]
                 if arm == "P21" else q.loc[q.archived_basic_quality & q.gt999])):
                stage_ledger.append(dict(arm=arm, law=law, stage=stage,
                                         pilot=256, valid_inference=len(pilot),
                                         classifier_accepted=int(pilot.gt999.sum()),
                                         native_success=len(q), selected=len(sub),
                                         a=float(sub.matched_filter.sum()),
                                         I=float(sub.information.sum()),
                                         G=float(sub.fixed_gain.sum()),
                                         CIDs="|".join(map(str, sorted(sub.CID)))))
    pd.DataFrame(stage_ledger).to_csv(dest / "stage-ledger.csv", index=False, float_format="%.17g")
    libs = np.array(sorted(set(scores.LIBID.astype(int))))
    rng = np.random.default_rng(26092691)
    draws = rng.integers(0, len(libs), size=(1000, len(libs)))
    mult = np.array([np.bincount(row, minlength=len(libs)) for row in draws], dtype=np.int16)
    libindex = {lib: i for i, lib in enumerate(libs)}
    out = []
    boots = []
    fields = []
    supports = []
    for law in ("approx_minus99", "exact_99"):
        for split in ("field_z", "field_z_snr15"):
            r = real.copy()
            r["transport_cell"] = r.cell if split == "field_z" else r.cell + "_" + r.snrbin
            subsets = {}
            for arm in ("P21", "G10"):
                base = scores.loc[(scores.arm == arm) & (scores.law == law) & scores.archived_basic_quality].copy()
                base["transport_cell"] = base.cell if split == "field_z" else base.cell + "_" + base.snrbin
                subsets[(arm, "pre_classifier_quality")] = base
                subsets[(arm, "post_classifier_quality")] = base.loc[base.gt999].copy()
                if arm == "P21":
                    subsets[(arm, "post_classifier_p21_nonnegative")] = base.loc[base.gt999 & (base.support_class == "nonnegative_declared_grid")].copy()
                else:
                    subsets[(arm, "post_classifier_p21_nonnegative")] = base.loc[base.gt999].copy()
            for support_source in ("own_stage", "post_classifier_common"):
                for stage in ("pre_classifier_quality", "post_classifier_quality", "post_classifier_p21_nonnegative"):
                    if support_source == "post_classifier_common" and stage == "post_classifier_p21_nonnegative":
                        continue
                    support_stage = stage if support_source == "own_stage" else "post_classifier_quality"
                    counts_p = subsets[("P21", support_stage)].transport_cell.value_counts()
                    counts_g = subsets[("G10", support_stage)].transport_cell.value_counts()
                    cells = sorted(set(counts_p[counts_p >= 2].index) & set(counts_g[counts_g >= 2].index)
                                   & set(r.transport_cell))
                    target = r.loc[r.transport_cell.isin(cells)]
                    assert len(target) and target.information.sum() > 0
                    for arm in ("P21", "G10"):
                        full = subsets[(arm, stage)]
                        q = full.loc[full.transport_cell.isin(cells)].copy()
                        assert len(q)
                        real_per = target.transport_cell.value_counts()
                        sim_per = q.transport_cell.value_counts()
                        # Fixed support is shared, but each evaluated stage has its own
                        # within-cell sample count and thus its own transport weights.
                        w = q.transport_cell.map(real_per).to_numpy(float) / q.transport_cell.map(sim_per).to_numpy(float)
                        assert abs(w.sum()-len(target)) < 1e-8
                        a = q.matched_filter.to_numpy(float)
                        information = q.information.to_numpy(float)
                        gain = q.fixed_gain.to_numpy(float)
                        suffix = dict(arm=arm, law=law, transport=split, stage=stage,
                                      support_source=support_source)
                        out.append(dict(**suffix, selected=len(full), supported_sim=len(q),
                                        supported_cells=len(cells), supported_real=len(target),
                                        real_fraction=len(target)/1006,
                                        raw_a=float(full.matched_filter.sum()),
                                        raw_I=float(full.information.sum()),
                                        raw_G=float(full.fixed_gain.sum()),
                                        weighted_a=float(w @ a), weighted_I=float(w @ information),
                                        weighted_G=float(w @ gain),
                                        transported_amplitude=float(w @ a / (w @ information)),
                                        transported_mean_gain=float(w @ gain / w.sum()),
                                        real_a=float(target.matched_filter.sum()),
                                        real_I=float(target.information.sum()),
                                        real_G=float(target.fixed_gain.sum()),
                                        real_amplitude=float(target.matched_filter.sum()/target.information.sum()),
                                        real_mean_gain=float(target.fixed_gain.mean()),
                                        cell_ids="|".join(cells)))
                        supports.extend(dict(**suffix, CID=int(x.CID), cell=x.transport_cell,
                                             base_weight=float(weight))
                                        for x, weight in zip(q.itertuples(), w))
                        for field, fq in q.groupby("field"):
                            fw = w[np.asarray(q.field == field)]
                            fields.append(dict(**suffix, field=field, n=len(fq),
                                               weighted_a=float(fw @ fq.matched_filter.to_numpy(float)),
                                               weighted_I=float(fw @ fq.information.to_numpy(float)),
                                               weighted_G=float(fw @ fq.fixed_gain.to_numpy(float))))
                        bweight = mult[:, [libindex[int(x)] for x in q.LIBID]] * w[None, :]
                        denom_i = bweight @ information
                        denom_w = bweight.sum(axis=1)
                        good = (denom_i > 0) & (denom_w > 0)
                        for j in range(1000):
                            boots.append(dict(**suffix, replicate=j, valid=bool(good[j]),
                                              amplitude=float((bweight[j] @ a)/denom_i[j]) if good[j] else np.nan,
                                              mean_gain=float((bweight[j] @ gain)/denom_w[j]) if good[j] else np.nan))
    result = pd.DataFrame(out)
    boot = pd.DataFrame(boots)
    for key, group in boot.groupby(["arm", "law", "transport", "stage", "support_source"]):
        columns = ("arm", "law", "transport", "stage", "support_source")
        mask = np.logical_and.reduce([result[col] == value for col, value in zip(columns, key)])
        valid = group.loc[group.valid]
        result.loc[mask, "bootstrap_failures"] = len(group)-len(valid)
        for name in ("amplitude", "mean_gain"):
            result.loc[mask, f"bootstrap_{name}_ci_low"] = valid[name].quantile(.025) if len(valid) else np.nan
            result.loc[mask, f"bootstrap_{name}_ci_high"] = valid[name].quantile(.975) if len(valid) else np.nan
    result.to_csv(dest / "transport-summary.csv", index=False, float_format="%.17g")
    pd.DataFrame(fields).to_csv(dest / "field-contributions.csv", index=False, float_format="%.17g")
    pd.DataFrame(supports).to_csv(dest / "supported-object-ledger.csv", index=False, float_format="%.17g")
    boot.to_csv(dest / "bootstrap.csv.gz", index=False, compression={"method":"gzip", "mtime":0}, float_format="%.17g")
    source = [MANIFEST, PROTO, Path(__file__), CONTROL/"score/object-scores.csv",
              inf/"manifest.json", inf/"simulation-probabilities.csv", inf/"real-probabilities.csv",
              VAL/"analysis/object-scores.csv", VAL/"fit.FITRES.TEXT",
              ROOT/"runs/research_2026_09_26/classifier_residual_sensitivity/membership_and_scores.csv"]
    (dest / "manifest.json").write_text(json.dumps({
        "inputs_sha256": {str(p.relative_to(ROOT)): digest(p) for p in source},
        "outputs_sha256": {str(p.relative_to(ROOT)): digest(p) for p in dest.iterdir() if p.name != "manifest.json"},
        "scope": "Post-primary common reconstructed-classifier selection sensitivity; no coefficient/template refit"}, indent=2) + "\n")
    return {"rows": len(result), "stage_rows": len(stage_ledger), "bootstrap_rows": len(boot),
            "max_bootstrap_failures": int(result.bootstrap_failures.max())}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("prepare_clump", "run_clump", "classify", "score"))
    parser.add_argument("arm", nargs="?", choices=("P21", "G10"))
    args = parser.parse_args()
    if args.action in ("prepare_clump", "run_clump"):
        assert args.arm
        result = prepare_clump(args.arm) if args.action == "prepare_clump" else run_clump(args.arm)
    elif args.action == "classify":
        result = classify()
    else:
        result = score()
    print(json.dumps(str(result) if isinstance(result, Path) else result, indent=2))


if __name__ == "__main__":
    main()
