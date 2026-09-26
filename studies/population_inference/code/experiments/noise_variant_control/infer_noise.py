"""Measured-peak clump and official SNN inference for frozen noise variants."""
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

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DESIGN = ROOT / "runs/research_2026_09_26/astra_design/noise_variant_design"
PROTO = ROOT / "runs/research_2026_09_26/astra_design/noise-variant-protocol.md"
HANDOFF = DESIGN / "handoff-manifest.json"
CLUMP = ROOT / "phase2/classification/reconstruction_20260926/clump_run"
COMMON = ROOT / "runs/research_2026_09_26/common_classifier_residual"
CLASSIFIER = ROOT / "scripts/research_2026_09_26/common_classifier_residual.py"
WRAPPER = ROOT / "scripts/phase2/classification/reconstruct.py"
ARMS = ("P21_rho000", "P21_rho090", "P21_noisetrue120")


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda: f.read(1048576), b""):
            h.update(part)
    return h.hexdigest()


def preflight():
    assert digest(PROTO) == "ca48b11823ac3d74f32b707421e95d698718c03111d8ffc8258a605ce8c91016"
    assert digest(HANDOFF) == "6670646abca020b65814c0d57c29e56c0226d1783a3a39ad3a7e454c1e6942d2"
    frozen = json.loads(HANDOFF.read_text())["sha256"]
    for rel in (f"runs/research_2026_09_26/astra_design/noise_variant_design/{arm}-cohort.csv" for arm in ARMS):
        assert digest(ROOT / rel) == frozen[rel]
    for rel in (f"runs/research_2026_09_26/astra_design/noise_variant_design/{arm}-cids.txt" for arm in ARMS):
        assert digest(ROOT / rel) == frozen[rel]


def module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(mod)
    return mod


def prepare(arm):
    preflight()
    assert arm in ARMS
    dest = HERE / "clump" / arm
    dest.mkdir(parents=True, exist_ok=False)
    version = f"PH2_pilot02_{arm}"
    base = (CLUMP / "clump.nml").read_text()
    old_path = f"PRIVATE_DATA_PATH = '{(ROOT/'phase2/official/inputs/SNDATA_ROOT/lcmerge/DES-SN5YR').resolve()}'"
    new_path = f"PRIVATE_DATA_PATH = '{(ROOT/'phase2/literature/simulations/outputs').resolve()}'"
    assert base.count(old_path) == 1
    base = base.replace(old_path, new_path)
    assert base.count("VERSION_PHOTOMETRY = 'DES-SN5YR_DES'") == 1
    base = base.replace("VERSION_PHOTOMETRY = 'DES-SN5YR_DES'", f"VERSION_PHOTOMETRY = '{version}'")
    assert base.count("TEXTFILE_PREFIX = 'DES-SN5YR_DES'") == 1
    base = base.replace("TEXTFILE_PREFIX = 'DES-SN5YR_DES'", f"TEXTFILE_PREFIX = '{version}'")
    assert "OPT_SETPKMJD = 16" in base and "PHOTFLAG_MSKREJ = 1016" in base
    (dest / "clump.nml").write_text(base)
    info = {"arm": arm, "version": version,
            "base_nml_sha256": digest(CLUMP / "clump.nml"),
            "nml_sha256": digest(dest / "clump.nml"),
            "binary_sha256": digest(CLUMP / "snana.exe"),
            "changed_fields": ["PRIVATE_DATA_PATH", "VERSION_PHOTOMETRY", "TEXTFILE_PREFIX"],
            "archived_fit_setpkmjd": 20, "common_clump_setpkmjd": 16,
            "truth_peak_used": False, "source_sha256": digest(Path(__file__))}
    (dest / "prepared.json").write_text(json.dumps(info, indent=2) + "\n")
    return info


def run_clump(arm):
    preflight()
    assert arm in ARMS
    path = HERE / "clump" / arm
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
    gate = {"arm": arm, "returncode": proc.returncode,
            "graceful": "ENDING PROGRAM GRACEFULLY." in (path / "clump.log").read_text(),
            "output_exists": output.exists(), "log_sha256": digest(path / "clump.log")}
    if output.exists():
        gate["output_sha256"] = digest(output)
    (path / "run-gate.json").write_text(json.dumps(gate, indent=2) + "\n")
    return gate


def infer():
    preflight()
    dest = HERE / "inference"
    dest.mkdir(exist_ok=False)
    cc = module(CLASSIFIER, "common_classifier")
    wrapped = module(WRAPPER, "official_reconstruction_wrapper")
    records, peaks_out, singles = [], [], []
    for arm in ARMS:
        clump = HERE / "clump" / arm
        gate = json.loads((clump / "run-gate.json").read_text())
        assert gate["returncode"] == 0 and gate["graceful"] and gate["output_exists"]
        ids = [s for s in (DESIGN / f"{arm}-cids.txt").read_text().split()]
        assert len(ids) == len(set(ids)) == 256
        peak = cc.read_peak(clump / f"PH2_pilot02_{arm}.SNANA.TEXT").set_index("CID").PKMJDINI
        assert set(ids) <= set(peak.index)
        q = pd.read_csv(DESIGN / f"{arm}-cohort.csv", dtype={"CID": str}).set_index("CID")
        for cid in ids:
            peaks_out.append(dict(arm=arm,CID=cid,common_peak=float(peak[cid]),
                                  archived_fit_peak=float(q.loc[cid,"PKMJDINI"]),
                                  common_minus_archived=float(peak[cid]-q.loc[cid,"PKMJDINI"])))
        stem = ROOT / f"phase2/literature/simulations/outputs/PH2_pilot02_{arm}/PH2_pilot02_{arm}"
        df, meta = wrapped.load_curves(stem, set(ids), peak.loc[ids].to_dict())
        assert len(meta) == 256 and set(meta.CID) == set(ids)
        check_ids = [ids[int(i)] for i in np.linspace(0,255,5,dtype=int)]
        alone = wrapped.infer(df.loc[df.SNID.isin(check_ids)],chunk_size=1).set_index("CID").pIa
        batched = wrapped.infer(df,chunk_size=128).set_index("CID").pIa
        delta = abs(alone-batched.loc[alone.index])
        assert delta.max() <= 1e-6
        assert ((alone>.999)==(batched.loc[alone.index]>.999)).all()
        singles.extend(dict(arm=arm,CID=cid,single=float(alone[cid]),batched=float(batched[cid]),
                            abs_delta=float(delta[cid])) for cid in alone.index)
        out = meta.merge(batched.rename("pIa"),left_on="CID",right_index=True,validate="one_to_one")
        out.insert(0,"arm",arm)
        out["gt999"] = out.pIa > .999
        records.append(out)
    allsim = pd.concat(records,ignore_index=True).sort_values(["arm","CID"])
    assert len(allsim) == 768 and not allsim.duplicated(["arm","CID"]).any()
    allsim.to_csv(dest/"variant-probabilities.csv",index=False,float_format="%.17g")
    pd.DataFrame(peaks_out).sort_values(["arm","CID"]).to_csv(dest/"peak-ledger.csv",index=False,float_format="%.17g")
    pd.DataFrame(singles).sort_values(["arm","CID"]).to_csv(dest/"single-batch-check.csv",index=False,float_format="%.17g")
    result = {"count":len(allsim),"by_arm":allsim.groupby("arm").size().astype(int).to_dict(),
              "gt999":allsim.groupby("arm").gt999.sum().astype(int).to_dict(),
              "max_single_batch_abs":float(pd.DataFrame(singles).abs_delta.max()),
              "threshold_within_1e6":allsim.loc[abs(allsim.pIa-.999)<=1e-6,["arm","CID","pIa"]].to_dict("records"),
              "raw_epochs":int(allsim.n_raw.sum()),"window_excluded":int(allsim.n_window_excluded.sum()),
              "flag_excluded":int(allsim.n_flag_excluded.sum()),"retained_epochs":int(allsim.n_retained.sum()),
              "truth_input_used":False}
    (dest/"inference-gate.json").write_text(json.dumps(result,indent=2)+"\n")
    sources = [Path(__file__),PROTO,HANDOFF,CLASSIFIER,WRAPPER,
               ROOT/"phase2/official/inputs/SNDATA_ROOT/models/classifiers/DES-SN5YR/SNNTRAINV19_z_TRAINDES_V19/model.pt"]
    sources += [HERE/"clump"/arm/f"PH2_pilot02_{arm}.SNANA.TEXT" for arm in ARMS]
    sources += [DESIGN/f"{arm}-cids.txt" for arm in ARMS]
    (dest/"manifest.json").write_text(json.dumps({
        "inputs_sha256":{str(p.relative_to(ROOT)):digest(p) for p in sources},
        "outputs_sha256":{str(p.relative_to(ROOT)):digest(p) for p in dest.iterdir() if p.name!="manifest.json"}},
        indent=2)+"\n")
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action",choices=("prepare","clump","infer"))
    ap.add_argument("arm",nargs="?",choices=ARMS)
    a=ap.parse_args()
    if a.action=="infer":
        result=infer()
    else:
        assert a.arm
        result=prepare(a.arm) if a.action=="prepare" else run_clump(a.arm)
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
