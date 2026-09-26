"""Execute the frozen disjoint P21/G10 same-catalogue holdout control."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "runs/research_2026_09_26/astra_design/simulation_holdout_design"
PROTO = ROOT / "runs/research_2026_09_26/astra_design/simulation-holdout-protocol.md"
OUT = ROOT / "runs/research_2026_09_26/simulation_holdout_control"
NATIVE_DESIGN = OUT / "native_design"
BASE_SOURCE = ROOT / "scripts/research_2026_09_26/simulation_residual_control.py"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for part in iter(lambda: stream.read(1048576), b""):
            h.update(part)
    return h.hexdigest()


def preflight():
    frozen = json.loads((DESIGN / "manifest.json").read_text())
    assert digest(PROTO) == "8f5cd945582f84486c380fc1d406c53fc4305f1278177fd30a8fa23b98d8a2a1"
    for rel, expected in frozen["sha256"].items():
        assert digest(ROOT / rel) == expected, rel
    common = ROOT / "runs/research_2026_09_26/common_classifier_residual"
    assert json.loads((common / "score-independent-verification.json").read_text())["result"] == "passed"
    assert (common / "inference-rerun-identity.json").exists()
    return frozen


def base_module():
    spec = importlib.util.spec_from_file_location("frozen_pilot_native_runner", BASE_SOURCE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    module.DESIGN = NATIVE_DESIGN
    module.OUT = OUT
    module.preflight = preflight
    return module


def prepare(arm):
    frozen = preflight()
    assert arm in ("P21", "G10")
    path = OUT / arm
    path.mkdir(parents=True, exist_ok=False)
    ids = [int(x) for x in (DESIGN / f"{arm}-cids.txt").read_text().split()]
    expected = frozen["counts"][arm]["holdout"]
    assert len(ids) == expected and len(set(ids)) == expected
    pilots = set(map(int, (ROOT / f"runs/research_2026_09_26/astra_design/simulation_design/{arm}-cids.txt").read_text().split()))
    assert not (set(ids) & pilots)
    cohort = pd.read_csv(DESIGN / f"{arm}-cohort.csv").set_index("CID").loc[ids]
    assert np.isfinite(cohort[["t0_double", "mB_double", "x1_double", "c_double"]].to_numpy(float)).all()
    assert np.all(np.isfinite(cohort.x0) & (cohort.x0 > 0))
    x0 = 10 ** ((10.635 - cohort.mB_double.to_numpy(float)) / 2.5)
    seed = path / "measured_initial_values.FITRES"
    seed.write_text("VARNAMES: CID PKMJD x0 x1 c\n" + "".join(
        f"SN: {cid} {t0:.17g} {amp:.17g} {x1:.17g} {color:.17g}\n"
        for cid, t0, amp, x1, color in zip(ids, cohort.t0_double, x0, cohort.x1_double, cohort.c_double)))
    base = ROOT / f"phase2/checkpoint/official-inputs/snana_forward_{arm.lower()}.nml"
    original = base.read_text()
    assert "OPT_MWCOLORLAW = 99" in original and "LFIXPAR_ALL" not in original
    assert "OPT_SNCID_LIST" not in original and "SNCID_LIST_FILE" not in original
    assert "SNTABLE_LIST      = 'FITRES(text:host)'" in original
    text = original.replace("&SNLCINP", "&SNLCINP\n"
                            "    MXLC_PLOT = 1000\n"
                            "    OPT_SNCID_LIST = 2\n"
                            f"    SNCID_LIST_FILE = '{seed.resolve()}'", 1)
    text = text.replace("SNTABLE_LIST      = 'FITRES(text:host)'",
                        "SNTABLE_LIST      = 'FITRES(text:host) LCPLOT(text:col)'", 1)
    text = re.sub(r"(?m)^    TEXTFILE_PREFIX\s*=.*$",
                  f"    TEXTFILE_PREFIX = '{(path / 'fit').resolve()}'", text, count=1)
    text = text.replace("OPT_MWCOLORLAW = 99", "OPT_MWCOLORLAW = -99", 1)
    nml = path / "fit.nml"
    nml.write_text(text)
    prepared = {"arm": arm, "law": "approx_minus99", "law_option": -99,
                "cohort": "holdout", "selected_ids": ids, "selected_count": len(ids),
                "opt_sncid_list": 2, "fitter_parameters_free": True,
                "truth_coordinates_used": False, "seed_sha256": digest(seed),
                "nml_sha256": digest(nml), "binary_sha256": digest(ROOT / "phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe"),
                "base_runner_sha256": digest(BASE_SOURCE), "prepared_source_sha256": digest(Path(__file__)),
                "source_manifest_sha256": digest(DESIGN / "manifest.json"),
                "base_nml_sha256": digest(base)}
    (path / "prepared.json").write_text(json.dumps(prepared, indent=2)+"\n")
    return {"arm": arm, "count": len(ids), "nml_sha256": digest(nml)}


def action(name, arm):
    preflight()
    path = OUT / arm
    native = base_module()
    return {"fit": native.run_fit, "export": native.export,
            "gate": native.gate, "tangent": native.tangent_gate,
            "mean_edges": native.mean_edges}[name](path)


def infer_holdout():
    preflight()
    dest = OUT / "inference"
    dest.mkdir(exist_ok=False)
    common = ROOT / "runs/research_2026_09_26/common_classifier_residual"
    common_source = ROOT / "scripts/research_2026_09_26/common_classifier_residual.py"
    spec = importlib.util.spec_from_file_location("closed_common_classifier", common_source)
    cc = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(cc)
    wrapped = cc.wrapper()
    records = []
    peaks_out = []
    single_rows = []
    for arm in ("P21", "G10"):
        ids = [s for s in (DESIGN / f"{arm}-cids.txt").read_text().split()]
        clump = common / f"clump/{arm}/PH2_pilot02_{arm}.SNANA.TEXT"
        peaks = cc.read_peak(clump).set_index("CID").PKMJDINI
        assert set(ids) <= set(peaks.index)
        cohort = pd.read_csv(DESIGN / f"{arm}-cohort.csv", dtype={"CID":str}).set_index("CID")
        for cid in ids:
            peaks_out.append(dict(arm=arm, CID=cid, common_peak=float(peaks[cid]),
                                  archived_fit_peak=float(cohort.loc[cid,"PKMJDINI"]),
                                  common_minus_archived=float(peaks[cid]-cohort.loc[cid,"PKMJDINI"])))
        stem = ROOT / f"phase2/literature/simulations/outputs/PH2_pilot02_{arm}/PH2_pilot02_{arm}"
        df, meta = wrapped.load_curves(stem, set(ids), peaks.loc[ids].to_dict())
        assert len(meta) == len(ids) and set(meta.CID) == set(ids)
        check_ids = [ids[int(i)] for i in np.linspace(0,len(ids)-1,5,dtype=int)]
        single = wrapped.infer(df.loc[df.SNID.isin(check_ids)], chunk_size=1).set_index("CID").pIa
        pred = wrapped.infer(df, chunk_size=128).set_index("CID").pIa
        delta = abs(single-pred.loc[single.index])
        assert delta.max() <= 1e-6
        assert ((single > .999)==(pred.loc[single.index] > .999)).all()
        single_rows.extend(dict(arm=arm,CID=cid,single=float(single[cid]),batched=float(pred[cid]),
                                abs_delta=float(delta[cid])) for cid in single.index)
        q = meta.merge(pred.rename("pIa"),left_on="CID",right_index=True,validate="one_to_one")
        q.insert(0,"arm",arm)
        q["gt999"] = q.pIa > .999
        records.append(q)
    probs = pd.concat(records,ignore_index=True).sort_values(["arm","CID"])
    assert len(probs) == 1569+1908
    probs.to_csv(dest/"holdout-probabilities.csv",index=False,float_format="%.17g")
    pd.DataFrame(peaks_out).sort_values(["arm","CID"]).to_csv(dest/"peak-ledger.csv",index=False,float_format="%.17g")
    pd.DataFrame(single_rows).sort_values(["arm","CID"]).to_csv(dest/"single-batch-check.csv",index=False,float_format="%.17g")
    gate = {"count":len(probs),"by_arm":probs.groupby("arm").size().astype(int).to_dict(),
            "gt999_by_arm":probs.groupby("arm").gt999.sum().astype(int).to_dict(),
            "max_single_batch_abs":float(pd.DataFrame(single_rows).abs_delta.max()),
            "threshold_within_1e6":probs.loc[abs(probs.pIa-.999)<=1e-6,["arm","CID","pIa"]].to_dict("records"),
            "raw_epochs":int(probs.n_raw.sum()),
            "window_excluded":int(probs.n_window_excluded.sum()),
            "flag_excluded":int(probs.n_flag_excluded.sum()),
            "retained_epochs":int(probs.n_retained.sum()),
            "truth_input_used":False}
    (dest/"inference-gate.json").write_text(json.dumps(gate,indent=2)+"\n")
    sources = [DESIGN/"manifest.json",PROTO,Path(__file__),common_source,
               common/"final_source_rerun/inference/manifest.json",
               common/"inference-rerun-identity.json"]
    sources += [common/f"clump/{arm}/PH2_pilot02_{arm}.SNANA.TEXT" for arm in ("P21","G10")]
    sources += [DESIGN/f"{arm}-cids.txt" for arm in ("P21","G10")]
    (dest/"manifest.json").write_text(json.dumps({
        "inputs_sha256":{str(p.relative_to(ROOT)):digest(p) for p in sources},
        "outputs_sha256":{str(p.relative_to(ROOT)):digest(p) for p in dest.iterdir() if p.name!="manifest.json"}},
        indent=2)+"\n")
    return gate


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=("prepare", "fit", "export", "gate", "tangent", "mean_edges", "infer"))
    ap.add_argument("arm", nargs="?", choices=("P21", "G10"))
    a = ap.parse_args()
    if a.action == "infer":
        result = infer_holdout()
    else:
        assert a.arm
        result = prepare(a.arm) if a.action == "prepare" else action(a.action, a.arm)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
