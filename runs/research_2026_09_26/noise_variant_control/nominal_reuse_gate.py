"""Close exact CID/attempt/settings provenance before reusing nominal exports."""
from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DESIGN = ROOT / "runs/research_2026_09_26/astra_design/noise_variant_design"
PILOT = ROOT / "runs/research_2026_09_26/simulation_residual_control/full256/P21/approx_minus99"
HOLD = ROOT / "runs/research_2026_09_26/simulation_holdout_control/P21"
COMMON = ROOT / "runs/research_2026_09_26/common_classifier_residual/inference/simulation-probabilities.csv"
HOLD_PROB = ROOT / "runs/research_2026_09_26/simulation_holdout_control/inference/holdout-probabilities.csv"


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1048576),b""):
            h.update(part)
    return h.hexdigest()


def normalize(nml):
    text=nml.read_text()
    text=re.sub(r"VERSION_PHOTOMETRY\s*=\s*'PH2_pilot02_[^']+'", "VERSION_PHOTOMETRY = '<version>'", text)
    text=re.sub(r"SNCID_LIST_FILE\s*=\s*'[^']+'", "SNCID_LIST_FILE = '<seed>'", text)
    text=re.sub(r"TEXTFILE_PREFIX\s*=\s*'[^']+'", "TEXTFILE_PREFIX = '<output>'", text)
    return text


def main():
    design=json.loads((DESIGN/"handoff-manifest.json").read_text())["sha256"]
    assert sha(DESIGN/"P21-cohort.csv")==design[str((DESIGN/"P21-cohort.csv").relative_to(ROOT))]
    ids=[int(x) for x in (DESIGN/"P21-cids.txt").read_text().split()]
    assert len(ids)==len(set(ids))==256
    q=pd.read_csv(DESIGN/"P21-cohort.csv").set_index("CID").loc[ids]
    refs={}
    records=[]
    max_seed_diff=0.0
    for label,path in (("pilot",PILOT),("holdout",HOLD)):
        prep=json.loads((path/"prepared.json").read_text())
        fit=json.loads((path/"fit-gate.json").read_text())
        prov=json.loads((path/"provenance-gate.json").read_text())
        tangent=json.loads((path/"tangent-gate.json").read_text())
        edge=json.loads((path/"mean-edge-gate.json").read_text())
        assert prep["law_option"]==-99 and prep["fitter_parameters_free"] and not prep["truth_coordinates_used"]
        binary=sha(ROOT/"phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe")
        assert (prep.get("binary_sha256") or prep["inputs_sha256"]["phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe"])==binary
        assert fit["returncode"]==0 and fit["graceful"] and prov["full_pilot_adequate"]
        assert tangent["all_rank4"] and edge["eligible_objects"]==prov["exported"]
        assert prov["max_raw_mjd_difference"]==prov["max_raw_flux_difference"]==prov["max_raw_error_difference"]==0
        assert sha(path/"provenance-ledger.csv")==prov["provenance_ledger_sha256"]
        assert sha(path/"mean-edge-ledger.csv")==edge["ledger_sha256"]
        assert sha(path/"objectives/manifest.json")
        source_ids=set(prep["selected_ids"])
        relevant=set(ids)&source_ids
        assert not (set(refs)&relevant)
        cohort_path=(ROOT/"runs/research_2026_09_26/astra_design/simulation_design/P21-cohort.csv") if label=="pilot" else (ROOT/"runs/research_2026_09_26/astra_design/simulation_holdout_design/P21-cohort.csv")
        old=pd.read_csv(cohort_path).set_index("CID").loc[sorted(relevant)]
        cols=["t0_double","mB_double","x1_double","c_double","PKMJDINI","zHEL","SIM_LIBID"]
        diffs={k:float(np.max(abs(old[k].to_numpy(float)-q.loc[sorted(relevant),k].to_numpy(float)))) for k in cols}
        assert all(v<1e-9 for v in diffs.values()),(label,diffs)
        max_seed_diff=max(max_seed_diff,*(diffs[k] for k in ("t0_double","mB_double","x1_double","c_double")))
        assert set(relevant)<=set(fit["fitres_cids"])&set(fit["objective_cids"])
        assert set(relevant)<=set(pd.read_csv(path/"provenance-ledger.csv").CID.astype(int))
        for cid in sorted(relevant):
            refs[cid]=(label,path)
            records.append(dict(CID=cid,generated_attempt_index=int(q.loc[cid,"generated_attempt_index"]),
                                original_scope=label,objective_path=str((path/"objectives"/f"objective_{cid}.npz").relative_to(ROOT)),
                                original_fitres_path=str((path/"fit.FITRES.TEXT").relative_to(ROOT)),
                                original_objective_sha256=sha(path/"objectives"/f"objective_{cid}.npz")))
    assert set(refs)==set(ids)
    # The source NMLs must be byte-identical after replacing the three
    # intentional path/version parameters; this includes clipping and priors.
    variants={arm:HERE/arm for arm in ("P21_rho000","P21_rho090","P21_noisetrue120")}
    expected=normalize(PILOT/"fit.nml")
    assert normalize(HOLD/"fit.nml")==expected
    for arm,path in variants.items():
        assert normalize(path/"fit.nml")==expected,arm
        p=json.loads((path/"prepared.json").read_text())
        assert p["binary_sha256"]==sha(ROOT/"phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe")
        assert p["base_nml_sha256"]==sha(ROOT/"phase2/checkpoint/official-inputs/snana_forward_p21.nml")
        other=pd.read_csv(DESIGN/f"{arm}-cohort.csv").set_index("CID").loc[ids]
        assert np.array_equal(other.generated_attempt_index.to_numpy(int),q.generated_attempt_index.to_numpy(int))
        assert np.array_equal(other.CID if "CID" in other else other.index.to_numpy(),q.index.to_numpy())
    pilot_prob=pd.read_csv(COMMON)
    hold_prob=pd.read_csv(HOLD_PROB)
    pilot_ids=set(pilot_prob.loc[pilot_prob.arm=="P21","CID"].astype(int))
    hold_ids=set(hold_prob.loc[hold_prob.arm=="P21","CID"].astype(int))
    assert {cid for cid,(label,_) in refs.items() if label=="pilot"}<=pilot_ids
    assert {cid for cid,(label,_) in refs.items() if label=="holdout"}<=hold_ids
    frame=pd.DataFrame(records).sort_values("CID")
    frame.to_csv(HERE/"nominal-reuse-ledger.csv",index=False)
    result={"nominal_ids":256,"pilot_reused":sum(x=="pilot" for x,_ in refs.values()),
            "holdout_reused":sum(x=="holdout" for x,_ in refs.values()),
            "all_successful_native_exports":True,"exact_normalized_fit_nml":True,
            "max_measurement_seed_difference":max_seed_diff,
            "nominal_ledger_sha256":sha(HERE/"nominal-reuse-ledger.csv"),
            "source_sha256":sha(Path(__file__)),
            "source_NML_sha256":sha(PILOT/"fit.nml"),
            "source_cohort_sha256":{label:sha(path) for label,path in (("pilot",ROOT/"runs/research_2026_09_26/astra_design/simulation_design/P21-cohort.csv"),("holdout",ROOT/"runs/research_2026_09_26/astra_design/simulation_holdout_design/P21-cohort.csv"))},
            "classifier_probabilities_sha256":{"pilot":sha(COMMON),"holdout":sha(HOLD_PROB)}}
    (HERE/"nominal-reuse-gate.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
