"""Independent arithmetic, membership and hash checks for holdout results."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def digest(p):
    h = hashlib.sha256()
    with p.open("rb") as stream:
        for block in iter(lambda:stream.read(1048576),b""):
            h.update(block)
    return h.hexdigest()


def check_manifest(path):
    m=json.loads(path.read_text())
    for rel,expected in {**m["inputs_sha256"],**m["outputs_sha256"]}.items():
        assert digest(ROOT/rel)==expected,rel
    return len(m["inputs_sha256"])+len(m["outputs_sha256"])


def main():
    manifests=[HERE/"object_scores/manifest.json",HERE/"transport/manifest.json"]
    count=sum(map(check_manifest,manifests))
    scores=pd.read_csv(HERE/"object_scores/holdout-object-scores.csv")
    with np.load(HERE/"object_scores/holdout-sufficient-arrays.npz",allow_pickle=False) as arr:
        assert np.array_equal(arr["CID"],scores.CID)
        assert np.array_equal(arr["arm"],scores.arm)
        c=np.asarray(json.loads(manifests[0].read_text())["frozen_coefficients"],float)
        a=arr["u"]@c;I=np.einsum("i,nij,j->n",c,arr["F"],c)
    assert np.max(abs(a-scores.matched_filter))<1e-12
    assert np.max(abs(I-scores.information))<1e-12
    assert np.max(abs(a-I/2-scores.fixed_gain))<1e-12
    assert not scores.duplicated(["arm","CID"]).any()
    pilot=pd.read_csv(ROOT/"runs/research_2026_09_26/simulation_residual_control/score/object-scores.csv")
    assert not (set(scores.loc[scores.arm=="P21","CID"]) & set(pilot.loc[pilot.arm=="P21","CID"]))
    assert not (set(scores.loc[scores.arm=="G10","CID"]) & set(pilot.loc[pilot.arm=="G10","CID"]))
    probs=pd.read_csv(HERE/"inference/holdout-probabilities.csv")[["arm","CID","pIa","gt999"]]
    joint=scores.merge(probs,on=["arm","CID"],validate="one_to_one",suffixes=("_score","_inference"))
    assert len(joint)==len(scores)
    assert np.max(abs(joint.pIa_score-joint.pIa_inference))<1e-12
    assert (joint.gt999_score==joint.gt999_inference).all()
    trans=pd.read_csv(HERE/"transport/transport-summary.csv")
    ledger=pd.read_csv(HERE/"transport/supported-object-ledger.csv")
    fields=pd.read_csv(HERE/"transport/field-contributions.csv")
    boot=pd.read_csv(HERE/"transport/bootstrap.csv.gz")
    keys=("scope","arm","transport","stage","support_source")
    for row in trans.itertuples():
        mask=np.logical_and.reduce([ledger[k]==getattr(row,k) for k in keys])
        weighted=ledger.loc[mask]
        assert len(weighted)==row.supported_sim and weighted.CID.is_unique
        source=scores.copy() if row.scope=="holdout" else pd.concat([scores,pilot.loc[pilot.law=="approx_minus99"]])
        q=source.loc[source.arm==row.arm].set_index("CID").loc[weighted.CID]
        w=weighted.base_weight.to_numpy(float)
        assert abs(w.sum()-row.supported_real)<1e-8
        for name,column in (("a","matched_filter"),("I","information"),("G","fixed_gain"),
                            ("chi2","projected_chi2"),("dimension","projected_dimension")):
            assert abs(w@q[column].to_numpy(float)-getattr(row,"weighted_"+name))<1e-9
        assert abs(row.weighted_a/row.weighted_I-row.transported_amplitude)<1e-12
        assert abs(row.weighted_chi2/row.weighted_dimension-row.pooled_Q)<1e-12
        fm=np.logical_and.reduce([fields[k]==getattr(row,k) for k in keys]);f=fields.loc[fm]
        assert f.n.sum()==row.supported_sim
        assert abs(f.weighted_a.sum()-row.weighted_a)<1e-9
        bm=np.logical_and.reduce([boot[k]==getattr(row,k) for k in keys]);b=boot.loc[bm]
        assert len(b)==1000 and (~b.valid).sum()==row.bootstrap_failures
        v=b.loc[b.valid]
        assert abs(v.amplitude.quantile(.025)-row.bootstrap_amplitude_ci_low)<1e-10
        assert abs(v.amplitude.quantile(.975)-row.bootstrap_amplitude_ci_high)<1e-10
    report={"manifest_hashes_verified":count,"object_scores":len(scores),"transport_rows":len(trans),
            "bootstrap_rows":len(boot),"max_score_arithmetic_error":float(max(np.max(abs(a-scores.matched_filter)),np.max(abs(I-scores.information)))),
            "status":"passed"}
    (HERE/"independent-verification.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
