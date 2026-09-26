"""Independent hash, sufficient-array, selection and transport verification."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ARMS=("P21","P21_rho000","P21_rho090","P21_noisetrue120")


def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for part in iter(lambda:f.read(1048576),b""):
            h.update(part)
    return h.hexdigest()


def verify_manifest(path):
    m=json.loads(path.read_text())
    for rel,wanted in {**m["inputs_sha256"],**m["outputs_sha256"]}.items():
        assert sha(ROOT/rel)==wanted,rel
    return len(m["inputs_sha256"])+len(m["outputs_sha256"])


def main():
    assert sha(HERE/"score-plan.md")=="92cc736b4a95033111655234c63328456320e77f885438b5ec17e0dea75268f0"
    manifests=[HERE/"inference/manifest.json",HERE/"object_scores/manifest.json",
               HERE/"transport/manifest.json",HERE/"common_mask/manifest.json"]
    count=sum(map(verify_manifest,manifests))
    assert sha(HERE/"inference/manifest.json")
    for source in ("native-runner.py","infer_noise.py","score_noise.py","transport_noise.py","common_mask.py","nominal_reuse_gate.py"):
        original=ROOT/"scripts/research_2026_09_26/noise_variant_control.py" if source=="native-runner.py" else HERE/source
        if source=="nominal_reuse_gate.py": original=HERE/source
        assert sha(original)==sha(HERE/"source_snapshots"/source),source
    score=pd.read_csv(HERE/"object_scores/object-scores.csv")
    select=pd.read_csv(HERE/"object_scores/selection-ledger.csv")
    assert len(score)==1023 and len(select)==1024
    assert not score.duplicated(["arm","CID"]).any()
    assert not select.duplicated(["arm","CID"]).any()
    assert set(score.arm)==set(ARMS)
    with np.load(HERE/"object_scores/sufficient-arrays.npz",allow_pickle=False) as arr:
        assert np.array_equal(arr["arm"],score.arm)
        assert np.array_equal(arr["CID"],score.CID)
        c=np.array(json.loads((HERE/"object_scores/score-gate.json").read_text())["frozen_coefficients"])
        M=arr["u"]@c;I=np.einsum("i,nij,j->n",c,arr["F"],c)
    max_arithmetic=float(max(np.max(abs(M-score.matched_filter)),np.max(abs(I-score.information)),
                             np.max(abs(M-I/2-score.fixed_gain))))
    assert max_arithmetic<1e-12
    missing=select.loc[~select.native_success,["arm","CID"]].to_dict("records")
    assert missing==[{"arm":"P21_noisetrue120","CID":10516}]
    pred=pd.read_csv(HERE/"inference/variant-probabilities.csv")
    for arm in ARMS[1:]:
        s=select.loc[select.arm==arm].set_index("CID")
        p=pred.loc[pred.arm==arm].set_index("CID")
        assert set(s.index)==set(p.index)
        p=p.loc[s.index]
        assert np.max(abs(s.pIa-p.pIa))<1e-12
        assert (s.gt999==p.gt999).all()
    summary=pd.read_csv(HERE/"transport/transport-summary.csv")
    support=pd.read_csv(HERE/"transport/supported-object-ledger.csv")
    fields=pd.read_csv(HERE/"transport/field-contributions.csv")
    boot=pd.read_csv(HERE/"transport/bootstrap.csv.gz")
    assert len(summary)==40 and len(boot)==40000
    keys=("stage","transport","arm")
    max_weighted=0.0
    for row in summary.itertuples():
        x=support.loc[(support.stage==row.stage)&(support.transport==row.transport)&(support.arm==row.arm)]
        assert len(x)==row.supported_sim and x.CID.is_unique
        q=score.loc[score.arm==row.arm].set_index("CID").loc[x.CID]
        w=x.base_weight.to_numpy(float)
        assert abs(w.sum()-row.supported_real)<1e-8
        for field,column in (("M","matched_filter"),("I","information"),("G","fixed_gain"),
                             ("Q","projected_chi2"),("dimension","projected_dimension")):
            diff=abs(w@q[column].to_numpy(float)-getattr(row,"weighted_"+field))
            max_weighted=max(max_weighted,diff)
            assert diff<1e-9,(row.stage,row.arm,field,diff)
        assert abs(row.weighted_M/row.weighted_I-row.amplitude)<1e-12
        assert abs(row.weighted_Q/row.weighted_dimension-row.pooled_Q)<1e-12
        assert abs(row.weighted_G/w.sum()-row.mean_gain)<1e-12
        f=fields.loc[(fields.stage==row.stage)&(fields.transport==row.transport)&(fields.arm==row.arm)]
        assert f.n.sum()==len(x)
        assert abs(f.weighted_M.sum()-row.weighted_M)<1e-9
        b=boot.loc[(boot.stage==row.stage)&(boot.transport==row.transport)&(boot.arm==row.arm)]
        assert len(b)==1000 and (~b.valid).sum()==row.bootstrap_failures
        for col in ("amplitude","mean_gain","pooled_Q","difference_amplitude","difference_mean_gain","difference_pooled_Q"):
            v=b.loc[b.valid if not col.startswith("difference") else b.paired_valid,col]
            assert abs(v.quantile(.025)-getattr(row,f"bootstrap_{col}_low"))<1e-10
            assert abs(v.quantile(.975)-getattr(row,f"bootstrap_{col}_high"))<1e-10
    for (stage,transport),q in summary.groupby(["stage","transport"]):
        nom=q.loc[q.arm=="P21"].iloc[0]
        assert (q.supported_cells==nom.supported_cells).all()
        assert (q.supported_real==nom.supported_real).all()
        for arm in ARMS[1:]:
            z=q.loc[q.arm==arm].iloc[0]
            b=boot.loc[(boot.stage==stage)&(boot.transport==transport)&(boot.arm==arm)]
            bn=boot.loc[(boot.stage==stage)&(boot.transport==transport)&(boot.arm=="P21")]
            assert np.array_equal(b.replicate,bn.replicate)
            assert np.max(abs(b.difference_amplitude-(b.amplitude.to_numpy()-bn.amplitude.to_numpy())))<1e-12
            if stage.startswith("common_native") or stage.startswith("all_four"):
                w1=support.loc[(support.stage==stage)&(support.transport==transport)&(support.arm=="P21")].sort_values("CID")
                w2=support.loc[(support.stage==stage)&(support.transport==transport)&(support.arm==arm)].sort_values("CID")
                assert np.array_equal(w1.CID,w2.CID)
                assert np.max(abs(w1.base_weight.to_numpy()-w2.base_weight.to_numpy()))<1e-12
    mask=pd.read_csv(HERE/"common_mask/common-mask-object-scores.csv")
    mg=json.loads((HERE/"common_mask/gate.json").read_text())
    assert len(mask)==mg["scored_rows"]
    if mg["failed_object_records"]==0:
        assert len(mask)==4*255
    report={"status":"passed","verified_manifest_hashes":count,
            "object_scores":len(score),"transport_rows":len(summary),"bootstrap_rows":len(boot),
            "max_sufficient_arithmetic_error":max_arithmetic,"max_transport_weighted_error":max_weighted,
            "mask_scored_rows":len(mask),"missing_fitres":missing}
    (HERE/"independent-verification.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
