"""Prespecified same-catalogue holdout/combined transport and block resampling."""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PILOT = ROOT / "runs/research_2026_09_26/simulation_residual_control"
COMMON = ROOT / "runs/research_2026_09_26/common_classifier_residual"
VAL = ROOT / "runs/research_2026_09_26/astra_design/validation1020"


def digest(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def main():
    assert digest(HERE/"score-plan.md") == "ee108850f9d64a3689f7f9c50e8667ee6d89cecb56547b1542b1893ee104288f"
    dest = HERE / "transport"
    dest.mkdir(exist_ok=False)
    hold = pd.read_csv(HERE / "object_scores/holdout-object-scores.csv")
    pilot = pd.read_csv(PILOT / "score/object-scores.csv")
    pilot = pilot.loc[pilot.law == "approx_minus99"].copy()
    pprob = pd.read_csv(COMMON / "inference/simulation-probabilities.csv")
    pilot = pilot.merge(pprob[["arm","CID","pIa","gt999"]],on=["arm","CID"],validate="one_to_one")
    assert len(pilot) == 510
    assert len(hold) >= .95*3477
    hold["cohort_scope"] = "holdout"
    pilot["cohort_scope"] = "pilot"
    combined = pd.concat([hold,pilot],ignore_index=True)
    assert not combined.duplicated(["arm","CID"]).any()
    assert combined.groupby("arm").CID.nunique().to_dict() == {
        "P21":len(hold.loc[hold.arm=="P21"])+255,
        "G10":len(hold.loc[hold.arm=="G10"])+255}
    combined["snrbin"] = np.where(combined.SNRMAX1_archived < 15,"lt15","ge15")
    realprob = pd.read_csv(COMMON / "inference/real-probabilities.csv")
    real_ids = set(realprob.loc[realprob.gt999,"CID"])
    assert len(real_ids) == 1006
    real = pd.read_csv(VAL / "analysis/object-scores.csv")
    real = real.loc[(real.arm == "published_mask") & real.CID.isin(real_ids)].copy()
    assert len(real) == 1006
    real["fixed_gain"] = real.fixed_prediction_gain
    edges = [.05,.2,.35,.5,.65,.8,1.2]
    real["zbin"] = pd.cut(real.zHEL,edges,include_lowest=True,labels=False)
    assert real.zbin.notna().all()
    real["zbin"] = real.zbin.astype(int)
    real["cell"] = real.field + "_" + real.zbin.astype(str)
    sys.path.insert(0,str(ROOT/"scripts/phase2/official"))
    from audit_fits import read_fit
    realfit = read_fit(VAL/"fit.FITRES.TEXT")[["CID","SNRMAX1"]]
    realfit.CID = realfit.CID.astype(int)
    real = real.merge(realfit,on="CID",validate="one_to_one")
    real["snrbin"] = np.where(real.SNRMAX1<15,"lt15","ge15")
    libs = np.array(sorted(set(combined.LIBID.astype(int))))
    rng = np.random.default_rng(26092691)
    draws = rng.integers(0,len(libs),size=(1000,len(libs)))
    multi = np.array([np.bincount(row,minlength=len(libs)) for row in draws],dtype=np.int16)
    lib_idx = {lib:i for i,lib in enumerate(libs)}
    summary=[];fields=[];supports=[];boots=[];stages=[]
    for scope in ("holdout","combined"):
        source = combined.loc[combined.cohort_scope == "holdout"].copy() if scope=="holdout" else combined.copy()
        for arm in ("P21","G10"):
            q = source.loc[source.arm==arm]
            for stage, sub in (
                ("all_success",q),
                ("archived_quality",q.loc[q.archived_basic_quality]),
                ("quality_classifier",q.loc[q.archived_basic_quality & q.gt999]),
                ("quality_classifier_p21_nonnegative",q.loc[q.archived_basic_quality & q.gt999 & (q.support_class=="nonnegative_declared_grid")]
                 if arm=="P21" else q.loc[q.archived_basic_quality & q.gt999])):
                stages.append(dict(scope=scope,arm=arm,stage=stage,selected=len(sub),
                                   raw_a=float(sub.matched_filter.sum()),raw_I=float(sub.information.sum()),
                                   raw_G=float(sub.fixed_gain.sum()),
                                   raw_chi2=float(sub.projected_chi2.sum()),
                                   raw_dimension=int(sub.projected_dimension.sum()),
                                   CIDs="|".join(map(str,sorted(sub.CID)))))
        for split in ("field_z","field_z_snr15"):
            r = real.copy()
            r["transport_cell"] = r.cell if split=="field_z" else r.cell+"_"+r.snrbin
            source["transport_cell"] = source.cell if split=="field_z" else source.cell+"_"+source.snrbin
            subsets={}
            for arm in ("P21","G10"):
                base = source.loc[source.arm==arm]
                subsets[(arm,"all_success")]=base
                subsets[(arm,"archived_quality")]=base.loc[base.archived_basic_quality]
                subsets[(arm,"quality_classifier")]=base.loc[base.archived_basic_quality & base.gt999]
                subsets[(arm,"quality_classifier_p21_nonnegative")]=(
                    base.loc[base.archived_basic_quality & base.gt999 & (base.support_class=="nonnegative_declared_grid")]
                    if arm=="P21" else base.loc[base.archived_basic_quality & base.gt999])
            for support_source in ("own_stage","main_stage_common"):
                for stage in ("all_success","archived_quality","quality_classifier","quality_classifier_p21_nonnegative"):
                    if support_source=="main_stage_common" and stage=="quality_classifier_p21_nonnegative":
                        continue
                    support_stage = stage if support_source=="own_stage" else "quality_classifier"
                    p = subsets[("P21",support_stage)].transport_cell.value_counts()
                    g = subsets[("G10",support_stage)].transport_cell.value_counts()
                    cells = sorted(set(p[p>=2].index)&set(g[g>=2].index)&set(r.transport_cell))
                    target = r.loc[r.transport_cell.isin(cells)]
                    assert len(target)>0 and target.information.sum()>0
                    for arm in ("P21","G10"):
                        full = subsets[(arm,stage)]
                        q = full.loc[full.transport_cell.isin(cells)].copy()
                        assert len(q)>0
                        np_cell = q.transport_cell.value_counts()
                        nr_cell = target.transport_cell.value_counts()
                        w = q.transport_cell.map(nr_cell).to_numpy(float)/q.transport_cell.map(np_cell).to_numpy(float)
                        assert abs(w.sum()-len(target))<1e-8
                        a=q.matched_filter.to_numpy(float);I=q.information.to_numpy(float);G=q.fixed_gain.to_numpy(float)
                        chi=q.projected_chi2.to_numpy(float);nu=q.projected_dimension.to_numpy(float)
                        suffix=dict(scope=scope,arm=arm,transport=split,stage=stage,support_source=support_source)
                        summary.append(dict(**suffix,selected=len(full),supported_sim=len(q),supported_cells=len(cells),
                            supported_real=len(target),real_fraction=len(target)/1006,
                            raw_a=float(full.matched_filter.sum()),raw_I=float(full.information.sum()),
                            raw_G=float(full.fixed_gain.sum()),weighted_a=float(w@a),weighted_I=float(w@I),
                            weighted_G=float(w@G),transported_amplitude=float((w@a)/(w@I)),
                            transported_mean_gain=float((w@G)/w.sum()),
                            weighted_chi2=float(w@chi),weighted_dimension=float(w@nu),
                            pooled_Q=float((w@chi)/(w@nu)),
                            effective_sample_size=float(w.sum()**2/(w@w)),
                            chi2_nu_q25=float(np.quantile(chi/nu,.25)),
                            chi2_nu_q50=float(np.quantile(chi/nu,.5)),
                            chi2_nu_q75=float(np.quantile(chi/nu,.75)),
                            real_a=float(target.matched_filter.sum()),real_I=float(target.information.sum()),
                            real_G=float(target.fixed_gain.sum()),
                            real_amplitude=float(target.matched_filter.sum()/target.information.sum()),
                            real_mean_gain=float(target.fixed_gain.mean()),
                            real_pooled_Q=float(target.projected_chi2.sum()/target.projected_dimension.sum()),
                            cell_ids="|".join(cells)))
                        supports.extend(dict(**suffix,CID=int(x.CID),cell=x.transport_cell,base_weight=float(ww))
                                        for x,ww in zip(q.itertuples(),w))
                        for field,fq in q.groupby("field"):
                            fw=w[np.asarray(q.field==field)]
                            fields.append(dict(**suffix,field=field,n=len(fq),
                                weighted_a=float(fw@fq.matched_filter.to_numpy(float)),
                                weighted_I=float(fw@fq.information.to_numpy(float)),
                                weighted_G=float(fw@fq.fixed_gain.to_numpy(float)),
                                weighted_chi2=float(fw@fq.projected_chi2.to_numpy(float)),
                                weighted_dimension=float(fw@fq.projected_dimension.to_numpy(float))))
                        bm = multi[:,[lib_idx[int(x)] for x in q.LIBID]]*w[None,:]
                        bi=bm@I;bw=bm.sum(axis=1)
                        good=(bi>0)&(bw>0)
                        for rep in range(1000):
                            boots.append(dict(**suffix,replicate=rep,valid=bool(good[rep]),
                                amplitude=float((bm[rep]@a)/bi[rep]) if good[rep] else np.nan,
                                mean_gain=float((bm[rep]@G)/bw[rep]) if good[rep] else np.nan,
                                pooled_Q=float((bm[rep]@chi)/(bm[rep]@nu)) if good[rep] else np.nan))
    frame=pd.DataFrame(summary)
    boot=pd.DataFrame(boots)
    keys=("scope","arm","transport","stage","support_source")
    for key,group in boot.groupby(list(keys)):
        mask=np.logical_and.reduce([frame[col]==value for col,value in zip(keys,key)])
        valid=group.loc[group.valid]
        frame.loc[mask,"bootstrap_failures"]=len(group)-len(valid)
        for col in ("amplitude","mean_gain","pooled_Q"):
            frame.loc[mask,f"bootstrap_{col}_ci_low"]=valid[col].quantile(.025) if len(valid) else np.nan
            frame.loc[mask,f"bootstrap_{col}_ci_high"]=valid[col].quantile(.975) if len(valid) else np.nan
    frame.to_csv(dest/"transport-summary.csv",index=False,float_format="%.17g")
    pd.DataFrame(stages).to_csv(dest/"stage-ledger.csv",index=False,float_format="%.17g")
    pd.DataFrame(fields).to_csv(dest/"field-contributions.csv",index=False,float_format="%.17g")
    pd.DataFrame(supports).to_csv(dest/"supported-object-ledger.csv",index=False,float_format="%.17g")
    boot.to_csv(dest/"bootstrap.csv.gz",index=False,compression={"method":"gzip","mtime":0},float_format="%.17g")
    sources=[Path(__file__),HERE/"score-plan.md",HERE/"object_scores/manifest.json",
             HERE/"object_scores/holdout-object-scores.csv",PILOT/"score/object-scores.csv",
             COMMON/"inference/simulation-probabilities.csv",COMMON/"inference/real-probabilities.csv",
             VAL/"analysis/object-scores.csv",VAL/"fit.FITRES.TEXT"]
    (dest/"manifest.json").write_text(json.dumps({
       "inputs_sha256":{str(p.relative_to(ROOT)):digest(p) for p in sources},
       "outputs_sha256":{str(p.relative_to(ROOT)):digest(p) for p in dest.iterdir() if p.name!="manifest.json"},
       "bootstrap_seed":26092691,"bootstrap_replicates":1000,
       "scope":"same-catalogue disjoint holdout first; combined labelled separately; no independent-seed null"},
       indent=2)+"\n")
    print(json.dumps({"transport_rows":len(frame),"bootstrap_rows":len(boot),
                      "max_bootstrap_failures":int(frame.bootstrap_failures.max())},indent=2))


if __name__=="__main__":
    main()
