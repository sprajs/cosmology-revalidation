"""Transport coupled-arm fixed scores on frozen cells with paired LIBID draws."""
from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
VAL=ROOT/"runs/research_2026_09_26/astra_design/validation1020"
COMMON=ROOT/"runs/research_2026_09_26/common_classifier_residual"
ARMS=("P21","P21_rho000","P21_rho090","P21_noisetrue120")
STAGES=("common_native_success","arm_specific_quality_classifier","all_four_quality_classifier_intersection",
        "arm_specific_quality_classifier_nonnegative_grid","all_four_intersection_nonnegative_grid")


def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for part in iter(lambda:f.read(1048576),b""):
            h.update(part)
    return h.hexdigest()


def main():
    assert sha(HERE/"score-plan.md")=="92cc736b4a95033111655234c63328456320e77f885438b5ec17e0dea75268f0"
    scoregate=json.loads((HERE/"object_scores/score-gate.json").read_text())
    assert scoregate["scored"]==1023 and scoregate["rank4_all"]
    dest=HERE/"transport";dest.mkdir(exist_ok=False)
    scores=pd.read_csv(HERE/"object_scores/object-scores.csv")
    select=pd.read_csv(HERE/"object_scores/selection-ledger.csv")
    assert len(scores)==1023 and len(select)==1024
    assert not scores.duplicated(["arm","CID"]).any()
    assert not select.duplicated(["arm","CID"]).any()
    common_ids=set.intersection(*(set(scores.loc[scores.arm==a,"CID"]) for a in ARMS))
    assert len(common_ids)==255
    both=select.loc[select.selected_archived_quality_classifier].groupby("CID").arm.nunique()
    intersect=set(both[both==4].index)&common_ids
    nonneg=set(select.loc[(select.arm=="P21")&(select.support_class=="nonnegative_declared_grid"),"CID"])
    q={a:scores.loc[scores.arm==a].copy() for a in ARMS}
    nominal_snr=q["P21"].set_index("CID").SNRMAX1_archived
    real_ids=set(pd.read_csv(COMMON/"inference/real-probabilities.csv").query("gt999").CID.astype(int))
    assert len(real_ids)==1006
    real=pd.read_csv(VAL/"analysis/object-scores.csv")
    real=real.loc[(real.arm=="published_mask")&real.CID.isin(real_ids)].copy()
    assert len(real)==1006
    real["fixed_gain"]=real.fixed_prediction_gain
    edges=[.05,.2,.35,.5,.65,.8,1.2]
    real["zbin"]=pd.cut(real.zHEL,edges,include_lowest=True,labels=False)
    assert real.zbin.notna().all()
    real["zbin"]=real.zbin.astype(int)
    real["cell"]=real.field+"_"+real.zbin.astype(str)
    sys.path.insert(0,str(ROOT/"scripts/phase2/official"))
    from audit_fits import read_fit
    rf=read_fit(VAL/"fit.FITRES.TEXT")[["CID","SNRMAX1"]]
    rf.CID=rf.CID.astype(int)
    real=real.merge(rf,on="CID",validate="one_to_one")
    real["snrbin"]=np.where(real.SNRMAX1<15,"lt15","ge15")
    libs=np.array(sorted(set(select.loc[select.arm=="P21","generated_attempt_index"].astype(int))))
    # SIM_LIBID is the cadence-block unit; several CIDs can share it.
    libs=np.array(sorted(set(scores.LIBID.astype(int))))
    rng=np.random.default_rng(26092694)
    draws=rng.integers(0,len(libs),size=(1000,len(libs)))
    multi=np.array([np.bincount(row,minlength=len(libs)) for row in draws],dtype=np.int16)
    libidx={lib:i for i,lib in enumerate(libs)}
    summaries=[];support=[];fields=[];boots=[];stages=[];flips=[]
    for arm in ARMS:
        s=select.loc[select.arm==arm]
        x=q[arm]
        for label,mask in (("attempted",np.ones(len(s),bool)),("native_success",s.native_success),
                           ("archived_quality",s.native_success&s.archived_basic_quality),
                           ("fresh_quality",s.native_success&s.fresh_basic_quality),
                           ("classifier",s.native_success&s.gt999),
                           ("archived_quality_classifier",s.selected_archived_quality_classifier),
                           ("all_four_quality_classifier_intersection",s.CID.isin(intersect))):
            ids=set(s.loc[mask,"CID"])
            y=x.loc[x.CID.isin(ids)]
            stages.append(dict(arm=arm,stage=label,attempts=len(ids),scored=len(y),
                               raw_M=float(y.matched_filter.sum()),raw_I=float(y.information.sum()),
                               raw_G=float(y.fixed_gain.sum()),raw_Q=float(y.projected_chi2.sum()),
                               raw_dimension=int(y.projected_dimension.sum()),
                               CIDs="|".join(map(str,sorted(ids)))))
        if arm!="P21":
            n=select.loc[select.arm=="P21"].set_index("CID").loc[s.CID]
            flips.append(dict(arm=arm,native_success_flips=int(np.sum(s.native_success.to_numpy()!=n.native_success.to_numpy())),
                              archived_quality_flips=int(np.sum(s.archived_basic_quality.to_numpy()!=n.archived_basic_quality.to_numpy())),
                              fresh_quality_flips=int(np.sum(s.fresh_basic_quality.to_numpy()!=n.fresh_basic_quality.to_numpy())),
                              classifier_flips=int(np.sum(s.gt999.to_numpy()!=n.gt999.to_numpy())),
                              combined_selection_flips=int(np.sum(s.selected_archived_quality_classifier.to_numpy()!=n.selected_archived_quality_classifier.to_numpy()))))
    for stage in STAGES:
        subsets={}
        for arm in ARMS:
            x=q[arm]
            if stage=="common_native_success":
                mask=x.CID.isin(common_ids)
            elif stage=="arm_specific_quality_classifier":
                mask=x.archived_basic_quality & x.gt999
            elif stage=="all_four_quality_classifier_intersection":
                mask=x.CID.isin(intersect)
            elif stage=="arm_specific_quality_classifier_nonnegative_grid":
                mask=x.archived_basic_quality & x.gt999 & x.CID.isin(nonneg)
            elif stage=="all_four_intersection_nonnegative_grid":
                mask=x.CID.isin(intersect & nonneg)
            else: raise AssertionError(stage)
            subsets[arm]=x.loc[mask].copy()
        for split in ("field_z","field_z_snr15"):
            r=real.copy()
            r["transport_cell"]=r.cell if split=="field_z" else r.cell+"_"+r.snrbin
            grouped={}
            for arm in ARMS:
                x=subsets[arm].copy()
                snr=x.SNRMAX1_archived if stage.startswith("arm_specific") else x.CID.map(nominal_snr)
                assert snr.notna().all()
                x["snrbin"]=np.where(snr<15,"lt15","ge15")
                x["transport_cell"]=x.cell if split=="field_z" else x.cell+"_"+x.snrbin
                grouped[arm]=x
            cellsets=[]
            for arm in ARMS:
                counts=grouped[arm].transport_cell.value_counts()
                cellsets.append(set(counts[counts>=2].index))
            cells=sorted(set.intersection(*cellsets)&set(r.transport_cell))
            assert cells,(stage,split)
            target=r.loc[r.transport_cell.isin(cells)]
            assert len(target)>0 and target.information.sum()>0
            boot_by_arm={}
            supported_by_arm={}
            for arm in ARMS:
                full=grouped[arm]
                x=full.loc[full.transport_cell.isin(cells)].copy()
                assert len(x)>0
                ns=x.transport_cell.value_counts();nr=target.transport_cell.value_counts()
                w=x.transport_cell.map(nr).to_numpy(float)/x.transport_cell.map(ns).to_numpy(float)
                assert abs(w.sum()-len(target))<1e-8
                M=x.matched_filter.to_numpy(float);I=x.information.to_numpy(float)
                G=x.fixed_gain.to_numpy(float);Q=x.projected_chi2.to_numpy(float)
                nu=x.projected_dimension.to_numpy(float)
                Cfrac=x.data_diagonal_trace_fraction.to_numpy(float)
                suffix=dict(stage=stage,transport=split,arm=arm)
                summaries.append(dict(**suffix,selected=len(full),supported_sim=len(x),
                                      supported_cells=len(cells),supported_real=len(target),
                                      real_fraction=len(target)/1006,
                                      weighted_M=float(w@M),weighted_I=float(w@I),
                                      weighted_G=float(w@G),weighted_Q=float(w@Q),
                                      weighted_dimension=float(w@nu),
                                      amplitude=float((w@M)/(w@I)),mean_gain=float((w@G)/w.sum()),
                                      pooled_Q=float((w@Q)/(w@nu)),
                                      effective_sample_size=float(w.sum()**2/(w@w)),
                                      weighted_data_trace_fraction=float((w@Cfrac)/w.sum()),
                                      weighted_epochs=float(w@x.epochs.to_numpy(float)),
                                      weighted_epochs_prepeak=float(w@x.epochs_prepeak.to_numpy(float)),
                                      weighted_epochs_postpeak=float(w@x.epochs_postpeak.to_numpy(float)),
                                      **{f"weighted_epochs_{band}":float(w@x[f"epochs_{band}"].to_numpy(float)) for band in "griz"},
                                      real_M=float(target.matched_filter.sum()),
                                      real_I=float(target.information.sum()),
                                      real_G=float(target.fixed_gain.sum()),
                                      real_amplitude=float(target.matched_filter.sum()/target.information.sum()),
                                      real_mean_gain=float(target.fixed_gain.mean()),
                                      real_pooled_Q=float(target.projected_chi2.sum()/target.projected_dimension.sum()),
                                      cell_ids="|".join(cells)))
                support.extend(dict(**suffix,CID=int(z.CID),LIBID=int(z.LIBID),
                                    cell=z.transport_cell,base_weight=float(ww)) for z,ww in zip(x.itertuples(),w))
                for field,fq in x.groupby("field"):
                    fw=w[np.asarray(x.field==field)]
                    fields.append(dict(**suffix,field=field,n=len(fq),weighted_M=float(fw@fq.matched_filter.to_numpy(float)),
                                       weighted_I=float(fw@fq.information.to_numpy(float)),
                                       weighted_G=float(fw@fq.fixed_gain.to_numpy(float)),
                                       weighted_Q=float(fw@fq.projected_chi2.to_numpy(float)),
                                       weighted_dimension=float(fw@fq.projected_dimension.to_numpy(float))))
                bm=multi[:,[libidx[int(z)] for z in x.LIBID]]*w[None,:]
                denI=bm@I;denW=bm.sum(axis=1);denNu=bm@nu
                good=(denI>0)&(denW>0)&(denNu>0)
                amp=np.divide(bm@M,denI,out=np.full(1000,np.nan),where=good)
                mg=np.divide(bm@G,denW,out=np.full(1000,np.nan),where=good)
                pq=np.divide(bm@Q,denNu,out=np.full(1000,np.nan),where=good)
                boot_by_arm[arm]=(good,amp,mg,pq)
                supported_by_arm[arm]=(set(x.CID),dict(zip(x.CID,w)))
            if stage in ("common_native_success","all_four_quality_classifier_intersection","all_four_intersection_nonnegative_grid"):
                baseline=supported_by_arm["P21"]
                for arm in ARMS[1:]:
                    other=supported_by_arm[arm]
                    assert baseline[0]==other[0]
                    assert all(abs(baseline[1][cid]-other[1][cid])<1e-12 for cid in baseline[0])
            nom=boot_by_arm["P21"]
            for arm in ARMS:
                good,amp,mg,pq=boot_by_arm[arm]
                both=good&nom[0]
                for j in range(1000):
                    boots.append(dict(stage=stage,transport=split,arm=arm,replicate=j,LIBID_seed=26092694,
                                      valid=bool(good[j]),paired_valid=bool(both[j]),
                                      amplitude=float(amp[j]),mean_gain=float(mg[j]),pooled_Q=float(pq[j]),
                                      difference_amplitude=float(amp[j]-nom[1][j]) if both[j] else np.nan,
                                      difference_mean_gain=float(mg[j]-nom[2][j]) if both[j] else np.nan,
                                      difference_pooled_Q=float(pq[j]-nom[3][j]) if both[j] else np.nan))
    summary=pd.DataFrame(summaries)
    boot=pd.DataFrame(boots)
    for row in summary.itertuples():
        mask=(boot.stage==row.stage)&(boot.transport==row.transport)&(boot.arm==row.arm)
        bb=boot.loc[mask]
        assert len(bb)==1000
        summary.loc[(summary.stage==row.stage)&(summary.transport==row.transport)&(summary.arm==row.arm),"bootstrap_failures"]=int((~bb.valid).sum())
        for column in ("amplitude","mean_gain","pooled_Q","difference_amplitude","difference_mean_gain","difference_pooled_Q"):
            v=bb.loc[bb.valid if not column.startswith("difference") else bb.paired_valid,column]
            summary.loc[(summary.stage==row.stage)&(summary.transport==row.transport)&(summary.arm==row.arm),f"bootstrap_{column}_low"]=float(v.quantile(.025)) if len(v) else np.nan
            summary.loc[(summary.stage==row.stage)&(summary.transport==row.transport)&(summary.arm==row.arm),f"bootstrap_{column}_high"]=float(v.quantile(.975)) if len(v) else np.nan
    summary.to_csv(dest/"transport-summary.csv",index=False,float_format="%.17g")
    pd.DataFrame(support).to_csv(dest/"supported-object-ledger.csv",index=False,float_format="%.17g")
    pd.DataFrame(fields).to_csv(dest/"field-contributions.csv",index=False,float_format="%.17g")
    pd.DataFrame(stages).to_csv(dest/"stage-ledger.csv",index=False,float_format="%.17g")
    pd.DataFrame(flips).to_csv(dest/"selection-flips.csv",index=False,float_format="%.17g")
    boot.to_csv(dest/"bootstrap.csv.gz",index=False,compression={"method":"gzip","mtime":0},float_format="%.17g")
    sources=[Path(__file__),HERE/"score-plan.md",HERE/"object_scores/manifest.json",
             HERE/"object_scores/object-scores.csv",HERE/"object_scores/selection-ledger.csv",
             COMMON/"inference/real-probabilities.csv",VAL/"analysis/object-scores.csv",VAL/"fit.FITRES.TEXT"]
    (dest/"manifest.json").write_text(json.dumps({"inputs_sha256":{str(x.relative_to(ROOT)):sha(x) for x in sources},
                                                   "outputs_sha256":{str(x.relative_to(ROOT)):sha(x) for x in dest.iterdir() if x.name!="manifest.json"},
                                                   "bootstrap_seed":26092694,"bootstrap_replicates":1000},indent=2)+"\n")
    print(json.dumps({"common_success":len(common_ids),"all_four_selected":len(intersect),
                      "transport_rows":len(summary),"bootstrap_rows":len(boot),
                      "maximum_bootstrap_failures":int(summary.bootstrap_failures.max())},indent=2))


if __name__=="__main__":
    main()
