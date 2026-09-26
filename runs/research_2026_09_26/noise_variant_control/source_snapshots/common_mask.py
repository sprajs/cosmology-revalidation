"""Intersect exact raw-PHOT accepted epochs across four coupled arms."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from collections import Counter, defaultdict, deque
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from astropy.io import fits

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ARMS=("P21","P21_rho000","P21_rho090","P21_noisetrue120")
VAL=ROOT/"runs/research_2026_09_26/astra_design/validation1020"


def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1048576),b""):
            h.update(part)
    return h.hexdigest()


def raw_data(arm,ids):
    version=f"PH2_pilot02_{arm}"
    base=ROOT/f"phase2/literature/simulations/outputs/{version}"
    with fits.open(base/f"{version}_HEAD.FITS",memmap=True) as hh, fits.open(base/f"{version}_PHOT.FITS",memmap=True) as pp:
        head=hh[1].data;phot=pp[1].data
        index={int(x["SNID"]):x for x in head}
        assert len(index)==len(head)
        out={}
        for cid in ids:
            h=index[cid];lo=int(h["PTROBS_MIN"]);hi=int(h["PTROBS_MAX"])
            q=phot[lo-1:hi].copy()
            assert len(q)==int(h["NOBS"])
            out[cid]=q
        return out


def native_object(arm,cid,nominal):
    if arm=="P21":
        path=ROOT/nominal.loc[cid,"objective_path"]
        assert sha(path)==nominal.loc[cid,"original_objective_sha256"]
    else:
        path=HERE/arm/"objectives"/f"objective_{cid}.npz"
    with np.load(path,allow_pickle=False) as f:
        return {k:f[k].copy() for k in f.files}


def row_key(mjd,band,flux,err):
    return (float(mjd),str(band).strip(),float(flux),float(err))


def accepted_raw_indices(obj,raw):
    pools=defaultdict(deque)
    for i,x in enumerate(raw):
        pools[row_key(x["MJD"],x["BAND"],x["FLUXCAL"],x["FLUXCALERR"])].append(i)
    requested=Counter(row_key(*x) for x in zip(obj["MJD"],obj["band"],obj["data_flux"],obj["data_fluxerr"]))
    ambiguous=sum(1 for key,count in requested.items() if len(pools[key])>1 and count<len(pools[key]))
    rows=[]
    for x in zip(obj["MJD"],obj["band"],obj["data_flux"],obj["data_fluxerr"]):
        key=row_key(*x)
        assert pools[key],key
        rows.append(pools[key].popleft())
    assert len(rows)==len(set(rows))
    return np.array(rows,int),ambiguous


def main():
    assert sha(HERE/"score-plan.md")=="92cc736b4a95033111655234c63328456320e77f885438b5ec17e0dea75268f0"
    gate=json.loads((HERE/"object_scores/score-gate.json").read_text())
    assert gate["scored"]==1023
    dest=HERE/"common_mask";dest.mkdir(exist_ok=False)
    primary=pd.read_csv(HERE/"object_scores/object-scores.csv").set_index(["arm","CID"])
    nominal=pd.read_csv(HERE/"nominal-reuse-ledger.csv").set_index("CID")
    ids=set.intersection(*(set(primary.loc[arm].index) for arm in ARMS))
    assert len(ids)==255
    raws={arm:raw_data(arm,ids) for arm in ARMS}
    spec=importlib.util.spec_from_file_location("noise_score",HERE/"score_noise.py")
    ns=importlib.util.module_from_spec(spec);assert spec.loader;spec.loader.exec_module(ns)
    frspec=importlib.util.spec_from_file_location("flux_response",ROOT/"scripts/salt_dust_audit/flux_response.py")
    fr=importlib.util.module_from_spec(frspec);assert frspec.loader;frspec.loader.exec_module(fr)
    model,bands,_,zp=fr.build_model()
    magoff={str(row["Filter Name"])[-1]:float(row["Primary Mag"]) for row in zp}
    with np.load(VAL/"frozen-discovery-coefficients.npz",allow_pickle=False) as f:
        coeff=f["basis_mean"].astype(float)
    rows=[];failures=[]
    for cid in sorted(ids):
        baseline=raws["P21"][cid]
        for arm in ARMS[1:]:
            other=raws[arm][cid]
            assert len(other)==len(baseline)
            assert np.array_equal(other["MJD"],baseline["MJD"])
            assert np.array_equal(other["BAND"],baseline["BAND"])
        objects={a:native_object(a,cid,nominal) for a in ARMS}
        mappings={a:accepted_raw_indices(objects[a],raws[a][cid]) for a in ARMS}
        ambiguous=sum(n for _,n in mappings.values())
        if ambiguous:
            failures.append(dict(CID=cid,reason="ambiguous duplicate raw-PHOT occurrence",ambiguous_groups=ambiguous))
            continue
        shared=set.intersection(*(set(x) for x,_ in mappings.values()))
        if len(shared)<5:
            failures.append(dict(CID=cid,reason="fewer than five common accepted raw epochs",shared_epochs=len(shared)))
            continue
        for arm in ARMS:
            obj=objects[arm]
            raw_idx=mappings[arm][0]
            pars=obj["parameters_x0_x1_c_t0"]
            row=SimpleNamespace(CID=str(cid),x0=float(pars[0]),x1=float(pars[1]),
                                c=float(pars[2]),PKMJD=float(pars[3]),zHEL=float(obj["zHEL"][0]))
            points=pd.DataFrame({"MJD":obj["MJD"],"BAND":obj["band"],
                                 "FLUXCAL":obj["data_flux"],"FLUXCALERR":obj["data_fluxerr"]})
            order=points.sort_values(["MJD","BAND"],kind="stable").index.to_numpy()
            _,audit,_,_=fr.analyze(row,points,{"MWEBV":obj["MWEBV"][0]},model,bands,magoff,False)
            selected=np.isin(raw_idx[order],list(shared))
            y=obj["data_flux"][order][selected];f=obj["model_flux"][order][selected]
            b=obj["band"][order][selected]
            C=obj["frozen_flux_covariance"][np.ix_(order,order)][np.ix_(selected,selected)]
            J=audit["jacobian_flux"].copy();J[:,0]=-fr.K*obj["model_flux"][order]
            J=J[selected]
            try:
                result,_,_=ns.projected(y,f,C,J,b,coeff)
            except (AssertionError,np.linalg.LinAlgError) as exc:
                failures.append(dict(CID=cid,arm=arm,reason="common-mask projection rank/information failure",detail=str(exc)))
                break
            p=primary.loc[(arm,cid)]
            original_recovery_error=0.0
            if len(shared)==len(raw_idx):
                original_recovery_error=max(abs(result[k]-float(p[k])) for k in
                                            ("projected_chi2","matched_filter","information","fixed_gain"))
                assert original_recovery_error<1e-9,(cid,arm,original_recovery_error)
            rows.append(dict(arm=arm,CID=cid,shared_epochs=len(shared),original_epochs=len(raw_idx),
                             dropped_epochs=len(raw_idx)-len(shared),
                             duplicate_ambiguous_groups=ambiguous,
                             original_recovery_error=original_recovery_error,
                             **result))
    frame=pd.DataFrame(rows)
    frame.to_csv(dest/"common-mask-object-scores.csv",index=False,float_format="%.17g")
    pd.DataFrame(failures).to_csv(dest/"common-mask-failures.csv",index=False)
    outcome={"common_native_success_objects":len(ids),"scored_rows":len(frame),
             "scored_complete_four_arm_objects":int(frame.groupby("CID").arm.nunique().eq(4).sum()) if len(frame) else 0,
             "failed_object_records":len(failures),
             "maximum_unmasked_recovery_error":float(frame.original_recovery_error.max()) if len(frame) else None,
             "total_epoch_drops_by_arm":frame.groupby("arm").dropped_epochs.sum().to_dict() if len(frame) else {}}
    (dest/"gate.json").write_text(json.dumps(outcome,indent=2)+"\n")
    # Keep the primary support and weights fixed for the selected-object mask comparison.
    supported=pd.read_csv(HERE/"transport/supported-object-ledger.csv")
    sens=[]
    for stage in ("common_native_success","all_four_quality_classifier_intersection",
                  "all_four_intersection_nonnegative_grid"):
        for split in ("field_z","field_z_snr15"):
            target=supported.loc[(supported.stage==stage)&(supported.transport==split)]
            for arm in ARMS:
                q=target.loc[target.arm==arm]
                s=frame.loc[frame.arm==arm].set_index("CID")
                missing=sorted(set(q.CID)-set(s.index))
                if missing:
                    sens.append(dict(stage=stage,transport=split,arm=arm,missing_count=len(missing),
                                     missing_CIDs="|".join(map(str,missing))))
                    continue
                x=s.loc[q.CID]
                w=q.base_weight.to_numpy(float)
                sens.append(dict(stage=stage,transport=split,arm=arm,missing_count=0,
                                 supported_sim=len(q),supported_real=float(w.sum()),
                                 weighted_M=float(w@x.matched_filter.to_numpy(float)),
                                 weighted_I=float(w@x.information.to_numpy(float)),
                                 amplitude=float((w@x.matched_filter.to_numpy(float))/(w@x.information.to_numpy(float))),
                                 mean_gain=float((w@x.fixed_gain.to_numpy(float))/w.sum()),
                                 pooled_Q=float((w@x.projected_chi2.to_numpy(float))/(w@x.projected_dimension.to_numpy(float))),
                                 weighted_dropped_epochs=float(w@x.dropped_epochs.to_numpy(float))))
    pd.DataFrame(sens).to_csv(dest/"fixed-support-sensitivity.csv",index=False,float_format="%.17g")
    sources=[Path(__file__),HERE/"score-plan.md",HERE/"object_scores/manifest.json",
             HERE/"transport/manifest.json",HERE/"nominal-reuse-ledger.csv",
             ROOT/"scripts/salt_dust_audit/flux_response.py"]
    (dest/"manifest.json").write_text(json.dumps({"inputs_sha256":{str(x.relative_to(ROOT)):sha(x) for x in sources},
                                                   "outputs_sha256":{str(x.relative_to(ROOT)):sha(x) for x in dest.iterdir() if x.name!="manifest.json"}},indent=2)+"\n")
    print(json.dumps(outcome,indent=2))


if __name__=="__main__":
    main()
