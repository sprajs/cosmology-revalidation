"""Score fixed observer direction on all successful coupled native exports."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DESIGN = ROOT / "runs/research_2026_09_26/astra_design/noise_variant_design"
PILOT = ROOT / "runs/research_2026_09_26/simulation_residual_control/full256/P21/approx_minus99"
HOLD = ROOT / "runs/research_2026_09_26/simulation_holdout_control/P21"
VAL = ROOT / "runs/research_2026_09_26/astra_design/validation1020"
COMMON = ROOT / "runs/research_2026_09_26/common_classifier_residual"
HOLDCLASS = ROOT / "runs/research_2026_09_26/simulation_holdout_control"
ARMS = ("P21", "P21_rho000", "P21_rho090", "P21_noisetrue120")
EXPECTED_PLAN = "92cc736b4a95033111655234c63328456320e77f885438b5ec17e0dea75268f0"


def digest(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for part in iter(lambda:f.read(1048576),b""):
            h.update(part)
    return h.hexdigest()


def gate():
    assert digest(HERE/"score-plan.md")==EXPECTED_PLAN
    nominal=json.loads((HERE/"nominal-reuse-gate.json").read_text())
    assert nominal["nominal_ids"]==256 and nominal["exact_normalized_fit_nml"]
    assert digest(HERE/"nominal-reuse-ledger.csv")==nominal["nominal_ledger_sha256"]
    inference=json.loads((HERE/"inference/inference-gate.json").read_text())
    assert inference["count"]==768 and inference["max_single_batch_abs"]<1e-6
    assert not inference["threshold_within_1e6"]
    assert digest(HERE/"inference/variant-probabilities.csv")==json.loads((HERE/"inference/manifest.json").read_text())["outputs_sha256"][str((HERE/"inference/variant-probabilities.csv").relative_to(ROOT))]
    for arm in ARMS[1:]:
        path=HERE/arm
        fit=json.loads((path/"fit-gate.json").read_text())
        prov=json.loads((path/"provenance-gate.json").read_text())
        tang=json.loads((path/"tangent-gate.json").read_text())
        edge=json.loads((path/"mean-edge-gate.json").read_text())
        assert fit["returncode"]==0 and fit["graceful"] and len(fit["fitres_cids"])>=.95*256
        assert prov["full_pilot_adequate"] and prov["exported"]==len(fit["fitres_cids"])
        assert tang["all_rank4"] and edge["eligible_objects"]==prov["exported"]
        assert prov["max_raw_mjd_difference"]==prov["max_raw_flux_difference"]==prov["max_raw_error_difference"]==0
        assert digest(path/"provenance-ledger.csv")==prov["provenance_ledger_sha256"]
        assert digest(path/"mean-edge-ledger.csv")==edge["ledger_sha256"]
    return nominal


def fit_table(path):
    sys.path.insert(0,str(ROOT/"scripts/phase2/official"))
    from audit_fits import read_fit
    x=read_fit(path/"fit.FITRES.TEXT")
    x.CID=x.CID.astype(int)
    assert x.CID.is_unique
    return x.set_index("CID")


def classifier_table():
    a=pd.read_csv(COMMON/"inference/simulation-probabilities.csv")
    b=pd.read_csv(HOLDCLASS/"inference/holdout-probabilities.csv")
    n=pd.concat([a.loc[a.arm=="P21"],b.loc[b.arm=="P21"]],ignore_index=True)
    assert n.CID.is_unique
    n["arm"]="P21"
    v=pd.read_csv(HERE/"inference/variant-probabilities.csv")
    p=pd.concat([n,v],ignore_index=True)
    assert not p.duplicated(["arm","CID"]).any()
    return p.set_index(["arm","CID"])


def support_table():
    p=pd.read_csv(ROOT/"runs/research_2026_09_26/simulation_residual_control/truth-support-ledger.csv")
    h=pd.read_csv(HOLDCLASS/"truth-support-ledger.csv")
    q=pd.concat([p.loc[p.arm=="P21"],h.loc[h.arm=="P21"]],ignore_index=True)
    assert q.CID.is_unique
    return q.set_index("CID")


def projected(y,f,C,J,b,coeff):
    assert len(y)==len(f)==len(b)==C.shape[0]==J.shape[0] and J.shape[1]==4
    L=np.linalg.cholesky(C)
    WJ=solve_triangular(L,J,lower=True)
    U,s,_=np.linalg.svd(WJ,full_matrices=True)
    rank=int(np.sum(s>s[0]*1e-10))
    assert rank==4
    Q=U[:,4:]
    res=Q.T@solve_triangular(L,y-f,lower=True)
    K=0.4*np.log(10)
    griz=np.column_stack([-K*f*(b==band) for band in "griz"])
    obs=np.column_stack([griz[:,0]-griz[:,1],griz[:,2]-griz[:,1],griz[:,3]-griz[:,1]])
    T=Q.T@solve_triangular(L,obs,lower=True)
    u=T.T@res
    F=T.T@T
    M=float(coeff@u)
    I=float(coeff@F@coeff)
    assert np.isfinite(M) and np.isfinite(I) and I>0
    return dict(projected_chi2=float(res@res),projected_dimension=len(res),
                matched_filter=M,information=I,fixed_gain=M-I/2,
                nuisance_condition=float(s[0]/s[-1])),u,F


def main():
    gate()
    dest=HERE/"object_scores"
    dest.mkdir(exist_ok=False)
    with np.load(VAL/"frozen-discovery-coefficients.npz",allow_pickle=False) as f:
        coeff=f["basis_mean"].astype(float)
    assert coeff.shape==(3,)
    spec=importlib.util.spec_from_file_location("flux_response",ROOT/"scripts/salt_dust_audit/flux_response.py")
    fr=importlib.util.module_from_spec(spec);assert spec.loader;spec.loader.exec_module(fr)
    model,bands,_,zp=fr.build_model()
    magoff={str(row["Filter Name"])[-1]:float(row["Primary Mag"]) for row in zp}
    qnom=pd.read_csv(HERE/"nominal-reuse-ledger.csv").set_index("CID")
    p=classifier_table()
    support=support_table()
    ids=[int(x) for x in (DESIGN/"P21-cids.txt").read_text().split()]
    assert len(ids)==256
    fits={"P21_pilot":fit_table(PILOT),"P21_holdout":fit_table(HOLD)}
    fits.update({arm:fit_table(HERE/arm) for arm in ARMS[1:]})
    rows=[]; ulist=[]; Flist=[]; selection=[]
    for arm in ARMS:
        cohort=pd.read_csv(DESIGN/f"{arm}-cohort.csv").set_index("CID").loc[ids]
        assert np.array_equal(cohort.generated_attempt_index.to_numpy(int),
                              pd.read_csv(DESIGN/"P21-cohort.csv").set_index("CID").loc[ids].generated_attempt_index.to_numpy(int))
        for j,cid in enumerate(ids):
            arch=cohort.loc[cid]
            fitkey=("P21_"+qnom.loc[cid,"original_scope"]) if arm=="P21" else arm
            fitrow=fits[fitkey].loc[cid] if cid in fits[fitkey].index else None
            if arm=="P21":
                path=ROOT/qnom.loc[cid,"objective_path"]
                assert digest(path)==qnom.loc[cid,"original_objective_sha256"]
            else:
                path=HERE/arm/"objectives"/f"objective_{cid}.npz"
            success=fitrow is not None and path.exists()
            prob=p.loc[(arm,cid)]
            assert bool(prob.gt999)==bool(float(prob.pIa)>.999)
            quality=(abs(float(fitrow.x1))<3 and abs(float(fitrow.c))<.3
                     and float(fitrow.x1ERR)<1 and float(fitrow.PKMJDERR)<2
                     and float(fitrow.cERR)<1.5 and float(fitrow.FITPROB)>.001
                     and .025<float(fitrow.zHD)<1.2) if success else False
            selection.append(dict(arm=arm,CID=cid,generated_attempt_index=int(arch.generated_attempt_index),
                                  native_success=success,archived_basic_quality=bool(arch.basic_quality_pass),
                                  fresh_basic_quality=quality,pIa=float(prob.pIa),gt999=bool(prob.gt999),
                                  selected_archived_quality_classifier=bool(arch.basic_quality_pass and prob.gt999 and success),
                                  support_class=str(support.loc[cid,"support_class"]),
                                  missing_fitres=fitrow is None))
            if not success:
                continue
            with np.load(path,allow_pickle=False) as obj:
                pars=obj["parameters_x0_x1_c_t0"]
                row=SimpleNamespace(CID=str(cid),x0=float(pars[0]),x1=float(pars[1]),
                                    c=float(pars[2]),PKMJD=float(pars[3]),zHEL=float(obj["zHEL"][0]))
                points=pd.DataFrame({"MJD":obj["MJD"],"BAND":obj["band"],
                                     "FLUXCAL":obj["data_flux"],"FLUXCALERR":obj["data_fluxerr"]})
                order=points.sort_values(["MJD","BAND"],kind="stable").index.to_numpy()
                _,audit,_,diag=fr.analyze(row,points,{"MWEBV":obj["MWEBV"][0]},model,bands,magoff,False)
                y=obj["data_flux"][order];f=obj["model_flux"][order]
                b=obj["band"][order];t=obj["rest_phase"][order]
                C=obj["frozen_flux_covariance"][np.ix_(order,order)]
                J=audit["jacobian_flux"].copy();J[:,0]=-fr.K*f
                result,u,F=projected(y,f,C,J,b,coeff)
                covset=obj["covariance_components"][order]
                assert covset.shape==(len(y),5)
                trace=float(np.trace(C));datatrace=float(np.sum(covset[:,3]**2))
                fudgetrace=float(np.sum(covset[:,4]**2))
                remainder=trace-datatrace-fudgetrace
                assert remainder>=-1e-6 and trace>0
                rows.append(dict(arm=arm,CID=cid,generated_attempt_index=int(arch.generated_attempt_index),
                                 field=str(arch.field),zHEL=float(arch.zHEL),zbin=int(arch.zbin),
                                 cell=str(arch.cell),LIBID=int(arch.SIM_LIBID),
                                 SNRMAX1_archived=float(arch.SNRMAX1),
                                 support_class=str(support.loc[cid,"support_class"]),
                                 archived_basic_quality=bool(arch.basic_quality_pass),
                                 fresh_basic_quality=bool(quality),pIa=float(prob.pIa),gt999=bool(prob.gt999),
                                 epochs=len(y),phase_min=float(np.min(t)),phase_max=float(np.max(t)),
                                 epochs_prepeak=int(np.sum(t<0)),epochs_postpeak=int(np.sum(t>=0)),
                                 **{f"epochs_{band}":int(np.sum(b==band)) for band in "griz"},
                                 covariance_trace=trace,data_diagonal_trace=datatrace,
                                 fudge_diagonal_trace=fudgetrace,remaining_trace=remainder,
                                 data_diagonal_trace_fraction=datatrace/trace,
                                 remaining_trace_fraction=remainder/trace,
                                 max_native_independent_mean_sigma=float(np.max(abs((audit["flux_model"]-f)/obj["data_fluxerr"][order]))),
                                 derivative_halfstep_error=float(diag["derivative_relative_error"]),
                                 **result))
                ulist.append(u);Flist.append(F)
        print(f"scored {arm}",flush=True)
    frame=pd.DataFrame(rows).sort_values(["arm","CID"])
    select=pd.DataFrame(selection).sort_values(["arm","CID"])
    assert len(select)==1024 and len(frame)==1023
    assert select.groupby("arm").native_success.sum().to_dict()=={"P21":256,"P21_rho000":256,"P21_rho090":256,"P21_noisetrue120":255}
    # Keep u/F aligned to the explicitly sorted object table.
    lookup={(r["arm"],r["CID"]):(u,F) for r,u,F in zip(rows,ulist,Flist)}
    frame.to_csv(dest/"object-scores.csv",index=False,float_format="%.17g")
    select.to_csv(dest/"selection-ledger.csv",index=False,float_format="%.17g")
    np.savez_compressed(dest/"sufficient-arrays.npz",
                        arm=frame.arm.to_numpy(dtype="U20"),CID=frame.CID.to_numpy(int),
                        u=np.array([lookup[(r.arm,r.CID)][0] for r in frame.itertuples()]),
                        F=np.array([lookup[(r.arm,r.CID)][1] for r in frame.itertuples()]))
    # Reused nominal scores must numerically reproduce immutable pilot/holdout products.
    pilot=pd.read_csv(ROOT/"runs/research_2026_09_26/simulation_residual_control/score/object-scores.csv")
    pilot=pilot.loc[(pilot.arm=="P21")&(pilot.law=="approx_minus99")]
    hold=pd.read_csv(HOLDCLASS/"object_scores/holdout-object-scores.csv")
    hold=hold.loc[hold.arm=="P21"]
    old=pd.concat([pilot,hold]).set_index("CID").loc[ids]
    new=frame.loc[frame.arm=="P21"].set_index("CID").loc[ids]
    columns=("projected_chi2","matched_filter","information","fixed_gain")
    maxdiff={col:float(np.max(abs(old[col]-new[col]))) for col in columns}
    assert all(v<1e-9 for v in maxdiff.values()),maxdiff
    report={"scored":len(frame),"by_arm":frame.groupby("arm").size().to_dict(),
            "frozen_coefficients":coeff.tolist(),"nominal_reuse_score_max_difference":maxdiff,
            "rank4_all":True,"positive_information_all":True,
            "maximum_halfstep":frame.groupby("arm").derivative_halfstep_error.max().to_dict(),
            "missing_fitres":select.loc[select.missing_fitres,["arm","CID"]].to_dict("records")}
    (dest/"score-gate.json").write_text(json.dumps(report,indent=2)+"\n")
    sources=[Path(__file__),HERE/"score-plan.md",HERE/"nominal-reuse-gate.json",
             HERE/"nominal-reuse-ledger.csv",HERE/"inference/manifest.json",
             ROOT/"scripts/salt_dust_audit/flux_response.py",VAL/"frozen-discovery-coefficients.npz"]
    sources += [HERE/arm/name for arm in ARMS[1:] for name in ("prepared.json","fit-gate.json","export-gate.json","provenance-gate.json","tangent-gate.json","mean-edge-gate.json","fit.FITRES.TEXT","objectives/manifest.json")]
    (dest/"manifest.json").write_text(json.dumps({"inputs_sha256":{str(x.relative_to(ROOT)):digest(x) for x in sources},
                                                   "outputs_sha256":{str(x.relative_to(ROOT)):digest(x) for x in dest.iterdir() if x.name!="manifest.json"}},indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
