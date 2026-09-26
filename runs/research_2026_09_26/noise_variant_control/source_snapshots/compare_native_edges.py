"""Native-J sensitivity with immutable primary masks, covariance and weights."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
VAL=ROOT/"runs/research_2026_09_26/astra_design/validation1020"
SELECT={"P21_rho000":(1294,7349,20101),"P21_rho090":(1294,7349),
        "P21_noisetrue120":(7349,20101)}
ARMS=("P21","P21_rho000","P21_rho090","P21_noisetrue120")
PLAN_SHA="583c2b783b66a4fff53b216670d3828f6b33e6875af5ab255b5d5eda57148dab"


def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for part in iter(lambda:f.read(1048576),b""):
            h.update(part)
    return h.hexdigest()


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);assert spec.loader;spec.loader.exec_module(mod)
    return mod


def npz(path):
    with np.load(path,allow_pickle=False) as f:
        return {k:f[k].copy() for k in f.files}


def probe(out,arm,cid,a,scale):
    steps={1:.001*scale,2:.0001*scale,3:.01*scale}
    J=np.zeros((len(a["MJD"]),4))
    J[:,0]=-(.4*np.log(10))*a["model_flux"]
    par_err=0.0
    for j,h in steps.items():
        pairs=[]
        for sign in (-1,1):
            b=npz(out/arm/f"p{j}_{'plus' if sign>0 else 'minus'}"/"objectives"/f"objective_{cid}.npz")
            for key in ("MJD","band","data_flux","data_fluxerr"):
                assert np.array_equal(a[key],b[key]),(arm,cid,j,sign,key)
            expect=a["parameters_x0_x1_c_t0"].copy();expect[j]+=sign*h
            err=float(np.max(abs(expect-b["parameters_x0_x1_c_t0"])))
            par_err=max(par_err,err)
            assert err<1e-12,(arm,cid,j,sign,err)
            pairs.append(b["model_flux"])
        J[:,j]=(pairs[1]-pairs[0])/(2*h)
    return J,par_err


def fit_arrays(a,fr,model,bands,magoff):
    pars=a["parameters_x0_x1_c_t0"]
    row=SimpleNamespace(CID="native_edge",x0=float(pars[0]),x1=float(pars[1]),
                        c=float(pars[2]),PKMJD=float(pars[3]),zHEL=float(a["zHEL"][0]))
    points=pd.DataFrame({"MJD":a["MJD"],"BAND":a["band"],
                         "FLUXCAL":a["data_flux"],"FLUXCALERR":a["data_fluxerr"]})
    order=points.sort_values(["MJD","BAND"],kind="stable").index.to_numpy()
    _,audit,_,_=fr.analyze(row,points,{"MWEBV":a["MWEBV"][0]},model,bands,magoff,False)
    y=a["data_flux"][order];f=a["model_flux"][order]
    b=a["band"][order]
    C=a["frozen_flux_covariance"][np.ix_(order,order)]
    J=audit["jacobian_flux"].copy();J[:,0]=-fr.K*f
    return order,y,f,b,C,J


def spectral(C,J1,J2):
    L=np.linalg.cholesky(C)
    a=solve_triangular(L,J1,lower=True)
    b=solve_triangular(L,J2,lower=True)
    u=np.linalg.svd(a,full_matrices=False)[0]
    v=np.linalg.svd(b,full_matrices=False)[0]
    return float(np.linalg.norm(u@u.T-v@v.T,2))


def transport_replacement(original,ledger,delta,key_cols):
    rows=[]
    for r in original.itertuples():
        q=ledger.loc[(ledger.stage==r.stage)&(ledger.transport==r.transport)&(ledger.arm==r.arm)]
        d=delta.loc[delta.arm==r.arm].set_index("CID")
        matched=q.loc[q.CID.isin(d.index)]
        D={col:float(np.sum(matched.base_weight.to_numpy(float)*d.loc[matched.CID,f"delta_{col}"].to_numpy(float)))
           for col in ("M","I","G","Q","dimension")}
        wM=r.weighted_M+D["M"];wI=r.weighted_I+D["I"]
        wG=r.weighted_G+D["G"];wQ=r.weighted_Q+D["Q"]
        wNu=r.weighted_dimension+D["dimension"]
        assert wI>0 and wNu>0
        rows.append(dict(stage=r.stage,transport=r.transport,arm=r.arm,
                         replaced_pairs=len(matched),replaced_CIDs="|".join(map(str,matched.CID)),
                         original_amplitude=r.amplitude,native_J_amplitude=wM/wI,
                         delta_amplitude=wM/wI-r.amplitude,
                         original_mean_gain=r.mean_gain,native_J_mean_gain=wG/r.supported_real,
                         delta_mean_gain=wG/r.supported_real-r.mean_gain,
                         original_pooled_Q=r.pooled_Q,native_J_pooled_Q=wQ/wNu,
                         delta_pooled_Q=wQ/wNu-r.pooled_Q,
                         weighted_delta_M=D["M"],weighted_delta_I=D["I"],
                         weighted_delta_G=D["G"],weighted_delta_Q=D["Q"],
                         supported_sim=r.supported_sim,supported_real=r.supported_real))
    return pd.DataFrame(rows)


def main():
    assert sha(HERE/"native-edge-plan.md")==PLAN_SHA
    full=HERE/"native_edge_full";half=HERE/"native_edge_half"
    for out in (full,half):
        assert json.loads((out/"run-gate.json").read_text())["all_exported"]
        assert sha(out/"mask-ledger.json")==sha(full/"mask-ledger.json")
    dest=HERE/"native_edge_sensitivity";dest.mkdir(exist_ok=False)
    ns=module(HERE/"score_noise.py","noise_score_native_comparison")
    cm=module(HERE/"common_mask.py","common_mask_native_comparison")
    fr=module(ROOT/"scripts/salt_dust_audit/flux_response.py","flux_response_native_comparison")
    model,bands,_,zp=fr.build_model()
    magoff={str(row["Filter Name"])[-1]:float(row["Primary Mag"]) for row in zp}
    coeff=npz(VAL/"frozen-discovery-coefficients.npz")["basis_mean"].astype(float)
    primary=pd.read_csv(HERE/"object_scores/object-scores.csv").set_index(["arm","CID"])
    mask_original=pd.read_csv(HERE/"common_mask/common-mask-object-scores.csv").set_index(["arm","CID"])
    ids=set(cid for values in SELECT.values() for cid in values)
    raw={arm:cm.raw_data(arm,ids) for arm in ARMS}
    results=[];masked=[];j_save={};max_py_recovery=0.;max_mask_py_recovery=0.;max_base_C_abs=0.;max_base_C_rel=0.
    for arm,selected in SELECT.items():
        for cid in selected:
            a=npz(HERE/arm/"objectives"/f"objective_{cid}.npz")
            b=npz(full/arm/"baseline"/"objectives"/f"objective_{cid}.npz")
            for key in ("MJD","band","data_flux","data_fluxerr","model_flux","parameters_x0_x1_c_t0"):
                assert np.array_equal(a[key],b[key]),(arm,cid,key)
            C=a["frozen_flux_covariance"];Cb=b["frozen_flux_covariance"]
            dC=Cb-C
            max_base_C_abs=max(max_base_C_abs,float(np.max(abs(dC))))
            max_base_C_rel=max(max_base_C_rel,float(np.linalg.norm(dC)/np.linalg.norm(C)))
            order,y,f,band,csort,Jpy=fit_arrays(a,fr,model,bands,magoff)
            pyscore,_,_=ns.projected(y,f,csort,Jpy,band,coeff)
            orig=primary.loc[(arm,cid)]
            py_recovery=max(abs(pyscore[k]-float(orig[k])) for k in
                            ("matched_filter","information","fixed_gain","projected_chi2"))
            max_py_recovery=max(max_py_recovery,py_recovery)
            assert py_recovery<1e-9,(arm,cid,py_recovery)
            Jfull,perr=probe(full,arm,cid,a,1.)
            Jhalf,herr=probe(half,arm,cid,a,.5)
            j_save[f"{arm}_{cid}_full"]=Jfull
            j_save[f"{arm}_{cid}_half"]=Jhalf
            fullscore,_,_=ns.projected(y,f,csort,Jfull[order],band,coeff)
            halfscore,_,_=ns.projected(y,f,csort,Jhalf[order],band,coeff)
            row=dict(arm=arm,CID=cid,epochs=len(y),baseline_coordinate_data_mean_exact=True,
                     baseline_C_max_abs=float(np.max(abs(dC))),
                     baseline_C_relative_frobenius=float(np.linalg.norm(dC)/np.linalg.norm(C)),
                     derivative_parameter_max_abs=max(perr,herr),
                     projector_spectral_difference=spectral(csort,Jpy,Jfull[order]),
                     original_python_recovery_error=py_recovery)
            for k,col in (("matched_filter","M"),("information","I"),("fixed_gain","G"),
                          ("projected_chi2","Q"),("projected_dimension","dimension")):
                row[f"original_{col}"]=float(orig[k])
                row[f"native_full_{col}"]=fullscore[k]
                row[f"native_half_{col}"]=halfscore[k]
                row[f"delta_{col}"]=fullscore[k]-float(orig[k])
                row[f"half_minus_full_{col}"]=halfscore[k]-fullscore[k]
            results.append(row)
            mappings={x:cm.accepted_raw_indices(cm.native_object(x,cid,pd.read_csv(HERE/"nominal-reuse-ledger.csv").set_index("CID")),raw[x][cid])[0] for x in ARMS}
            shared=set.intersection(*(set(v) for v in mappings.values()))
            keep=np.isin(mappings[arm][order],list(shared))
            assert keep.sum()>=5
            yk=y[keep];fk=f[keep];bk=band[keep]
            Ck=csort[np.ix_(keep,keep)]
            oldmask=mask_original.loc[(arm,cid)]
            py_mask,_,_=ns.projected(yk,fk,Ck,Jpy[keep],bk,coeff)
            mask_recovery=max(abs(py_mask[k]-float(oldmask[k])) for k in
                              ("matched_filter","information","fixed_gain","projected_chi2"))
            max_mask_py_recovery=max(max_mask_py_recovery,mask_recovery)
            assert mask_recovery<1e-9,(arm,cid,mask_recovery)
            nat_mask,_,_=ns.projected(yk,fk,Ck,Jfull[order][keep],bk,coeff)
            half_mask,_,_=ns.projected(yk,fk,Ck,Jhalf[order][keep],bk,coeff)
            mr=dict(arm=arm,CID=cid,shared_epochs=int(keep.sum()),
                    original_python_mask_recovery_error=mask_recovery)
            for k,col in (("matched_filter","M"),("information","I"),("fixed_gain","G"),
                          ("projected_chi2","Q"),("projected_dimension","dimension")):
                mr[f"original_{col}"]=float(oldmask[k])
                mr[f"native_full_{col}"]=nat_mask[k]
                mr[f"native_half_{col}"]=half_mask[k]
                mr[f"delta_{col}"]=nat_mask[k]-float(oldmask[k])
                mr[f"half_minus_full_{col}"]=half_mask[k]-nat_mask[k]
            masked.append(mr)
    frame=pd.DataFrame(results).sort_values(["arm","CID"])
    maskframe=pd.DataFrame(masked).sort_values(["arm","CID"])
    assert len(frame)==len(maskframe)==7
    frame.to_csv(dest/"native-edge-object-scores.csv",index=False,float_format="%.17g")
    maskframe.to_csv(dest/"native-edge-common-mask-scores.csv",index=False,float_format="%.17g")
    np.savez_compressed(dest/"native-jacobians.npz",**j_save)
    primary_sum=pd.read_csv(HERE/"transport/transport-summary.csv")
    support=pd.read_csv(HERE/"transport/supported-object-ledger.csv")
    fixed=transport_replacement(primary_sum,support,frame,("stage","transport","arm"))
    fixed.to_csv(dest/"fixed-transport-sensitivity.csv",index=False,float_format="%.17g")
    old_mask=pd.read_csv(HERE/"common_mask/fixed-support-sensitivity.csv")
    mask_rows=[]
    for r in old_mask.itertuples():
        q=support.loc[(support.stage==r.stage)&(support.transport==r.transport)&(support.arm==r.arm)]
        d=maskframe.loc[maskframe.arm==r.arm].set_index("CID")
        m=q.loc[q.CID.isin(d.index)]
        dM=float(np.sum(m.base_weight.to_numpy(float)*d.loc[m.CID,"delta_M"].to_numpy(float)))
        dI=float(np.sum(m.base_weight.to_numpy(float)*d.loc[m.CID,"delta_I"].to_numpy(float)))
        dG=float(np.sum(m.base_weight.to_numpy(float)*d.loc[m.CID,"delta_G"].to_numpy(float)))
        dQ=float(np.sum(m.base_weight.to_numpy(float)*d.loc[m.CID,"delta_Q"].to_numpy(float)))
        dNu=float(np.sum(m.base_weight.to_numpy(float)*d.loc[m.CID,"delta_dimension"].to_numpy(float)))
        # Original mask table gives weighted M/I but not weighted chi2/dim;
        # derive the latter from its per-object ledger on the same fixed support.
        mm=mask_original.loc[[(r.arm,int(cid)) for cid in q.CID]]
        w=q.base_weight.to_numpy(float)
        M0=float(w@mm.matched_filter.to_numpy(float));I0=float(w@mm.information.to_numpy(float))
        G0=float(w@mm.fixed_gain.to_numpy(float));Q0=float(w@mm.projected_chi2.to_numpy(float))
        Nu0=float(w@mm.projected_dimension.to_numpy(float))
        assert abs(M0-r.weighted_M)<1e-8 and abs(I0-r.weighted_I)<1e-8
        mask_rows.append(dict(stage=r.stage,transport=r.transport,arm=r.arm,
                              replaced_pairs=len(m),replaced_CIDs="|".join(map(str,m.CID)),
                              original_amplitude=r.amplitude,native_J_amplitude=(M0+dM)/(I0+dI),
                              delta_amplitude=(M0+dM)/(I0+dI)-r.amplitude,
                              original_mean_gain=r.mean_gain,native_J_mean_gain=(G0+dG)/w.sum(),
                              original_pooled_Q=r.pooled_Q,native_J_pooled_Q=(Q0+dQ)/(Nu0+dNu),
                              delta_pooled_Q=(Q0+dQ)/(Nu0+dNu)-r.pooled_Q,
                              supported_sim=r.supported_sim,supported_real=r.supported_real))
    pd.DataFrame(mask_rows).to_csv(dest/"fixed-common-mask-transport-sensitivity.csv",index=False,float_format="%.17g")
    outcome={"selected_pairs":7,"all_coordinate_data_mean_exact":True,
             "maximum_baseline_C_absolute_change":max_base_C_abs,
             "maximum_baseline_C_relative_frobenius_change":max_base_C_rel,
             "maximum_original_python_score_recovery_error":max_py_recovery,
             "maximum_original_python_mask_recovery_error":max_mask_py_recovery,
             "maximum_native_half_minus_full_gain":float(np.max(abs(frame.half_minus_full_G))),
             "maximum_native_half_minus_full_mask_gain":float(np.max(abs(maskframe.half_minus_full_G))),
             "scope":"Seven selected edge arm-object pairs only; original free-fit C/means/masks/weights retained."}
    (dest/"result.json").write_text(json.dumps(outcome,indent=2)+"\n")
    sources=[Path(__file__),HERE/"native-edge-plan.md",full/"run-gate.json",half/"run-gate.json",
             full/"mask-ledger.json",HERE/"object_scores/manifest.json",HERE/"transport/manifest.json",
             HERE/"common_mask/manifest.json",ROOT/"scripts/salt_dust_audit/flux_response.py"]
    (dest/"manifest.json").write_text(json.dumps({"inputs_sha256":{str(x.relative_to(ROOT)):sha(x) for x in sources},
                                                   "outputs_sha256":{str(x.relative_to(ROOT)):sha(x) for x in dest.iterdir() if x.name!="manifest.json"}},indent=2)+"\n")
    print(json.dumps(outcome,indent=2))


if __name__=="__main__":
    main()
