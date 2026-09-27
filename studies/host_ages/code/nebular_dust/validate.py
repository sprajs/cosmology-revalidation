#!/usr/bin/env python3
"""Independent source, selection, gradient and heldout-arithmetic validation."""
from pathlib import Path
import ast,datetime,hashlib,json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).parent
WORK=ROOT/".work/nebular-dust"
OUT=ROOT/"studies/host_ages/results/nebular_dust"

def sha(path):
    with Path(path).open("rb") as stream:return hashlib.file_digest(stream,"sha256").hexdigest()

def main():
    summary=json.loads((OUT/"summary.json").read_text())
    design=json.loads((HERE/"design.json").read_text())
    assert sha(HERE/"design.json")==summary["design_sha256"]
    assert sha(HERE/"analyze.py")==summary["code_sha256"]
    assert sha(OUT/"sources.json")==summary["source_context_sha256"]
    for kind in ["input_sha256","dependency_sha256","output_sha256"]:
        for path,digest in summary[kind].items():assert sha(ROOT/path)==digest,path
    for source in json.loads((OUT/"sources.json").read_text())["files"]:
        assert sha(ROOT/source["path"])==source["sha256"]
    original=pd.read_csv(ROOT/".work/galaxy-validation/host-spectrum-age-ledger.csv",dtype={"specObjID":str})
    base=pd.read_csv(ROOT/".work/galaxy-validation/brightness-cohort.csv",dtype={"specObjID":str})
    ledger=pd.read_csv(WORK/"cohort-ledger.csv",dtype={"specObjID":str})
    assert ledger.source_id.tolist()==base.source_id.tolist()
    q=ledger.merge(original,on=["source_id","specObjID","physical_host"],validate="one_to_one",suffixes=("_saved","_source"))
    for line in ["h_alpha","h_beta","oiii_5007","nii_6584"]:
        for suffix in ["_flux","_flux_err","_eqw"]:
            name=line+suffix
            assert np.allclose(q[name+"_saved"],q[name+"_source"],rtol=0,atol=0,equal_nan=True)
    assert ledger.fold.tolist()==base.fold.tolist()
    assert ledger.groupby("physical_host").fold.nunique().max()==1
    a,b=ledger.h_alpha_flux.to_numpy(),ledger.h_beta_flux.to_numpy()
    sa,sb=ledger.h_alpha_flux_err.to_numpy(),ledger.h_beta_flux_err.to_numpy()
    valid=np.isfinite(np.c_[a,b,sa,sb]).all(axis=1)&(sa>0)&(sb>0)&~np.isin(np.c_[a,b],[-9999.,9999.]).any(axis=1)
    assert np.array_equal(valid,ledger.all_valid_lines)
    assert np.array_equal(valid&(a>3*sa)&(b>3*sb),ledger.balmer_snr3)
    assert np.array_equal(valid&(a>5*sa)&(b>5*sb),ledger.balmer_snr5)
    assert {name:int(ledger[name].sum()) for name in summary["selection_counts"]}==summary["selection_counts"]
    # Contrast is invariant to flux units and does not discard signed values.
    a,b,sa,sb=[v[valid] for v in [a,b,sa,sb]]
    r=2.86
    def contrast(x,y):return (x-r*y)/np.sqrt(x*x+(r*y)**2+sa*sa+(r*sb)**2)
    den=np.sqrt(a*a+(r*b)**2+sa*sa+(r*sb)**2)
    da=1/den-(a-r*b)*a/den**3
    db=-r/den-(a-r*b)*r*r*b/den**3
    stepsa=1e-5*np.maximum(abs(a),sa);stepsb=1e-5*np.maximum(abs(b),sb)
    finite_a=(contrast(a+stepsa,b)-contrast(a-stepsa,b))/(2*stepsa)
    finite_b=(contrast(a,b+stepsb)-contrast(a,b-stepsb))/(2*stepsb)
    gradient_error=max(float(np.max(abs(finite_a-da)/np.maximum(abs(da),1e-12))),
                       float(np.max(abs(finite_b-db)/np.maximum(abs(db),1e-12))))
    assert gradient_error<1e-7
    scale=1e-17
    unit_error=float(np.max(abs(contrast(a,b)-(a*scale-r*b*scale)/np.sqrt((a*scale)**2+(r*b*scale)**2+(sa*scale)**2+(r*sb*scale)**2))))
    assert unit_error<1e-14
    # This is a local propagation diagnostic, not a latent-line posterior:
    # Gaussian draws are centred at measured values; no S/N cut on draws.
    rng=np.random.default_rng(design["seed"]+81)
    aa=a+rng.normal(size=(10000,len(a)))*sa
    bb=b+rng.normal(size=(10000,len(b)))*sb
    values=(aa-r*bb)/np.sqrt(aa*aa+(r*bb)**2+sa*sa+(r*sb)**2)
    delta_variance=(da*sa)**2+(db*sb)**2
    mc_variance=values.var(axis=0,ddof=1)
    ratios=mc_variance/delta_variance
    low=(a/sa<=3)|(b/sb<=3)
    mc={"scope":"Gaussian perturbations around measurements; tests first-order propagation only, not a flux-population posterior or repeated-sample coverage.",
        "draws":10000,"ratio_MC_variance_to_delta_variance_all_quantiles":np.quantile(ratios,[.05,.5,.95]).tolist(),
        "weak_or_nondetected_objects":int(low.sum()),"weak_ratio_quantiles":np.quantile(ratios[low],[.05,.5,.95]).tolist(),
        "contrast_bias_max_abs":float(abs(values.mean(axis=0)-contrast(a,b)).max())}
    heldout=pd.read_csv(WORK/"heldout.csv")
    error=float(np.max(abs(heldout.score+.5*(np.log(2*np.pi*heldout.variance)+(heldout.y-heldout.prediction)**2/heldout.variance))))
    assert error<1e-12 and (heldout.variance>0).all()
    score_count=0
    for run in summary["runs"]:
        assert not run["gates"],run["gates"]
        for name,fit in run["fits"].items():
            assert fit["valid_minimum"]
            cov=np.array(fit["coefficient_covariance"])
            assert np.allclose(cov,cov.T,atol=1e-12) and np.linalg.eigvalsh(cov).min()>0
            if "heldout_log_score" not in fit:continue
            h=heldout[(heldout.cohort==run["cohort"])&(heldout.SED_controls==run["SED_controls"])&(heldout.setting==run["setting"])&(heldout.model==name)]
            assert len(h)==run["n"] and h.source_id.is_unique
            assert h.groupby("physical_host").fold.nunique().max()==1
            assert abs(h.score.sum()-fit["heldout_log_score"])<1e-9
            assert abs(np.sqrt(np.mean((h.y-h.prediction)**2))-fit["heldout_RMSE_mag"])<1e-12
            score_count+=1
    checks=summary["independent_fit_checks"]
    assert all(c["passed"] and c["independent_optimizer_success"] for c in checks)
    for path in HERE.glob("*.py"):ast.parse(path.read_text(),filename=str(path))
    result={"created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"passed":True,
            "exact_source_associations":len(ledger),"physical_hosts":int(ledger.physical_host.nunique()),
            "signed_Balmer_objects":int(valid.sum()),"selected_nonpositive_Balmer_objects":int(((a<=0)|(b<=0)).sum()),
            "contrast_gradient_max_relative_error":gradient_error,"contrast_flux_unit_identity_error":unit_error,
            "nonlinear_uncertainty_diagnostic":mc,"heldout_score_max_identity_error":error,"heldout_model_checks":score_count,
            "independent_full_fit_checks":len(checks),"independent_objective_max_difference":max(abs(c["independent_minus2loglike_difference"]) for c in checks),
            "failed_fit_or_heldout_gates":sum(len(r["gates"]) for r in summary["runs"]),
            "input_sha256":summary["input_sha256"],"summary_sha256":sha(OUT/"summary.json"),
            "code_sha256":{str(p.relative_to(ROOT)):sha(p) for p in HERE.glob("*") if p.is_file()}}
    (OUT/"validation.json").write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":main()
