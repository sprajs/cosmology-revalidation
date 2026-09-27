#!/usr/bin/env python3
"""Common-host nebular-line and age predictive comparisons; no dust equivalence."""
import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
from pathlib import Path
import datetime, hashlib, json, sys, time
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import norm, spearmanr

ROOT = Path(__file__).resolve().parents[4]
HERE = Path(__file__).parent
WORK = ROOT / ".work/nebular-dust"
OUT = ROOT / "studies/host_ages/results/nebular_dust"
ENV = ROOT / "studies/host_ages/code/environment_validation"
sys.path.insert(0, str(ENV))
import ztf, titan
DESIGN = json.loads((HERE / "design.json").read_text())
SCALES = dict(h_alpha=2.473, h_beta=1.882, oiii_5007=1.566, nii_6584=2.039)
LINE_COLUMNS = [line + suffix for line in SCALES for suffix in ("_flux", "_flux_err", "_eqw")]
AUDIT = []
PREDICTIONS = []


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def safe(value):
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return safe(value.tolist())
    if isinstance(value, np.generic):
        return safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def load():
    for path, digest in DESIGN["inputs"].items():
        assert sha(ROOT / path) == digest, path
    b = pd.read_csv(ROOT / ".work/galaxy-validation/brightness-cohort.csv", dtype={"specObjID": str})
    l = pd.read_csv(ROOT / ".work/galaxy-validation/host-spectrum-age-ledger.csv", dtype={"specObjID": str})
    keep = ["source_id", "specObjID", "physical_host", "SN_to_fibre_arcsec", "fibre_diameter_kpc"] + LINE_COLUMNS
    d = b.merge(l[keep], on=["source_id", "specObjID", "physical_host"], how="left", validate="one_to_one", indicator=True)
    assert d._merge.eq("both").all() and len(d) == len(b)
    fold = d.fold.copy()
    saved_y = d.y.to_numpy()
    d = ztf.prepare(d)
    assert np.max(abs(saved_y - d.y)) < 1e-12
    d["fold"] = fold
    assert d.groupby("physical_host").fold.nunique().max() == 1
    return d


def line_values(d, scaled=False, groves=False):
    a, b = d.h_alpha_flux.to_numpy(), d.h_beta_flux.to_numpy()
    sa, sb = d.h_alpha_flux_err.to_numpy(), d.h_beta_flux_err.to_numpy()
    if scaled:
        sa, sb = sa * SCALES["h_alpha"], sb * SCALES["h_beta"]
    if groves:
        factor = 1 + .35 / abs(d.h_beta_eqw.to_numpy())
        b, sb = b * factor, sb * factor
    return a, b, sa, sb


def proxy(d, kind, scaled=False, rho=0., groves=False):
    a, b, sa, sb = line_values(d, scaled, groves)
    if kind == "log":
        value = np.log(a / b / 2.86)
        da, db = 1 / a, -1 / b
    else:
        r = 2.86
        den = np.sqrt(a*a + (r*b)**2 + sa*sa + (r*sb)**2)
        value = (a-r*b) / den
        da = 1/den - (a-r*b)*a/den**3
        db = -r/den - (a-r*b)*r*r*b/den**3
    variance = (da*sa)**2 + (db*sb)**2 + 2*rho*da*db*sa*sb
    assert np.isfinite(value).all() and np.isfinite(variance).all() and (variance >= 0).all()
    q = d.copy()
    q["nebular"] = value
    q["nebular_variance"] = variance
    return q


def valid_lines(d):
    x = d[["h_alpha_flux", "h_beta_flux", "h_alpha_flux_err", "h_beta_flux_err"]].to_numpy()
    return np.isfinite(x).all(axis=1) & (x[:, 2:] > 0).all(axis=1) & ~np.isin(x[:, :2], [-9999., 9999.]).any(axis=1)


def masks(d):
    valid = valid_lines(d)
    ha, hb = d.h_alpha_flux / d.h_alpha_flux_err, d.h_beta_flux / d.h_beta_flux_err
    sf = valid.copy()
    sf_scaled = valid.copy()
    for line, scale in SCALES.items():
        snr = d[line + "_flux"] / d[line + "_flux_err"]
        sf &= (snr > 3) & (d[line + "_flux_err"] > 0)
        sf_scaled &= (snr > 3 * scale) & (d[line + "_flux_err"] > 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        x = np.log10(d.nii_6584_flux / d.h_alpha_flux)
        y = np.log10(d.oiii_5007_flux / d.h_beta_flux)
    bpt = (x < .05) & (y < .61/(x-.05)+1.3)
    return {"balmer_snr3": valid & (ha>3) & (hb>3),
            "balmer_snr5": valid & (ha>5) & (hb>5),
            "starforming_bpt": sf & bpt,
            "all_valid_lines": valid,
            "scaled_error_snr3": valid & (ha>3*SCALES["h_alpha"]) & (hb>3*SCALES["h_beta"]),
            "scaled_error_starforming_bpt": sf_scaled & bpt}


def matrix(d, names):
    return np.column_stack([d.nebular.to_numpy() if name == "nebular" else titan.matrix(d, [name])[:, 0] for name in names])


def variance(d, b, names, age_error=False):
    v = ztf.variance(d, b, names)
    if "nebular" in names:
        v += b[names.index("nebular")]**2 * d.nebular_variance.to_numpy()
    if age_error and "host_age" in names:
        v += b[names.index("host_age")]**2 * ((d.mass_weighted_age_84 - d.mass_weighted_age_16).to_numpy()/2)**2
    return v


def independent_variance(d, b, names, age_error=False):
    """Expand the light-curve covariance directly from native fitted columns."""
    p = dict(zip(names, b))
    dx = p["x1"] + 2*p.get("x12", 0)*d.x1.to_numpy()
    dc = p["c"] + 2*p.get("c2", 0)*d.c.to_numpy()
    dm = -2.5 / np.log(10) / d.x0.to_numpy()
    v = (dm*d.x0_err.to_numpy())**2 + (dx*d.x1_err.to_numpy())**2 + (dc*d.c_err.to_numpy())**2
    v -= 2*dm*dx*d.cov_x0_x1.to_numpy() + 2*dm*dc*d.cov_x0_c.to_numpy()
    v += 2*dx*dc*d.cov_x1_c.to_numpy() + d.mu_z_variance.to_numpy()
    if "localcolour" in names:
        v += p["localcolour"]**2*d.restframe_gz_err_local.to_numpy()**2
    if "mass_step" in names:
        probability = norm.cdf((d.mass_global.to_numpy()-10)/np.maximum(d.mass_err_global.to_numpy(),1e-6))
        v += p["mass_step"]**2*probability*(1-probability)
    if "nebular" in names:
        v += p["nebular"]**2*d.nebular_variance.to_numpy()
    if age_error and "host_age" in names:
        v += p["host_age"]**2*((d.mass_weighted_age_84-d.mass_weighted_age_16).to_numpy()/2)**2
    return v


def fit_checked(d, names, age_error, label, audit=True):
    x = matrix(d, names)
    if len(d) < 2*len(names)+5 or np.linalg.matrix_rank(x) != len(names):
        raise ValueError(f"Insufficient rank/support: n={len(d)}, p={len(names)}, rank={np.linalg.matrix_rank(x)}")
    vf = lambda dd, bb, nn: variance(dd, bb, nn, age_error)
    result, b, cov, sig = ztf.fit(d, names, matrix, vf)
    if audit:
        y = d.y.to_numpy()
        def objective(coeff, sigma):
            v = independent_variance(d, coeff, names, age_error) + sigma*sigma
            return float(np.sum(np.log(2*np.pi*v)+(y-x@coeff)**2/v))
        direct = objective(b, sig)
        verror = float(np.max(abs(independent_variance(d,b,names,age_error)-vf(d,b,names))))
        assert verror < 1e-12 and abs(direct-result["minus2loglike"]) < 1e-9
        # Condition parameters using data-column scales; independent optimizer
        # receives a perturbed solution, not an implementation-specific Hessian.
        scale = 1/np.maximum(np.std(x, axis=0), .001)
        scale[0] = 1.
        point = np.r_[b/scale, np.log(sig)]
        rng = np.random.default_rng(DESIGN["seed"] + len(AUDIT))
        initial = point + rng.normal(size=len(point))*.002
        cost = lambda v: objective(v[:-1]*scale, np.exp(v[-1]))
        other = minimize(cost, initial, method="SLSQP", bounds=[(None,None)]*len(b)+[(np.log(1e-5),0)],
                         options={"ftol":1e-10,"maxiter":2000})
        delta = float(other.fun-direct)
        verified = bool(other.success and abs(delta)<1e-5)
        AUDIT.append({"model":label,"n":len(d),"variance_max_error":verror,"direct_minus2loglike_difference":direct-result["minus2loglike"],
                      "independent_optimizer_success":bool(other.success),"independent_optimizer_message":str(other.message),
                      "independent_minus2loglike_difference":delta,"passed":verified})
        if not verified:
            raise RuntimeError("Independent optimizer check failed: "+str(AUDIT[-1]))
    if "host_age" in names:
        value,se = result["parameters"]["host_age"],result["standard_errors"]["host_age"]
        result["conditional_age_slope_95"] = [value-1.96*se,value+1.96*se]
        delta=.03/se
        result["normal_power_for_abs003_conditional"] = float(norm.cdf(-1.96-delta)+norm.sf(1.96-delta))
    return result,b,cov,sig


def run_comparison(d, cohort, controls, setting="formal", cv=True, age_error=False):
    base = DESIGN["base"] + (["titan_mass","titan_AV","titan_metallicity"] if controls else [])
    models = {"base":base,"age":base+["host_age"],"nebular":base+["nebular"],"nebular_age":base+["nebular","host_age"]}
    output = {"n":len(d),"hosts":int(d.physical_host.nunique()),"fold_counts":d.groupby("fold").size().to_dict(),
              "cohort":cohort,"SED_controls":controls,"setting":setting,"models":models,"fits":{},"gates":[]}
    scores={}
    for name,names in models.items():
        label=f"{cohort}:{controls}:{setting}:{name}"
        try:
            result,b,cov,sig=fit_checked(d,names,age_error,label)
        except Exception as exc:
            output["gates"].append({"model":name,"stage":"fullfit","error":str(exc)})
            continue
        output["fits"][name]=result
        if cv:
            predictions=np.full(len(d),np.nan);vpred=np.full(len(d),np.nan)
            try:
                for fold in range(5):
                    train=d.fold!=fold;test=~train
                    ff,bb,cc,ss=fit_checked(d[train],names,age_error,label+f":fold{fold}",audit=False)
                    xx=matrix(d[test],names)
                    predictions[test]=xx@bb
                    vpred[test]=variance(d[test],bb,names,age_error)+ss*ss+np.einsum("ni,ij,nj->n",xx,cc,xx)
                score=-.5*(np.log(2*np.pi*vpred)+(d.y.to_numpy()-predictions)**2/vpred)
                assert np.isfinite(score).all() and (vpred>0).all()
                scores[name]=score
                result["heldout_log_score"]=float(score.sum())
                result["heldout_RMSE_mag"]=float(np.sqrt(np.mean((d.y.to_numpy()-predictions)**2)))
                for i,(_,row) in enumerate(d.iterrows()):
                    PREDICTIONS.append(dict(cohort=cohort,SED_controls=controls,setting=setting,model=name,source_id=row.source_id,
                                            physical_host=row.physical_host,fold=int(row.fold),prediction=predictions[i],variance=vpred[i],score=score[i],y=float(row.y)))
            except Exception as exc:
                output["gates"].append({"model":name,"stage":"heldout","error":str(exc)})
    output["comparisons"]={}
    rng=np.random.default_rng(DESIGN["seed"])
    for added,old in [("age","base"),("nebular","base"),("nebular_age","nebular"),("nebular_age","age")]:
        pair={}
        if added in output["fits"] and old in output["fits"]:
            pair["minus2loglike_improvement"]=output["fits"][old]["minus2loglike"]-output["fits"][added]["minus2loglike"]
        if added in scores and old in scores:
            delta=scores[added]-scores[old]
            units=np.array([delta[d.physical_host.eq(host)].sum() for host in sorted(set(d.physical_host))])
            boot=rng.choice(units,(2000,len(units)),replace=True).sum(axis=1)
            pair.update(delta_heldout_log_score=float(delta.sum()),host_bootstrap95_conditional_on_folds=np.quantile(boot,[.025,.975]).tolist())
        output["comparisons"][added+"_minus_"+old]=pair
    x=matrix(d,base)
    ages=d.mass_weighted_age_50.to_numpy()
    ar=ages-x@np.linalg.lstsq(x,ages,rcond=None)[0]
    xr=matrix(d,base+["nebular"])
    arr=ages-xr@np.linalg.lstsq(xr,ages,rcond=None)[0]
    output["age_support"]={"sd_after_base_Gyr":float(np.std(ar)),"sd_after_base_and_nebular_Gyr":float(np.std(arr)),
                           "median_marginal_age_error_Gyr":float(np.median((d.mass_weighted_age_84-d.mass_weighted_age_16)/2)),
                           "nebular_age_spearman":float(spearmanr(d.nebular,ages).statistic),
                           "outside_fibre_radius":int((d.SN_to_fibre_arcsec>1.5).sum())}
    return output


def main():
    start=time.perf_counter()
    WORK.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    d=load();selection=masks(d)
    manifest=d[["source_id","specObjID","physical_host","fold"]+LINE_COLUMNS].copy()
    for name,mask in selection.items():manifest[name]=mask
    manifest.to_csv(WORK/"cohort-ledger.csv",index=False)
    runs=[]
    for name in ["balmer_snr3","balmer_snr5","starforming_bpt","all_valid_lines"]:
        dd=proxy(d[selection[name]].copy(),"contrast" if name=="all_valid_lines" else "log")
        for controls in [False,True]:runs.append(run_comparison(dd,name,controls))
    # Fixed primary common cohort for changes to measurement-error assumptions.
    for setting,rho,scaled,ageerr,groves in [
        ("rho_minus05",-.5,False,False,False),("rho_plus05",.5,False,False,False),
        ("historical_error_scale",0,True,False,False),("independent_age_error",0,False,True,False),
        ("conditional_Hbeta_EW035",0,False,False,True)]:
        mask=selection["balmer_snr3"].copy()
        if groves:mask &= d.h_beta_eqw<0
        dd=proxy(d[mask].copy(),"log",scaled,rho,groves)
        for controls in [False,True]:runs.append(run_comparison(dd,"balmer_snr3",controls,setting,cv=False,age_error=ageerr))
    for name in ["scaled_error_snr3","scaled_error_starforming_bpt"]:
        dd=proxy(d[selection[name]].copy(),"log",scaled=True)
        for controls in [False,True]:runs.append(run_comparison(dd,name,controls,setting="historical_error_scale",cv=True))
    pd.DataFrame(PREDICTIONS).to_csv(WORK/"heldout.csv",index=False)
    files=[WORK/"cohort-ledger.csv",WORK/"heldout.csv"]
    summary={"created_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"design_sha256":sha(HERE/"design.json"),
             "input_sha256":DESIGN["inputs"],"source_context_sha256":sha(OUT/"sources.json"),
             "code_sha256":sha(Path(__file__)),"dependency_sha256":{str(p.relative_to(ROOT)):sha(p) for p in [ENV/"ztf.py",ENV/"titan.py",ENV/"analyse.py"]},
             "parent_objects":len(d),"selection_counts":{k:int(v.sum()) for k,v in selection.items()},
             "invalid_line_rows":d.loc[~selection["all_valid_lines"],["source_id","h_alpha_flux","h_alpha_flux_err","h_beta_flux","h_beta_flux_err"]].to_dict("records"),
             "runs":runs,"independent_fit_checks":AUDIT,"output_sha256":{str(p.relative_to(ROOT)):sha(p) for p in files},
             "seconds":time.perf_counter()-start,
             "limits":["Nebular line fluxes already use a model for continuum subtraction and a foreground correction; no second foreground correction applied.",
                       "Diagonal line-error propagation, unknown line covariance and continuum systematics; low-SNR contrast has no physical attenuation interpretation.",
                       "SED age/dust/mass/metallicity joint likelihood absent; median conditioning and independent-age-error sensitivity do not solve shared errors.",
                       "Central fibres and selected emission-line galaxies are not all SN environments; possible mediator conditioning is not a causal test.",
                       "Exact individual ZTF x0 blinding transform unknown. These conditional predictive associations do not establish an unblinded cosmological correction."]}
    (OUT/"summary.json").write_text(json.dumps(safe(summary),indent=2,allow_nan=False)+"\n")
    print(json.dumps({"counts":summary["selection_counts"],"runs":len(runs),"fit_checks":len(AUDIT),"failed_checks":sum(not a["passed"] for a in AUDIT),"gates":sum(len(r["gates"]) for r in runs),"seconds":summary["seconds"]},indent=2))


if __name__=="__main__":main()
