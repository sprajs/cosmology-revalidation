#!/usr/bin/env python3
"""Independent age/SFH semantics audit; never changes published age results.

Requires the pinned historical FSPS source and original R19 circle archive in
runs/assumption_audit/pantheon/sources. Algebra below integrates *formed* mass,
the same estimand as MC-Age. It does not refit photometry or surviving mass.
"""
from pathlib import Path
import ast
import hashlib
import json
import logging
import tarfile

import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import numpy as np
import pandas as pd
from scipy import integrate
from scipy.special import gammainc
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/assumption_audit/pantheon"
SRC = OUT / "sources"
MC = ROOT / "sources/repos/benjaminrose__mc-age/calculateAge.py"


def moments(tau, start, transition, slope, cosmic_age, relative):
    """Exact piecewise moments; relative=True implements historical FSPS.

    The early SFR is t exp(-t/tau). FSPS uses (t/tau) exp(-t/tau);
    multiplying *both* its phases by tau cancels in the age moment ratio.
    Its late SFR must therefore be K*(1+slope*u), K=T exp(-T/tau).
    MC-Age instead uses K+slope*u, u=t-T.
    """
    tau, start, transition, slope = np.broadcast_arrays(tau, start, transition, slope)
    age = cosmic_age - start
    turn = transition - start
    k = turn * np.exp(-turn / tau)
    gradient = slope * k if relative else slope
    duration = np.maximum(0., age - turn)
    stop = np.full_like(duration, np.inf, dtype=float)
    np.divide(-k, gradient, out=stop, where=gradient < 0)
    duration = np.where(gradient < 0, np.minimum(duration, stop), duration)
    upper = np.minimum(age, turn)
    early_mass = tau**2 * gammainc(2, upper / tau)
    early_moment = 2 * tau**3 * gammainc(3, upper / tau)
    late_mass = k*duration + gradient*duration**2/2
    late_moment = turn*late_mass + k*duration**2/2 + gradient*duration**3/3
    total = early_mass + late_mass
    return age-(early_moment+late_moment)/total, late_mass/total


def source_functions():
    """Execute only inspected pure functions from the preserved MC-Age file."""
    names = {"integrate_age", "star_formation_gupta", "t_star_formation_gupta", "lnprior", "runFSPS"}
    tree = ast.parse(MC.read_text())
    selected = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    ns = dict(np=np, integrate=integrate, logging=logging, u=u,
              cosmo=FlatLambdaCDM(H0=70, Om0=.27), ramp=lambda x: max(0., x))
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(MC), "exec"), ns)
    return ns


def verify_quadrature(tau, start, transition, slope, cosmic_age, relative):
    """Independent numerical integration with all breaks explicitly supplied."""
    age, turn = cosmic_age-start, transition-start
    k = turn*np.exp(-turn/tau)
    g = slope*k if relative else slope
    def sfr(t):
        return t*np.exp(-t/tau) if t <= turn else max(0., k+g*(t-turn))
    points = [turn] if 0 < turn < age else []
    if g < 0 and turn < turn-k/g < age:
        points.append(turn-k/g)
    mass = integrate.quad(sfr, 0, age, points=points, epsabs=1e-12, epsrel=1e-11)[0]
    first = integrate.quad(lambda t:t*sfr(t), 0, age, points=points, epsabs=1e-12, epsrel=1e-11)[0]
    return age-first/mass


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    ns = source_functions()
    cosmic_age = float(ns["cosmo"].age(.05).value)
    # R19 Table 5, with a separate interior-prior burst case; model 8's tau=.1
    # and slope=20 are at the stated prior boundary and are not the main proof.
    cases = [(1,.5,1.5,9.,-1.), (2,.5,1.5,9.,15.), (3,7.,3.,10.,15.),
             (4,7.,3.,13.,0.), (5,.5,1.5,9.,-1.), (6,7.,3.,10.,15.),
             (7,.5,1.5,6.,15.), (8,.1,8.,12.,20.), (9,.2,7.5,11.,10.)]
    rows=[]
    for ident, tau, start, transition, slope in cases:
        old, oldfraction = moments(tau,start,transition,slope,cosmic_age,False)
        new, newfraction = moments(tau,start,transition,slope,cosmic_age,True)
        original = ns["integrate_age"](tau,start,transition,slope,.05)[0]
        prior = ns["lnprior"]([-.5,.3,tau,start,transition,np.arctan(slope),-25],.05)
        spy=type("SP",(),{"params":{},"get_mags":lambda self,**kw:kw})()
        ns["runFSPS"](spy,.05,-.5,.3,tau,start,transition,np.arctan(slope))
        assert abs(spy.params["sf_slope"]-slope) < 1e-12
        qold=verify_quadrature(tau,start,transition,slope,cosmic_age,False)
        qnew=verify_quadrature(tau,start,transition,slope,cosmic_age,True)
        assert abs(qold-old)<1e-8 and abs(qnew-new)<1e-8
        # The original does not split quadrature at the transition; its late
        # zero-slope test differs by 4.3e-5 Gyr. Preserve this diagnostic.
        assert abs(original-old)<1e-4
        rows.append(dict(case=ident,tau=tau,start=start,transition=transition,slope=slope,
                         valid_under_pinned_prior=bool(np.isfinite(prior)),
                         mc_age_gyr=float(old),fsps_sfh_age_gyr=float(new),
                         fsps_minus_mc_age_gyr=float(new-old),
                         mc_linear_formed_mass_fraction=float(oldfraction),
                         fsps_linear_formed_mass_fraction=float(newfraction),
                         actual_source_minus_analytic_gyr=float(original-old)))
    pd.DataFrame(rows).to_csv(OUT/"sfh-semantics-cases.csv",index=False)
    # Check original authors' archived posterior age columns, not synthetic
    # chains created here. Historic Astropy defaults included Tcmb0=2.725 K.
    old_cosmic_age=float(FlatLambdaCDM(H0=70,Om0=.27,Tcmb0=2.725).age(.05).value)
    archive=SRC/"rose2019-circle.tar.gz"
    chains=[]
    with tarfile.open(archive) as bundle:
        for member in bundle:
            if not member.isfile() or not member.name.endswith("_chain.tsv"):
                continue
            d=pd.read_csv(bundle.extractfile(member),sep="\t",comment="#",header=None).to_numpy()
            # Every stored row; no burn-in/quality selection is inferred here.
            tau,start,transition,phi=d[:,2:6].T
            old,_=moments(tau,start,transition,np.tan(phi),old_cosmic_age,False)
            new,_=moments(tau,start,transition,np.tan(phi),old_cosmic_age,True)
            chains.append(dict(member=member.name,rows=len(d),
                               archived_age_quantiles=np.quantile(d[:,7],[.16,.5,.84]).tolist(),
                               code_age_quantiles=np.quantile(old,[.16,.5,.84]).tolist(),
                               fsps_consistent_age_quantiles=np.quantile(new,[.16,.5,.84]).tolist(),
                               corrected_minus_archived_median=float(np.median(new)-np.median(d[:,7])),
                               abs_archived_minus_code_p50_p99_max=np.quantile(abs(d[:,7]-old),[.5,.99,1]).tolist()))
    # Arbitrary first regular file of the observed global-host archive, streamed
    # without selecting a host by its effect size. Unlike circle, this archive's
    # age column agrees with the later zero-radiation Astropy default.
    observed_path=SRC/"SN5916_campbellG_chain.tsv"
    observed=None
    if observed_path.exists():
        photometry=ROOT/"sources/repos/benjaminrose__mc-age/data/campbell_global.tsv"
        z=float(pd.read_csv(photometry,sep="\t").set_index("SNID").loc[5916,"redshift"])
        d=pd.read_csv(observed_path,sep="\t",comment="#",header=None).to_numpy()
        tau,start,transition,phi=d[:,2:6].T
        ca=float(ns["cosmo"].age(z).value)
        old,_=moments(tau,start,transition,np.tan(phi),ca,False)
        new,_=moments(tau,start,transition,np.tan(phi),ca,True)
        observed=dict(member="campbellG/SN5916_campbellG_chain.tsv",rows=len(d),z=z,
            cosmology_Tcmb0_K=0.,archived_age_quantiles=np.quantile(d[:,7],[.16,.5,.84]).tolist(),
            fsps_consistent_age_quantiles=np.quantile(new,[.16,.5,.84]).tolist(),
            corrected_minus_archived_median=float(np.median(new)-np.median(d[:,7])),
            abs_archived_minus_code_p50_p99_max=np.quantile(abs(d[:,7]-old),[.5,.99,1]).tolist(),
            age_logmetallicity_dust2_correlation=np.corrcoef(d[:,[7,0,1]].T).tolist(),
            limitation="One arbitrarily encountered original R19 host; not the revised C25 age posterior or a population effect estimate.")
    # Quantile convention and Gaussian support are approximation diagnostics.
    age_support=[]
    for policy in ["G11_first","R19_first"]:
        d=pd.read_csv(ROOT/f"data/derived/age_signal/matched_{policy}.csv")
        d=d[(d.zHD>.06)&(d.zHD<.42)]
        maxage=FlatLambdaCDM(H0=70,Om0=.27).age(d.zHD.to_numpy()).value
        negative=norm.cdf(-d.age/d.age_err)
        too_old=norm.sf((maxage-d.age)/d.age_err)
        age_support.append(dict(policy=policy,n=len(d),
            summed_gaussian_negative_probability=float(negative.sum()),
            summed_gaussian_older_than_cosmos_probability=float(too_old.sum()),
            n_more_than_5pct_unphysical=int(((negative+too_old)>.05).sum()),
            gaussian_sigma_for_claimed_86p4_quantile=float(norm.ppf(.864))))
    result=dict(
        status="Confirmed source-code mismatch; C25 deployment not established",
        no_photometry_refit=True,
        current_astropy_cosmic_age_zp05_gyr=cosmic_age,
        historical_default_approx_cosmic_age_zp05_gyr=old_cosmic_age,
        cases=rows,original_archived_circle_chains=chains,original_observed_host=observed,
        gaussian_age_summary_diagnostics=age_support,
        interpretation="FSPS slope is fractional per Gyr; MC-Age integrates it as an additive SFR slope. Age postprocessing differs from the fitted SED model. Archived circle chains establish use in released R19 validation outputs, not by themselves in C25 tables.")
    (OUT/"age-semantics.json").write_text(json.dumps(result,indent=2)+"\n")
    inputs=[MC,archive,ROOT/"papers/text/1902.01433v1.txt",ROOT/"papers/text/chung2025-published.txt"]+list(SRC.glob("fsps-*.f90"))
    if observed_path.exists(): inputs.extend([observed_path,photometry])
    manifest=dict(script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
                  versions={"numpy":np.__version__,"pandas":pd.__version__},
                  fsps_revision="ae31b2f63d865354ce944e5c22eba6e93e01e67d",
                  mc_age_revision="92713be96a89da991fe53bffcc596a5c0942fc37")
    (OUT/"age-semantics-manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
