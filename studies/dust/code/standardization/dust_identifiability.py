#!/usr/bin/env python3
"""STD-02 analytic dust-geometry counterexample and latent-moment degeneracy.

These constructed populations are not fits to the Universe.
"""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, subprocess, platform
import numpy as np
import pandas as pd
from scipy.integrate import quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/standardization/dust_identifiability'; OUT.mkdir(parents=True,exist_ok=True)
K=2.5/np.log(10)

def slab(tau):
    """Absorption-only uniform-slab attenuation, including the transparent limit."""
    tau=np.asarray(tau,dtype=float)
    if np.any(~np.isfinite(tau)) or np.any(tau<0):
        raise ValueError('Slab optical depth must be finite and nonnegative')
    transmission=np.ones_like(tau)
    np.divide(-np.expm1(-tau),tau,out=transmission,where=tau!=0)
    attenuation=np.asarray(-2.5*np.log10(transmission))
    # The logarithm loses relative accuracy when transmission approaches one.
    # -ln[(1-exp(-tau))/tau] = tau/2 - tau^2/24 + tau^4/2880 + O(tau^6).
    thin=tau<1e-4
    t=tau[thin]
    attenuation[thin]=K*(t/2-t*t/24+t**4/2880)
    return attenuation
def result(rv,tau,limit=.5):
    av=float(slab(tau));ab=float(slab(tau*(1+1/rv)))
    f=min(1,limit/(K*tau*(1+1/rv)))
    # Uniform SN depth u~U[0,1], extinction=K*tau*u; deterministic A_B<limit.
    return {'microscopic_RV':rv,'tau_V':tau,'galaxy_A_V':av,'galaxy_A_B':ab,'galaxy_effective_RV':av/(ab-av),
            'SN_screen_RV':rv,'selected_fraction_AB_lt_0p5':f,'unselected_mean_SN_E':K*tau/(2*rv),
            'selected_mean_SN_E':K*tau*f/(2*rv)}

def main():
    rows=[result(rv,tau) for rv in [2.,3.1,4.] for tau in [.01,.1,.3,1,3,10]]
    pd.DataFrame(rows).to_csv(OUT/'geometry_grid.csv',index=False)
    checks=[]
    for rv in [2.,3.1,4.]:
        thin=result(rv,1e-6)
        checks.append({'RV':rv,'thin_limit_abs_error':abs(thin['galaxy_effective_RV']-rv)})
        assert abs(thin['galaxy_effective_RV']-rv)<2e-6
        curve=np.array([result(rv,t)['galaxy_effective_RV'] for t in np.geomspace(.001,100,100)])
        assert np.all(np.diff(curve)>0)
    integration_errors=[abs(quad(lambda u:np.exp(-tau*u),0,1)[0]-(-np.expm1(-tau)/tau)) for tau in [.01,.1,1,3,10]]
    assert max(integration_errors)<1e-12
    # Same two-observable covariance, different intrinsic/extinction decomposition.
    mixtures=[]
    for name,beta,rb,vi,ve,vnoise in [('A',2.,4.,.005,.005,.0225),('B',2.5,4.5,.0075,.0025,.025)]:
        covariance=np.array([[vi+ve,beta*vi+rb*ve],[beta*vi+rb*ve,beta**2*vi+rb**2*ve+vnoise]])
        mixtures.append({'name':name,'beta_intrinsic':beta,'R_B':rb,'variance_intrinsic_colour':vi,'variance_reddening':ve,
                         'variance_grey_noise':vnoise,'observable_covariance_c_m':covariance.tolist(),
                         'beta_effective':covariance[0,1]/covariance[0,0]})
    assert np.allclose(mixtures[0]['observable_covariance_c_m'],mixtures[1]['observable_covariance_c_m'])
    # Exact latent reparameterization, no random generation necessary.
    age=np.linspace(0,10,101); e=.03+.004*age; ci=-.02+.001*age
    beta=2.;rb=4.;b=-.03;k=.005
    c=ci+e; m=beta*ci+rb*e+b*age
    ep=e+k*age; cip=ci-k*age; bp=b-(rb-beta)*k
    cp=cip+ep;mp=beta*cip+rb*ep+bp*age
    exact={'original_age_slope_mag_per_Gyr':b,'alternative_age_slope_mag_per_Gyr':bp,'reddening_shift_mag_per_Gyr':k,
           'max_abs_observed_colour_difference':float(np.max(abs(c-cp))),
           'max_abs_observed_brightness_difference':float(np.max(abs(m-mp))),
           'min_alternative_reddening':float(ep.min()),
           'meaning':'Exact degeneracy if intrinsic colour and reddening population means may vary with age. Priors or extra data can break it.'}
    assert np.max(abs(m-mp))<1e-12
    paired=[dict(population='illustrative lower-opacity',**result(3.1,.1)),dict(population='illustrative higher-opacity',**result(2.,3.))]
    res={'status':'Analytic counterexample, not empirical population inference','geometry_assumptions':'Absorption-only mixed homogeneous plane-parallel slab; equal B/V spatial emissivity distribution; no scattering; SNe uniform in depth.',
         'selection_assumptions':'Detection only if SN extinction A_B<0.5 mag; intrinsic brightness, redshift, cadence held fixed. Real surveys need richer selection.',
         'grid':rows,'illustrative_opposite_trend_pair':paired,'analytic_checks':checks,'integration_max_abs_error':max(integration_errors),
         'same_observable_moments_different_dust_models':mixtures,
         'moment_warning':'Equal covariance is not equal full distributions for exponential reddening: skewness/multiband information and specified latent distributions can distinguish models.',
         'exact_age_dust_colour_reparameterization':exact}
    (OUT/'results.json').write_text(json.dumps(res,indent=2)+'\n')
    fig,ax=plt.subplots(1,2,figsize=(11,4))
    t=np.geomspace(.001,30,300)
    for rv in [2.,3.1,4.]:
        ax[0].semilogx(t,[result(rv,x)['galaxy_effective_RV'] for x in t],label=f'Sightline RV={rv}')
        ax[1].semilogx(t,[result(rv,x)['selected_fraction_AB_lt_0p5'] for x in t],label=f'Sightline RV={rv}')
    ax[0].set(xlabel='Slab optical depth tau_V',ylabel='Integrated-light effective attenuation RV',ylim=(1.5,14))
    ax[1].set(xlabel='Slab optical depth tau_V',ylabel='SN fraction passing illustrative extinction limit',ylim=(0,1.05))
    for a in ax:a.legend(fontsize=8)
    fig.suptitle('Constructed absorption-only geometry: extinction and attenuation differ')
    fig.tight_layout();fig.savefig(OUT/'dust_geometry_counterexample.png',dpi=170);plt.close(fig)
    def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
    manifest={'experiment':'STD-02','utc':datetime.now(timezone.utc).isoformat(),'code_revision':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'script_sha256':sha(Path(__file__)),'plan_sha256':sha(ROOT/'docs/experiments/standardization-plan.md'),'inputs':'Analytic constructed model specified fully in code; no observational measurements used.',
              'configuration':{'RV':[2.,3.1,4.],'tau_V':[.01,.1,.3,1,3,10],'SN_extinction_cut_A_B':.5,'seed':None},
              'environment':{'python':platform.python_version(),'numpy':np.__version__},'outputs':[str(p.relative_to(ROOT)) for p in OUT.iterdir() if p.name!='manifest.json']}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(res,indent=2))

if __name__=='__main__':main()
