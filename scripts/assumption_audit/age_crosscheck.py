"""Independent native-FSPS-normalization quadrature on observed R19 draws."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.integrate import quad
from astropy.cosmology import FlatLambdaCDM
from pantheon_age_semantics import moments, ROOT


def main():
    p=ROOT/'runs/assumption_audit/pantheon/sources/SN5916_campbellG_chain.tsv'
    a=pd.read_csv(p,sep='\t',comment='#',header=None,nrows=10000).to_numpy()
    rng=np.random.default_rng(192106)
    indices=rng.choice(len(a),24,replace=False)
    cosmic=float(FlatLambdaCDM(H0=70,Om0=.27).age(.1740905).value)
    rows=[]
    for i in indices:
        tau,start,trans,phi=a[i,2:6]
        s=np.tan(phi);duration=cosmic-start;turn=trans-start
        # Native FSPS early normalization; no analytic moment decomposition.
        k=turn/tau*np.exp(-turn/tau)
        def sfr(t):
            return t/tau*np.exp(-t/tau) if t<=turn else max(k*(1+s*(t-turn)),0)
        breaks=[turn] if 0<turn<duration else []
        if s<0 and turn<turn-1/s<duration:breaks.append(turn-1/s)
        mass=quad(sfr,0,duration,points=breaks,epsabs=1e-11,epsrel=1e-10)[0]
        numerator=quad(lambda t:(duration-t)*sfr(t),0,duration,points=breaks,epsabs=1e-11,epsrel=1e-10)[0]
        direct=numerator/mass
        analytic=float(moments(tau,start,trans,s,cosmic,True)[0])
        rows.append({'row':int(i),'native_FSPS_SFH_direct_mean_age':direct,
                     'analytic_mean_age':analytic,'archived_age':float(a[i,7]),
                     'difference_direct_minus_analytic':direct-analytic})
    error=max(abs(r['difference_direct_minus_analytic']) for r in rows)
    assert error<1e-8
    files=[p,Path(__file__),Path(__file__).with_name('pantheon_age_semantics.py')]
    result={'method':'Native FSPS SFR normalization and direct lookback-age quadrature on 24 fixed-seed observed posterior rows; no SED refit',
            'seed':192106,'rows':rows,'max_abs_difference_Gyr':error,
            'sources_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
            'limitations':'Bounded numerical/source review; does not estimate sample-wide cosmological effect.'}
    target=ROOT/'runs/assumption_audit/age-independent-crosscheck.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(error)


if __name__=='__main__':main()
