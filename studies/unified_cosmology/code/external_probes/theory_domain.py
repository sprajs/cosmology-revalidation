"""Expose CAMB's analytic CPL boundary; explore, never adopt, bypass routes.

Separate processes preserve actual solver failures. Direct-field assignment is
an unsupported bypass of the Python validation and is diagnostic only. The
documented tabulated-w interface has different spline/extrapolation semantics.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
WORK=ROOT/'.work/unified-cosmology/external-probes'
RESULTS=HERE.parents[1]/'results/external_probes'


def probe(mode,w,wa,n):
    import camb
    p=camb.set_params(H0=67.36,ombh2=.02237,omch2=.12,As=2.0989e-9,
                      ns=.9649,tau=.0544,mnu=.06,nnu=3.044,
                      dark_energy_model='ppf',lmax=2500,
                      lens_potential_accuracy=1,halofit_version='mead')
    start=time.monotonic()
    out={'mode':mode,'w':w,'wa':wa,'grid_n':n}
    try:
        if mode=='analytic': p.DarkEnergy.set_params(w=w,wa=wa)
        elif mode=='direct_fields': p.DarkEnergy.w=w; p.DarkEnergy.wa=wa
        elif mode=='table':
            a=np.geomspace(1e-9,1,n)
            p.DarkEnergy.set_w_a_table(a,w+wa*(1-a))
        r=camb.get_results(p)
        cls=r.get_cmb_power_spectra(p,CMB_unit='muK')['total']
        scale=np.array([1/1101.,.1,.5,1.])
        rho,actualw=r.get_dark_energy_rho_w(scale)
        expected=scale**(-3*(1+w+wa))*np.exp(-3*wa*(1-scale))
        out.update(status='finite',rho_max_relative_error=float(np.max(abs(rho/expected-1))),
            w_max_error=float(np.max(abs(actualw-(w+wa*(1-scale))))),
            rdrag=float(r.get_derived_params()['rdrag']),
            omega_de_at_z1100=float(r.get_Omega('de',1100.)),
            TT_TE_EE_at_ells=cls[[10,100,500,1000,2000]][:,[0,3,1]].tolist())
        assert np.isfinite(cls).all()
    except Exception as e:
        out.update(status='error',exception=type(e).__name__,message=str(e))
    out['seconds']=time.monotonic()-start
    print('DOMAIN_RESULT '+json.dumps(out),flush=True)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--probe',nargs=4)
    args=parser.parse_args()
    if args.probe:
        mode,w,wa,n=args.probe
        return probe(mode,float(w),float(wa),int(n))
    rows=[]
    for w,wa in [(-.8,-.6),(-.5,.7)]:
        for mode,n in [('analytic',0),('direct_fields',0),('table',2000),('table',4000)]:
            cmd=[sys.executable,__file__,'--probe',mode,str(w),str(wa),str(n)]
            try:
                r=subprocess.run(cmd,text=True,capture_output=True,timeout=90)
                log=WORK/f'domain-{mode}-{w}-{wa}-{n}.log'
                log.write_text(r.stdout+r.stderr)
                found=[x.removeprefix('DOMAIN_RESULT ') for x in r.stdout.splitlines()
                       if x.startswith('DOMAIN_RESULT ')]
                row=json.loads(found[-1]) if found else {'status':'process_failure',
                    'mode':mode,'w':w,'wa':wa,'grid_n':n,'returncode':r.returncode}
                row['log']=str(log.relative_to(ROOT)); rows.append(row)
            except subprocess.TimeoutExpired:
                rows.append({'status':'timeout','mode':mode,'w':w,'wa':wa,'grid_n':n})
    import camb
    source=Path(camb.dark_energy.__file__)
    out={'effective_primary_CPL_domain':'w+wa <= 0 from analytic CAMB1.6.6 Python validation',
         'no_scientific_adapter_changed':True,'probes':rows,
         'scope':'Numerical feasibility only. Bypassing validation is not an established physical/theoretical validity extension, and table extrapolation and PPF perturbation assumptions need separate validation. Late-only unconstrained CPL sensitivity is separate.',
         'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
         'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (RESULTS/'theory-domain.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__': main()
