"""Nonlinear response to a frozen, hypothetical four-band recalibration.

This tests a numerical response, not the physical origin or truth of a correction.
Observed flux and its full frozen covariance transform together.
"""
from pathlib import Path
import sys
import json
import hashlib
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.salt_dust_audit.flux_response import build_model,K


def main():
    base=ROOT/'runs/research_2026_09_26/astra_design'
    out=ROOT/'runs/research_2026_09_26/calibration_pattern_refit'
    out.mkdir(exist_ok=False)
    source=Path(__file__).read_bytes();(out/'executed_source.py').write_bytes(source)
    coefpath=base/'validation1020/frozen-discovery-coefficients.npz'
    coef=np.load(coefpath);dmag=-coef['gauge_griz']@coef['basis_mean']
    meta=pd.read_csv(ROOT/'runs/salt_dust_audit/flux_response/selected_objects.csv',dtype={'CID':str})
    meta=meta[meta.PROB_SNNV19>.999].sort_values(['zHEL','CID'])
    chosen=meta.iloc[np.rint(np.linspace(0,len(meta)-1,6)).astype(int)]
    chosen[['CID','zHEL']].to_csv(out/'cohort.csv',index=False)
    (out/'protocol.txt').write_text('Six fixed redshift ranks of discovery43; frozen exactC ridge.02 observer coefficients; apply opposite signed magnitude vector to observed flux, scaling full covariance consistently; same SALT surface and no BBC/selection. Compare native nonlinear pair to a fresh local Fisher response at its baseline optimum. No validation1020 outcomes used.\n')
    model,bands,paths,zp=build_model();offsets={str(r['Filter Name'])[-1]:float(r['Primary Mag']) for r in zp}
    D=np.array([1,.16087,-3.1178,0]);steps=np.array([1e-4,1e-3,1e-4,.01]);rows=[];inputs=[coefpath]
    for row in chosen.itertuples():
        p=base/f'exact43/objectives/objective_{row.CID}.npz';inputs.append(p);o=np.load(p)
        t=o['MJD'];band=o['band'];bs=np.array([bands[b] for b in band],dtype=object)
        fc=np.array([10**(-.4*(.27+offsets[b])) for b in band])
        correction=np.array([dmag['griz'.index(b)] for b in band]);scale=np.exp(-K*correction)
        x0,x1,c,t0=o['parameters_x0_x1_c_t0'];z=float(o['zHEL'][0]);ebv=float(o['MWEBV'][0])
        y=o['data_flux'];L=np.linalg.cholesky(o['frozen_flux_covariance'])
        def flux(theta):
            model.set(z=z,t0=t0+theta[3],x0=x0*np.exp(-K*theta[0]),x1=x1+theta[1],c=c+theta[2],
                mwebv=ebv,mwrv=3.1,hostebv=0,hostrv=3.1)
            return fc*model.bandflux(bs,t,zp=27.5,zpsys='ab')
        def fit(s,start):
            # Cholesky(S C S) = S L for diagonal positive S.
            fun=lambda theta:solve_triangular(L,flux(theta)/s-y,lower=True)
            return least_squares(fun,start,xtol=1e-12,ftol=1e-12,gtol=1e-10,max_nfev=400,
                bounds=([-2,-8,-1,-30],[2,8,1,30]))
        old=fit(np.ones(len(y)),np.zeros(4));f=flux(old.x)
        J=np.column_stack([(flux(old.x+np.eye(4)[j]*h)-flux(old.x-np.eye(4)[j]*h))/(2*h) for j,h in enumerate(steps)])
        jw=solve_triangular(L,J,lower=True)
        # At fixed transformed noise the model change is +K*f*dmag;
        # the compensating fitted change is -J^+*(+K*f*dmag).
        linear=np.linalg.lstsq(jw,solve_triangular(L,-K*f*correction,lower=True),rcond=None)[0]
        new=fit(scale,old.x+linear);shift=new.x-old.x
        rows.append(dict(CID=row.CID,zHEL=z,epochs=len(y),baseline_success=bool(old.success),alternative_success=bool(new.success),
            baseline_optimality=float(old.optimality),alternative_optimality=float(new.optimality),
            baseline_chi2=float(old.fun@old.fun),alternative_chi2=float(new.fun@new.fun),
            old_theta=old.x.tolist(),new_theta=new.x.tolist(),linear_theta=linear.tolist(),
            nonlinear_delta_standardized=float(D@shift),linear_delta_standardized=float(D@linear),
            difference_mag=float(D@(shift-linear))))
        print(row.CID,rows[-1]['nonlinear_delta_standardized'],rows[-1]['difference_mag'],flush=True)
    pd.DataFrame(rows).to_csv(out/'paired_refits.csv',index=False)
    report={'scope':'Hypothetical removal of fixed exact43 calibration-shaped pattern; native SALT refits with exact exported SNANA reference covariance scaled with data flux. No proof of true calibration bias, no retraining, alpha/beta fit or BBC/selection.',
        'delta_observed_magnitude_griz':dmag.tolist(),'all_success':all(r['baseline_success'] and r['alternative_success'] for r in rows),
        'maximum_absolute_nonlinear_minus_linear_mag':max(abs(r['difference_mag']) for r in rows),'records':rows}
    (out/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    inputs.extend(paths+[ROOT/'scripts/salt_dust_audit/flux_response.py',ROOT/'runs/salt_dust_audit/flux_response/selected_objects.csv'])
    outputs=list(out.glob('*'))
    (out/'manifest.json').write_text(json.dumps({'source_sha256':hashlib.sha256(source).hexdigest(),
        'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
        'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in outputs}},indent=2)+'\n')
    assert report['all_success']


if __name__=='__main__':main()
