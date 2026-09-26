"""Paired nonlinear foreground refits on fixed observed flux and covariance."""
from pathlib import Path
import sys
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import least_squares

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
from scripts.salt_dust_audit.flux_response import build_model, K


def main():
    out=ROOT/'runs/research_2026_09_26/foreground_response/nonlinear_check'
    out.mkdir(exist_ok=False)
    src=Path(__file__).read_bytes();(out/'executed_source.py').write_bytes(src)
    base=ROOT/'runs/research_2026_09_26/foreground_response'
    maps=pd.read_csv(base/'map_samples.csv',dtype={'CID':str}).sort_values('zHEL')
    ids=list(dict.fromkeys([maps.iloc[int(i)].CID for i in np.rint(np.linspace(0,len(maps)-1,5))]+[maps.loc[maps.delta_ebv_cleaning.abs().idxmax(),'CID']]))
    arrays=np.load(ROOT/'runs/salt_dust_audit/flux_response/matrices.npz')
    model,bands,paths,zp=build_model()
    offsets={str(r['Filter Name'])[-1]:float(r['Primary Mag']) for r in zp}
    D=np.array([1,.16087,-3.11780,0]);rows=[]
    for row in maps[maps.CID.isin(ids)].itertuples():
        prefix=row.CID+'__';t=arrays[prefix+'mjd'];band=arrays[prefix+'band']
        bs=np.array([bands[b] for b in band],dtype=object)
        scale=np.array([10**(-.4*(.27+offsets[b])) for b in band])
        L=np.linalg.cholesky(arrays[prefix+'measurement_plus_model_covariance'])
        observed=arrays[prefix+'flux_observed']
        def flux(theta,delta):
            model.set(z=row.zHEL,t0=row.PKMJD+theta[3],x0=row.x0*np.exp(-K*theta[0]),
                x1=row.x1+theta[1],c=row.c+theta[2],mwebv=row.HEAD_MWEBV+delta,
                mwrv=3.1,hostebv=0,hostrv=3.1)
            return scale*model.bandflux(bs,t,zp=27.5,zpsys='ab')
        assert np.allclose(flux(np.zeros(4),0),arrays[prefix+'flux_model'],rtol=1e-10,atol=1e-9)
        def fit(target,delta,start):
            fun=lambda theta:solve_triangular(L,flux(theta,delta)-target,lower=True)
            return least_squares(fun,start,xtol=1e-12,ftol=1e-12,gtol=1e-10,max_nfev=300,
                                 bounds=([-2,-8,-1,-30],[2,8,1,30]))
        linear=-arrays[prefix+'measurement_plus_model_response'][:,4]*row.delta_ebv_cleaning
        for kind,target in [('noiseless_reference',arrays[prefix+'flux_model']),('actual_observed',observed)]:
            old=fit(target,0,np.zeros(4));new=fit(target,row.delta_ebv_cleaning,old.x+linear)
            shift=new.x-old.x
            rows.append(dict(CID=row.CID,zHEL=row.zHEL,target=kind,delta_EBV=row.delta_ebv_cleaning,
                baseline_success=bool(old.success),alternative_success=bool(new.success),
                baseline_chi2=float(old.fun@old.fun),alternative_chi2=float(new.fun@new.fun),
                linear_delta_standardized=float(D@linear),nonlinear_delta_standardized=float(D@shift),
                difference_mag=float(D@(shift-linear)),old_theta=old.x.tolist(),new_theta=new.x.tolist(),
                baseline_optimality=float(old.optimality),alternative_optimality=float(new.optimality)))
        print(row.CID,flush=True)
    pd.DataFrame(rows).to_csv(out/'paired_refits.csv',index=False)
    summary={'selection':'Five fixed redshift ranks plus largest absolute map delta; deduplicated',
        'scope':'Same trained SALT3, fixed accepted epochs and reference native model covariance. No SNANA/BBC/selection closure.',
        'all_success':all(r['baseline_success'] and r['alternative_success'] for r in rows),
        'max_absolute_linear_error_by_target':{k:max(abs(r['difference_mag']) for r in rows if r['target']==k) for k in ['noiseless_reference','actual_observed']},
        'records':rows}
    (out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    inputs=[base/'map_samples.csv',ROOT/'runs/salt_dust_audit/flux_response/matrices.npz',*paths,
            ROOT/'scripts/salt_dust_audit/flux_response.py']
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps({'source_sha256':hashlib.sha256(src).hexdigest(),
        'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs}},indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items() if k!='records'},indent=2))
    assert summary['all_success']


if __name__=='__main__':main()
