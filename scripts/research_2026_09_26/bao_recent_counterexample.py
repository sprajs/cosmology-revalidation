"""Exploratory constructive test: past acceleration is not present acceleration.

Fits the same 13 compressed observations with continuous H and bounded piecewise
q. These fits construct compatible histories; they are not calibrated model
selection, posterior probabilities, or a survey-level likelihood reconstruction.
"""
from pathlib import Path
import json
import hashlib
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import minimize
from scipy.integrate import quad

from bao_shape import ROOT, OUT, INPUT

SOURCE=Path(__file__).read_bytes()


def geometry(z, logscale, q, edges):
    t=np.log1p(np.asarray(z));te=np.log1p(edges);dt=np.diff(te)
    scale=np.exp(logscale)
    gs=scale*np.exp(-np.r_[0,np.cumsum(q*dt)[:-1]])
    integ=np.empty(len(q));np.divide(-np.expm1(-q*dt),q,out=integ,where=np.abs(q)>1e-9)
    integ=np.where(np.abs(q)>1e-9,integ,dt)
    dstart=np.r_[0,np.cumsum(gs*integ)[:-1]]
    k=np.minimum(np.searchsorted(te[1:],t,side='right'),len(q)-1)
    width=t-te[k];qq=q[k]
    loc=np.empty_like(t);np.divide(-np.expm1(-qq*width),qq,out=loc,where=np.abs(qq)>1e-9)
    loc=np.where(np.abs(qq)>1e-9,loc,width)
    d=dstart[k]+gs[k]*loc;g=gs[k]*np.exp(-qq*width)
    return d,g/(1+np.asarray(z))


def main():
    out=OUT/'recent_counterexample'
    if out.exists():raise FileExistsError(out)
    out.mkdir();(out/'executed_source.py').write_bytes(SOURCE)
    inp=[INPUT/'desi_gaussian_bao_ALL_GCcomb_mean.txt',INPUT/'desi_gaussian_bao_ALL_GCcomb_cov.txt']
    df=pd.read_csv(inp[0],sep=r'\s+',comment='#',names=['z','value','kind'])
    cov=np.loadtxt(inp[1]);L=np.linalg.cholesky(cov);z=df.z.to_numpy();kind=df.kind.to_numpy();y=df.value.to_numpy()
    rng=np.random.default_rng(260929);allresults={};tables=[]
    for label,edges in [('main',np.array([0,.1,.3,.51,.706,.934,1.321,1.484,2.33])),
                        ('coarser',np.array([0,.3,.51,.934,1.484,2.33]))]:
        n=len(edges)-1
        def predict(par):
            dm,dh=geometry(z,par[0],par[1:],edges)
            return np.where(kind=='DM_over_rs',dm,np.where(kind=='DH_over_rs',dh,(z*dm**2*dh)**(1/3)))
        def objective(par):
            return float(np.sum(solve_triangular(L,y-predict(par),lower=True)**2))
        for branch in ['free','recent_nonnegative','recent_q_plus_half','always_nonnegative']:
            bounds=[(np.log(10),np.log(100))]+[(-3,3)]*n
            for i in range(n):
                if branch=='always_nonnegative' or (branch=='recent_nonnegative' and edges[i+1]<=.3):bounds[i+1]=(0,3)
                if branch=='recent_q_plus_half' and edges[i+1]<=.3:bounds[i+1]=(.5,.5)
            fits=[]
            for k in range(10):
                x=np.r_[np.log(30),rng.uniform(-.7,.7,n)]
                x=np.array([np.clip(v,*b) for v,b in zip(x,bounds)])
                f=minimize(objective,x,bounds=bounds,method='L-BFGS-B',
                           options={'maxiter':4000,'ftol':1e-12,'gtol':1e-7})
                fits.append({'chi2':float(f.fun),'x':f.x.tolist(),'success':bool(f.success)})
            best=min(fits,key=lambda a:a['chi2'])
            # Verify exact integral with independent quadrature of E^-1.
            par=np.array(best['x']);scale=np.exp(par[0]);qq=par[1:];te=np.log1p(edges)
            def inverse_e(zz):
                tz=np.log1p(zz);v=0.
                for q,lo,hi in zip(qq,te[:-1],te[1:]):v+=q*max(0,min(tz,hi)-lo)
                return np.exp(-v)/(1+zz)
            dm,dh=geometry(z,par[0],qq,edges)
            direct=np.array([scale*quad(inverse_e,0,zz,points=[e for e in edges[1:-1] if e<zz],epsabs=1e-11,epsrel=1e-11)[0] for zz in z])
            err=float(np.max(np.abs(dm-direct)));assert err<1e-8
            key=label+'__'+branch
            pred=predict(par)
            for i,row in df.iterrows():tables.append({'fit':key,'z':row.z,'kind':row.kind,'observed':row.value,'predicted':pred[i],'marginal_sigma':np.sqrt(cov[i,i])})
            allresults[key]={'best':best,'all_starts':fits,'edges':edges.tolist(),
                'c_over_H0rd':float(scale),'quadrature_max_error':err,
                'parameters_on_bounds':[i for i,(v,b) in enumerate(zip(par,bounds)) if min(abs(v-b[0]),abs(v-b[1]))<1e-4],
                'q_is_binwise_not_point_measurement':True}
            print(key,best['chi2'],best['x'],flush=True)
    (out/'results.json').write_text(json.dumps({'status':'exploratory constructive histories; no calibrated significance or posterior',
        'assumptions':'flat constant-ruler effective-redshift Gaussian BAO, free scale, finite q bins',
        'seed':260929,'fits':allresults},indent=2)+'\n')
    pd.DataFrame(tables).to_csv(out/'predictions.csv',index=False)
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps({'utc':datetime.now(timezone.utc).isoformat(),
        'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inp},
        'source_sha256':hashlib.sha256(SOURCE).hexdigest(),
        'outputs_sha256':{p.name:sha(p) for p in [out/'results.json',out/'predictions.csv']}},indent=2)+'\n')


if __name__=='__main__':main()
