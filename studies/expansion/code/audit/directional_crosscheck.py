"""Independent direct-matrix and exact-distance checks of the directional branch."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import numpy as np
import pandas as pd
import json, hashlib, datetime
from scipy.linalg import cholesky,solve_triangular
from scipy.optimize import minimize
from independent import ode_distance

ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'runs/audit'
BASE=ROOT/'sources/repos/Shin107__Anisotropy-in-Pantheon-Plus/Analysis C1'
data=pd.read_csv(BASE/'Z_mbcorr.csv')
raw=np.load(BASE/'statsys_mbcorr.npy',allow_pickle=False)
# Reproduce source's lower-triangle covariance convention, explicitly.
cov=np.tril(raw)+np.tril(raw,-1).T
old=json.loads((ROOT/'runs/directional/c1-fits.json').read_text())
rec=next(r for r in old if r['frame']=='zHEL' and r['reverse_bias'] and not r['age_correction'])
d=data[(data.zHEL>.00937)&(data.zHEL<.8)].sort_values('zHEL');ids=d.index.to_numpy();z=d.zHEL.to_numpy()
c=cov[np.ix_(ids,ids)];y=(d.m_b_corr+d.biasCor_m_b).to_numpy()
def cosmographic_mu(z,zh,q,j):
    # Cubic FLRW luminosity-distance polynomial, independent scalar implementation.
    dl=z+(1-q)/2*z*z-(1-q-3*q*q+j)/6*z**3
    return 25+5*np.log10(299792.458/70*dl*(1+zh)/(1+z))
results={}
for label in ['isotropic','dipole']:
    fit=rec[label];q=fit['q']
    if label=='dipole':
        ra,dec=np.deg2rad(d.RA),np.deg2rad(d.DEC);r0=np.deg2rad(168.);d0=np.deg2rad(-7.)
        cosine=np.sin(dec)*np.sin(d0)+np.cos(dec)*np.cos(d0)*np.cos(ra-r0)
        q=q+fit['qd']*cosine*np.exp(-z/fit['S'])
    residual=y-cosmographic_mu(z,z,q,fit['j'])
    L=cholesky(c+fit['s']**2*np.eye(len(c)),lower=True)
    w=solve_triangular(L,residual,lower=True);u=solve_triangular(L,np.ones(len(c)),lower=True)
    M=w@u/(u@u);w-=M*u;chi=w@w;logdet=2*np.log(np.diag(L)).sum();nll=chi+logdet+len(c)*np.log(2*np.pi)
    results[label]={'rows':len(c),'chi2':float(chi),'M':float(M),'minus2logL':float(nll),
           'minus2logL_difference_from_eigen_solver':float(nll-fit['nll2'])}
# Noiseless exact-LCDM apparent magnitudes at real Pantheon sky/redshift covariance.
d=data[(data.zHEL>.01)&(data.zHEL<.8)].sort_values('zHEL');ids=d.index.to_numpy();z=d.zHD.to_numpy();zh=d.zHEL.to_numpy()
L=cholesky(cov[np.ix_(ids,ids)],lower=True);u=solve_triangular(L,np.ones(len(d)),lower=True)
chi,_=ode_distance([.3,-1,0],z);y=25+5*np.log10(299792.458/70*(1+zh)*chi)-19.3
def objective(t):
    r=solve_triangular(L,y-cosmographic_mu(z,zh,*t),lower=True,check_finite=False)
    r-=u*(u@r)/(u@u)
    return r@r
fit=minimize(objective,[-.55,1.],method='Nelder-Mead',options={'xatol':1e-10,'fatol':1e-12,'maxiter':2000})
syn=json.loads((ROOT/'runs/directional/synthetic-summary.json').read_text())
ref=next(r for r in syn if r['model']=='lcdm' and r['zmax']==.8)
results['noiseless_exact_lcdm_fit']={'rows':len(d),'q_truth':-.55,'q_fit':float(fit.x[0]),'j_fit':float(fit.x[1]),
       'q_bias':float(fit.x[0]+.55),'q_difference_from_directional_fit':float(fit.x[0]-ref['noiseless']['q']),
       'chi2':float(fit.fun),'optimizer_success':bool(fit.success)}
grid=np.array([.1,.4,.8]);chi,_=ode_distance([.3,-1,0],grid)
err=cosmographic_mu(grid,grid,-.55,1)-(25+5*np.log10(299792.458/70*(1+grid)*chi))
results['lcdm_taylor_error']={'z':grid.tolist(),'mu_error_mag':err.tolist()}
path=OUT/'directional-crosscheck.json';path.write_text(json.dumps(results,indent=2)+'\n')
files=[BASE/'Z_mbcorr.csv',BASE/'statsys_mbcorr.npy',ROOT/'runs/directional/c1-fits.json',ROOT/'runs/directional/synthetic-summary.json',Path(__file__),ROOT/'scripts/audit/independent.py']
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest={'time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Direct Cholesky cross-check of directional eigenvalue likelihood and cubic bias',
          'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in files},'outputs_sha256':{str(path.relative_to(ROOT)):sha(path)},'random_seed':None}
(OUT/'directional-manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(results,indent=2))
