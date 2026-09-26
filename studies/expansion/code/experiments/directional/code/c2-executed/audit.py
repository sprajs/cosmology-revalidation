"""Independent C1 cosmography reproduction and preregistered synthetic checks.

Run with OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/audit.py.
No author Python code is executed and no pickle objects are loaded.
"""
import json
import hashlib
import pathlib
import sys
import datetime
import numpy as np
import pandas as pd
from scipy import linalg, optimize, integrate
from scipy.interpolate import CubicSpline
from astropy.cosmology import FlatLambdaCDM

ROOT = pathlib.Path(__file__).resolve().parents[2]
SRC = ROOT / 'sources/repos/Shin107__Anisotropy-in-Pantheon-Plus'
OUT = ROOT / 'runs/directional'
CONST = 25 + 5*np.log10(299792.458/70)

def cosine(ra,dec,ra0=168,dec0=-7):
    r,d,r0,d0=np.deg2rad([ra,dec,np.full_like(ra,ra0),np.full_like(dec,dec0)])
    return np.sin(d)*np.sin(d0)+np.cos(d)*np.cos(d0)*np.cos(r-r0)

def load_c1(zmin=.00937,zmax=.8,frame='zHEL'):
    z=pd.read_csv(SRC/'Analysis C1/Z_mbcorr.csv')
    p=pd.read_csv(SRC/'Analysis C1/Pantheon+SH0ES.dat',sep=r'\s+')
    for col in z: assert np.array_equal(z[col],p[col]),col
    c=np.load(SRC/'Analysis C1/statsys_mbcorr.npy',allow_pickle=False)
    # Match author's lower=True Cholesky convention, despite tiny rounding asymmetry.
    c=np.tril(c)+np.tril(c,-1).T
    z['zLG']=(1+z.zHEL)*np.sqrt((1+299*cosine(z.RA.values,z.DEC.values,333.53277784,49.33111122)/299792.458)/(1-299*cosine(z.RA.values,z.DEC.values,333.53277784,49.33111122)/299792.458))-1
    z=z[(z.zHEL>zmin)&(z.zHEL<zmax)].sort_values('zHEL')
    ix=z.index.to_numpy(); c=c[np.ix_(ix,ix)]
    zc=z[frame].to_numpy(); zh=z.zHEL.to_numpy()
    cd=cosine(z.RA.to_numpy(),z.DEC.to_numpy(),162.95389715,-25.96734154) if frame=='zLG' else cosine(z.RA.to_numpy(),z.DEC.to_numpy())
    return z,zc,zh,cd,c

def cubic_mu(z,zh,q,j):
    b=1+.5*(1-q)*z-(1-q-3*q*q+j)*z*z/6
    if np.any(b<=0): return None,None
    mu=CONST+5*np.log10(z*b*(1+zh)/(1+z))
    dq=5/np.log(10)*(-z/2+(1+6*q)*z*z/6)/b
    dj=5/np.log(10)*(-z*z/6)/b
    return mu,(dq,dj)

class Fit:
    def __init__(self,z,zh,cd,c):
        self.z,self.zh,self.cd=z,zh,cd
        self.e,self.u=linalg.eigh(c,check_finite=False)
        assert self.e.min()>0
        self.one=self.u.T@np.ones(len(z))
    def objective(self,p,y,dip=False,scatter=False,detail=False):
        q,j=p[:2]; n=2
        if dip:
            d,logS=p[2:4]; S=np.exp(logS); ex=self.cd*np.exp(-self.z/S); q=q+d*ex;n=4
        mu,der=cubic_mu(self.z,self.zh,q,j)
        if mu is None: return (1e40,np.zeros(len(p)))
        v=p[n] if scatter else 0
        den=self.e+v
        r=self.u.T@(y-mu)
        m=np.sum(self.one*r/den)/np.sum(self.one**2/den)
        r-=m*self.one
        chi=np.sum(r*r/den)
        val=chi+np.log(den).sum()+len(y)*np.log(2*np.pi)
        gmu=[der[0],der[1]]
        if dip:gmu.extend([der[0]*ex,der[0]*d*ex*self.z/S])
        wr=self.u@(r/den)
        grad=[-2*np.dot(x,wr) for x in gmu]
        if scatter:grad.append(np.sum(1/den-r*r/den**2))
        if detail:return dict(q=float(p[0]),j=float(p[1]),M=float(m),s=float(np.sqrt(v)),qd=float(p[2]) if dip else 0,S=float(np.exp(p[3])) if dip else None,nll2=float(val),chi2=float(chi))
        return val,np.array(grad)
    def fit(self,y,dip=False,scatter=False,start=None,multistart=True):
        p=[-.4,.5]+([-6,np.log(.025)] if dip else [])+([.001] if scatter else [])
        if start is not None:p=list(start)
        bounds=[(-3,3),(-20,20)]+([(-150,150),(np.log(self.z[0]),np.log(1.))] if dip else [])+([(0,.25)] if scatter else [])
        starts=[p]
        if dip and multistart:
            for d,s in [(-30,.01),(20,.012),(-6,.07),(0,.04)]:
                a=p.copy();a[2]=d;a[3]=np.log(max(s,self.z[0]));starts.append(a)
        fits=[optimize.minimize(self.objective,a,args=(y,dip,scatter),jac=True,bounds=bounds,method='L-BFGS-B',options={'ftol':1e-13,'gtol':1e-7,'maxiter':1500,'maxls':40}) for a in starts]
        best=min(fits,key=lambda r:r.fun)
        out=self.objective(best.x,y,dip,scatter,True)
        out.update(success=bool(best.success),message=str(best.message),parameters=best.x.tolist(),start_objectives=[float(r.fun) for r in fits],gradient_max=float(abs(best.jac).max()))
        return out
    def sigmaq(self,q,j):
        _,der=cubic_mu(self.z,self.zh,q,j)
        a=np.array([np.ones(len(self.z)),*der]).T
        rot=self.u.T@a
        return float(np.sqrt(np.linalg.inv(rot.T@(rot/self.e[:,None]))[1,1]))

def exact(z,kind):
    if kind=='lcdm':
        v=np.array([integrate.quad(lambda t:1/np.sqrt(.3*(1+t)**3+.7),0,float(x),epsabs=1e-12,epsrel=1e-12)[0] for x in z])
    elif kind=='eds':v=2*(1-1/np.sqrt(1+z))
    elif kind=='coasting':v=np.log1p(z)
    return (1+z)*v

def write(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')

def data_fits():
    results=[]
    age=pd.read_csv(SRC/'median_deltaage.csv',header=None).sort_values(0)
    spline=CubicSpline(age[0],age[1])
    for frame in ['zHEL','zCMB','zHD','zLG']:
        tab,z,zh,cd,c=load_c1(frame=frame);f=Fit(z,zh,cd,c)
        for rev in [True,False]:
            for correction in [False,True]:
                y=tab.m_b_corr.to_numpy()+(tab.biasCor_m_b.to_numpy() if rev else 0)-(.03*spline(zh) if correction else 0)
                iso=f.fit(y,scatter=True)
                di=f.fit(y,dip=True,scatter=True)
                rec=dict(frame=frame,reverse_bias=rev,age_correction=correction,n=len(z),isotropic=iso,dipole=di,delta_nll2=iso['nll2']-di['nll2'])
                results.append(rec);print(json.dumps(rec),flush=True);write('c1-fits.json',results)

def synthetic():
    models={'lcdm':(-.55,1),'eds':(.5,1),'coasting':(0,0)}
    grid=np.array([.01,.1,.2,.4,.5,.8])
    errors=[]
    for kind,(q,j) in models.items():
        mu,_=cubic_mu(grid,grid,q,j);dm=mu-CONST-5*np.log10(exact(grid,kind))
        errors.append(dict(model=kind,z=grid.tolist(),mu_error=dm.tolist(),relative_distance_error=(10**(dm/5)-1).tolist()))
    ast=FlatLambdaCDM(H0=70,Om0=.3,Tcmb0=0)
    assert np.max(np.abs(ast.luminosity_distance(grid).value/(299792.458/70)/exact(grid,'lcdm')-1))<1e-12
    write('taylor-errors.json',errors)
    rng=np.random.default_rng(260609650);results=[];mockrows=[]
    for zmax in [.1,.2,.4,.8]:
        tab,z,zh,cd,c=load_c1(zmin=.01,zmax=zmax,frame='zHD');f=Fit(z,zh,cd,c)
        # Same noise draws across models for paired comparison.
        noise=(f.u*np.sqrt(f.e))@rng.normal(size=(len(z),100))
        for kind,(q,j) in models.items():
            y=CONST+5*np.log10(exact(z,kind)*(1+zh)/(1+z))-19.3
            fit=f.fit(y,start=[q,j]);sq=f.sigmaq(fit['q'],fit['j'])
            qmock=[]
            for k in range(100):
                a=f.fit(y+noise[:,k],start=[fit['q'],fit['j']],multistart=False);qmock.append(a['q']);mockrows.append(dict(model=kind,zmax=zmax,realization=k,q=a['q'],j=a['j'],success=a['success']))
            rec=dict(model=kind,zmax=zmax,n=len(z),truth_q=q,truth_j=j,noiseless=fit,forecast_sigma_q=sq,bias_q=fit['q']-q,bias_in_sigma=(fit['q']-q)/sq,noise_mean_q=float(np.mean(qmock)),noise_sd_q=float(np.std(qmock,ddof=1)),noise_q_quantiles=np.quantile(qmock,[.025,.5,.975]).tolist())
            if zmax==.8 and kind=='lcdm':
                rec['dipole_exact_injection']=f.fit(y,dip=True,start=[q,j,0,np.log(.025)])
                ty=cubic_mu(z,zh,q,j)[0]-19.3
                rec['dipole_cubic_control']=f.fit(ty,dip=True,start=[q,j,0,np.log(.025)])
                rec['cubic_control']=f.fit(ty,start=[q,j])
            results.append(rec);print(json.dumps(rec),flush=True);write('synthetic-summary.json',results)
            pd.DataFrame(mockrows).to_csv(OUT/'synthetic-realizations.csv',index=False)

if __name__=='__main__':
    OUT.mkdir(exist_ok=True,parents=True)
    mode=sys.argv[1] if len(sys.argv)>1 else 'all'
    files=[SRC/'Analysis C1/Z_mbcorr.csv',SRC/'Analysis C1/Pantheon+SH0ES.dat',SRC/'Analysis C1/statsys_mbcorr.npy',SRC/'median_deltaage.csv',pathlib.Path(__file__)]
    write('manifest-'+mode+'.json',dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),command='OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/audit.py '+mode,seed=260609650,inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},versions={'numpy':np.__version__,'pandas':pd.__version__},plan='docs/experiments/directional-plan.md'))
    if mode in ['all','data']:data_fits()
    if mode in ['all','synthetic']:synthetic()
