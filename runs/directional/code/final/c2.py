"""Faithful released-input C2 likelihood with exact analytic covariance gradients."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import sys
import time
import json
import hashlib
import pathlib
import numpy as np
import pandas as pd
from scipy.linalg import cholesky,cho_solve
from scipy.linalg.lapack import dpotri
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
from audit import SRC,OUT,ROOT,cosine,cubic_mu,write

class C2:
    def __init__(self,frame='zHEL'):
        p=pd.read_csv(SRC/'Analysis C2/Pantheon+SH0ES.dat',sep=r'\s+')
        ix=np.load(SRC/'index_sorted_lane.npy',allow_pickle=False)
        tab=p.iloc[ix].reset_index(drop=True)
        good=(tab.zHEL>.00937)&(tab.zHEL<.8)
        block=np.flatnonzero(good)
        ind=(3*block[:,None]+np.arange(3)).ravel()
        cov=np.load(OUT/'sources/cov_final.npy',allow_pickle=False)
        self.c=np.array(cov[np.ix_(ind,ind)],order='F');self.t=tab.loc[good]
        self.n=len(self.t);self.dim=3*self.n
        self.z=self.t[frame].to_numpy();self.zh=self.t.zHEL.to_numpy()
        self.cd=cosine(self.t.RA.to_numpy(),self.t.DEC.to_numpy())
        self.y=self.t[['mB','x1','c']].to_numpy().copy()
        self.y[:,0]+=self.t.biasCor_m_b.to_numpy()
        self.design=np.tile(np.eye(3),(self.n,1))
        self.ii=np.arange(self.n)*3
        self.ncall=0
    def calc(self,p,y,dip,gradient=True):
        q,j,a,b,lm,lx,lc=p[:7]
        vm,vx,vc=np.exp(2*np.array([lm,lx,lc]))
        xvec=np.array([-a,1,0]);cvec=np.array([b,0,1]);em=np.diag([1.,0,0])
        B=vm*em+vx*np.outer(xvec,xvec)+vc*np.outer(cvec,cvec)
        if dip:
            d,ls=p[7:9];S=np.exp(ls);ex=self.cd*np.exp(-self.z/S);q=q+d*ex
        mu,der=cubic_mu(self.z,self.zh,q,j)
        if mu is None:return 1e40,np.zeros(len(p))
        c=self.c.copy(order='F')
        for k in range(3):
            for l in range(3):c[self.ii+k,self.ii+l]+=B[k,l]
        L=cholesky(c,lower=True,overwrite_a=True,check_finite=False)
        res=y.copy();res[:,0]-=mu;res=res.ravel()
        sol=cho_solve((L,True),np.column_stack((res,self.design)),check_finite=False)
        mean=np.linalg.solve(self.design.T@sol[:,1:],self.design.T@sol[:,0])
        res-=self.design@mean
        wr=sol[:,0]-sol[:,1:]@mean
        chi=np.dot(res,wr);nll=chi+2*np.log(np.diag(L)).sum()+self.dim*np.log(2*np.pi)
        self.last=dict(q=float(p[0]),j=float(p[1]),alpha=float(a),beta=float(b),sigmaM=float(np.sqrt(vm)),sigmaX=float(np.sqrt(vx)),sigmaC=float(np.sqrt(vc)),M=float(mean[0]+a*mean[1]-b*mean[2]),X=float(mean[1]),C=float(mean[2]),qd=float(p[7]) if dip else 0.,S=float(np.exp(p[8])) if dip else None,nll2=float(nll),chi2=float(chi),parameters=p.tolist())
        if not gradient:return float(nll)
        inv,info=dpotri(L,lower=1,overwrite_c=1)
        assert info==0
        G=np.empty((3,3))
        r3=wr.reshape(-1,3)
        for k in range(3):
            for l in range(k+1):G[k,l]=G[l,k]=np.sum(inv[self.ii+k,self.ii+l])-np.dot(r3[:,k],r3[:,l])
        dBa=vx*np.array([[2*a,-1,0],[-1,0,0],[0,0,0]])
        dBb=vc*np.array([[2*b,0,1],[0,0,0],[1,0,0]])
        grad=[-2*np.dot(der[0],r3[:,0]),-2*np.dot(der[1],r3[:,0])]+[np.sum(G*x) for x in [dBa,dBb,2*vm*em,2*vx*np.outer(xvec,xvec),2*vc*np.outer(cvec,cvec)]]
        if dip:grad.extend([-2*np.dot(der[0]*ex,r3[:,0]),-2*np.dot(der[0]*d*ex*self.z/S,r3[:,0])])
        self.ncall+=1
        if self.ncall%10==0:print('iteration',self.ncall,'nll2',nll,'q',p[0],'qd',p[7] if dip else 0,flush=True)
        return float(nll),np.array(grad)
    def run(self,y,dip,start):
        bounds=[(-3,3),(-20,20),(.01,.5),(.1,6),(-8,np.log(.5)),(np.log(.2),np.log(2)),(np.log(.01),np.log(.2))]+([(-150,150),(np.log(self.z[0]),0)] if dip else [])
        t=time.time();r=minimize(self.calc,start,args=(y,dip),jac=True,method='L-BFGS-B',bounds=bounds,options={'ftol':1e-12,'gtol':1e-6,'maxiter':500,'maxls':30})
        self.calc(r.x,y,dip,False)
        return dict(**self.last,success=bool(r.success),message=str(r.message),gradient=r.jac.tolist(),seconds=time.time()-t,nfev=int(r.nfev))

if __name__=='__main__':
    f=C2();base=np.array([.01,-.65,.175,3.84,np.log(.121),np.log(.966),np.log(.058),-31.8,np.log(.0094)])
    age=pd.read_csv(SRC/'median_deltaage.csv',header=None).sort_values(0);sp=CubicSpline(age[0],age[1])
    mode=sys.argv[1] if len(sys.argv)>1 else 'run'
    if mode=='check':
        p=base.copy();p[-1]=np.log(.015);v,g=f.calc(p,f.y,True);fd=[]
        for i in range(len(p)):
            step=1e-5*max(abs(p[i]),1);dp=np.zeros(len(p));dp[i]=step
            fd.append((f.calc(p+dp,f.y,True,False)-f.calc(p-dp,f.y,True,False))/(2*step))
        write('c2-gradient-check.json',dict(analytic=g.tolist(),finite_difference=fd,max_scaled_error=float(np.max(abs(g-fd)/np.maximum(1,abs(g))))));print('gradient check',g,fd,flush=True)
    else:
        results=[]
        for corr in [False,True]:
            y=f.y.copy();y[:,0]-=.03*sp(f.zh) if corr else 0
            p=base.copy();p[0]=.35 if corr else .01
            d=f.run(y,True,p);print('DIPOLE',corr,json.dumps(d),flush=True)
            p2=np.array(d['parameters']);p2[0]+=.08;p2[1]+=.25;p2[7]*=.85;p2[8]=np.log(.014)
            d2=f.run(y,True,p2);print('DIPOLE_SECOND',corr,json.dumps(d2),flush=True)
            if d2['nll2']<d['nll2']:d,d2=d2,d
            iso=f.run(y,False,np.array(d['parameters'][:7]));print('ISOTROPIC',corr,json.dumps(iso),flush=True)
            results.append(dict(frame='zHEL',n=f.n,reverse_bias=True,age_correction=corr,dipole=d,second_start=d2,isotropic=iso,delta_nll2=iso['nll2']-d['nll2']))
            write('c2-fits.json',results)
        write('manifest-c2.json',dict(command='OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/c2.py',inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [pathlib.Path(__file__),ROOT/'scripts/directional/audit.py',OUT/'sources/cov_final.npy',SRC/'index_sorted_lane.npy',SRC/'Analysis C2/Pantheon+SH0ES.dat',SRC/'median_deltaage.csv']},configuration={'cut':'.00937 < zHEL < .8','frame':'zHEL','reverse_bias':True,'means':'profiled exactly','covariance_derivatives':'exact analytic','bounds':'in c2.py'}))
