#!/usr/bin/env python3
"""Gaussian latent-age population sensitivity, full distance covariance.
This distributional approximation is distinct from original LINMIX.
"""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
from scipy.linalg import cholesky,solve_triangular
from scipy.optimize import minimize
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyse import ROOT,OUT,DER,PP,save_manifest

def run():
    p=PP/'Pantheon+SH0ES_STAT+SYS.cov'; flat=np.loadtxt(p);n=int(flat[0]); C=flat[1:].reshape(n,n)
    d=pd.read_csv(DER/'matched_G11_first.csv');d=d[(d.zHD>.06)&(d.zHD<.42)];ids=d.pp_row.to_numpy(dtype=int);C=C[np.ix_(ids,ids)];C=(C+C.T)/2
    age=d.age.to_numpy();sa=d.age_err.to_numpy();z=d.zHD.to_numpy()-d.zHD.mean();mass=d.HOST_LOGMASS.to_numpy()-d.HOST_LOGMASS.mean();colour=d.c.to_numpy()-d.c.mean()
    summaries=[]; curves=[]
    for ycol in ['hr_corrected','hr_no_bias']:
      y=d[ycol].to_numpy()
      for covname,joint in [('full',False),('diagonal',False),('full',True)]:
        name=f'{ycol}_{covname}_{"joint" if joint else "base"}';co=C if covname=='full' else np.diag(np.diag(C))
        # theta=m,logtau,alpha,b,logscatter[,age_z,hr_z,hr_mass,hr_colour]
        bounds=[(0,14),(-4,2.6),(-3,3),(-.2,.15),(-9,0)] + ([(-30,30),(-5,5),(-.5,.5),(-3,3)] if joint else [])
        def nll(theta):
            m,lt,alpha,b,ls=theta[:5]; tau2=np.exp(2*lt);s2=np.exp(2*ls);mi=np.full(len(age),m);extra=np.zeros(len(age))
            if joint: mi=mi+theta[5]*z;extra=theta[6]*z+theta[7]*mass+theta[8]*colour
            av=tau2+sa**2;am=mi+tau2/av*(age-mi);apv=tau2*sa**2/av
            cov=co+np.diag(s2+b*b*apv)
            try:L=cholesky(cov,lower=True,check_finite=False)
            except np.linalg.LinAlgError:return 1e20
            residual=y-alpha-b*am-extra;rv=solve_triangular(L,residual,lower=True,check_finite=False)
            return float(.5*np.sum(np.log(av)+(age-mi)**2/av)+np.sum(np.log(np.diag(L)))+.5*(rv@rv))
        starts=[np.array([age.mean(),np.log(t),y.mean()-b*age.mean(),b,np.log(.07)]+([0,0,0,0] if joint else [])) for t,b in [(1,-.03),(2,-.01),(.3,-.06)]]
        fits=[minimize(nll,s,method='L-BFGS-B',bounds=bounds,options={'maxiter':2000,'ftol':1e-12,'gtol':1e-6}) for s in starts]; best=min(fits,key=lambda f:f.fun)
        grid=np.arange(-.10,.040001,.002);profile=[];theta=best.x.copy()
        keep=np.array([i for i in range(len(theta)) if i!=3]); nb=[bounds[i] for i in keep]
        def optimise_at(b,t):
            def f(par):q=t.copy();q[keep]=par;q[3]=b;return nll(q)
            return minimize(f,t[keep],method='L-BFGS-B',bounds=nb,options={'maxiter':1000,'ftol':1e-11,'gtol':1e-6})
        for b in grid:
            fits=[optimise_at(b,theta),optimise_at(b,best.x)];rr=min(fits,key=lambda f:f.fun);theta[keep]=rr.x;theta[3]=b
            profile.append(rr.fun);curves.append({'name':name,'slope':b,'nll':rr.fun,'optim_success':bool(rr.success),'nuisance_boundary':any(abs(theta[i]-bounds[i][j])<1e-4 for i in keep for j in [0,1])})
        base=min(best.fun,min(profile));delta=2*(np.asarray(profile)-base);interval={}
        for thresh,label in [(1.,'delta_chi2_1'),(3.84,'delta_chi2_3p84')]:
            allowed=grid[delta<=thresh];interval[label]=[float(allowed.min()),float(allowed.max())] if len(allowed) else None
        null=min(profile,key=lambda v:abs(v)) if False else profile[int(np.argmin(abs(grid)))]
        r={'name':name,'n':len(d),'optimizer_success':bool(best.success),'best_theta':best.x.tolist(),'slope_mle':float(best.x[3]),'sigma_intrinsic_mle':float(np.exp(best.x[4])),'tau_age_mle':float(np.exp(best.x[1])),'minus2loglike_null_delta':float(2*(null-base)),'profile_intervals':interval,'nll':float(base),'profile_grid_bounds':[-.1,.04],'note':'likelihood profile, no posterior probability; Gaussian age-summary assumption'}
        summaries.append(r);print(json.dumps(r),flush=True)
        # incremental durable output so interrupted fits are not silently lost
        (OUT/'fullcov-eiv-summary.json').write_text(json.dumps(summaries,indent=2));pd.DataFrame(curves).to_csv(OUT/'fullcov-eiv-profiles.csv',index=False)
    fig,ax=plt.subplots(figsize=(7,4),layout='constrained');df=pd.DataFrame(curves)
    for name,q in df.groupby('name'):ax.plot(q.slope,2*(q.nll-q.nll.min()),label=name)
    ax.axhline(3.84,color='black',ls=':');ax.axvline(0,color='black',lw=.6);ax.set(ylim=(0,12),xlabel='Host-age slope (mag/Gyr)',ylabel='Profile −2 Δ log likelihood',title='Gaussian latent-age distribution sensitivity');ax.legend(fontsize=7);fig.savefig(OUT/'fullcov-eiv-profiles.png',dpi=170);plt.close(fig)
    save_manifest('A2-fullcov-eiv',[p,DER/'matched_G11_first.csv',Path(__file__)],{'covariance':'STAT+SYS original row subset, symmetrized','age_population':'normal, mean constant or linear z','intrinsic_scatter':'free','grid':[-.1,.04,.002],'optimiser':'L-BFGS-B multi-start','mass_colour_errors':'not modelled','profile':'not Bayesian posterior'},[OUT/'fullcov-eiv-summary.json',OUT/'fullcov-eiv-profiles.csv',OUT/'fullcov-eiv-profiles.png'])
if __name__=='__main__':run()
