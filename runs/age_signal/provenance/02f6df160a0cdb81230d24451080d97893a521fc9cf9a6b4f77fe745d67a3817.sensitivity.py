#!/usr/bin/env python3
"""Checks for input algebra, full covariance, redshift/cosmology, and age leverage."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyse import ROOT,OUT,DER,PP,AGES,CUTS,mu,wls,save_manifest

def run():
    inputs=[PP/'Pantheon+SH0ES_STAT+SYS.cov',PP/'Pantheon+SH0ES_STATONLY.cov']
    covs={}
    for p in inputs:
        v=np.loadtxt(p); n=int(v[0]); assert len(v)-1==n*n; covs[p.stem.split('_')[-1]]=v[1:].reshape(n,n)
    summary={}; rows=[]
    for pol in ['G11_first','R19_first']:
        d=pd.read_csv(DER/f'matched_{pol}.csv'); d=d[(d.zHD>.06)&(d.zHD<.42)].copy(); idx=d.pp_row.to_numpy(dtype=int)
        sub={k:v[np.ix_(idx,idx)] for k,v in covs.items()}
        # weighted removal/projection relation for a fitted reference cosmology
        x=d.age.to_numpy(); y=d.hr_corrected.to_numpy(); e=d.error_diag.to_numpy(); z=d.zHD.to_numpy(); w=1/e**2
        x0=x-np.average(x,weights=w); g=mu(z,d.zHEL,.301)-mu(z,d.zHEL,.299);g=g-np.average(g,weights=w)
        fraction=(np.sum(w*x0*g)**2)/(np.sum(w*x0*x0)*np.sum(w*g*g))
        bins=np.floor(z/.05).astype(int); within=np.empty_like(x)
        for k in np.unique(bins):
            m=bins==k;within[m]=x[m]-x[m].mean()
        raw=np.sum(within**2)/np.sum((x-x.mean())**2)
        corrected_var=np.var(x,ddof=1)-np.mean(d.age_err**2)
        summary[pol]={'age_std':float(x.std()),'mean_age_error':float(d.age_err.mean()),'observed_age_z_r':float(np.corrcoef(x,z)[0,1]),'bin_centering_attenuation_observed_age_OLS':float(raw),'cosmology_projection_loss_weighted_age_fraction':float(fraction),'moment_corrected_age_var':float(corrected_var),'latent_age_z_r2_moment_approx':float(np.cov(x,z)[0,1]**2/(corrected_var*np.var(z,ddof=1))),'covariances':{k:{'max_asymmetry':float(np.max(abs(v-v.T))),'min_eigenvalue':float(np.linalg.eigvalsh(v).min()),'diag_error_rms_delta':float(np.sqrt(np.mean((np.sqrt(np.diag(v))-e)**2))),'offdiag_max':float(np.max(abs(v-np.diag(np.diag(v)))))} for k,v in sub.items()}}
        for cut in [.2,.42]:
            m=z<cut; q=d[m]; ii=np.where(m)[0]
            for correction in ['hr_corrected','hr_no_bias']:
                for om in [.2,.3,.4]:
                    for zframe in ['zHD_zHEL','zHD_zHD']:
                        yy=q[correction].to_numpy()+q.mu_model.to_numpy()-mu(q.zHD,q.zHEL if zframe=='zHD_zHEL' else q.zHD,om)
                        for covname,c in [('DIAG',np.diag(q.error_diag**2))]+[(k,v[np.ix_(ii,ii)]) for k,v in sub.items()]:
                            b,C,chi=wls(q.age.to_numpy(),yy,q.error_diag.to_numpy(),cov=c)
                            rows.append({'policy':pol,'zcut':cut,'correction':correction,'Omega_m':om,'redshift_frame':zframe,'covariance':covname,'slope':b[1],'slope_se':np.sqrt(C[1,1]),'chi2':chi})
    pd.DataFrame(rows).to_csv(OUT/'sensitivity-wls.csv',index=False)
    (OUT/'sensitivity-audit.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
    # Verify effect direction and identify released standardization coefficients algebraically.
    q=pd.read_csv(DER/'matched_G11_first.csv'); X=np.c_[np.ones(len(q)),q.x1,q.c,(q.HOST_LOGMASS>10)*1]; target=q.m_b_corr-q.mB+q.biasCor_m_b
    beta=np.linalg.lstsq(X,target,rcond=None)[0]
    audit={'release_identity':'m_b_corr = mB + alpha*x1 - beta*c - biasCor_m_b + constant + mass_step*I(logM>10)','coefficient_labels':['constant','alpha','minus_beta','mass_step'],'coefficients':beta.tolist(),'max_abs_remainder_mag':float(np.max(np.abs(target-X@beta))),'note':'Empirical algebra check, not re-estimation from light curves; host dependence is inside simulation bias term.'}
    (OUT/'correction-algebra.json').write_text(json.dumps(audit,indent=2))
    save_manifest('A2-sensitivity',inputs+[DER/f'matched_{p}.csv' for p in ['G11_first','R19_first']]+[Path(__file__)],{'Omega_m':[.2,.3,.4],'age_errors':'ignored in WLS sensitivity, not final physical slope','cosmology_projection':'local single-parameter derivative at .3'},[OUT/'sensitivity-wls.csv',OUT/'sensitivity-audit.json',OUT/'correction-algebra.json'])
if __name__=='__main__':run()
