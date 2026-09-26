#!/usr/bin/env python3
"""Bounded age slope sensitivity to independently reconstructed low-z revision."""
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyse import ROOT,OUT,DER,PP,save_manifest,rhat_four,wls
VENDOR=ROOT/'sources/external/age_signal/linmix-933dbb1359dcb5404cf881bfdd5cf433b0195152';sys.path.insert(0,str(VENDOR));import linmix

def main():
    original=PP/'Pantheon+SH0ES.dat';revision=ROOT/'data/derived/standardization/Pantheon_W26_lowz_mass_revision.dat'
    old=pd.read_csv(original,sep=r'\s+',dtype={'CID':str});new=pd.read_csv(revision,sep=r'\s+',dtype={'CID':str})
    assert np.array_equal(old.CID,new.CID);assert np.array_equal(old.IDSURVEY,new.IDSURVEY)
    delta=(new.m_b_corr-old.m_b_corr).to_numpy();p=PP/'Pantheon+SH0ES_STAT+SYS.cov';v=np.loadtxt(p);n=int(v[0]);C=v[1:].reshape(n,n);summ=[]
    for k,pol in enumerate(['G11_first','R19_first']):
        d=pd.read_csv(DER/f'matched_{pol}.csv');d=d[(d.zHD>.06)&(d.zHD<.42)];idx=d.pp_row.to_numpy(dtype=int);sC=C[np.ix_(idx,idx)];sC=(sC+sC.T)/2
        y=d.hr_corrected.to_numpy()+delta[idx];age=d.age.to_numpy();e=np.sqrt(np.diag(sC));name=f'{pol}_lowz_revision_covdiag';seed=20260970+k
        lm=linmix.LinMix(age,y,xsig=d.age_err.to_numpy(),ysig=e,nchains=4,parallelize=True,seed=seed,K=3);lm.run_mcmc(miniter=20000,maxiter=100000,silent=False);b=lm.chain['beta']
        target=OUT/'mass-revision';target.mkdir(exist_ok=True);np.savez_compressed(target/f'{name}.npz',chain=lm.chain)
        coef,cov,chi=wls(age,y,e,cov=sC)
        r={'name':name,'n':len(d),'changed_rows':int(np.sum(abs(delta[idx])>1e-8)),'mean_delta_all':float(delta[idx].mean()),'mean':float(b.mean()),'sd':float(b.std()),'q025_975':np.quantile(b,[.025,.975]).tolist(),'p_negative':float(np.mean(b<0)),'seed':seed,'retained_draws':len(b),'split_rhat':{key:rhat_four(lm.chain[key]) for key in ['alpha','beta','sigsqr']},'fullcov_gls_slope':float(coef[1]),'fullcov_gls_se':float(np.sqrt(cov[1,1])),'limitation':'Original covariance retained; LINMIX uses diagonal, GLS full covariance treats measured age fixed'}
        summ.append(r);print(json.dumps(r),flush=True)
        d[['CID','age_source','pp_row','age']].assign(delta_hr=delta[idx],hr_revised=y).to_csv(DER/f'mass_revision_{pol}.csv',index=False)
    result=OUT/'mass-revision-summary.json';result.write_text(json.dumps(summ,indent=2));save_manifest('A2-mass-revision',[original,revision,p,Path(__file__)]+[DER/f'matched_{p}.csv' for p in ['G11_first','R19_first']],{'seed_base':20260970,'nchains':4,'miniter':20000,'maxiter':100000,'delta':'new.m_b_corr-old.m_b_corr'},[result]+list(DER.glob('mass_revision_*.csv')))
if __name__=='__main__':main()
