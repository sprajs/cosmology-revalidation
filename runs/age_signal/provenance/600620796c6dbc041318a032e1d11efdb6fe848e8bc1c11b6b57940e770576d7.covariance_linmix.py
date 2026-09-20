#!/usr/bin/env python3
from pathlib import Path
import sys,json
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyse import ROOT,OUT,DER,PP,save_manifest,rhat_four
VENDOR=ROOT/'sources/external/age_signal/linmix-933dbb1359dcb5404cf881bfdd5cf433b0195152'
sys.path.insert(0,str(VENDOR));import linmix

def main():
    path=PP/'Pantheon+SH0ES_STAT+SYS.cov';v=np.loadtxt(path);n=int(v[0]);C=v[1:].reshape(n,n); results=[]
    for i,pol in enumerate(['G11_first','R19_first']):
        d=pd.read_csv(DER/f'matched_{pol}.csv');d=d[(d.zHD>.06)&(d.zHD<.42)]
        for j,y in enumerate(['hr_corrected','hr_no_bias']):
            for k,cut in enumerate([.2,.42]):
                q=d[d.zHD<cut]; seed=20261020+i*4+j*2+k;name=f'{pol}_{y}_covdiag_{cut}'
                target=OUT/'covariance-linmix'/f'{name}.json';target.parent.mkdir(exist_ok=True)
                if target.exists(): results.append(json.loads(target.read_text()));continue
                lm=linmix.LinMix(q.age.to_numpy(),q[y].to_numpy(),xsig=q.age_err.to_numpy(),ysig=np.sqrt(np.diag(C)[q.pp_row.to_numpy(dtype=int)]),nchains=4,parallelize=True,seed=seed,K=3)
                print('START',name,flush=True);lm.run_mcmc(miniter=20000,maxiter=100000,silent=False);b=lm.chain['beta']
                r={'name':name,'n':len(q),'seed':seed,'mean':float(b.mean()),'sd':float(b.std()),'q025_975':np.quantile(b,[.025,.975]).tolist(),'p_negative':float(np.mean(b<0)),'scatter_mean':float(np.mean(np.sqrt(lm.chain['sigsqr']))),'retained_draws':len(b),'split_rhat':{key:rhat_four(lm.chain[key]) for key in ['alpha','beta','sigsqr']}}
                target.write_text(json.dumps(r,indent=2));np.savez_compressed(target.with_suffix('.npz'),chain=lm.chain);results.append(r);print(json.dumps(r),flush=True)
    summary=OUT/'covariance-linmix-summary.json';summary.write_text(json.dumps(results,indent=2))
    save_manifest('A2-actual-covdiag-linmix',[path,Path(__file__),VENDOR/'linmix/linmix.py']+[DER/f'matched_{p}.csv' for p in ['G11_first','R19_first']],{'nchains':4,'miniter':20000,'maxiter':100000,'seed_base':20261020,'approximation':'actual covariance diagonal used, offdiagonal neglected by LINMIX; full covariance GLS separate'},[summary]+list((OUT/'covariance-linmix').glob('*.json')))
if __name__=='__main__':main()
