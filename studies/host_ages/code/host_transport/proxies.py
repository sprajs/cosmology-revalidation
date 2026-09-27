#!/usr/bin/env python3
"""Execute cross-survey selected-host support with explicitly model-derived proxies."""
from pathlib import Path
import hashlib,json,datetime
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.optimize import minimize, linprog
ROOT=Path(__file__).resolve().parents[4];W=ROOT/'.work/host-transport';O=ROOT/'studies/host_ages/results/host_transport';HERE=Path(__file__).parent

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def balance(x,t):
    mean=x.mean(0);sd=x.std(0);a=(x-mean)/sd;b=(t-mean)/sd
    feasible=linprog(np.zeros(len(a)),A_eq=np.vstack([np.ones(len(a)),a.T]),b_eq=np.r_[1,b],bounds=(0,None),method='highs')
    if not feasible.success:raise ValueError('Target moments outside sampled convex hull')
    def fun(k):
        q=a@k;l=logsumexp(q);w=np.exp(q-l)
        return l-b@k,(w@a-b)
    opt=minimize(fun,np.zeros(x.shape[1]),jac=True,method='BFGS',options={'gtol':1e-9,'maxiter':2000})
    w=np.exp(a@opt.x-logsumexp(a@opt.x));err=float(np.max(abs(w@a-b)))
    if not np.isfinite(err) or not np.isfinite(w).all() or err>1e-6:raise ValueError(f'Unbalanced or unsupported target: {err}')
    return w,{'standardized_moment_residual':err,'ess':float(1/(w@w)),'largest_weight':float(w.max()),'parameters':opt.x.tolist()}
def main():
    plan=json.loads((HERE/'proxy-design.json').read_text());rng=np.random.default_rng(plan['seed']);p=W/'roman-properties.csv';d=pd.read_csv(p)
    valid=np.isfinite(d[['uv_local','uv_global','euv_local','euv_global','stellar_mass','redshift']]).all(axis=1)&(d.euv_local>0)&(d.euv_global>0)&(d.uv_local<90)&(d.uv_global<90)&(d.stellar_mass>0)
    d=d[valid].copy();d['local_minus_global']=d.uv_local-d.uv_global;out=[]
    for survey in ['all','SDSS','SNLS5','CfAIII','CfAIV','CSP']:
        for lo,hi in plan['redshift_bins']:
            t=d[(d.redshift>=lo)&(d.redshift<hi)&((d.survey==survey) if survey!='all' else True)]
            if len(t)==0:continue
            y=t.local_minus_global.to_numpy();n=len(y)
            boot=y[rng.integers(0,n,(plan['bootstrap_draws'],n))].mean(1)
            meas={str(rho):float(np.sqrt(np.sum(t.euv_local**2+t.euv_global**2-2*rho*t.euv_local*t.euv_global))/n) for rho in [-.5,0,.5,.9]}
            out.append({'survey':survey,'z_bin':[lo,hi],'n':n,'mean_local_uv':float(t.uv_local.mean()),'mean_global_uv':float(t.uv_global.mean()),'mean_local_minus_global':float(y.mean()),'median_local_minus_global':float(np.median(y)),'SN_sampling_bootstrap95':np.quantile(boot,[.025,.975]).tolist(),'measurement_only_SE_by_assumed_rho':meas})
    common=d[d.redshift.between(.2,.4,inclusive='left')&d.survey.isin(['SDSS','SNLS5'])].copy();tests=[]
    for source,target in [('SDSS','SNLS5'),('SNLS5','SDSS')]:
        a=common[common.survey==source];b=common[common.survey==target];cols=['redshift','uv_global','stellar_mass'];x=a[cols].to_numpy();t=b[cols].to_numpy()
        w,diag=balance(x,t.mean(0));delta=a.local_minus_global.to_numpy();target_delta=b.local_minus_global.to_numpy()
        pred=np.column_stack([np.ones(len(t)),t])@np.linalg.lstsq(np.column_stack([np.ones(len(x)),x]),a.uv_local,rcond=None)[0]
        err=pred-b.uv_local.to_numpy();base=b.uv_global.to_numpy()-b.uv_local.to_numpy()
        sampling=[];failed=0
        for _ in range(plan['bootstrap_draws']):
            ai=rng.integers(0,len(a),len(a));bi=rng.integers(0,len(b),len(b))
            try:
                wb,_=balance(x[ai],t[bi].mean(0));sampling.append(float(wb@delta[ai]-target_delta[bi].mean()))
            except ValueError:failed+=1
        tests.append({'source':source,'target':target,'n_source':len(a),'n_target':len(b),'source_z_range':[float(a.redshift.min()),float(a.redshift.max())],'target_z_range':[float(b.redshift.min()),float(b.redshift.max())],**diag,'unweighted_source_local_global_offset':float(delta.mean()),'weighted_source_local_global_offset':float(w@delta),'target_local_global_offset':float(target_delta.mean()),'transport_minus_target':float(w@delta-target_delta.mean()),'transport_SN_sampling_bootstrap95':np.quantile(sampling,[.025,.975]).tolist(),'bootstrap_failed_balance':failed,'heldout_localUV_RMSE':float(np.sqrt(np.mean(err**2))),'heldout_mean_prediction_error':float(err.mean()),'global_as_local_RMSE':float(np.sqrt(np.mean(base**2))),'global_as_local_mean_prediction_error':float(base.mean())})
        a.assign(transport_weight=w).to_csv(W/f'roman-transport-{source}-to-{target}.csv',index=False)
    result={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'method':plan,'eligible_events':len(d),'descriptive_bins':out,'cross_survey_tests':tests,'interpretation':'These are observed selected-host proxy distributions conditioned on a published SED procedure. Neither color nor its between-bin change is converted to Gyr or a distance correction; age/dust/SFH degeneracy and survey selection remain.', 'input_sha256':{str(p.relative_to(ROOT)):sha(p)},'code_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'proxy-design.json']}}
    (O/'proxy-transport.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(tests,indent=2))
if __name__=='__main__':main()
