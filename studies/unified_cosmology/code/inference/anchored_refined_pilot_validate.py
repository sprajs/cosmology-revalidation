"""Replay every toy pilot decision with independent normalized SciPy densities."""
import argparse
from contextlib import ExitStack
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal,multivariate_t

import anchored_refined_pilot as pilot
from independence_proposal import FrozenMixture


def validate():
    import camb,cobaya.model,independence_runtime
    def forbidden(*args,**kwargs):raise RuntimeError('No physical calls permitted.')
    names=['H0','x'];mean=np.array([68.,0.]);cov=np.array([[.3,.04],[.04,.7]])
    mixture=FrozenMixture(mean,cov,names)
    target_mean=np.array([68.2,.25]);target_cov=np.array([[.4,-.1],[-.1,.5]])
    info={'params':{n:{'prior':{'min':-100,'max':100}}for n in names},'likelihood':{'synthetic':{}}}
    target={'identity':'synthetic-pilot-only','configuration':info}
    provenance={'proposal_kind':pilot.proposal.KIND,'synthetic_only':True}
    available=os.sched_getaffinity(0);base=pilot.ROOT/'.work/unified-cosmology/anchored-local-validation';base.mkdir(exist_ok=True,parents=True)
    controls=[];density_error=0.;decisions=0
    with tempfile.TemporaryDirectory(dir=base,prefix='pilot-') as temporary,ExitStack() as stack:
        root=Path(temporary)
        for name in ['get_results','get_background','get_transfer_functions']:
            stack.enter_context(patch.object(camb,name,forbidden))
        stack.enter_context(patch.object(pilot.s,'runtime',return_value={'synthetic_only':True}))
        stack.enter_context(patch.object(pilot.s,'validation_guard',return_value={'synthetic_only':True}))
        stack.enter_context(patch.object(pilot.proposal,'validation_guard',return_value={'synthetic_only':True}))
        stack.enter_context(patch.object(pilot,'validation_guard',return_value={'synthetic_only':True}))
        stack.enter_context(patch.object(pilot.s,'configuration',return_value=info))
        stack.enter_context(patch.object(pilot.s,'identify',return_value=target))
        stack.enter_context(patch.object(pilot.proposal,'load',return_value=(mixture,provenance)))
        stack.enter_context(patch.object(independence_runtime,'native_initialize',return_value={'synthetic_only':True,'actual_calls':0}))
        for kind in ['Gaussian_with_prior_rejections','caught_CAMB_failure','successful_fallback_cap']:
            failure=root/(kind+'-failure.jsonl')
            class FakeModel:
                likelihood={'synthetic':{}}
                parameterization=SimpleNamespace(sampled_params=lambda:names)
                def __init__(self):self.theory={'spectral_surrogate':SimpleNamespace(exact_calls=0,surrogate_calls=0)}
                def __enter__(self):return self
                def __exit__(self,*args):return False
                def logposterior(self,point,cached=False):
                    assert cached is False
                    x=np.array([point[n]for n in names]);t=self.theory['spectral_surrogate'];t.surrogate_calls+=1
                    if kind=='successful_fallback_cap':t.exact_calls+=1
                    logp=float(multivariate_normal.logpdf(x,mean=target_mean,cov=target_cov))
                    if abs(x[1])>2:logp=-np.inf
                    if kind=='caught_CAMB_failure':
                        failure.write_text(json.dumps({'parameters':point,'error':'synthetic caught failure'})+'\n');logp=-np.inf
                    return SimpleNamespace(logpost=logp,logprior=0.,loglikes=[logp])
            with patch.object(cobaya.model,'get_model',side_effect=lambda configuration:FakeModel()),patch.object(pilot.proposal,'failure_log_path',return_value=failure):
                folder=root/kind;result=pilot.pilot({'model':'cpl'},root/'unused',folder,min(available))
            rows=[json.loads(line)for line in(folder/'records.jsonl').read_text().splitlines()]
            if kind=='Gaussian_with_prior_rejections':
                assert result['status']=='completed_declared_request_budget'and len(rows)==400
                with np.load(folder/'requests.npz',allow_pickle=False)as request:
                    points=request['points'];logq=np.logaddexp(np.log(.9)+multivariate_normal.logpdf(points,mean=mean,cov=1.05**2*cov),
                        np.log(.1)+multivariate_t.logpdf(points,loc=mean,shape=2.4*cov,df=5))
                    density_error=float(np.max(abs(logq-request['logq'])));assert density_error<1e-10
                    current=None;held=[];weights=[]
                    logtarget=multivariate_normal.logpdf(points,mean=target_mean,cov=target_cov)
                    logtarget[np.abs(points[:,1])>2]=-np.inf
                    for i,row in enumerate(rows):
                        accept=bool(np.isfinite(logtarget[i])and(current is None or request['logu'][i]<min(0.,logtarget[i]-logtarget[current]+logq[current]-logq[i])))
                        if accept:current=i
                        assert row['accepted']==accept and row['occupied_candidate']==current
                        if current is not None:held.append(current)
                        if np.isfinite(logtarget[i]):weights.append(logtarget[i]-logq[i])
                    w=np.exp(np.array(weights)-logsumexp(weights))
                    assert abs(result['candidate_raw_importance_ESS']-1/(w@w))<1e-8
                    assert result['occupied_states_recorded_including_rejections']==len(held)
                    assert result['qualified_scientific_measurement']is False
                    assert len(set(held))<len(held)and any(not r['finite_target']for r in rows)
                    decisions=len(rows)
            elif kind=='caught_CAMB_failure':
                assert result['status']=='failed_evaluation'and len(rows)==1 and result['successful_native_fallbacks']==0
                assert 'synthetic caught failure'in rows[0]['CAMB_failure_log_new_text']
            else:
                assert result['status']=='stopped_native_fallback_budget'and len(rows)==8
            controls.append({'case':kind,'status':result['status'],'synthetic_calls':len(rows),
                             'measurement_qualified':result['qualified_scientific_measurement']})
    os.sched_setaffinity(0,available)
    return {'status':'passed_synthetic_refined_anchored_pilot_controls','physical_calls':0,
        'independent_MH_decisions':decisions,'maximum_independent_logq_error':density_error,
        'controls':controls,'dependency_sha256':pilot.dependencies(),
        'source_sha256':{pilot.s.relative(p):pilot.s.digest(p)for p in [Path(__file__),Path(pilot.__file__),pilot.DESIGN]},
        'scope':'Physical APIs forbidden. Factory, qualification and initialization explicitly mocked for consumer mechanics only. No observational posterior or real efficiency claim.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    value=validate();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':value['status'],'physical_calls':0}))


if __name__=='__main__':main()
