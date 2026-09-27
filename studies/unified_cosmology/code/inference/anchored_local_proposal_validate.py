"""Independent Gaussian algebra and mocked-consumer controls; no physical calls."""
import argparse
from contextlib import ExitStack
import copy
import json
import os
from pathlib import Path
import tempfile
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
from scipy.stats import multivariate_normal,multivariate_t

import anchored_local_proposal as producer
from independence_proposal import FrozenMixture


def validate():
    rng=np.random.default_rng(8274301);design=json.loads(producer.DESIGN.read_text())
    errors=[];covariance_errors=[];density_errors=[];fallbacks=[]
    for d in [3,7,13]:
        L=np.tril(rng.normal(size=(d,d))*.08)+np.eye(d)
        reference=FrozenMixture(rng.normal(size=d),L@L.T,['H0']+[f'x{i}'for i in range(1,d)])
        Q=rng.normal(size=(d,d));precision=Q@Q.T/d+np.eye(d)*.5
        centre=rng.normal(size=d)*.2
        fun=lambda x:.5*(x-centre)@precision@(x-centre)+13.4
        hs=[producer.finite_hessian(fun,centre,s)for s in design['Hessian']['central_difference_steps']]
        gs=[producer.finite_gradient(fun,centre,s)for s in design['Hessian']['central_difference_steps']]
        errors.append(max(np.max(abs(h-precision))for h in hs))
        assert max(np.max(abs(g))for g in gs)<1e-10
        covariance,details=producer.covariance_from_hessians(hs,reference,design)
        expected=reference.L@np.linalg.inv(precision)@reference.L.T
        covariance_errors.append(np.max(abs(covariance-expected)))
        assert details['status']=='two_step_curvature_proposal'
        candidate=FrozenMixture(reference.mean+reference.L@centre,covariance,reference.names)
        points=rng.normal(size=(64,d))+candidate.mean
        g=multivariate_normal.logpdf(points,mean=candidate.mean,cov=1.05**2*covariance)
        t=multivariate_t.logpdf(points,loc=candidate.mean,shape=2.4*covariance,df=5)
        independent=np.logaddexp(np.log(.9)+g,np.log(.1)+t)
        density_errors.append(np.max(abs(candidate.logpdf(points)-independent)))
        for case in [('nonfinite',[], 'nonfinite_stencil'),('step_disagreement',[precision,precision*2],None),
                     ('too_many_flat',[np.zeros((d,d))]*2,None)]:
            cov,diagnosis=producer.covariance_from_hessians(case[1],reference,design,case[2])
            assert np.allclose(cov,4*reference.covariance)and diagnosis['status']=='declared_broad_reference_covariance_fallback'
            fallbacks.append({'dimension':d,'case':case[0]})
    assert max(errors)<1e-8 and max(covariance_errors)<1e-8 and max(density_errors)<1e-10
    # A zero curvature in one auxiliary direction is explicitly regularized,
    # never interpreted as a measured posterior uncertainty.
    reference=FrozenMixture(np.array([68.,1.,1.]),np.diag([.5,.1,.2])**2,['H0','x1','x2'])
    hs=[np.diag([2.,1.,0.])]*2
    cov,details=producer.covariance_from_hessians(hs,reference,design)
    assert details['regularized_dimensions']==1 and np.linalg.eigvalsh(cov).min()>0
    starts=producer.starts(reference)
    actual=reference.mean+starts@reference.L.T
    assert actual[1,0]-actual[0,0]==1. and actual[2,0]-actual[0,0]==-1.
    contracts=[];guard_refusals=0
    import camb,cobaya.model,independence_runtime
    def forbidden(*args,**kwargs):raise RuntimeError('No physical call allowed by this validator.')
    available=os.sched_getaffinity(0)
    base=producer.ROOT/'.work/unified-cosmology/anchored-local-validation';base.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=base,prefix='synthetic-') as scratch,ExitStack() as stack:
        scratch=Path(scratch)
        for name in ['get_results','get_transfer_functions','get_background']:
            stack.enter_context(patch.object(camb,name,forbidden))
        names=['H0','x1','x2'];mean=np.array([68.,.1,-.2]);C=np.diag([.5,.2,.3])**2
        reference=FrozenMixture(mean,C,names);mode=np.array([.4,-.3,.2]);precision=np.diag([1.3,.8,1.1])
        info={'params':{n:{'prior':{'min':-100,'max':100}}for n in names},'likelihood':{'synthetic':{}}}
        target={'identity':'synthetic-only-target','configuration':info}
        reference_folder=scratch/'reference';reference_folder.mkdir()
        (reference_folder/'proposal.json').write_text('{}\n');np.savez(reference_folder/'proposal.npz',mean=mean,cov=C,names=names)
        stack.enter_context(patch.object(producer,'read_reference',return_value=(reference,{})))
        stack.enter_context(patch.object(producer.sampling,'runtime',return_value={'synthetic_only':True}))
        stack.enter_context(patch.object(producer.sampling,'validation_guard',return_value={'synthetic_only':True}))
        stack.enter_context(patch.object(producer,'validation_guard',return_value={'synthetic_only':True}))
        stack.enter_context(patch.object(producer.sampling,'configuration',return_value=info))
        stack.enter_context(patch.object(producer.sampling,'identify',return_value=target))
        stack.enter_context(patch.object(producer,'preflight',return_value={'rows':[{'status':'evaluated_background_polynomial'}]*3,'synthetic_only':True}))
        stack.enter_context(patch.object(independence_runtime,'native_initialize',return_value={'synthetic_only':True,'actual_calls':0}))
        for kind in ['gaussian','native_budget','NaN','caught_CAMB_failure']:
            fake_failure_path=scratch/(kind+'-CAMB-failure.jsonl')
            class FakeModel:
                likelihood={'synthetic':{}}
                parameterization=SimpleNamespace(sampled_params=lambda:names)
                def __init__(self):self.theory={'spectral_surrogate':SimpleNamespace(exact_calls=0,surrogate_calls=0)}
                def __enter__(self):return self
                def __exit__(self,*args):return False
                def logposterior(self,point,cached=False):
                    assert cached is False
                    t=self.theory['spectral_surrogate'];t.surrogate_calls+=1
                    if kind=='native_budget':t.exact_calls+=1
                    z=np.linalg.solve(reference.L,np.array([point[n]for n in names])-mean)
                    logp=float(-.5*(z-mode)@precision@(z-mode))
                    if kind=='NaN':logp=float('nan')
                    if kind=='caught_CAMB_failure':
                        fake_failure_path.write_text(json.dumps({'parameters':point,'error':'synthetic recorded CAMB error'})+'\n')
                        logp=-np.inf
                    return SimpleNamespace(logpost=logp,logprior=0.,loglikes=[logp])
            with patch.object(cobaya.model,'get_model',side_effect=lambda config:FakeModel()),patch.object(producer,'failure_log_path',return_value=fake_failure_path):
                output=scratch/kind
                result=producer.optimize({'model':'synthetic','surrogate':str(producer.ROOT/design['surrogate']['path'])},reference_folder,output,min(available))
            if kind=='gaussian':
                assert result['proposal_frozen']
                proposal,record=producer.load(output,info)
                assert np.max(abs(proposal.mean-(mean+reference.L@mode)))<1e-5
                assert np.max(abs(proposal.covariance-reference.L@np.linalg.inv(precision)@reference.L.T))<1e-7
                for tamper in ['configuration','proposal_bytes']:
                    rejected=False
                    if tamper=='configuration':
                        changed=copy.deepcopy(info);changed['params']['H0']['prior']['max']=99
                        try:producer.load(output,changed)
                        except AssertionError:rejected=True
                    else:
                        path=output/'proposal.npz';original=path.read_bytes();path.write_bytes(original+b'changed')
                        try:producer.load(output,info)
                        except AssertionError:rejected=True
                        finally:path.write_bytes(original)
                    assert rejected;guard_refusals+=1
            elif kind=='native_budget':
                assert result['status']=='stopped_declared_budget:maximum_successful_native_fallbacks'
                assert result['target_calls']==8 and not result['proposal_frozen']
            else:
                assert result['status'].startswith('failed_evaluation:')and not result['proposal_frozen']
                row=json.loads((output/'evaluations.jsonl').read_text())
                assert row['status']=='failed_target_evaluation'and result['target_calls']==1
                if kind=='caught_CAMB_failure':
                    assert 'synthetic recorded CAMB error'in row['CAMB_failure_log_new_text']
                    assert result['successful_native_fallbacks']==0 and 'CAMB_failure_log_snapshot'in result
            contracts.append({'case':kind,'status':result['status'],'synthetic_target_calls':result['target_calls'],
                              'proposal_frozen':result['proposal_frozen']})
    os.sched_setaffinity(0,available)
    return {'status':'passed_synthetic_anchored_direct_proposal_controls_no_physical_calls',
        'Gaussian_dimensions':[3,7,13],'maximum_Hessian_absolute_error':float(max(errors)),
        'maximum_covariance_absolute_error':float(max(covariance_errors)),
        'maximum_independent_logmixture_error':float(max(density_errors)),
        'explicit_fallback_controls':fallbacks,'synthetic_consumer_contracts':contracts,'tampering_refusals':guard_refusals,
        'dependency_sha256':producer.dependencies(),
        'source_sha256':{producer.sampling.relative(p):producer.sampling.digest(p)for p in [Path(__file__),Path(producer.__file__),producer.DESIGN]},
        'physical_calls':0,'scope':'All physical entrypoints forbidden; target/factory/initialization/qualification mocked only for labelled synthetic consumer controls. No actual proposal, optimization, timing or cosmology is certified.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=validate();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'physical_calls':0}))


if __name__=='__main__':main()
