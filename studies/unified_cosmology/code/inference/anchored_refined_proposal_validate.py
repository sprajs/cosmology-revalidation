"""Independent refinement/lineage controls with every physical call forbidden."""
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

import anchored_refined_proposal as producer
from independence_proposal import FrozenMixture


def validate():
    import camb,cobaya.model,independence_runtime
    def forbidden(*a,**k):raise RuntimeError('This validation forbids physical/model calls.')
    original_affinity=os.sched_getaffinity(0);contracts=[];refusals=[]
    design=json.loads(producer.DESIGN.read_text())
    assert design['maximum_target_calls']==1200 and design['maximum_successful_native_fallbacks']==16
    assert design['maximum_evaluation_wall_seconds']==1800 and design['optimizer']['maximum_starts']==1
    assert producer.finite_hessian is producer.original.finite_hessian
    assert producer.finite_gradient is producer.original.finite_gradient
    assert producer.covariance_from_hessians is producer.original.covariance_from_hessians
    old_design=json.loads(producer.original.DESIGN.read_text())
    assert design['Hessian']==old_design['Hessian']
    for key in ['method','xtol','ftol','white_search_box','maximum_function_evaluations_per_start']:
        assert design['optimizer'][key]==old_design['optimizer'][key]
    proof=producer.original.validation_guard()
    base=producer.ROOT/'.work/unified-cosmology/anchored-local-validation';base.mkdir(exist_ok=True,parents=True)
    with ExitStack()as stack,tempfile.TemporaryDirectory(dir=base,prefix='refinement-test-')as scratch:
        scratch=Path(scratch)
        for name in ['get_results','get_transfer_functions','get_background']:
            stack.enter_context(patch.object(camb,name,forbidden))
        stack.enter_context(patch.object(cobaya.model,'get_model',side_effect=forbidden))
        # Replay the real stopped-ledger selection, source/asset identity only.
        parent_folder=producer.ROOT/next(iter(design['stopped_parent']));parent_folder=parent_folder.parent
        actual_manifest=json.loads((parent_folder/'manifest.json').read_text())
        ref_file=[producer.ROOT/p for p in actual_manifest['input_sha256']if p.endswith('proposal.npz')][0]
        actual_ref=FrozenMixture.read(ref_file,actual_manifest['input_sha256'][producer.sampling.relative(ref_file)])
        _,selected,pins=producer.stopped_parent(parent_folder,actual_manifest['arguments'],actual_ref)
        raw=[json.loads(x)for x in(parent_folder/'evaluations.jsonl').read_text().splitlines()]
        expected=sorted((x for x in raw if x['finite_target']),key=lambda x:(-x['target_logpost'],x['index']))[0]
        assert selected==expected
        changed=copy.deepcopy(actual_manifest['arguments']);changed['native_accuracy']=1
        try:producer.stopped_parent(parent_folder,changed,actual_ref)
        except AssertionError:refusals.append('changed_parent_target')
        else:raise AssertionError('Changed target accepted.')
        changed=copy.deepcopy(selected);changed['loglikes'][next(iter(changed['loglikes']))]+=.01
        try:producer.verify_density_replay(selected,changed,1e-7)
        except AssertionError:refusals.append('component_density_change_despite_total_unchanged')
        else:raise AssertionError('Changed component accepted.')
        changed=copy.deepcopy(selected);changed['target_logpost']=float('nan')
        try:producer.verify_density_replay(selected,changed,1e-7)
        except AssertionError:refusals.append('NaN_replay')
        else:raise AssertionError('NaN accepted.')
        names=['H0','x1','x2'];mean=np.array([68.,.1,-.2]);C=np.diag([.5,.2,.3])**2
        reference=FrozenMixture(mean,C,names);mode=np.array([.4,-.3,.2]);precision=np.diag([1.3,.8,1.1])
        origin=np.array([.2,-.15,.1]);point=dict(zip(names,map(float,mean+reference.L@origin)))
        ll=float(-.5*(origin-mode)@precision@(origin-mode))
        parent_best={'index':19,'point':point,'white':origin.tolist(),'finite_target':True,
                     'target_logpost':ll,'logprior':0.,'loglikes':{'synthetic':ll}}
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
        stack.enter_context(patch.object(producer,'stopped_parent',return_value=({},parent_best,{})))
        stack.enter_context(patch.object(independence_runtime,'native_initialize',return_value={'synthetic_only':True,'actual_calls':0}))
        mixture_error=None
        for kind in ['gaussian','native_budget','NaN','caught_CAMB_failure','changed_replay']:
            failure=scratch/(kind+'-CAMB-failure.jsonl')
            class FakeModel:
                likelihood={'synthetic':{}}
                parameterization=SimpleNamespace(sampled_params=lambda:names)
                def __init__(self):self.theory={'spectral_surrogate':SimpleNamespace(exact_calls=0,surrogate_calls=0)}
                def __enter__(self):return self
                def __exit__(self,*args):return False
                def logposterior(self,p,cached=False):
                    assert cached is False
                    t=self.theory['spectral_surrogate'];t.surrogate_calls+=1
                    if kind=='native_budget':t.exact_calls+=1
                    z=np.linalg.solve(reference.L,np.array([p[n]for n in names])-mean)
                    logp=float(-.5*(z-mode)@precision@(z-mode))
                    if kind=='NaN':logp=float('nan')
                    if kind=='changed_replay':logp+=.01
                    if kind=='caught_CAMB_failure':
                        failure.write_text(json.dumps({'parameters':p,'error':'synthetic recorded CAMB error'})+'\n');logp=-np.inf
                    return SimpleNamespace(logpost=logp,logprior=0.,loglikes=[logp])
            with patch.object(cobaya.model,'get_model',side_effect=lambda config:FakeModel()),patch.object(producer,'failure_log_path',return_value=failure):
                output=scratch/kind
                result=producer.optimize({'model':'synthetic','surrogate':str(producer.ROOT/design['surrogate']['path'])},reference_folder,parent_folder,output,min(original_affinity))
            rows=[json.loads(x)for x in(output/'evaluations.jsonl').read_text().splitlines()]
            assert rows[0]['point']==parent_best['point']
            if kind=='gaussian':
                assert result['proposal_frozen']and len(result['optimizer_runs'])==1
                candidate,record=producer.load(output,info)
                assert np.max(abs(candidate.mean-(mean+reference.L@mode)))<1e-5
                expected=reference.L@np.linalg.inv(precision)@reference.L.T
                assert np.max(abs(candidate.covariance-expected))<1e-7
                rng=np.random.default_rng(8274701);x=rng.normal(size=(128,3))+candidate.mean
                independent=np.logaddexp(np.log(.9)+multivariate_normal.logpdf(x,mean=candidate.mean,cov=1.05**2*candidate.covariance),
                    np.log(.1)+multivariate_t.logpdf(x,loc=candidate.mean,shape=2.4*candidate.covariance,df=5))
                mixture_error=float(np.max(abs(independent-candidate.logpdf(x))));assert mixture_error<1e-10
                changed_info=copy.deepcopy(info);changed_info['params']['H0']['prior']['max']=99
                try:producer.load(output,changed_info)
                except AssertionError:refusals.append('changed_candidate_target')
                else:raise AssertionError('Changed target accepted.')
                path=output/'proposal.npz';original=path.read_bytes();path.write_bytes(original+b'changed')
                try:producer.load(output,info)
                except AssertionError:refusals.append('changed_candidate_bytes')
                else:raise AssertionError('Changed NPZ accepted.')
                path.write_bytes(original)
            elif kind=='native_budget':
                assert result['status']=='stopped_declared_budget:maximum_successful_native_fallbacks'
                assert result['target_calls']==16 and not result['proposal_frozen']
            else:
                assert result['status'].startswith('failed_evaluation:')and not result['proposal_frozen']
                assert result['target_calls']==1
                if kind=='caught_CAMB_failure':
                    assert result['successful_native_fallbacks']==0
                    assert 'synthetic recorded CAMB error'in rows[0]['CAMB_failure_log_new_text']
            contracts.append({'case':kind,'status':result['status'],'synthetic_target_calls':result['target_calls'],
                              'proposal_frozen':result['proposal_frozen']})
    os.sched_setaffinity(0,original_affinity)
    return {'status':'passed_synthetic_anchored_refinement_controls_no_physical_calls',
        'actual_stopped_parent_rows':len(raw),'actual_best_saved_index':selected['index'],
        'original_hessian_gradient_covariance_kernel_identity':True,'original_validation_sha256':proof,
        'maximum_independent_logmixture_error':mixture_error,'synthetic_consumer_contracts':contracts,
        'refusal_controls':refusals,'dependency_sha256':producer.dependencies(),
        'source_sha256':{producer.sampling.relative(p):producer.sampling.digest(p)for p in[Path(__file__),Path(producer.__file__),producer.DESIGN]},
        'input_sha256':pins,'physical_calls':0,
        'scope':'One-start refinement and stopped-parent contracts tested. Physical APIs forbidden, factory and initialization mocked only for labelled Gaussian cases. No observational proposal or posterior certified.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    result=validate();a.output.parent.mkdir(exist_ok=True,parents=True);a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'physical_calls':0}))


if __name__=='__main__':main()
