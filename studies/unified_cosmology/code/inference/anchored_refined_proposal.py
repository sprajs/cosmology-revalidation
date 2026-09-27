"""Separate one-start anchored proposal refinement; preserved original search."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import time

import numpy as np
from scipy.optimize import minimize

import anchored_local_proposal as original
from anchored_proposal import read_reference
from independence_proposal import FrozenMixture
from target_identity import canonical

sampling=original.sampling
ROOT=sampling.ROOT;HERE=sampling.HERE
DESIGN=HERE/'anchored-refined-proposal-design.json'
VALIDATION=ROOT/'studies/unified_cosmology/results/inference/anchored-refined-proposal-validation.json'
SCHEMA='anchored-refined-density-proposal-v1'
KIND='anchored_refined_density_v1'
SCOPE='Separate numerical refinement only; no posterior or unseen-mode qualification.'
finite_gradient=original.finite_gradient
finite_hessian=original.finite_hessian
covariance_from_hessians=original.covariance_from_hessians
failure_log_path=original.failure_log_path
BudgetStop=original.BudgetStop


def dependencies():
    bound=original.dependencies()
    for path in [Path(__file__),DESIGN]:bound[sampling.relative(path)]=sampling.digest(path)
    return bound


def validation_guard():
    value=json.loads(VALIDATION.read_text())
    assert value['status']=='passed_synthetic_anchored_refinement_controls_no_physical_calls'
    assert value['dependency_sha256']==dependencies()
    sampling.verify_hashes(value['source_sha256']);assert value['physical_calls']==0
    return {sampling.relative(VALIDATION):sampling.digest(VALIDATION)}


def stopped_parent(folder,arguments,reference):
    folder=Path(folder).resolve();design=json.loads(DESIGN.read_text())
    pins=design['stopped_parent'];sampling.verify_hashes(pins)
    assert all((ROOT/path).parent==folder for path in pins)
    manifest=json.loads((folder/'manifest.json').read_text())
    result=json.loads((folder/'result.json').read_text())
    rows=[json.loads(line)for line in(folder/'evaluations.jsonl').read_text().splitlines()]
    assert manifest['schema']==original.SCHEMA
    assert manifest['arguments']==arguments
    assert manifest['dependency_sha256']==original.dependencies()
    sampling.verify_hashes(manifest['dependency_sha256']);sampling.verify_hashes(manifest['input_sha256'])
    assert manifest['target_identity']==sampling.identify(arguments)
    assert result['status']=='stopped_declared_budget:maximum_successful_native_fallbacks'
    assert result['proposal_frozen']is False and result['target_calls']==len(rows)
    assert not (folder/'proposal.json').exists() and not (folder/'proposal.npz').exists()
    for i,row in enumerate(rows):
        assert row['index']==i and row.get('status')!='failed_target_evaluation'
        reconstructed=dict(zip(reference.names,map(float,reference.mean+reference.L@np.asarray(row['white']))))
        assert row['point']==reconstructed
        if row['finite_target']:
            assert np.isfinite([row['target_logpost'],row['logprior'],*row['loglikes'].values()]).all()
            assert abs(row['target_logpost']-row['logprior']-sum(row['loglikes'].values()))<1e-9
    best=max((row for row in rows if row['finite_target']),key=lambda row:row['target_logpost'])
    return manifest,best,pins


def verify_density_replay(old,new,tolerance):
    assert old['point']==new['point'] and old['finite_target'] and new['finite_target']
    assert old['loglikes'].keys()==new['loglikes'].keys()
    a=np.array([old['target_logpost'],old['logprior'],*old['loglikes'].values()])
    b=np.array([new['target_logpost'],new['logprior'],*new['loglikes'].values()])
    assert np.isfinite(a).all() and np.isfinite(b).all()
    delta=float(np.max(abs(a-b)));assert delta<=tolerance,'Fresh starting density changed.'
    return {'status':'passed_single_fresh_parent_density_replay','maximum_absolute_density_difference':delta,
            'absolute_tolerance':tolerance,'parent_index':old['index'],'new_index':new['index']}


def optimize(arguments,reference_folder,parent_folder,output,cpu):
    environment=sampling.runtime();sampling.validation_guard();validation=validation_guard();bound=dependencies()
    original.verify_surrogate(arguments)
    design=json.loads(DESIGN.read_text());reference,reference_record=read_reference(reference_folder)
    info=sampling.configuration(arguments);target=sampling.identify(arguments)
    parent,best_parent,parent_bindings=stopped_parent(parent_folder,arguments,reference)
    names=[k for k,v in info['params'].items() if isinstance(v,dict) and 'prior'in v]
    assert names==reference.names
    output=Path(output).resolve();assert output.is_relative_to(ROOT/'.work') and not output.exists()
    assert cpu in os.sched_getaffinity(0);os.sched_setaffinity(0,{cpu})
    output.mkdir(parents=True)
    inputs={sampling.relative(Path(reference_folder)/name):sampling.digest(Path(reference_folder)/name) for name in ['proposal.json','proposal.npz']}
    inputs.update(parent_bindings)
    sampling.dump_new(output/'manifest.json',{'schema':SCHEMA,'arguments':arguments,'target_identity':target,
        'dependency_sha256':bound,'input_sha256':inputs,'validation':validation,'runtime':environment,
        'reference_role':design['reference_role'],'CPU_affinity':[cpu],'scope':SCOPE,
        'stopped_parent_folder':sampling.relative(Path(parent_folder).resolve()),'parent_best_index':best_parent['index']})
    from independence_runtime import native_initialize
    initialization=native_initialize(info);sampling.dump_new(output/'native-initialization.json',initialization)
    from cobaya.model import get_model
    rows=[];runs=[];best=None;hessians=[];gradients=[];status='completed_declared_refinement';stencil_error=None;replay=None
    started=time.monotonic()
    with get_model(info) as model:
        assert list(model.parameterization.sampled_params())==reference.names
        theory=model.theory['spectral_surrogate']
        failure_path=failure_log_path(arguments)
        def evaluate(white):
            nonlocal best
            if len(rows)>=design['maximum_target_calls']:raise BudgetStop('maximum_target_calls')
            if time.monotonic()-started>=design['maximum_evaluation_wall_seconds']:raise BudgetStop('maximum_evaluation_wall_seconds')
            if theory.exact_calls>=design['maximum_successful_native_fallbacks']:raise BudgetStop('maximum_successful_native_fallbacks')
            point=dict(zip(reference.names,map(float,reference.mean+reference.L@np.asarray(white))))
            tick=time.monotonic()
            failure_before=failure_path.read_bytes() if failure_path.exists() else b''
            failure_delta=b''
            try:
                result=model.logposterior(point,cached=False);value=float(result.logpost)
                failure_after=failure_path.read_bytes() if failure_path.exists() else b''
                assert failure_after.startswith(failure_before),'CAMB failure log was changed/truncated.'
                failure_delta=failure_after[len(failure_before):]
                if failure_delta:
                    raise RuntimeError('Frozen surrogate recorded a CAMB failure; construction stops even if Cobaya returned -inf.')
                assert not np.isnan(value) and value!=np.inf
                finite=bool(np.isfinite(value));closure=None
                if finite:
                    closure=value-float(result.logprior)-sum(result.loglikes);assert abs(closure)<1e-9
                row={'index':len(rows),'point':point,'white':np.asarray(white).tolist(),'finite_target':finite,
                    'target_logpost':value if finite else None,'logprior':float(result.logprior) if np.isfinite(result.logprior) else None,
                    'loglikes':{k:float(v) if np.isfinite(v) else None for k,v in zip(model.likelihood,result.loglikes)},
                    'density_accounting_closure':closure,'successful_native_fallbacks':int(theory.exact_calls),
                    'surrogate_calls':int(theory.surrogate_calls),'seconds':time.monotonic()-tick}
            except Exception as error:
                row={'index':len(rows),'point':point,'white':np.asarray(white).tolist(),
                     'status':'failed_target_evaluation','error':repr(error),'seconds':time.monotonic()-tick,
                     'CAMB_failure_log_new_text':failure_delta.decode(errors='replace')}
                rows.append(row)
                with (output/'evaluations.jsonl').open('a') as stream:stream.write(json.dumps(row,allow_nan=False)+'\n')
                raise
            rows.append(row)
            with (output/'evaluations.jsonl').open('a') as stream:stream.write(json.dumps(row,allow_nan=False)+'\n')
            if finite and (best is None or value>best['target_logpost']):best=row
            return -value if finite else np.inf
        try:
            start=np.asarray(best_parent['white'])
            evaluate(start)
            replay=verify_density_replay(best_parent,rows[-1],design['density_replay_absolute_tolerance'])
            options=design['optimizer']
            result=minimize(lambda x:min(evaluate(x),1e100),start,method=options['method'],
                bounds=[options['white_search_box']]*reference.d,
                options={'maxfev':options['maximum_function_evaluations_per_start'],'xtol':options['xtol'],'ftol':options['ftol']})
            runs.append({'start_index':0,'parent_best_index':best_parent['index'],'success':bool(result.success),
                         'message':str(result.message),'nfev':int(result.nfev),'white':result.x.tolist(),'objective':float(result.fun)})
            assert best is not None,'No finite anchored target point found.'
            # Hold the chosen center fixed: stencil evaluations must not move it.
            center=np.asarray(best['white']);mode_record_index=best['index']
            try:
                for step in design['Hessian']['central_difference_steps']:
                    gradients.append(finite_gradient(evaluate,center,step))
                    hessians.append(finite_hessian(evaluate,center,step))
                if max(np.max(abs(g)) for g in gradients)>design['Hessian']['gradient_maximum_absolute']:
                    stencil_error='declared_stationarity_check_failed'
            except ValueError as error:stencil_error=repr(error)
        except BudgetStop as error:status='stopped_declared_budget:'+str(error)
        except Exception as error:status='failed_evaluation:'+repr(error)
        fallback_count=int(theory.exact_calls);surrogate_count=int(theory.surrogate_calls)
    result={'status':status,'proposal_frozen':False,'target_calls':len(rows),'seconds':time.monotonic()-started,
        'successful_native_fallbacks':fallback_count,'surrogate_calls':surrogate_count,'explicit_native_initialization_calls':1,
        'density_replay':replay,'parent_best_index':best_parent['index'],'optimizer_runs':runs,'gradient_checks':[g.tolist() for g in gradients],'dependency_sha256':bound,'scope':SCOPE}
    if failure_path.exists():
        snapshot=output/'observed-CAMB-failure-log.jsonl';snapshot.write_bytes(failure_path.read_bytes())
        result['CAMB_failure_log_snapshot']={'path':sampling.relative(snapshot),'sha256':sampling.digest(snapshot),
                                            'original_path':sampling.relative(failure_path)}
    sampling.verify_hashes(bound);assert sampling.identify(arguments)==target
    if status=='completed_declared_refinement':
        covariance,details=covariance_from_hessians(hessians,reference,design,stencil_error)
        mean=reference.mean+reference.L@center;mixture=FrozenMixture(mean,covariance,reference.names)
        np.savez(output/'proposal.npz',mean=mixture.mean,cov=mixture.covariance,names=mixture.names)
        np.savez(output/'curvature.npz',center_white=center,hessians=np.asarray(hessians),reference_mean=reference.mean,
                 reference_covariance=reference.covariance,names=reference.names)
        bindings=dict(inputs)
        for path in output.iterdir():
            if path.is_file():bindings[sampling.relative(path)]=sampling.digest(path)
        record={'schema':SCHEMA,'proposal_kind':KIND,'created_utc':datetime.now(timezone.utc).isoformat(),
            'proposal_sha256':sampling.digest(output/'proposal.npz'),'source_sha256':bound,'input_sha256':bindings,
            'target_identity':target,'training_settings':arguments,'stopped_parent_folder':sampling.relative(Path(parent_folder).resolve()),
            'parent_best_index':best_parent['index'],'density_replay':replay,'names':mixture.names,'mode_record_index':mode_record_index,
            'curvature_diagnosis':details,'design':design,'scope':SCOPE,'posterior_qualification':False}
        sampling.dump_new(output/'proposal.json',record)
        result.update(proposal_frozen=True,proposal_sha256=record['proposal_sha256'],curvature_diagnosis=details)
    sampling.dump_new(output/'result.json',result)
    return result


def load(folder,info):
    folder=Path(folder).resolve();record=json.loads((folder/'proposal.json').read_text())
    assert record['schema']==SCHEMA and record['proposal_kind']==KIND and record['scope']==SCOPE
    original.verify_surrogate(record['training_settings'])
    assert record['source_sha256']==dependencies() and record['design']==json.loads(DESIGN.read_text())
    sampling.verify_hashes(record['source_sha256']);sampling.verify_hashes(record['input_sha256'])
    assert canonical(info)==record['target_identity']['configuration']
    assert sampling.identify(record['training_settings'])==record['target_identity']
    candidate=FrozenMixture.read(folder/'proposal.npz',record['proposal_sha256'])
    with np.load(folder/'curvature.npz',allow_pickle=False)as archive:
        reference=FrozenMixture(archive['reference_mean'],archive['reference_covariance'],archive['names'].tolist())
        covariance,details=covariance_from_hessians(list(archive['hessians']),reference,record['design'],record['curvature_diagnosis']['stencil_error'])
        assert np.array_equal(candidate.mean,reference.mean+reference.L@archive['center_white'])
    _,parent_best,_=stopped_parent(ROOT/record['stopped_parent_folder'],record['training_settings'],reference)
    rows=[json.loads(line)for line in(folder/'evaluations.jsonl').read_text().splitlines()]
    assert record['parent_best_index']==parent_best['index']
    assert record['density_replay']==verify_density_replay(parent_best,rows[0],record['design']['density_replay_absolute_tolerance'])
    assert details==record['curvature_diagnosis'] and np.array_equal(candidate.covariance,covariance)
    assert candidate.names==record['names']
    return candidate,{'proposal_kind':KIND,'proposal_file':str(folder/'proposal.npz'),
        'proposal_sha256':record['proposal_sha256'],'proposal_record_path':sampling.relative(folder/'proposal.json'),
        'proposal_record_sha256':sampling.digest(folder/'proposal.json'),
        'training_target_identity':record['target_identity']['identity'],'training_settings':record['training_settings'],'scope':SCOPE}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',choices=['lcdm','cpl'],required=True);p.add_argument('--surrogate',type=Path,required=True)
    p.add_argument('--native-accuracy',type=int,choices=[1,2],default=2)
    p.add_argument('--reference-proposal',type=Path,required=True);p.add_argument('--parent-folder',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--cpu',type=int,required=True)
    a=p.parse_args();result=optimize(sampling.settings(a.model,a.surrogate,a.native_accuracy),a.reference_proposal,a.parent_folder,a.output,a.cpu)
    print(json.dumps({k:result[k]for k in ['status','proposal_frozen','target_calls','successful_native_fallbacks']}))


if __name__=='__main__':main()
