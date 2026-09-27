"""Direct anchored-density mode/curvature proposal; never a measurement."""
import argparse
from contextlib import ExitStack
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import time
from unittest.mock import patch

import numpy as np
from scipy.linalg import solve_triangular
from scipy.optimize import minimize

import anchored_sampling as sampling
from anchored_proposal import read_reference
from independence_proposal import FrozenMixture
from target_identity import canonical

ROOT=sampling.ROOT; HERE=sampling.HERE
DESIGN=HERE/'anchored-local-proposal-design.json'
VALIDATION=ROOT/'studies/unified_cosmology/results/inference/anchored-local-proposal-validation.json'
SCHEMA='anchored-direct-density-proposal-v1'
KIND='anchored_direct_density_v1'
SCOPE='Numerical proposal only; no scientific posterior qualification follows from optimization.'


def dependencies():
    bound=sampling.dependencies()
    for path in [Path(__file__),DESIGN,HERE/'anchored_adapter.py',HERE/'anchored_proposal.py']:
        bound[sampling.relative(path)]=sampling.digest(path)
    return bound


def starts(reference):
    j=reference.names.index('H0')
    shift=2*reference.covariance[:,j]/np.sqrt(reference.covariance[j,j])
    physical=np.array([reference.mean,reference.mean+shift,reference.mean-shift])
    return solve_triangular(reference.L,(physical-reference.mean).T,lower=True).T


def verify_surrogate(arguments):
    expected=json.loads(DESIGN.read_text())['surrogate'];path=Path(arguments['surrogate']).resolve()
    assert path==(ROOT/expected['path']).resolve()
    assert sampling.digest(path)==expected['sha256']


def failure_log_path(arguments):
    # Exact path used by the frozen SpectralSurrogate.record_failure. Observe
    # it without modifying the model, its logging method or its return value.
    return Path(arguments['surrogate']).parent/f'surrogate-camb-failures-{os.getpid()}.jsonl'


def finite_hessian(function,point,step):
    """Central stencil, independently recomputed at each declared step size."""
    point=np.asarray(point,float);d=len(point);result=np.zeros((d,d));center=float(function(point))
    if not np.isfinite(center):raise ValueError('Nonfinite Hessian center.')
    directions=np.eye(d)*step
    for i in range(d):
        values=[float(function(point+directions[i])),float(function(point-directions[i]))]
        if not np.isfinite(values).all():raise ValueError('Nonfinite diagonal Hessian stencil.')
        result[i,i]=(sum(values)-2*center)/step**2
        for j in range(i):
            values=[float(function(point+si*directions[i]+sj*directions[j])) for si,sj in [(1,1),(1,-1),(-1,1),(-1,-1)]]
            if not np.isfinite(values).all():raise ValueError('Nonfinite cross Hessian stencil.')
            result[i,j]=result[j,i]=(values[0]-values[1]-values[2]+values[3])/(4*step**2)
    return result


def finite_gradient(function,point,step):
    point=np.asarray(point,float);directions=np.eye(len(point))*step
    values=np.array([(function(point+v)-function(point-v))/(2*step) for v in directions])
    if not np.isfinite(values).all():raise ValueError('Nonfinite gradient stencil.')
    return values


def covariance_from_hessians(hessians,reference,design,stencil_error=None):
    spec=design['Hessian'];details={'stencil_error':stencil_error,'scope':'Proposal covariance only.'}
    fallback=stencil_error is not None or len(hessians)!=2
    if not fallback:
        first,second=map(np.asarray,hessians)
        assert first.shape==second.shape==(reference.d,reference.d)
        assert np.isfinite(first).all() and np.isfinite(second).all()
        delta=float(np.linalg.norm(first-second)/max(np.linalg.norm(first),np.linalg.norm(second),1e-12))
        combined=(first+second)/2;values,vectors=np.linalg.eigh(combined)
        clipped=values<spec['minimum_curvature_eigenvalue']
        details.update(relative_Frobenius_step_disagreement=delta,curvature_eigenvalues=values.tolist(),
                       regularized_dimensions=int(clipped.sum()))
        fallback=(delta>spec['relative_Frobenius_agreement_maximum'] or int(clipped.sum())>spec['maximum_clipped_dimensions'])
        if not fallback:
            kept=np.maximum(values,spec['minimum_curvature_eigenvalue'])
            white=(vectors*(1/kept))@vectors.T
            covariance=reference.L@white@reference.L.T*spec['proposal_covariance_multiplier']
            details.update(status='two_step_curvature_proposal',regularized_eigenvalues=kept.tolist())
    if fallback:
        covariance=reference.covariance*spec['fallback_reference_covariance_multiplier']
        details.update(status='declared_broad_reference_covariance_fallback',
                       covariance_multiplier=spec['fallback_reference_covariance_multiplier'])
    covariance=(covariance+covariance.T)/2
    np.linalg.cholesky(covariance)
    return covariance,details


def validation_guard():
    value=json.loads(VALIDATION.read_text())
    assert value['status']=='passed_synthetic_anchored_direct_proposal_controls_no_physical_calls'
    assert value['dependency_sha256']==dependencies()
    sampling.verify_hashes(value['source_sha256'])
    assert value['physical_calls']==0
    return {sampling.relative(VALIDATION):sampling.digest(VALIDATION)}


def preflight(reference,arguments):
    """Three backgrounds/spectrum polynomials; explicitly forbid native spectra."""
    import camb
    from spectral_surrogate import coordinates,polynomial,SPECTRA
    verify_surrogate(arguments);info=sampling.configuration(arguments)
    extra=info['theory']['spectral_surrogate']['extra_args']
    with np.load(arguments['surrogate'],allow_pickle=False) as archive:
        centre=archive['centre'];chol=archive['coordinate_cholesky'];exponents=archive['exponents']
        coefficients=archive['coefficients'];scale=archive['output_scale'];length=int(archive['length'])
    def forbidden(*args,**kwargs):raise RuntimeError('Preflight cannot evaluate native spectra.')
    rows=[]
    with ExitStack() as stack:
        for name in ['get_results','get_transfer_functions']:stack.enter_context(patch.object(camb,name,forbidden))
        for white in starts(reference):
            point=dict(zip(reference.names,map(float,reference.mean+reference.L@white)))
            started=time.monotonic()
            try:
                pars=camb.set_params(H0=point['H0'],ombh2=point['ombh2'],omch2=point['omch2'],
                    As=1e-10*np.exp(point['logA']),ns=point['ns'],tau=point['tau'],
                    w=point.get('w',-1.),wa=point.get('wa',0.),**extra)
                background=camb.get_background(pars)
                coordinate=np.linalg.solve(chol,coordinates(pars,background)-centre)
                over=bool(np.max(abs(coordinate))>4.)
                invalid=[]
                if not over:
                    spectra=dict(zip(SPECTRA,((polynomial(coordinate,exponents)@coefficients)*scale).reshape(len(SPECTRA),length)))
                    invalid=[name for name in SPECTRA if not np.isfinite(spectra[name]).all()]
                    invalid+= [name+'_nonpositive' for name in ['tt','ee','pp'] if np.any(spectra[name][2:]<=0)]
                row={'status':'evaluated_background_polynomial','point':point,'spectral_coordinates':coordinate.tolist(),
                     'outside_envelope':over,'invalid_spectra':invalid,'would_request_native_fallback':over or bool(invalid)}
            except Exception as error:row={'status':'failed_background_preflight','point':point,'error':repr(error)}
            row['seconds']=time.monotonic()-started;rows.append(row)
    return {'status':'background_only_fallback_preflight','rows':rows,'native_spectrum_calls':0,
            'warm_state':'No native warmup in this preflight; optimization performs the separately recorded unchanged native initialization.'}


class BudgetStop(RuntimeError):pass


def optimize(arguments,reference_folder,output,cpu):
    environment=sampling.runtime();sampling.validation_guard();validation=validation_guard();bound=dependencies()
    verify_surrogate(arguments)
    design=json.loads(DESIGN.read_text());reference,reference_record=read_reference(reference_folder)
    info=sampling.configuration(arguments);target=sampling.identify(arguments)
    names=[k for k,v in info['params'].items() if isinstance(v,dict) and 'prior'in v]
    assert names==reference.names
    output=Path(output).resolve();assert output.is_relative_to(ROOT/'.work') and not output.exists()
    assert cpu in os.sched_getaffinity(0);os.sched_setaffinity(0,{cpu})
    output.mkdir(parents=True)
    inputs={sampling.relative(Path(reference_folder)/name):sampling.digest(Path(reference_folder)/name) for name in ['proposal.json','proposal.npz']}
    sampling.dump_new(output/'manifest.json',{'schema':SCHEMA,'arguments':arguments,'target_identity':target,
        'dependency_sha256':bound,'input_sha256':inputs,'validation':validation,'runtime':environment,
        'reference_role':design['reference_role'],'CPU_affinity':[cpu],'scope':SCOPE})
    pre=preflight(reference,arguments);sampling.dump_new(output/'fallback-preflight.json',pre)
    if any(row['status']!='evaluated_background_polynomial' for row in pre['rows']):
        result={'status':'failed_preflight_no_optimization','proposal_frozen':False,'preflight':pre}
        sampling.dump_new(output/'result.json',result);return result
    from independence_runtime import native_initialize
    initialization=native_initialize(info);sampling.dump_new(output/'native-initialization.json',initialization)
    from cobaya.model import get_model
    rows=[];runs=[];best=None;hessians=[];gradients=[];status='completed_declared_construction';stencil_error=None
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
            for index,start in enumerate(starts(reference)):
                evaluate(start)
                options=design['optimizer']
                result=minimize(lambda x:min(evaluate(x),1e100),start,method=options['method'],
                    bounds=[options['white_search_box']]*reference.d,
                    options={'maxfev':options['maximum_function_evaluations_per_start'],'xtol':options['xtol'],'ftol':options['ftol']})
                runs.append({'start_index':index,'success':bool(result.success),'message':str(result.message),
                             'nfev':int(result.nfev),'white':result.x.tolist(),'objective':float(result.fun)})
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
        'optimizer_runs':runs,'gradient_checks':[g.tolist() for g in gradients],'dependency_sha256':bound,'scope':SCOPE}
    if failure_path.exists():
        snapshot=output/'observed-CAMB-failure-log.jsonl';snapshot.write_bytes(failure_path.read_bytes())
        result['CAMB_failure_log_snapshot']={'path':sampling.relative(snapshot),'sha256':sampling.digest(snapshot),
                                            'original_path':sampling.relative(failure_path)}
    sampling.verify_hashes(bound);assert sampling.identify(arguments)==target
    if status=='completed_declared_construction':
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
            'target_identity':target,'training_settings':arguments,'names':mixture.names,'mode_record_index':mode_record_index,
            'curvature_diagnosis':details,'design':design,'scope':SCOPE,'posterior_qualification':False}
        sampling.dump_new(output/'proposal.json',record)
        result.update(proposal_frozen=True,proposal_sha256=record['proposal_sha256'],curvature_diagnosis=details)
    sampling.dump_new(output/'result.json',result)
    return result


def load(folder,info):
    folder=Path(folder).resolve();record=json.loads((folder/'proposal.json').read_text())
    assert record['schema']==SCHEMA and record['proposal_kind']==KIND and record['scope']==SCOPE
    verify_surrogate(record['training_settings'])
    assert record['source_sha256']==dependencies() and record['design']==json.loads(DESIGN.read_text())
    sampling.verify_hashes(record['source_sha256']);sampling.verify_hashes(record['input_sha256'])
    assert canonical(info)==record['target_identity']['configuration']
    assert sampling.identify(record['training_settings'])==record['target_identity']
    proposal=FrozenMixture.read(folder/'proposal.npz',record['proposal_sha256'])
    with np.load(folder/'curvature.npz',allow_pickle=False) as archive:
        reference=FrozenMixture(archive['reference_mean'],archive['reference_covariance'],archive['names'].tolist())
        covariance,details=covariance_from_hessians(list(archive['hessians']),reference,record['design'],record['curvature_diagnosis']['stencil_error'])
        assert np.array_equal(proposal.mean,reference.mean+reference.L@archive['center_white'])
    assert details==record['curvature_diagnosis'] and np.array_equal(proposal.covariance,covariance)
    assert proposal.names==record['names']
    return proposal,{'proposal_kind':KIND,'proposal_file':str(folder/'proposal.npz'),
        'proposal_sha256':record['proposal_sha256'],'proposal_record_path':sampling.relative(folder/'proposal.json'),
        'proposal_record_sha256':sampling.digest(folder/'proposal.json'),
        'training_target_identity':record['target_identity']['identity'],'training_settings':record['training_settings'],'scope':SCOPE}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=['preflight','optimize'])
    p.add_argument('--model',choices=['lcdm','cpl'],required=True);p.add_argument('--surrogate',type=Path,required=True)
    p.add_argument('--native-accuracy',type=int,choices=[1,2],default=2);p.add_argument('--reference-proposal',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--cpu',type=int)
    a=p.parse_args();arguments=sampling.settings(a.model,a.surrogate,a.native_accuracy)
    if a.action=='preflight':
        sampling.runtime();sampling.validation_guard();validation_guard()
        reference,_=read_reference(a.reference_proposal);result=preflight(reference,arguments)
        result.update(dependency_sha256=dependencies(),target_identity=sampling.identify(arguments))
        sampling.dump_new(a.output,result)
    else:
        assert a.cpu is not None;result=optimize(arguments,a.reference_proposal,a.output,a.cpu)
    print(json.dumps({k:result[k] for k in ['status','proposal_frozen','target_calls','successful_native_fallbacks'] if k in result}))


if __name__=='__main__':main()
