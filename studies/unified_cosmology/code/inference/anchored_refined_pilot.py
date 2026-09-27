"""Bounded MH diagnostic for the separately identified anchored refinement."""
import argparse
import json
import os
from pathlib import Path
import time
import numpy as np

import anchored_refined_proposal as proposal
from anchored_pilot import requests,transition,metrics

s=proposal.sampling;HERE=s.HERE;ROOT=s.ROOT
DESIGN=HERE/'anchored-refined-pilot-design.json'
SCHEMA='anchored-refined-proposal-pilot-v1'
VALIDATION=ROOT/'studies/unified_cosmology/results/inference/anchored-refined-pilot-validation.json'


def dependencies():
    bound=proposal.dependencies()
    for path in [Path(__file__),DESIGN,HERE/'anchored_pilot.py']:
        bound[s.relative(path)]=s.digest(path)
    return bound


def validation_guard():
    checked=json.loads(VALIDATION.read_text())
    assert checked['status']=='passed_synthetic_refined_anchored_pilot_controls'
    assert checked['dependency_sha256']==dependencies()
    s.verify_hashes(checked['source_sha256'])
    assert checked['physical_calls']==0
    return {s.relative(VALIDATION):s.digest(VALIDATION)}


def pilot(arguments,proposal_folder,output,cpu):
    s.runtime();s.validation_guard();proposal.validation_guard();checked=validation_guard()
    bound=dependencies();design=json.loads(DESIGN.read_text())
    info=s.configuration(arguments);target=s.identify(arguments)
    mixture,frozen=proposal.load(proposal_folder,info)
    assert frozen['proposal_kind']==proposal.KIND
    output=Path(output).resolve();assert output.is_relative_to(ROOT/'.work') and not output.exists()
    assert cpu in os.sched_getaffinity(0);os.sched_setaffinity(0,{cpu});output.mkdir(parents=True)
    request=requests(mixture,arguments['model'],design)
    path=output/'requests.npz';np.savez_compressed(path,**request)
    s.dump_new(output/'manifest.json',{'schema':SCHEMA,'arguments':arguments,
        'target_identity':target,'proposal':frozen,'dependency_sha256':bound,'validation':checked,
        'requests_sha256':s.digest(path),'CPU_affinity':[cpu],'scope':design['scope']})
    from independence_runtime import native_initialize
    warmup=native_initialize(info);s.dump_new(output/'native-initialization.json',warmup)
    from cobaya.model import get_model
    rows=[];current=None;current_logp=None;current_logq=None
    status='completed_declared_request_budget';started=time.monotonic()
    failure_path=proposal.failure_log_path(arguments)
    with get_model(info) as model:
        assert list(model.parameterization.sampled_params())==mixture.names
        theory=model.theory['spectral_surrogate']
        for index,x in enumerate(request['points']):
            if theory.exact_calls>=design['maximum_successful_native_fallbacks']:
                status='stopped_native_fallback_budget';break
            if time.monotonic()-started>=design['maximum_evaluation_wall_seconds']:
                status='stopped_wall_budget_after_completed_call';break
            s.verify_hashes(bound);point=dict(zip(mixture.names,map(float,x)));tick=time.monotonic()
            before=failure_path.read_bytes() if failure_path.exists() else b'';delta=b''
            try:
                value=model.logposterior(point,cached=False);logp=float(value.logpost)
                after=failure_path.read_bytes() if failure_path.exists() else b''
                assert after.startswith(before),'Recorded CAMB failure log changed/truncated.'
                delta=after[len(before):]
                if delta:raise RuntimeError('Frozen surrogate recorded CAMB failure; pilot stops.')
                assert not np.isnan(logp) and logp!=np.inf
                closure=logp-float(value.logprior)-sum(value.loglikes) if np.isfinite(logp) else None
                if closure is not None:assert abs(closure)<1e-9
                accept,alpha=transition(logp,current_logp,current_logq,float(request['logq'][index]),float(request['logu'][index]))
                if accept:current=index;current_logp=logp;current_logq=float(request['logq'][index])
                row={'index':index,'point':point,'component':str(request['components'][index]),
                    'logq':float(request['logq'][index]),'log_uniform':float(request['logu'][index]),
                    'target_logpost':logp if np.isfinite(logp) else None,'finite_target':bool(np.isfinite(logp)),
                    'loglikes':{k:float(v) if np.isfinite(v) else None for k,v in zip(model.likelihood,value.loglikes)},
                    'logprior':float(value.logprior) if np.isfinite(value.logprior) else None,
                    'density_accounting_closure':closure,'accepted':accept,'log_acceptance':alpha,'occupied_candidate':current,
                    'successful_native_fallbacks':int(theory.exact_calls),'surrogate_calls':int(theory.surrogate_calls)}
            except Exception as error:
                row={'index':index,'point':point,'status':'failed_evaluation','error':repr(error),
                     'CAMB_failure_log_new_text':delta.decode(errors='replace'),'occupied_candidate':None}
                status='failed_evaluation'
            row['seconds']=time.monotonic()-tick;rows.append(row)
            with (output/'records.jsonl').open('a') as stream:stream.write(json.dumps(row,allow_nan=False)+'\n')
            if status=='failed_evaluation':break
            if(index+1)%25==0:print(json.dumps({'target_calls':index+1,'seconds':time.monotonic()-started}),flush=True)
        fallback_count=int(theory.exact_calls);surrogate_count=int(theory.surrogate_calls)
    if failure_path.exists():(output/'observed-CAMB-failure-log.jsonl').write_bytes(failure_path.read_bytes())
    s.verify_hashes(bound);assert dependencies()==bound and s.identify(arguments)==target
    _,after=proposal.load(proposal_folder,info);assert after==frozen
    result=metrics(rows,request)
    result.update(status=status,seconds=time.monotonic()-started,target_identity=target['identity'],
        successful_native_fallbacks=fallback_count,surrogate_calls=surrogate_count,
        explicit_native_initialization_calls=1,dependency_sha256=bound,
        input_sha256={s.relative(p):s.digest(p) for p in output.iterdir() if p.is_file()},limits=design['limits'])
    s.dump_new(output/'result.json',result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',choices=['lcdm','cpl'],required=True);p.add_argument('--surrogate',type=Path,required=True)
    p.add_argument('--native-accuracy',type=int,choices=[1,2],default=2)
    p.add_argument('--proposal-folder',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cpu',type=int,required=True);a=p.parse_args()
    result=pilot(s.settings(a.model,a.surrogate,a.native_accuracy),a.proposal_folder,a.output,a.cpu)
    print(json.dumps({k:result[k] for k in ['status','target_calls_recorded','finite_target_calls','successful_native_fallbacks']}))


if __name__=='__main__':main()
