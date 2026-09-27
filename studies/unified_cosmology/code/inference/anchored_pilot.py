"""Bounded anchored proposal-efficiency pilot; never produces a measurement."""
import argparse
import json
import os
from pathlib import Path
import time
import numpy as np
from scipy.special import logsumexp

import anchored_sampling as s

HERE=s.HERE;ROOT=s.ROOT
DESIGN=HERE/'anchored-pilot-design.json'
VALIDATION=ROOT/'studies/unified_cosmology/results/inference/anchored-pilot-validation.json'


def dependencies():
    bound=s.dependencies()
    for path in [Path(__file__),DESIGN,HERE/'anchored_pilot_validate.py']:
        bound[s.relative(path)]=s.digest(path)
    return bound


def transition(logp,current_logp,current_logq,trial_logq,log_uniform):
    assert not np.isnan(logp) and logp!=np.inf
    assert np.isfinite(trial_logq) and np.isfinite(log_uniform) and log_uniform<0
    if current_logp is None:return bool(np.isfinite(logp)),None
    assert np.isfinite(current_logp) and np.isfinite(current_logq)
    alpha=min(0.,logp-current_logp+current_logq-trial_logq)
    return bool(np.isfinite(logp) and log_uniform<alpha),float(alpha) if np.isfinite(alpha) else None


def requests(mixture,model,design):
    seeds=design['seed'][model]
    rng=np.random.default_rng(seeds['proposals'])
    candidates=[mixture.draw(rng) for _ in range(design['maximum_target_calls'])]
    points=np.array([p for p,_ in candidates]);components=np.array([c for _,c in candidates])
    logq=np.asarray(mixture.logpdf(points));logu=np.log(np.random.default_rng(seeds['acceptance']).random(len(points)))
    assert np.isfinite(points).all() and np.isfinite(logq).all() and np.isfinite(logu).all()
    return {'points':points,'components':components,'logq':logq,'logu':logu,'names':np.asarray(mixture.names)}


def metrics(rows,request):
    finite=[r for r in rows if r.get('finite_target')]
    occupied=np.array([r['occupied_candidate'] for r in rows if r.get('occupied_candidate') is not None],dtype=int)
    result={'scope':'Proposal-efficiency diagnostic only; no posterior or convergence qualification.',
            'target_calls_recorded':len(rows),'finite_target_calls':len(finite),
            'occupied_states_recorded_including_rejections':len(occupied),
            'qualified_scientific_measurement':False}
    if len(occupied):
        _,counts=np.unique(occupied,return_counts=True)
        moves=sum(bool(r.get('accepted')) for r in rows)-1
        result.update(accepted_transitions_after_initialization=int(moves),
                      transition_acceptance_fraction=float(moves/(len(occupied)-1)) if len(occupied)>1 else None,
                      longest_hold=int(counts.max()),largest_occupied_state_fraction=float(counts.max()/len(occupied)))
        import arviz as az
        ess={}
        for j,name in enumerate(request['names']):
            values=request['points'][occupied,j]
            if len(values)<20 or np.ptp(values)==0:
                ess[str(name)]={'status':'insufficient_or_constant','bulk_ess':None,'identity_ess':None};continue
            a=float(az.ess(values,method='bulk'));b=float(az.ess(values,method='identity'))
            ess[str(name)]={'status':'short_single_chain_diagnostic','bulk_ess':a if np.isfinite(a) else None,
                            'identity_ess':b if np.isfinite(b) else None}
        result['fixed_coordinate_ESS_diagnostic']=ess
    if finite:
        lw=np.array([r['target_logpost']-r['logq'] for r in finite]);w=np.exp(lw-logsumexp(lw))
        result.update(candidate_raw_importance_ESS=float(1/(w@w)),
                      candidate_largest_normalized_weight=float(w.max()),
                      importance_scope='Independent q candidates from the finite recorded request subset, not posterior qualification.')
    return result


def pilot(arguments,proposal_folder,output,cpu):
    s.runtime();s.validation_guard()
    checked=json.loads(VALIDATION.read_text());bound=dependencies()
    assert checked['status']=='passed_synthetic_pilot_density_and_retention_checks'
    assert checked['dependency_sha256']==bound
    design=json.loads(DESIGN.read_text())
    info=s.configuration(arguments);target=s.identify(arguments)
    mixture,proposal=s.load_proposal(proposal_folder,info)
    assert proposal.get('proposal_kind')=='anchored_bridge_weighted_v1','Require genuinely audited anchored bridge training.'
    output=Path(output).resolve();assert output.is_relative_to(ROOT/'.work') and not output.exists()
    available=os.sched_getaffinity(0);assert cpu in available
    os.sched_setaffinity(0,{cpu});assert os.sched_getaffinity(0)=={cpu}
    output.mkdir(parents=True)
    request=requests(mixture,arguments['model'],design)
    path=output/'requests.npz';np.savez_compressed(path,**request)
    s.dump_new(output/'manifest.json',{'schema':'anchored-proposal-pilot-v1','arguments':arguments,
        'target_identity':target,'proposal':proposal,'dependency_sha256':bound,
        'pilot_validation_sha256':s.digest(VALIDATION),'requests_sha256':s.digest(path),
        'CPU_affinity':[cpu],'scope':design['scope']})
    from independence_runtime import native_initialize
    warmup=native_initialize(info)
    s.dump_new(output/'native-initialization.json',warmup)
    from cobaya.model import get_model
    rows=[];current=None;current_logp=None;current_logq=None
    status='completed_declared_request_budget'
    start=time.monotonic()
    with get_model(info) as model:
        assert list(model.parameterization.sampled_params())==mixture.names
        theory=model.theory['spectral_surrogate']
        for i,x in enumerate(request['points']):
            s.verify_hashes(bound)
            point=dict(zip(mixture.names,map(float,x)));call_start=time.monotonic()
            try:
                result=model.logposterior(point);logp=float(result.logpost)
                closure=float(logp-float(result.logprior)-sum(result.loglikes)) if np.isfinite(logp) else None
                if closure is not None:assert abs(closure)<1e-9
                accept,alpha=transition(logp,current_logp,current_logq,float(request['logq'][i]),float(request['logu'][i]))
                if accept:current=i;current_logp=logp;current_logq=float(request['logq'][i])
                row={'index':i,'point':point,'component':str(request['components'][i]),
                     'logq':float(request['logq'][i]),'log_uniform':float(request['logu'][i]),
                     'target_logpost':logp if np.isfinite(logp) else None,'finite_target':bool(np.isfinite(logp)),
                     'loglikes':{k:float(v) if np.isfinite(v) else None for k,v in zip(model.likelihood,result.loglikes)},
                     'logprior':float(result.logprior) if np.isfinite(result.logprior) else None,
                     'density_accounting_closure':closure,
                     'accepted':accept,'log_acceptance':alpha,'occupied_candidate':current,
                     'successful_native_fallbacks':int(theory.exact_calls),'surrogate_calls':int(theory.surrogate_calls)}
            except Exception as error:
                row={'index':i,'point':point,'status':'failed_evaluation','error':repr(error),
                     'occupied_candidate':None};status='failed_evaluation'
            row['seconds']=time.monotonic()-call_start;rows.append(row)
            with (output/'records.jsonl').open('a') as stream:stream.write(json.dumps(row,allow_nan=False)+'\n')
            if status=='failed_evaluation':break
            if theory.exact_calls>=design['maximum_successful_native_fallbacks']:
                status='stopped_native_fallback_budget';break
            if time.monotonic()-start>=design['maximum_evaluation_wall_seconds']:
                status='stopped_wall_budget_after_completed_call';break
            if (i+1)%25==0:print(json.dumps({'target_calls':i+1,'finite':sum(r.get('finite_target',False) for r in rows),'seconds':time.monotonic()-start}),flush=True)
        fallback_count=int(theory.exact_calls);surrogate_count=int(theory.surrogate_calls)
    assert dependencies()==bound and s.identify(arguments)==target
    _,after=s.load_proposal(proposal_folder,info);assert after==proposal
    result=metrics(rows,request)
    result.update(status=status,seconds=time.monotonic()-start,target_identity=target['identity'],
        successful_native_fallbacks=fallback_count,surrogate_calls=surrogate_count,
        explicit_native_initialization_calls=1,
        native_count_limits='Fallback count records successful native spectra; failed internal native attempts are not instrumented. The fixed initialization is separate.',
        dependency_sha256=bound,input_sha256={s.relative(p):s.digest(p) for p in output.iterdir() if p.is_file()},
        limits=design['limits'])
    s.dump_new(output/'result.json',result)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model',choices=['lcdm','cpl'],required=True);p.add_argument('--surrogate',type=Path,required=True)
    p.add_argument('--native-accuracy',type=int,choices=[1,2],default=1)
    p.add_argument('--proposal-folder',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--cpu',type=int,required=True)
    a=p.parse_args();result=pilot(s.settings(a.model,a.surrogate,a.native_accuracy),a.proposal_folder,a.output,a.cpu)
    print(json.dumps({k:result[k] for k in ['status','target_calls_recorded','finite_target_calls','successful_native_fallbacks']}))

if __name__=='__main__':main()
