"""Native-CAMB correction of a converged numerical-proposal posterior.

Every selected point, failure and exact/proposal likelihood contribution is kept.
Pareto smoothing is a diagnostic only: reported weights are never clipped.
"""
import argparse
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import re
import time
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
_exact = _proposal = None


def record_digest(record):
    """Seal numerical payloads, including retained failed/nonfinite outcomes."""
    payload = {key: value for key, value in record.items() if key != 'payload_sha256'}
    return hashlib.sha256(json.dumps(payload, sort_keys=True,
                                    separators=(',', ':')).encode()).hexdigest()


def verify_record(record):
    assert record.get('payload_sha256') == record_digest(record), 'Native cached payload changed.'


def correction_dependencies():
    """Bind evaluation and every preregistered weighted-stability decision."""
    paths = [Path(__file__),HERE/'luminosity_sensitivity.py',
             HERE/'luminosity-sensitivity-design.json']
    return {str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths}


def initialize(settings):
    global _exact,_proposal
    from cobaya.model import get_model
    if settings.get('gpu'):
        from modern_gpu import configuration
    elif settings.get('fast_lensing'):
        from modern_fast import configuration
    else:
        from modern_run import configuration
    arguments = {k:settings[k] for k in ['model','evolution','sample','calibration']}
    _exact = get_model(configuration(**arguments))
    _proposal = get_model(configuration(**arguments,surrogate=Path(settings['surrogate'])))


def evaluate(task):
    index,point,destination,identity,stored_logpost = task
    path = Path(destination)/f'{index:05d}.json'
    if path.exists():
        old = json.loads(path.read_text())
        verify_record(old)
        assert old['index']==index and old['point']==point and old['target_identity']==identity
        return old
    started = time.monotonic()
    try:
        exact = _exact.logposterior(point)
        approx = _proposal.logposterior(point)
        record = {'index':index,'point':point,'status':'finite',
            'exact_logpost':float(exact.logpost),'proposal_logpost':float(approx.logpost),
            'log_weight':float(exact.logpost-approx.logpost),
            'exact_loglikes':dict(zip(_exact.likelihood,map(float,exact.loglikes))),
            'proposal_loglikes':dict(zip(_proposal.likelihood,map(float,approx.loglikes))),
            'derived':dict(zip(_exact.parameterization.derived_params(),map(float,exact.derived)))}
        if not np.isfinite(record['log_weight']):record['status']='nonfinite'
        if not all(np.isfinite(v) for v in record['derived'].values()):record['status']='nonfinite_derived'
        record['proposal_chain_logpost_difference']=float(approx.logpost-stored_logpost)
        if abs(record['proposal_chain_logpost_difference'])>.01:record['status']='proposal_density_mismatch'
    except Exception as error:
        record = {'index':index,'point':point,'status':'exception','exception':repr(error)}
    record['seconds'] = time.monotonic()-started
    record['target_identity']=identity
    record['payload_sha256']=record_digest(record)
    temporary = path.with_suffix('.part')
    temporary.write_text(json.dumps(record,indent=2)+'\n');temporary.replace(path)
    return record


def weighted_summary(values,weights):
    order = np.argsort(values);x = np.asarray(values)[order];w = weights[order]
    positions = np.cumsum(w)-w/2
    mean = float(np.sum(x*w));variance = float(np.sum(w*(x-mean)**2))
    return {'mean':mean,'sd':float(np.sqrt(variance)),
            'quantiles_025_16_50_84_975':np.interp([.025,.16,.5,.84,.975],positions,x).tolist()}


def summarize(records,groups):
    from luminosity_sensitivity import weight_diagnostics, DESIGN_PATH
    if any(row['status']!='finite' for row in records):
        # Individual records retain the native outcome. A failed summary still
        # needs valid JSON: nonfinite numerical values are explicitly null.
        def finite_json(value):
            if isinstance(value,dict):return {k:finite_json(v) for k,v in value.items()}
            if isinstance(value,list):return [finite_json(v) for v in value]
            if isinstance(value,(float,np.floating)) and not np.isfinite(value):return None
            return value
        return {'status':'failed_exact_evaluation',
            'failures':finite_json([r for r in records if r['status']!='finite']),
            'nonfinite_failure_numbers':'null; individual evaluation records preserve the original outcome'}
    gates=json.loads(DESIGN_PATH.read_text())['overlap_gates']
    groups=np.asarray(groups)
    assert len(groups)==len(records)
    lw=np.array([r['log_weight'] for r in records])
    names = sorted(set(records[0]['point'])|set(records[0]['derived']))
    values={name:np.array([dict(r['point'],**r['derived'])[name] for r in records])
            for name in names}
    if not np.isfinite(lw).all() or not all(np.isfinite(x).all() for x in values.values()):
        return {'status':'failed_exact_evaluation','failures':['Nonfinite correction weight or parameter value.']}
    check=weight_diagnostics(lw,values,groups,gates)
    failed=list(check['failed_gates'])
    if len(records)<gates['minimum_exact_points']:
        failed.append('minimum_exact_points')
    posterior = {};per_chain = {}
    for name in names:
        posterior[name]=dict(check['weighted_summaries_for_diagnostics'][name])
        stability=check['weighted_chain_stability'][name]
        per_chain[name]=[stability['independent_chain_weighted_means'][str(group)]
                         for group in np.unique(groups)]
        # Retain legacy field names, with exactly the shared tested arithmetic.
        batch_checks={key:{'mean_mcse':value['mean_mcse'],
            'minimum_points_per_batch':value['minimum_points_per_batch'],
            'minimum_batch_weight':value['minimum_weight_per_batch'],
            'maximum_batch_weight':value['maximum_weight_per_batch']}
            for key,value in stability['batch_sensitivity'].items()}
        posterior[name]['weighted_batch_mean_mcse']=batch_checks['20']['mean_mcse']
        posterior[name]['batch_sensitivity']=batch_checks
    status = 'passed_importance_weight_gates' if not failed else 'insufficient_importance_overlap'
    return {'status':status,'exact_points':len(records),
            'raw_weight_ess':check['raw_weight_ESS'],'pareto_k':check['Pareto_k'],
            'pareto_k_status':check['Pareto_k_status'],
            'largest_normalized_weight':check['largest_normalized_weight'],
            'log_weight_quantiles':check['log_weight_quantiles'],
            'posterior':posterior,'independent_chain_weighted_means':per_chain,
            'per_chain_weight_diagnostics':check['per_chain'],
            'weighted_chain_stability':check['weighted_chain_stability'],
            'weighted_stability_gates':gates,'failed_gates':failed,
            'qualification':'All registered weight, chain-agreement and batch-stability gates are required. Parent-chain convergence is separately enforced at selection. Weight ESS is not autocorrelation-adjusted ESS; finite diagnostics do not prove global mode coverage.'}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('folder',type=Path)
    p.add_argument('--diagnostics',type=Path,required=True)
    p.add_argument('--points',type=int,default=2000)
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--name',default='exact-correction',
                   help='Separate preserved correction folder for a new point selection.')
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();folder=a.folder.resolve()
    assert re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',a.name),'Use a simple distinct correction name.'
    diagnostics_path=str(a.diagnostics.resolve().relative_to(ROOT))
    diagnostics_sha256=hashlib.sha256(a.diagnostics.read_bytes()).hexdigest()
    check=json.loads(a.diagnostics.read_text());assert check['status']=='passed'
    actual_paths=sorted(folder.glob('chain.[0-9]*.txt'))
    assert len(actual_paths)==4
    assert set(check['input_sha256'])=={str(p.relative_to(ROOT)) for p in actual_paths},'Diagnostics and selected chain cohort differ.'
    assert a.points>=2000 and a.points%4==0
    manifest=json.loads((folder/'run-0.json').read_text());settings=manifest['arguments']
    assert settings['surrogate']
    assert check['discard_fraction_after_sampler_burnin']==.3
    for path,digest in check['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'Diagnostics refer to different chain bytes.'
    from late_geometry import sample_path
    if settings.get('gpu'):
        from modern_gpu import configuration, identify
    elif settings.get('fast_lensing'):
        from modern_fast import configuration, identify
    else:
        from modern_run import configuration
        from target_identity import identify
    proposal_info=configuration(**{k:settings[k] for k in ['model','evolution','sample','calibration']},surrogate=Path(settings['surrogate']))
    identity=identify(proposal_info,sample_path(settings['sample']),Path(settings['surrogate']))
    assert identity==manifest['target_identity'],'Sampled proposal source, data, prior or spectrum identity changed.'
    for rank in range(4):
        assert json.loads((folder/f'run-{rank}.json').read_text())['target_identity']==identity
    dependencies=correction_dependencies()
    correction_identity=hashlib.sha256(json.dumps(
        {'target_identity':identity['identity'],'correction_dependencies':dependencies},
        sort_keys=True,separators=(',',':')).encode()).hexdigest()
    out=folder/a.name;out.mkdir(exist_ok=True)
    design=out/'selection.json'
    if a.output.exists():
        previous = json.loads(a.output.read_text())
        assert previous['selection_path'] == str(design.relative_to(ROOT))
        for name, expected in previous['native_record_sha256'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == expected, 'Previously summarized native file changed.'
    if design.exists():
        selection=json.loads(design.read_text());assert len(selection['points'])==a.points
        assert selection['correction_identity']==correction_identity
        assert selection['correction_dependency_sha256']==dependencies
        assert selection['proposal_target_identity']==identity['identity']
        assert selection['diagnostics_path']==diagnostics_path
        assert selection['diagnostics_sha256']==diagnostics_sha256
        assert selection['surrogate_sha256']==hashlib.sha256(Path(settings['surrogate']).read_bytes()).hexdigest()
        for path,digest in selection['chain_sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'Parent chains changed after point selection.'
    else:
        # Already rebuilt and verified against every immutable rank manifest.
        # Do not infer parameter names from concurrently written shared YAML.
        config=proposal_info
        names=[k for k,v in config['params'].items() if isinstance(v,dict) and 'prior' in v]
        paths=actual_paths
        selected=[];groups=[];locations=[];stored_logposts=[]
        rng=np.random.default_rng(272639)
        for group,path in enumerate(paths):
            columns=path.open().readline().lstrip('#').split();raw=np.loadtxt(path,ndmin=2)
            assert np.allclose(raw[:,0],np.round(raw[:,0]))
            origin=np.repeat(np.arange(len(raw)),raw[:,0].astype(int))
            common=check['equal_chain_length_for_diagnostics']
            assert len(origin)-int(.3*len(origin))>=common
            discarded=len(origin)-common;origin=origin[discarded:]
            expanded=raw[origin]
            n=a.points//4;assert len(expanded)>=n
            indices=np.floor((np.arange(n)+rng.random(n))*len(expanded)/n).astype(int)
            selected.extend([{name:float(expanded[i,columns.index(name)]) for name in names} for i in indices])
            groups.extend([group]*n)
            locations.extend([{'chain':path.name,'row':int(origin[i]),'expanded_index':int(discarded+i)} for i in indices])
            stored_logposts.extend([-float(expanded[i,columns.index('minuslogpost')]) for i in indices])
        selection={'points':selected,'groups':groups,'settings':settings,
            'locations':locations,'stored_logposts':stored_logposts,'correction_identity':correction_identity,
            'correction_dependency_sha256':dependencies,
            'proposal_target_identity':identity['identity'],
            'chain_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            'diagnostics_path':diagnostics_path,'diagnostics_sha256':diagnostics_sha256,
            'surrogate_sha256':hashlib.sha256(Path(settings['surrogate']).read_bytes()).hexdigest()}
        design.write_text(json.dumps(selection,indent=2)+'\n')
    tasks=[(i,point,str(out),correction_identity,selection['stored_logposts'][i]) for i,point in enumerate(selection['points'])]
    with mp.get_context('spawn').Pool(a.workers,initializer=initialize,initargs=(settings,)) as pool:
        records=[]
        for row in pool.imap(evaluate,tasks,chunksize=1):
            records.append(row)
            if len(records)%50==0:print(json.dumps({'completed':len(records),'total':len(tasks)}),flush=True)
    result=summarize(records,np.array(selection['groups']))
    assert correction_dependencies()==dependencies,'Correction or gate implementation changed during evaluation.'
    assert hashlib.sha256(a.diagnostics.read_bytes()).hexdigest()==diagnostics_sha256,'Parent diagnostics changed during evaluation.'
    result['selection_path']=str(design.relative_to(ROOT))
    result['selection_sha256']=hashlib.sha256(design.read_bytes()).hexdigest()
    result['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['correction_dependency_sha256']=dependencies
    result['correction_identity']=correction_identity
    result['proposal_target_identity']=identity['identity']
    result['diagnostics_path']=diagnostics_path
    result['diagnostics_sha256']=diagnostics_sha256
    result['native_record_sha256']={str((out/f'{i:05d}.json').relative_to(ROOT)):
        hashlib.sha256((out/f'{i:05d}.json').read_bytes()).hexdigest() for i in range(len(records))}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['posterior','independent_chain_weighted_means']},indent=2))


if __name__=='__main__':main()
