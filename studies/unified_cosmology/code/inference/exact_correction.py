"""Native-CAMB correction of a converged numerical-proposal posterior.

Every selected point, failure and exact/proposal likelihood contribution is kept.
Pareto smoothing is a diagnostic only: reported weights are never clipped.
"""
import argparse
import hashlib
import json
import multiprocessing as mp
from pathlib import Path
import time
import numpy as np
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
_exact = _proposal = None


def initialize(settings):
    global _exact,_proposal
    from cobaya.model import get_model
    if settings.get('fast_lensing'):
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
        assert old['point']==point and old['target_identity']==identity
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
    import arviz as az
    if any(row['status']!='finite' for row in records):
        return {'status':'failed_exact_evaluation','failures':[r for r in records if r['status']!='finite']}
    lw = np.array([r['log_weight'] for r in records]);weights = np.exp(lw-logsumexp(lw))
    _,k = az.psislw(lw)
    ess = float(1/np.sum(weights**2))
    names = sorted(set(records[0]['point'])|set(records[0]['derived']))
    posterior = {};per_chain = {}
    for name in names:
        values = np.array([dict(r['point'],**r['derived'])[name] for r in records])
        posterior[name] = weighted_summary(values,weights)
        means = []
        for group in np.unique(groups):
            keep = groups==group;cw = weights[keep]/weights[keep].sum()
            means.append(float(np.sum(values[keep]*cw)))
        per_chain[name] = means
        # Sensitivity to batch length is visible. This is not guaranteed
        # conservative for long autocorrelation or rare, concentrated weights.
        batch_checks={}
        for batches in [10,20,40]:
            numerators=[];denominators=[];sizes=[]
            for group in np.unique(groups):
                indices=np.flatnonzero(groups==group)
                for block in np.array_split(indices,batches):
                    numerators.append(np.sum(weights[block]*(values[block]-posterior[name]['mean'])))
                    denominators.append(np.sum(weights[block]));sizes.append(len(block))
            n=len(numerators)
            mcse=float(np.sqrt(np.var(numerators,ddof=1)/n)/np.mean(denominators))
            batch_checks[str(batches)]={'mean_mcse':mcse,'minimum_points_per_batch':min(sizes),
                'minimum_batch_weight':float(min(denominators)),'maximum_batch_weight':float(max(denominators))}
        posterior[name]['weighted_batch_mean_mcse']=batch_checks['20']['mean_mcse']
        posterior[name]['batch_sensitivity']=batch_checks
    status = 'passed_importance_weight_gates' if ess>=400 and float(k)<.7 else 'insufficient_importance_overlap'
    return {'status':status,'exact_points':len(records),'raw_weight_ess':ess,'pareto_k':float(k),
            'largest_normalized_weight':float(weights.max()),
            'log_weight_quantiles':np.quantile(lw,[0,.025,.5,.975,1]).tolist(),
            'posterior':posterior,'independent_chain_weighted_means':per_chain,
            'qualification':'Must also pass parent-chain convergence and inspect weighted Monte Carlo stability. Finite diagnostics do not prove global mode coverage.'}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('folder',type=Path)
    p.add_argument('--diagnostics',type=Path,required=True)
    p.add_argument('--points',type=int,default=2000)
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();folder=a.folder.resolve()
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
    if settings.get('fast_lensing'):
        from modern_fast import configuration, identify
    else:
        from modern_run import configuration
        from target_identity import identify
    proposal_info=configuration(**{k:settings[k] for k in ['model','evolution','sample','calibration']},surrogate=Path(settings['surrogate']))
    identity=identify(proposal_info,sample_path(settings['sample']),Path(settings['surrogate']))
    assert identity==manifest['target_identity'],'Sampled proposal source, data, prior or spectrum identity changed.'
    for rank in range(4):
        assert json.loads((folder/f'run-{rank}.json').read_text())['target_identity']==identity
    correction_identity=hashlib.sha256((identity['identity']+hashlib.sha256(Path(__file__).read_bytes()).hexdigest()).encode()).hexdigest()
    out=folder/'exact-correction';out.mkdir(exist_ok=True)
    design=out/'selection.json'
    if design.exists():
        selection=json.loads(design.read_text());assert len(selection['points'])==a.points
        assert selection['correction_identity']==correction_identity
        assert selection['diagnostics_sha256']==hashlib.sha256(a.diagnostics.read_bytes()).hexdigest()
        assert selection['surrogate_sha256']==hashlib.sha256(Path(settings['surrogate']).read_bytes()).hexdigest()
        for path,digest in selection['chain_sha256'].items():
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'Parent chains changed after point selection.'
    else:
        from cobaya.yaml import yaml_load_file
        config=yaml_load_file(str(folder/'chain.updated.yaml'))
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
            'chain_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            'diagnostics_sha256':hashlib.sha256(a.diagnostics.read_bytes()).hexdigest(),
            'surrogate_sha256':hashlib.sha256(Path(settings['surrogate']).read_bytes()).hexdigest()}
        design.write_text(json.dumps(selection,indent=2)+'\n')
    tasks=[(i,point,str(out),correction_identity,selection['stored_logposts'][i]) for i,point in enumerate(selection['points'])]
    with mp.get_context('spawn').Pool(a.workers,initializer=initialize,initargs=(settings,)) as pool:
        records=[]
        for row in pool.imap(evaluate,tasks,chunksize=1):
            records.append(row)
            if len(records)%50==0:print(json.dumps({'completed':len(records),'total':len(tasks)}),flush=True)
    result=summarize(records,np.array(selection['groups']))
    result['selection_sha256']=hashlib.sha256(design.read_bytes()).hexdigest()
    result['code_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['posterior','independent_chain_weighted_means']},indent=2))


if __name__=='__main__':main()
