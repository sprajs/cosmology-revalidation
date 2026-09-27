"""Native correction for a separately identified anchored proposal cohort."""
import argparse
import json
import multiprocessing as mp
import os
from pathlib import Path
import re
import time
import numpy as np

import anchored_sampling as sampling
from diagnostics import check_mpi
from exact_correction import record_digest,verify_record,summarize

HERE=sampling.HERE;ROOT=sampling.ROOT
_exact=_proposal=None
_bound=None
_initialization=None


def dependencies():
    names=['anchored_correction.py','anchored_sampling.py','anchored-sampling-design.json',
           'exact_correction.py','luminosity_sensitivity.py','luminosity-sensitivity-design.json',
           'diagnostics.py','independence_runtime.py']
    return {sampling.relative(HERE/name):sampling.digest(HERE/name) for name in names}


def initialize(manifest,bound,destination,identity):
    global _exact,_proposal,_bound,_initialization
    sampling.runtime();sampling.verify_hashes(bound)
    assert dependencies()==bound
    sampling.verify_manifest(manifest)
    info=sampling.configuration(manifest['arguments'])
    from independence_runtime import native_initialize
    initial=native_initialize(info)
    path=Path(destination)/f'native-initialization-{os.getpid()}.json'
    sampling.dump_new(path,{'correction_identity':identity,'source_sha256':bound,
                           'explicit_native_initialization_calls':1,'initialization':initial})
    _initialization={sampling.relative(path):sampling.digest(path)}
    from cobaya.model import get_model
    _exact=get_model(sampling.configuration(manifest['arguments'],native=True))
    _proposal=get_model(info)
    sampling.verify_hashes(bound);_bound=bound


def read_cached(task):
    index,point,destination,identity,stored_logpost=task
    path=Path(destination)/f'{index:05d}.json'
    if not path.exists():return None
    record=json.loads(path.read_text());verify_record(record)
    assert record['schema']=='dedicated-anchored-native-point-v1'
    assert record['index']==index and record['point']==point and record['target_identity']==identity
    assert record['stored_chain_logpost']==stored_logpost
    sampling.verify_hashes(record['initialization_record_sha256'])
    return record


def evaluate(task):
    index,point,destination,identity,stored_logpost=task
    sampling.verify_hashes(_bound)
    old=read_cached(task)
    if old is not None:return old
    path=Path(destination)/f'{index:05d}.json' 
    started=time.monotonic()
    try:
        exact=_exact.logposterior(point);proposal=_proposal.logposterior(point)
        record={'index':index,'point':point,'status':'finite',
            'exact_logpost':float(exact.logpost),'proposal_logpost':float(proposal.logpost),
            'log_weight':float(exact.logpost-proposal.logpost),
            'exact_loglikes':dict(zip(_exact.likelihood,map(float,exact.loglikes))),
            'proposal_loglikes':dict(zip(_proposal.likelihood,map(float,proposal.loglikes))),
            'derived':dict(zip(_exact.parameterization.derived_params(),map(float,exact.derived))),
            'proposal_chain_logpost_difference':float(proposal.logpost-stored_logpost)}
        if not np.isfinite(record['log_weight']):record['status']='nonfinite'
        if not all(np.isfinite(v) for v in record['derived'].values()):record['status']='nonfinite_derived'
        if abs(record['proposal_chain_logpost_difference'])>.01:record['status']='proposal_density_mismatch'
        ep=record['exact_logpost']-sum(record['exact_loglikes'].values())
        pp=record['proposal_logpost']-sum(record['proposal_loglikes'].values())
        record['native_proposal_logprior_difference']=float(ep-pp)
        if not np.isfinite(ep-pp) or abs(ep-pp)>1e-9:record['status']='prior_density_mismatch'
    except Exception as error:
        record={'index':index,'point':point,'status':'exception','exception':repr(error)}
    sampling.verify_hashes(_bound)
    record.update(schema='dedicated-anchored-native-point-v1',seconds=time.monotonic()-started,
                  target_identity=identity,stored_chain_logpost=stored_logpost,initialization_record_sha256=_initialization)
    record['payload_sha256']=record_digest(record)
    temporary=path.with_suffix('.part')
    temporary.write_text(json.dumps(record,indent=2)+'\n');temporary.replace(path)
    return record


def completion_evidence(folder,manifest):
    """Reject an active, finite-budget or incompletely stopped sampler cohort."""
    from cobaya.yaml import yaml_load_file
    checkpoint=folder/'chain.checkpoint'
    saved=yaml_load_file(str(checkpoint))['sampler']
    assert list(saved)==['independence_sampler.IndependenceMCMC']
    state=next(iter(saved.values()))
    assert state['converged'] is True and state['mpi_size']==4 and state['burn_in']==0
    convergence=folder/'chain_independence-convergence.jsonl'
    with convergence.open() as stream:
        last=None
        for line in stream:
            if line.strip():last=json.loads(line)
    assert last is not None and last['converged'] is True and last['original_convergence_passed'] is True
    assert last['proposal_sha256']==manifest['proposal']['proposal_sha256']
    for key in ['original_multivariate_previous','original_multivariate_current']:
        assert last[key] is not None and np.isfinite(last[key]) and last[key]<.005
    assert last['original_bound_statistic'] is not None and np.isfinite(last['original_bound_statistic']) and last['original_bound_statistic']<.05
    assert state['Rminus1_last']==last['original_multivariate_current']
    assert last['rank_gate']['passed'] is True and not last['rank_gate']['failed']
    assert last['rank_gate']['discard_fraction_after_sampler_burnin']==.3
    paths=[checkpoint,convergence]
    for rank in range(4):
        path=folder/f'chain_independence-terminal-{rank}.json'
        terminal=json.loads(path.read_text())
        ledger=folder/f'chain_independence-candidates-{rank}.jsonl'
        assert terminal['converged'] is True and terminal['terminal_weight']==1
        assert terminal['terminal_state_in_collection'] is False
        assert terminal['ledger_sha256']==sampling.digest(ledger)
        paths.extend([path,ledger])
    return {'status':'all_original_sampler_and_terminal_gates_passed','final_convergence':last,
            'input_sha256':{sampling.relative(p):sampling.digest(p) for p in paths}}


def qualified_chain_inputs(folder,diagnostics_path):
    folder=Path(folder).resolve();diagnostics_path=Path(diagnostics_path).resolve()
    manifests=[json.loads((folder/f'run-{rank}.json').read_text()) for rank in range(4)]
    first=manifests[0]
    for rank,manifest in enumerate(manifests):
        assert manifest['schema']=='dedicated-anchored-chain-v1' and manifest['rank']==rank and manifest['MPI_size']==4
        assert manifest['arguments']==first['arguments'] and manifest['target_identity']==first['target_identity']
        assert manifest['sampler_source_sha256']==first['sampler_source_sha256']
    sampling.verify_manifest(first)
    completed=completion_evidence(folder,first)
    check=json.loads(diagnostics_path.read_text())
    assert check['status']=='passed' and not check['failed_gates']
    assert check['discard_fraction_after_sampler_burnin']==.3
    assert check['independent_chains']==4
    assert check['code_sha256']==sampling.digest(HERE/'diagnostics.py')
    sampling.verify_hashes(check['input_sha256'])
    assert check_mpi(folder,.3)==check,'Recomputed parent diagnostics disagree.'
    return first,check,completed


def select_points(paths,names,check,points):
    """Same stratification and RLE expansion as the original correction."""
    assert len(paths)==4 and points>=2000 and points%4==0
    rng=np.random.default_rng(json.loads(sampling.DESIGN.read_text())['selection_seed'])
    selected=[];groups=[];locations=[];stored=[]
    for group,path in enumerate(paths):
        columns=path.open().readline().lstrip('#').split();raw=np.loadtxt(path,ndmin=2)
        assert np.allclose(raw[:,0],np.round(raw[:,0])) and np.all(raw[:,0]>0)
        origin=np.repeat(np.arange(len(raw)),raw[:,0].astype(int))
        common=check['equal_chain_length_for_diagnostics']
        assert len(origin)-int(.3*len(origin))>=common
        discarded=len(origin)-common;origin=origin[discarded:];expanded=raw[origin]
        count=points//4;assert len(expanded)>=count
        indices=np.floor((np.arange(count)+rng.random(count))*len(expanded)/count).astype(int)
        selected.extend([{name:float(expanded[i,columns.index(name)]) for name in names} for i in indices])
        groups.extend([group]*count)
        locations.extend([{'chain':path.name,'row':int(origin[i]),'expanded_index':int(discarded+i)} for i in indices])
        stored.extend([-float(expanded[i,columns.index('minuslogpost')]) for i in indices])
    return dict(points=selected,groups=groups,locations=locations,stored_logposts=stored)


def run(folder,diagnostics_path,output,points=2000,workers=4,name='anchored-exact-correction'):
    sampling.runtime();sampling.validation_guard()
    assert re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',name)
    assert 1<=workers<=4
    folder=Path(folder).resolve();output=Path(output).resolve();diagnostics_path=Path(diagnostics_path).resolve()
    manifest,check,completed=qualified_chain_inputs(folder,diagnostics_path)
    bound=dependencies();target=manifest['target_identity']
    correction_id=sampling.identity({'target_identity':target['identity'],'correction_dependencies':bound})
    paths=sorted(folder.glob('chain.[0-9]*.txt'))
    assert len(paths)==4
    assert set(check['input_sha256'])=={sampling.relative(p) for p in paths}
    info=sampling.configuration(manifest['arguments'])
    names=[k for k,v in info['params'].items() if isinstance(v,dict) and 'prior'in v]
    selection=select_points(paths,names,check,points)
    selection.update(schema='dedicated-anchored-native-selection-v1',settings=manifest['arguments'],
        correction_identity=correction_id,correction_dependency_sha256=bound,
        proposal_target_identity=target['identity'],native_accuracy=target['native_accuracy'],
        sampler_completion=completed,chain_sha256=check['input_sha256'],diagnostics_path=sampling.relative(diagnostics_path),
        diagnostics_sha256=sampling.digest(diagnostics_path),
        rank_manifest_sha256={sampling.relative(folder/f'run-{rank}.json'):sampling.digest(folder/f'run-{rank}.json') for rank in range(4)})
    destination=folder/name;design=destination/'selection.json'
    sampling.dump_new(design,selection)
    if output.exists():
        previous=json.loads(output.read_text())
        assert previous['selection_sha256']==sampling.digest(design)
        sampling.verify_hashes(previous['native_record_sha256'])
    tasks=[(i,p,str(destination),correction_id,selection['stored_logposts'][i]) for i,p in enumerate(selection['points'])]
    records=[read_cached(task) for task in tasks]
    pending=[task for task,record in zip(tasks,records) if record is None]
    if pending:
        with mp.get_context('spawn').Pool(min(workers,len(pending)),initializer=initialize,
                initargs=(manifest,bound,str(destination),correction_id)) as pool:
            for row in pool.imap(evaluate,pending,chunksize=1):
                records[row['index']]=row
                completed_count=sum(r is not None for r in records)
                if completed_count%50==0:print(json.dumps({'anchored_native_points':completed_count,'total':points}),flush=True)
    result=summarize(records,np.array(selection['groups']))
    assert dependencies()==bound
    sampling.verify_manifest(manifest)
    sampling.verify_hashes(selection['chain_sha256']);sampling.verify_hashes(selection['rank_manifest_sha256'])
    sampling.verify_hashes(completed['input_sha256'])
    assert sampling.digest(diagnostics_path)==selection['diagnostics_sha256']
    result.update(schema='dedicated-anchored-native-summary-v1',selection_path=sampling.relative(design),
        selection_sha256=sampling.digest(design),code_sha256=sampling.digest(__file__),
        correction_dependency_sha256=bound,correction_identity=correction_id,
        proposal_target_identity=target['identity'],native_accuracy=target['native_accuracy'],
        diagnostics_path=selection['diagnostics_path'],diagnostics_sha256=selection['diagnostics_sha256'],
        initialization_record_sha256={path:hashvalue for row in records for path,hashvalue in row['initialization_record_sha256'].items()},
        native_record_sha256={sampling.relative(destination/f'{i:05d}.json'):sampling.digest(destination/f'{i:05d}.json') for i in range(points)})
    sampling.dump_new(output,result)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder',type=Path);p.add_argument('--diagnostics',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--points',type=int,default=2000)
    p.add_argument('--workers',type=int,default=4);p.add_argument('--name',default='anchored-exact-correction')
    a=p.parse_args();r=run(a.folder,a.diagnostics,a.output,a.points,a.workers,a.name)
    print(json.dumps({'status':r['status'],'native_accuracy':r['native_accuracy']}))

if __name__=='__main__':main()
