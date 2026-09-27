"""Shared strict runtime for separate native2 pilot and full-cohort execution."""
import copy
from contextlib import ExitStack
import importlib.metadata
import json
import os
from pathlib import Path
import time
from unittest.mock import patch

import numpy as np

import native_accuracy_correction as contracts
import native_posterior_precision as original
import native_precision_review_consumer as consumer
from target_identity import canonical

ROOT=contracts.ROOT
HERE=Path(__file__).resolve().parent
DESIGN=HERE/'native-accuracy-execution-design-v2.json'
VALIDATION=ROOT/'studies/unified_cosmology/results/inference/native-accuracy-execution-validation.json'
GRID=np.array([0.,.5,1.,2.33,1100.])


def relative(path):return str(Path(path).resolve().relative_to(ROOT))
def digest(path):return contracts.digest(path)
def identity(value):return contracts.identity({k:v for k,v in value.items() if k!='payload_sha256'})


def sources():
    out=contracts.sources()
    for name in ['native_accuracy_runtime.py','native_accuracy_pilot.py','native_accuracy_execute.py',
                 'native_accuracy_qualify.py','native_accuracy_execution_validate.py',DESIGN.name]:
        p=HERE/name;out[relative(p)]=digest(p)
    return out


def merge(*items):return consumer.merge(*items)
def verify_hashes(value):consumer.hashes(value)


def write_new(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    sealed=dict(value,payload_sha256=identity(value))
    part=path.with_name(path.name+'.part')
    with part.open('x') as stream:
        json.dump(sealed,stream,indent=2,allow_nan=False);stream.write('\n');stream.flush();os.fsync(stream.fileno())
    os.link(part,path);part.unlink() # Atomic exclusive publication, no overwrite.


def read(path):
    value=json.loads(Path(path).read_text());consumer.finite_tree(value)
    assert value['payload_sha256']==identity(value),'Changed sealed native2 payload.'
    return value


def plain(value):return {k:v for k,v in value.items() if k!='payload_sha256'}


def validation_guard():
    value=json.loads(VALIDATION.read_text())
    assert value['status']=='passed_synthetic_native_accuracy_execution_validation'
    assert value['source_sha256']==sources() and value['physical_calls']==0
    return {relative(VALIDATION):digest(VALIDATION)}


def verify_environment():
    for name,value in json.loads(DESIGN.read_text())['environment'].items():
        assert os.environ.get(name)==value,(name,'must equal',value)


def evidence(screen_path,review_path):
    """No theory calculations: freshly qualify and reconstruct original evidence."""
    receipt=consumer.verify(Path(screen_path),Path(review_path))
    assert receipt['reviewed_screen_status']=='numerical_sensitivity_requires_followup'
    flags=receipt['reviewed_diagnostic_flags']
    assert flags and all(x in ('total_centered_RMS','total_centered_max_absolute') or x.startswith('component_variation:') for x in flags)
    assert not receipt['reviewed_failures']
    screen=json.loads(Path(screen_path).read_text())
    old_path=ROOT/screen['plan_path'];old=original.sealed_read(old_path)
    one,two=contracts.configuration_pair(old['settings'])
    assert canonical(one)==old['nominal_native_configuration']
    assert canonical(two)==old['doubled_native_configuration']
    refs=[]
    for i,(selected,row,status) in enumerate(zip(old['selected'],screen['points'],receipt['point_statuses'])):
        assert status['reviewed_status']=='finite_native_precision' and not status['reviewed_failed_checks']
        assert row['point']==selected['point'] and row['audit_index']==i
        refs.append(dict(selected,audit_index=i,record_path=row['record_path'],record_sha256=row['record_sha256'],
                         spectrum_path=row['spectrum_path'],spectrum_sha256=row['spectrum_sha256'],
                         original_status=row['status'],original_failed_checks=row['failed_checks']))
    assert len(refs)==32 and refs[0]['point']!=refs[31]['point']
    bound=merge(receipt['input_sha256'],{relative(Path(screen_path)):digest(screen_path),
                 relative(Path(review_path)):digest(review_path),relative(old_path):digest(old_path)})
    return dict(receipt_identity=identity(receipt),parent_proposal_identity=receipt['parent_target_identity'],
                screen_path=relative(screen_path),review_path=relative(review_path),
                original_plan_path=relative(old_path),settings=old['settings'],references=refs,
                native_configuration=canonical(two),versions=old['frozen_target']['versions'],
                input_sha256=bound,original_screen_status=receipt['original_screen_status'],
                thermal_reviewed_status=receipt['reviewed_screen_status']),two


def base_plan(screen,review,kind):
    evidence_value,_=evidence(screen,review)
    return dict(schema='native-accuracy-execution-plan-v2',kind=kind,
                design=json.loads(DESIGN.read_text()),evidence=evidence_value,
                source_sha256=sources(),validation_sha256=validation_guard(),
                converter_contract=contracts.converter_contract())


def finalize_plan(plan,work):
    work=Path(work).resolve();assert work.is_relative_to(ROOT/'.work') and not work.exists()
    plan['identity']=identity(plan)
    work.mkdir(parents=True)
    write_new(work/'plan.json',plan)
    return read(work/'plan.json')


def verify_base(plan_path):
    plan=read(plan_path)
    assert plan['identity']==identity({k:v for k,v in plan.items() if k not in ('identity','payload_sha256')})
    assert plan['design']==json.loads(DESIGN.read_text()) and plan['source_sha256']==sources()
    assert plan['validation_sha256']==validation_guard()
    assert plan['converter_contract']==contracts.converter_contract()
    rebuilt,two=evidence(ROOT/plan['evidence']['screen_path'],ROOT/plan['evidence']['review_path'])
    assert rebuilt==plan['evidence']
    verify_hashes(plan['source_sha256']);verify_hashes(rebuilt['input_sha256'])
    assert all(importlib.metadata.version(k)==v for k,v in rebuilt['versions'].items())
    return plan,two


def load_spectrum(path):
    with np.load(path,allow_pickle=False) as archive:
        value={k:archive[k].copy() for k in archive.files}
    assert set(value)=={'ell',*contracts.SPECTRA}
    assert np.array_equal(value['ell'],np.arange(len(value['ell'])))
    for k in contracts.SPECTRA:assert value[k].shape==value['ell'].shape and np.isfinite(value[k]).all()
    return value


def save_spectrum(path,value):
    with Path(path).open('xb') as stream:np.savez_compressed(stream,**value)


def old_reference(ref):
    assert digest(ROOT/ref['record_path'])==ref['record_sha256']
    row=original.sealed_read(ROOT/ref['record_path'])
    assert row['point']==ref['point'] and row['audit_index']==ref['audit_index']
    assert row['status']==ref['original_status'] and row['failed_checks']==ref['original_failed_checks']
    assert digest(ROOT/ref['spectrum_path'])==ref['spectrum_sha256']
    raw=contracts.raw_from_dl(load_spectrum(ROOT/ref['spectrum_path']))
    return dict(loglikes=row['high_loglikes'],logpriors=row['high_logpriors'],logpost=row['high_logpost'],
                derived=row['high_derived'],H=row['background']['high_H_km_s_Mpc'],
                metadata={'finalized_theory_extra_args':row['finalized_theory_extra_args'],
                          'CAMB_Params_max_l':row['CAMB_Params_max_l'],'provider_length':row['provider_Dl_length'],
                          'units':'FIRASmuK2','ell_factor':False,'pp':'dimensionless C_phi_phi'}),raw


def compare(value,raw,reference,reference_raw,limits):
    assert set(value['loglikes'])==set(reference['loglikes']) and set(value['derived'])==set(reference['derived'])
    assert canonical(value['metadata'])==canonical(reference['metadata'])
    assert np.array_equal(raw['ell'],reference_raw['ell'])
    errors={'component_loglike_absolute_max':max(abs(value['loglikes'][k]-v) for k,v in reference['loglikes'].items()),
            'prior_sum_absolute_max':abs(sum(value['logpriors'])-sum(reference['logpriors'])),
            'logpost_absolute_max':abs(value['logpost']-reference['logpost']),
            'derived_scaled_max':max(abs(value['derived'][k]-v)/(1+abs(v)) for k,v in reference['derived'].items()),
            'background_H_relative_max':float(np.max(abs(np.array(value['H'])/reference['H']-1))),
            'spectrum_max_scale_absolute_max':max(float(np.max(abs(raw[k]-reference_raw[k])))/max(float(np.max(abs(reference_raw[k]))),1e-300) for k in contracts.SPECTRA)}
    assert np.isfinite(list(errors.values())).all()
    return {'errors':errors,'failed_checks':[k for k,v in errors.items() if v>limits[k]]}


def native_value(model,point):
    """Exactly one requested density; model already belongs to this worker."""
    result=model.logposterior(point)
    assert np.isfinite(result.logpost)
    value={'loglikes':dict(zip(model.likelihood,map(float,result.loglikes))),
           'logpriors':list(map(float,result.logpriors)),'logpost':float(result.logpost),
           'derived':dict(zip(model.parameterization.derived_params(),map(float,result.derived)))}
    background=model.provider.get_CAMBdata()
    raw={k:np.array(v,copy=True) for k,v in model.provider.get_Cl(ell_factor=False,units='FIRASmuK2').items() if k in ('ell',*contracts.SPECTRA)}
    dl={k:np.array(v,copy=True) for k,v in model.provider.get_Cl(ell_factor=True,units='FIRASmuK2').items() if k in ('ell',*contracts.SPECTRA)}
    converted=contracts.raw_from_dl(dl)
    assert np.array_equal(raw['ell'],converted['ell'])
    assert max(float(np.max(abs(raw[k]-converted[k])))/max(float(np.max(abs(raw[k]))),1e-300) for k in contracts.SPECTRA)<1e-13
    value['H']=background.hubble_parameter(GRID).tolist()
    value['metadata']={'finalized_theory_extra_args':dict(model.theory['camb'].extra_args),
                       'CAMB_Params_max_l':int(background.Params.max_l),'provider_length':len(raw['ell']),
                       'units':'FIRASmuK2','ell_factor':False,'pp':'dimensionless C_phi_phi'}
    consumer.finite_tree(value)
    return value,raw


def request_paths(work,index):
    base=Path(work)/f'{index:05d}'
    return {kind:base.with_suffix(suffix) for kind,suffix in [('record','.json'),('ticket','.ticket.json'),
            ('claim','.claim.json'),('spectrum','.raw-cl.npz')]}


def request_identity(plan,request):
    return {'plan_identity':plan['identity'],'kind':plan['kind'],'request':request}


def completed(work,plan,request):
    paths=request_paths(work,request['index']);row=read(paths['record'])
    assert row['binding']==request_identity(plan,request)
    assert row['native_logposterior_invocations'] in (0,1)
    if row['origin']=='persistent_native':
        for kind in ('ticket','claim'):
            ticket=read(paths[kind]);assert ticket['binding']==request_identity(plan,request)
            assert row[kind+'_sha256']==digest(paths[kind])
    if row['status']=='finite_native_accuracy2':
        assert row['native_logposterior_invocations']==(0 if row['origin']=='reused_fixed32' else 1)
        assert not row['failed_checks']
        assert row['spectrum_path']==relative(paths['spectrum']) and digest(paths['spectrum'])==row['spectrum_sha256']
        load_spectrum(paths['spectrum'])
    row.update(record_path=relative(paths['record']),record_sha256=digest(paths['record']))
    return row


def incomplete(work,plan,request):
    paths=request_paths(work,request['index'])
    assert not paths['record'].exists()
    claimed=paths['claim'].exists()
    return {'status':'incomplete_attempt_no_retry','binding':request_identity(plan,request),
            'unknown_native_invocation_count':claimed,'native_logposterior_invocations':0,
            'known_unstarted':not claimed,'existing_files':{relative(p):digest(p) for p in paths.values() if p.exists()}}


def collect(work,plan):
    ledger=Path(work)/'record-hashes.json'
    if ledger.exists():
        old=read(ledger);assert old['plan_identity']==plan['identity'];verify_hashes(old['files'])
    rows=[]
    for request in plan['requests']:
        paths=request_paths(work,request['index'])
        rows.append(completed(work,plan,request) if paths['record'].exists() else incomplete(work,plan,request))
    files={}
    for request in plan['requests']:
        files.update({relative(p):digest(p) for p in request_paths(work,request['index']).values() if p.exists()})
    for p in sorted(Path(work).glob('worker-*')):
        if p.is_file():files[relative(p)]=digest(p)
    return rows,files


def guards_without_physics():
    import camb
    import cobaya.model
    def forbidden(*args,**kwargs):raise RuntimeError('Read-only native2 qualification forbids physical calculations.')
    stack=ExitStack()
    for owner,name in [(camb,'get_results'),(camb,'get_background'),(camb,'get_transfer_functions'),
                       (cobaya.model,'get_model'),(cobaya.model.Model,'__init__')]:
        stack.enter_context(patch.object(owner,name,forbidden))
    return stack


def construct_model(info):
    from cobaya.model import get_model
    model=get_model(copy.deepcopy(info))
    assert set(model.theory)=={'camb','camb.transfers'}
    model.add_requirements({'CAMBdata':None})
    return model


def persistent_worker(work,worker,verify_plan,build_record):
    """One process, one model, a fixed sequence; never replaces/retries a worker."""
    work=Path(work).resolve();plan,info=verify_plan(work/'plan.json');verify_environment()
    requests=[r for r in plan['requests'] if r.get('worker')==worker and r['origin']=='persistent_native']
    assert requests and [r['ordinal'] for r in requests]==list(range(len(requests)))
    receipt=work/f'worker-{worker}.json';assert not receipt.exists()
    started=time.monotonic();model=construct_model(info);construction=time.monotonic()-started
    seen={};completed_indices=[];failure=None
    try:
        for request in requests:
            paths=request_paths(work,request['index']);binding=request_identity(plan,request)
            assert read(paths['ticket'])['binding']==binding
            assert not paths['claim'].exists() and not paths['record'].exists() and not paths['spectrum'].exists()
            write_new(paths['claim'],{'binding':binding,'pid':os.getpid(),'worker':worker,'ordinal':request['ordinal']})
            start=time.monotonic();invocations=0
            row={'binding':binding,'origin':'persistent_native','point':request['point'],
                 'pid':os.getpid(),'worker':worker,'ordinal':request['ordinal'],
                 'model_construction_seconds':construction,
                 'ticket_sha256':digest(paths['ticket']),'claim_sha256':digest(paths['claim'])}
            try:
                invocations=1
                value,raw=native_value(model,request['point'])
                key=identity(request['point'])
                repeated=compare(value,raw,*seen[key],plan['design']['pilot_tolerances']) if key in seen else None
                row.update(build_record(plan,request,value,raw))
                row['repeat_replay']=repeated
                if repeated and repeated['failed_checks']:row['failed_checks']+=['repeat:'+x for x in repeated['failed_checks']]
                seen[key]=(copy.deepcopy(value),{k:v.copy() for k,v in raw.items()})
                save_spectrum(paths['spectrum'],raw)
                row.update(value=value,spectrum_path=relative(paths['spectrum']),spectrum_sha256=digest(paths['spectrum']))
                row['status']='failed_native_accuracy2_checks' if row['failed_checks'] else 'finite_native_accuracy2'
            except Exception as exc:
                row.update(status='failed_native_accuracy2_evaluation',failed_checks=['evaluation_exception'],
                           error_type=type(exc).__name__,error=str(exc))
            row.update(native_logposterior_invocations=invocations,request_seconds=time.monotonic()-start)
            write_new(paths['record'],row);completed_indices.append(request['index'])
            if row['status']!='finite_native_accuracy2':failure=request['index'];break
    finally:
        if hasattr(model,'close'):model.close()
        # Fail closed even after completed requests if code/input bytes changed.
        verify_plan(work/'plan.json')
        write_new(receipt,{'plan_identity':plan['identity'],'worker':worker,'pid':os.getpid(),
                          'completed_indices':completed_indices,'failed_index':failure,
                          'model_constructions':1,'model_construction_seconds':construction,
                          'elapsed_seconds':time.monotonic()-started})
    return failure is None


def execute_workers(work,module,before_workers=None):
    """Explicit fixed subprocess pool; no automatic resubmission or replacement."""
    import subprocess
    import sys
    from concurrent.futures import ThreadPoolExecutor
    work=Path(work).resolve();plan=read(work/'plan.json');verify_environment()
    execution=work/'execution.json'
    assert not execution.exists() and not any(work.glob('worker-*'))
    assert not any(request_paths(work,r['index'])['claim'].exists() for r in plan['requests'])
    write_new(execution,{'plan_identity':plan['identity'],'coordinator_pid':os.getpid(),
                         'workers':plan['workers'],'retry_permitted':False})
    if before_workers is not None:before_workers()
    for request in plan['requests']:
        if request['origin']=='persistent_native':
            write_new(request_paths(work,request['index'])['ticket'],{'binding':request_identity(plan,request)})
    def launch(worker):
        log=work/f'worker-{worker}.log'
        with log.open('x') as stream:
            child=subprocess.Popen([sys.executable,str(HERE/module),'worker','--work',str(work),'--worker',str(worker)],
                                   cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,env=os.environ.copy())
            code=child.wait()
        write_new(work/f'worker-{worker}-process.json',{'plan_identity':plan['identity'],
                  'worker':worker,'pid':child.pid,'returncode':code,'log_path':relative(log),'log_sha256':digest(log)})
        return code
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
        codes=list(pool.map(launch,range(plan['workers'])))
    return codes


def seal_summary(work,plan,report,files,output=None):
    work=Path(work)
    write_new(work/'record-hashes.json',{'plan_identity':plan['identity'],'files':files})
    report=dict(report,record_manifest_path=relative(work/'record-hashes.json'),
                record_manifest_sha256=digest(work/'record-hashes.json'))
    write_new(work/'summary.json',report)
    if output is not None:write_new(output,report)
    return report


def verify_saved_summary(work,report,files,summary_path=None):
    work=Path(work);ledger=read(work/'record-hashes.json')
    assert ledger==dict(plan_identity=report['plan_identity'],files=files,payload_sha256=identity({'plan_identity':report['plan_identity'],'files':files}))
    saved=read(work/'summary.json')
    expected=dict(report,record_manifest_path=relative(work/'record-hashes.json'),record_manifest_sha256=digest(work/'record-hashes.json'))
    assert plain(saved)==expected
    if summary_path is not None:assert plain(read(summary_path))==expected
    return saved


def process_integrity(work,plan,rows):
    """All successful rows must originate in their single fixed worker process."""
    pids=[];failures=[];bindings={}
    for worker in range(plan['workers']):
        paths=[Path(work)/f'worker-{worker}.json',Path(work)/f'worker-{worker}-process.json']
        if not all(p.exists() for p in paths):failures.append(f'worker:{worker}:incomplete');continue
        receipt,process=map(read,paths)
        expected=[r['index'] for r in plan['requests'] if r.get('worker')==worker and r['origin']=='persistent_native']
        assert receipt['plan_identity']==process['plan_identity']==plan['identity']
        assert receipt['worker']==process['worker']==worker and receipt['pid']==process['pid']
        assert process['log_path']==relative(Path(work)/f'worker-{worker}.log')
        assert digest(ROOT/process['log_path'])==process['log_sha256']
        bindings.update({relative(p):digest(p) for p in paths})
        if process['returncode']!=0 or receipt['failed_index'] is not None or receipt['completed_indices']!=expected:
            failures.append(f'worker:{worker}:failed_or_incomplete')
        assert receipt['model_constructions']==1;pids.append(receipt['pid'])
        for index in receipt['completed_indices']:
            row=rows[index];assert row['worker']==worker and row['pid']==receipt['pid']
            assert read(request_paths(work,index)['claim'])['pid']==receipt['pid']
    assert len(pids)==len(set(pids))
    return failures
