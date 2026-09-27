"""Separate gated accuracy2 replay and same32point accuracy3 screen."""
import os
for _key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[_key]='1'
import argparse
import copy
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

import native_posterior_precision as original
import native_precision_review_consumer as consumer
from target_identity import canonical

ROOT=original.ROOT
HERE=Path(__file__).resolve().parent
DESIGN=HERE/'native_precision_refinement_design.json'
VALIDATION=ROOT/'studies/unified_cosmology/results/inference/native_precision_refinement_validation.json'


def dependencies():
    paths=[Path(__file__),DESIGN,HERE/'native_precision_refinement_validate.py',
           Path(original.__file__),original.DESIGN,Path(consumer.__file__),consumer.DESIGN,
           Path(consumer.thermal.__file__),consumer.thermal.DESIGN,
           HERE/'measurement_summary.py',HERE/'expansion_history.py',HERE/'likelihood.py',
           HERE/'exact_correction.py',HERE/'target_identity.py']
    return {original.relative(p):original.digest(p) for p in paths}


def write_once(path,value):
    path=Path(path)
    if path.exists():
        assert {k:v for k,v in original.sealed_read(path).items() if k!='payload_sha256'}==value, 'Immutable record already differs.'
    else:
        path.parent.mkdir(parents=True,exist_ok=True)
        original.sealed_write(path,value)


def validation_guard():
    value=json.loads(VALIDATION.read_text())
    assert value['status']=='passed_synthetic_native_refinement_validation'
    assert value['source_sha256']==dependencies() and value['physical_calls']==0
    return original.digest(VALIDATION)


def configurations(settings,old_design,new_design):
    _,two=original.native_configurations(settings,old_design)
    assert new_design['numerical_controls']==old_design['numerical_controls']
    assert new_design['reference_value']==old_design['comparison_value']==2
    three=copy.deepcopy(two)
    for key in new_design['numerical_controls']:
        assert two['theory']['camb']['extra_args'][key]==2
        three['theory']['camb']['extra_args'][key]=new_design['refinement_value']
    restored=copy.deepcopy(three)
    for key in new_design['numerical_controls']:
        restored['theory']['camb']['extra_args'][key]=2
    assert canonical(restored)==canonical(two)
    return two,three


def prepare(screen_path,review_path,work):
    source=dependencies()
    receipt=consumer.verify(screen_path,review_path) # Always before cache/model work.
    assert receipt['reviewed_screen_status']=='numerical_sensitivity_requires_followup', 'Require genuine surviving numerical density variation.'
    flags=receipt['reviewed_diagnostic_flags']
    assert flags and all(x in ('total_centered_RMS','total_centered_max_absolute') or x.startswith('component_variation:') for x in flags)
    validation=validation_guard()
    screen=json.loads(Path(screen_path).read_text())
    old_plan_path=ROOT/screen['plan_path'];old_plan=original.sealed_read(old_plan_path)
    design=json.loads(DESIGN.read_text())
    two,three=configurations(old_plan['settings'],old_plan['design'],design)
    assert canonical(two)==old_plan['doubled_native_configuration']
    points=[]
    for index,(item,row) in enumerate(zip(old_plan['selected'],screen['points'])):
        assert item['point']==row['point'] and row['audit_index']==index
        points.append(dict(item,audit_index=index,accuracy2_record_path=row['record_path'],
                           accuracy2_record_sha256=row['record_sha256'],
                           accuracy2_spectrum_path=row['spectrum_path'],accuracy2_spectrum_sha256=row['spectrum_sha256']))
    assert len(points)==32
    plan=dict(schema='native-precision-refinement-plan-v1',design=design,points=points,
              source_sha256=source,validation_sha256=validation,
              screen_path=original.relative(screen_path),screen_sha256=original.digest(screen_path),
              review_path=original.relative(review_path),review_sha256=original.digest(review_path),
              original_plan_path=original.relative(old_plan_path),original_plan_sha256=original.digest(old_plan_path),
              verified_original_receipt_identity=original.identity(receipt),
              original_repaired_one_to_two=receipt['reviewed_screen'],
              original_unrepaired_screen=receipt['original_screen'],
              settings=old_plan['settings'],original_design=old_plan['design'],
              accuracy2_configuration=canonical(two),accuracy3_configuration=canonical(three),
              parent_target_identity=receipt['parent_target_identity'])
    assert dependencies()==source
    plan['identity']=original.identity(plan)
    work=Path(work).resolve();assert work.is_relative_to(ROOT/'.work')
    write_once(work/'plan.json',plan)
    return original.sealed_read(work/'plan.json')


def verify_plan(plan_path):
    plan=original.sealed_read(plan_path)
    assert plan['identity']==original.identity({k:v for k,v in plan.items() if k not in ('identity','payload_sha256')})
    assert plan['source_sha256']==dependencies() and plan['validation_sha256']==validation_guard()
    assert plan['design']==json.loads(DESIGN.read_text())
    for prefix in ('screen','review','original_plan'):
        assert original.digest(ROOT/plan[prefix+'_path'])==plan[prefix+'_sha256']
    receipt=consumer.verify(ROOT/plan['screen_path'],ROOT/plan['review_path'])
    assert original.identity(receipt)==plan['verified_original_receipt_identity']
    assert receipt['reviewed_screen_status']=='numerical_sensitivity_requires_followup'
    assert plan['parent_target_identity']==receipt['parent_target_identity']
    screen=json.loads((ROOT/plan['screen_path']).read_text())
    assert screen['plan_path']==plan['original_plan_path']
    old_plan=original.sealed_read(ROOT/plan['original_plan_path'])
    assert plan['settings']==old_plan['settings'] and plan['original_design']==old_plan['design']
    assert plan['original_repaired_one_to_two']==receipt['reviewed_screen']
    assert plan['original_unrepaired_screen']==receipt['original_screen']
    assert len(plan['points'])==32
    for i,(item,old) in enumerate(zip(plan['points'],old_plan['selected'])):
        assert item['audit_index']==i and all(item[k]==v for k,v in old.items())
        for label in ('record','spectrum'):
            assert item['accuracy2_'+label+'_path']==screen['points'][i][label+'_path']
            assert item['accuracy2_'+label+'_sha256']==screen['points'][i][label+'_sha256']
        for label in ('accuracy2_record','accuracy2_spectrum'):
            assert original.digest(ROOT/item[label+'_path'])==item[label+'_sha256']
    two,three=configurations(plan['settings'],plan['original_design'],plan['design'])
    assert canonical(two)==plan['accuracy2_configuration'] and canonical(three)==plan['accuracy3_configuration']
    return plan


def reference(plan,index):
    item=plan['points'][index]
    row=original.sealed_read(ROOT/item['accuracy2_record_path'])
    assert row['audit_index']==index and row['point']==item['point']
    assert row['status'] in ('finite_native_precision','failed_native_precision_checks')
    return row


def replay_checks(reference_row,row,spectra,old_spectra,design):
    limits=design['replay_tolerances']
    assert set(reference_row['high_loglikes'])==set(row['loglikes'])
    assert set(reference_row['high_derived'])==set(row['derived'])
    errors={'component_loglike_absolute_max':max(abs(row['loglikes'][k]-v) for k,v in reference_row['high_loglikes'].items()),
            'prior_sum_absolute_max':abs(sum(row['logpriors'])-sum(reference_row['high_logpriors'])),
            'logpost_absolute_max':abs(row['logpost']-reference_row['high_logpost']),
            'derived_scaled_max':max(abs(row['derived'][k]-v)/(1+abs(v)) for k,v in reference_row['high_derived'].items())}
    assert set(spectra)==set(old_spectra)=={'ell','tt','ee','bb','te','pp'}
    assert np.array_equal(spectra['ell'],old_spectra['ell'])
    spectral={}
    for key in ('tt','ee','bb','te','pp'):
        x,y=np.asarray(spectra[key]),np.asarray(old_spectra[key])
        assert x.shape==y.shape and np.isfinite(x).all() and np.isfinite(y).all()
        spectral[key]=float(np.max(abs(x-y))/max(float(np.max(abs(y))),1e-300))
    errors['spectrum_max_scale_absolute_max']=max(spectral.values())
    assert np.isfinite(list(errors.values())).all()
    return dict(errors={k:float(v) for k,v in errors.items()},spectral_errors=spectral,
                failed_checks=[k for k,v in errors.items() if v>limits[k]])


def background_checks(reference_row,row,values,old_design):
    limits=old_design['diagnostic_thresholds']
    high=reference_row['high_derived']
    assert np.isfinite([*high.values(),values['rdrag'],values['omegam'],*values['H'],
                        *values['expansion'].values(),*reference_row['background']['high_H_km_s_Mpc']]).all()
    errors={'rdrag_relative':abs(values['rdrag']/high['rdrag']-1),
            'omegam_absolute':abs(values['omegam']-high['omegam']),
            'q_scaled':max(abs(values['expansion'][k]-high[k])/(1+abs(high[k])) for k in ('q0','q05','q1')),
            'j_scaled':max(abs(values['expansion'][k]-high[k])/(1+abs(high[k])) for k in ('j0','j05','j1')),
            'H_relative':float(np.max(abs(np.asarray(values['H'])/reference_row['background']['high_H_km_s_Mpc']-1)))}
    thresholds={k:limits['nominal_background_'+k+'_closure_max'] for k in errors if k!='H_relative'}
    thresholds['H_relative']=json.loads(consumer.thermal.DESIGN.read_text())['high_background_H_relative_max']
    assert np.isfinite(list(errors.values())).all()
    return dict(errors=errors,failed_checks=['reference_accuracy2_background_'+k for k,v in errors.items() if v>thresholds[k]])


def file_paths(work,stage,index):
    stem='replay-00' if stage=='replay' else f'refine-{index:02d}'
    return Path(work)/(stem+'.json'),Path(work)/(stem+'.log')


def attempt_paths(work,stage,index):
    path,_=file_paths(work,stage,index)
    return path.with_name(path.stem+'-ticket.json'),path.with_name(path.stem+'-claimed.json')


def exclusive_payload(path,value):
    """The process claim, unlike an idempotent result replay, is one-use."""
    with Path(path).open('x') as stream:
        json.dump(dict(value,payload_sha256=original.identity(value)),stream,indent=2,allow_nan=False)
        stream.write('\n')


def finite_payload(value):
    if isinstance(value,dict):return {k:finite_payload(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [finite_payload(v) for v in value]
    if isinstance(value,(float,np.floating)) and not np.isfinite(value):return None
    return value


def verify_ledger(work,stage):
    ledger=Path(work)/(stage+'-record-hashes.json')
    if ledger.exists():
        mapping=original.sealed_read(ledger)['files']
        consumer.hashes(mapping)
    return ledger


def require_replay(work,plan):
    ledger=verify_ledger(work,'replay');assert ledger.exists(), 'Replay hash ledger is required.'
    path=Path(work)/'replay-summary.json';summary=original.sealed_read(path)
    assert summary['plan_identity']==plan['identity']
    assert summary['status']=='accuracy2_replay_passed', 'Accuracy2 replay has not passed; refinement prohibited.'
    row=read_completed(work,'replay',0,plan)
    assert row['status']=='accuracy2_replay_passed' and summary['point_record_sha256']==row['record_sha256']
    assert summary['known_native_logposterior_invocations']==1 and summary['unknown_count_attempts']==0
    bound=original.sealed_read(ledger)
    assert bound['plan_identity']==plan['identity']
    assert bound['files']=={row[k+'_path']:row[k+'_sha256'] for k in ('record','log','spectrum','ticket','claim')}
    return original.digest(path)


def read_completed(work,stage,index,plan):
    path,log=file_paths(work,stage,index)
    row=original.sealed_read(path)
    assert row['plan_identity']==plan['identity'] and row['stage']==stage and row['audit_index']==index
    assert row['point']==plan['points'][index]['point']
    assert row['chain']==plan['points'][index]['chain']
    assert row['accuracy2_reference_sha256']==plan['points'][index]['accuracy2_record_sha256']
    assert row['plan_sha256']==original.digest(Path(work)/'plan.json')
    assert row['native_logposterior_invocations'] in (0,1)
    if row['status'] in ('finite_native_precision','accuracy2_replay_passed'):
        assert row['native_logposterior_invocations']==1 and not row['failed_checks']
    if stage=='refine':assert row['replay_gate_sha256']==require_replay(work,plan)
    else:assert row['replay_gate_sha256'] is None
    for kind,p in zip(('ticket','claim'),attempt_paths(work,stage,index)):
        value=original.sealed_read(p)
        assert value['plan_identity']==plan['identity'] and value['stage']==stage and value['index']==index
        row[kind+'_path']=original.relative(p);row[kind+'_sha256']=original.digest(p)
    if 'spectrum_path' in row:
        assert original.digest(ROOT/row['spectrum_path'])==row['spectrum_sha256']
    row.update(record_path=original.relative(path),record_sha256=original.digest(path),
               log_path=original.relative(log),log_sha256=original.digest(log))
    return row


def worker(plan_path,stage,index):
    plan_path=Path(plan_path).resolve();work=plan_path.parent
    plan=verify_plan(plan_path)
    assert stage in ('replay','refine') and 0<=index<32
    if stage=='replay':assert index==0
    gate=require_replay(work,plan) if stage=='refine' else None
    path,_=file_paths(work,stage,index);assert not path.exists()
    ticket,claim=attempt_paths(work,stage,index)
    reserved=original.sealed_read(ticket)
    assert {k:v for k,v in reserved.items() if k!='payload_sha256'}==dict(plan_identity=plan['identity'],stage=stage,index=index)
    exclusive_payload(claim,dict(plan_identity=plan['identity'],stage=stage,index=index))
    item=plan['points'][index];old=reference(plan,index)
    row=dict(stage=stage,audit_index=index,chain=item['chain'],point=item['point'],
             plan_identity=plan['identity'],plan_sha256=original.digest(plan_path),
             accuracy2_reference_sha256=item['accuracy2_record_sha256'],
             native_logposterior_invocations=0,auxiliary_background_calls=0,replay_gate_sha256=gate)
    started=time.monotonic()
    try:
        from cobaya.model import get_model
        from expansion_history import native_thermal_parameters
        from likelihood import expansion_diagnostics
        import camb
        two,three=configurations(plan['settings'],plan['original_design'],plan['design'])
        with get_model(two if stage=='replay' else three) as model:
            model.add_requirements({'CAMBdata':None})
            assert 'camb' in model.theory and set(model.theory)<={'camb','camb.transfers'}
            row['native_logposterior_invocations']=1
            evaluated=model.logposterior(item['point'])
            assert np.isfinite(evaluated.logpost)
            row.update(loglikes=dict(zip(model.likelihood,map(float,evaluated.loglikes))),
                       logpriors=list(map(float,evaluated.logpriors)),logpost=float(evaluated.logpost),
                       derived=dict(zip(model.parameterization.derived_params(),map(float,evaluated.derived))))
            full=model.provider.get_CAMBdata();raw=model.provider.get_Cl(ell_factor=True)
            spectra={k:np.asarray(raw[k]) for k in ('ell','tt','ee','bb','te','pp')}
            assert all(np.isfinite(x).all() for x in spectra.values())
            spectrum=path.with_name(path.stem+'-spectra.npz');assert not spectrum.exists()
            np.savez_compressed(spectrum,**spectra)
            args=dict(model.theory['camb'].extra_args)
            restored=copy.deepcopy(args)
            for key in plan['design']['numerical_controls']:restored[key]=2
            assert canonical(restored)==canonical(old['finalized_theory_extra_args'])
            row.update(finalized_theory_extra_args=args,CAMB_Params_max_l=int(full.Params.max_l),
                       provider_Dl_length=len(raw['ell']),spectrum_path=original.relative(spectrum),
                       spectrum_sha256=original.digest(spectrum))
            low_parameters=full.Params.copy()
            for key in plan['design']['numerical_controls']:setattr(low_parameters.Accuracy,key,2)
            row['auxiliary_background_calls']=1
            low=camb.get_background(native_thermal_parameters(low_parameters))
            z=np.array([0.,.5,1.,2.33,1100.]);grid=np.concatenate([np.arange(5)*.001+v for v in (0.,.5,1.)])
            values=dict(rdrag=float(low.get_derived_params()['rdrag']),omegam=float(low_parameters.omegam),
                        H=low.hubble_parameter(z).tolist(),expansion=expansion_diagnostics(grid,low.hubble_parameter(grid)))
            row['reference_accuracy2_background']=values
            row['background_checks']=background_checks(old,row,values,plan['original_design'])
            row['high_H_km_s_Mpc']=full.hubble_parameter(z).tolist()
            low_density={'exact_loglikes':old['high_loglikes'],'exact_logpost':old['high_logpost']}
            row['comparison']=original.density_comparison(low_density,row['loglikes'],row['logpriors'],row['logpost'],plan['original_design']['diagnostic_thresholds'])
            failed=row['comparison']['failed_density_checks']+row['background_checks']['failed_checks']
            if stage=='replay':
                with np.load(ROOT/item['accuracy2_spectrum_path']) as saved:
                    row['replay_checks']=replay_checks(old,row,spectra,dict(saved),plan['design'])
                failed+=row['replay_checks']['failed_checks']
            row['failed_checks']=failed
            row['status']=('accuracy2_replay_passed' if not failed else 'accuracy2_replay_failed') if stage=='replay' else ('finite_native_precision' if not failed else 'failed_native_precision_checks')
            json.dumps(row,allow_nan=False)
    except Exception as error:
        # Preserve every available native field; nonfinite quantities are null
        # only in an explicitly failed record, never an accepted comparison.
        row=finite_payload(row)
        row.update(status='exception',error=repr(error))
        spectrum=path.with_name(path.stem+'-spectra.npz')
        if spectrum.exists():row.update(spectrum_path=original.relative(spectrum),spectrum_sha256=original.digest(spectrum))
    row['seconds']=time.monotonic()-started
    verify_plan(plan_path)
    if stage=='refine':assert require_replay(work,plan)==gate
    write_once(path,row)


def collect_or_execute(work,stage,index,plan,execute):
    path,log=file_paths(work,stage,index)
    ticket,claim=attempt_paths(work,stage,index)
    if not path.exists():
        if any(p.exists() for p in (log,ticket,claim)):
            value=dict(audit_index=index,status='incomplete_attempt_no_retry',unknown_native_invocation_count=True)
            for kind,p in [('log',log),('ticket',ticket),('claim',claim)]:
                if p.exists():value[kind+'_path']=original.relative(p);value[kind+'_sha256']=original.digest(p)
            return value
        if not execute:return dict(audit_index=index,status='not_attempted',unknown_native_invocation_count=False)
        exclusive_payload(ticket,dict(plan_identity=plan['identity'],stage=stage,index=index))
        with log.open('x') as output:
            run=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--plan',str(Path(work)/'plan.json'),
                                '--stage',stage,'--index',str(index)],stdout=output,stderr=subprocess.STDOUT)
        if not path.exists():
            value=collect_or_execute(work,stage,index,plan,False)
            value.update(status='worker_failed_before_record',returncode=run.returncode)
            return value
    return read_completed(work,stage,index,plan)


def run_stage(work,stage,execute,workers):
    work=Path(work).resolve();plan=verify_plan(work/'plan.json')
    assert stage in ('replay','refine') and 1<=workers<=plan['design']['maximum_workers']
    if stage=='refine':gate=require_replay(work,plan)
    else:gate=None
    ledger=verify_ledger(work,stage)
    ids=[0] if stage=='replay' else list(range(32))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        rows=list(pool.map(lambda i:collect_or_execute(work,stage,i,plan,execute),ids))
    files={}
    for row in rows:
        for kind in ('record','log','spectrum','ticket','claim'):
            if kind+'_path' in row:files[row[kind+'_path']]=row[kind+'_sha256']
    if execute:write_once(ledger,dict(plan_identity=plan['identity'],files=files))
    known=sum(r.get('native_logposterior_invocations',0) for r in rows)
    unknown=sum(bool(r.get('unknown_native_invocation_count')) for r in rows)
    if stage=='replay':
        result=dict(status=rows[0]['status'],point_record_sha256=rows[0].get('record_sha256'))
    else:
        design=copy.deepcopy(plan['original_design'])
        design.update(threshold_interpretation=plan['design']['interpretation'],limitations=plan['design']['interpretation'])
        result=original.summarize(rows,design)
        result['comparison']='accuracy3 minus accuracy2 at the unchanged original32 points'
    result.update(stage=stage,plan_identity=plan['identity'],plan_path=original.relative(work/'plan.json'),
                  plan_sha256=original.digest(work/'plan.json'),points=rows,
                  original_repaired_one_to_two=plan['original_repaired_one_to_two'],
                  original_unrepaired_screen=plan['original_unrepaired_screen'],
                  known_native_logposterior_invocations=known,unknown_count_attempts=unknown,
                  cumulative_known_native_logposterior_invocations=known+(1 if stage=='refine' else 0),
                  known_auxiliary_background_calls=sum(r.get('auxiliary_background_calls',0) for r in rows),
                  posterior_qualification=False,posterior_reweighting_performed=False,
                  source_sha256=plan['source_sha256'],replay_gate_sha256=gate)
    verify_plan(work/'plan.json')
    if stage=='refine':assert require_replay(work,plan)==gate
    if execute:write_once(work/(stage+'-summary.json'),result)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=('prepare','replay','refine','report','worker'))
    p.add_argument('--screen',type=Path);p.add_argument('--thermal-review',type=Path);p.add_argument('--work',type=Path)
    p.add_argument('--plan',type=Path);p.add_argument('--stage',choices=('replay','refine'));p.add_argument('--index',type=int)
    p.add_argument('--workers',type=int,default=1);p.add_argument('--output',type=Path)
    a=p.parse_args()
    if a.action=='worker':worker(a.plan,a.stage,a.index);return
    if a.action=='prepare':
        assert a.screen and a.thermal_review and a.work
        value=prepare(a.screen.resolve(),a.thermal_review.resolve(),a.work.resolve())
        print(json.dumps({'status':'prepared_only','identity':value['identity'],'native_logposterior_invocations':0}));return
    assert a.work
    stage=a.stage if a.action=='report' else a.action
    assert stage
    value=run_stage(a.work,stage,a.action!='report',a.workers)
    if a.output:write_once(a.output,value)
    print(json.dumps({k:value[k] for k in ('status','known_native_logposterior_invocations','unknown_count_attempts','posterior_qualification')}))


if __name__=='__main__':main()
