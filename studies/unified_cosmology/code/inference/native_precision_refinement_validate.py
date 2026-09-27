"""Synthetic-only refinement contracts; no model, spectra or background call."""
import argparse
from contextlib import ExitStack
import copy
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

import numpy as np

import native_precision_refinement as new
from native_precision_thermal_review_validate import fixture

ROOT=new.ROOT


def refuses(operation):
    try:operation()
    except (AssertionError,KeyError,ValueError,FileExistsError,FileNotFoundError):return
    raise AssertionError('Invalid refinement fixture accepted.')


def arithmetic():
    design=json.loads(new.DESIGN.read_text());old_design=json.loads(new.original.DESIGN.read_text())
    reference,_,values=fixture(old_design)
    row=dict(loglikes=copy.deepcopy(reference['high_loglikes']),logpriors=reference['high_logpriors'],
             logpost=reference['high_logpost'],derived=copy.deepcopy(reference['high_derived']))
    spectra={k:np.arange(8,dtype=float)+1 for k in ('ell','tt','ee','bb','te','pp')}
    assert not new.replay_checks(reference,row,spectra,spectra,design)['failed_checks']
    failures=[]
    for name in ('component','prior','logpost','derived','spectrum'):
        r=copy.deepcopy(row);s=copy.deepcopy(spectra)
        if name=='component':r['loglikes']['cmb']+=1e-5
        if name=='prior':r['logpriors']=[3.00001]
        if name=='logpost':r['logpost']+=1e-5
        if name=='derived':r['derived']['rdrag']+=.01
        if name=='spectrum':s['tt'][3]+=.01
        failed=new.replay_checks(reference,r,s,spectra,design)['failed_checks'];assert failed
        failures.append({'mutation':name,'failed_checks':failed})
    background=dict(rdrag=values['nominal']['rdrag'],omegam=values['nominal']['omegam'],
                    expansion=values['nominal']['expansion'],H=values['nominal']['H_km_s_Mpc'])
    assert not new.background_checks(reference,row,background,old_design)['failed_checks']
    for key in ('rdrag','omegam','H','expansion'):
        b=copy.deepcopy(background)
        if key=='H':b[key][3]=np.nan
        elif key=='expansion':b[key]['q1']=np.nan
        else:b[key]=np.nan
        refuses(lambda:new.background_checks(reference,row,b,old_design))
    summaries={}
    for mode in ('constant','cancelling','total'):
        rows=[]
        for i in range(32):
            r,stored,_=fixture(old_design,i)
            delta=20. if mode=='constant' else (-1 if i%2 else 1)*.3
            r['high_loglikes']['cmb']+=delta
            if mode=='cancelling':r['high_loglikes']['bao']-=delta
            r['high_logpost']=sum(r['high_loglikes'].values())+3
            r['comparison']=new.original.density_comparison(stored,r['high_loglikes'],[3.],r['high_logpost'],old_design['diagnostic_thresholds'])
            r.update(status='finite_native_precision',failed_checks=[]);rows.append(r)
        summaries[mode]=new.original.summarize(rows,old_design)
    assert summaries['constant']['status']=='no_large_variation_detected_on_fixed32'
    assert {'component_variation:cmb','component_variation:bao'}<=set(summaries['cancelling']['diagnostic_flags'])
    assert 'total_centered_RMS' in summaries['total']['diagnostic_flags']
    return {'replay_mutations':failures,'nonfirst_background_NaN_refused':True,
            'spread_flags':{k:v['diagnostic_flags'] for k,v in summaries.items()}}


def make_fixture(directory):
    old_design=json.loads(new.original.DESIGN.read_text())
    selected=[];screen_rows=[]
    spectra={k:np.arange(8,dtype=float)+1 for k in ('ell','tt','ee','bb','te','pp')}
    for i in range(32):
        row,_,_=fixture(old_design,i)
        path=directory/f'old-{i:02d}.json';new.write_once(path,row)
        spectrum=directory/f'old-{i:02d}.npz';np.savez(spectrum,**spectra)
        item=dict(point=row['point'],chain=i//8,stratum=i%8,parent_index=i*8,expanded_index=i*100)
        selected.append(item)
        screen_rows.append(dict(row,record_path=new.original.relative(path),record_sha256=new.original.digest(path),
                                spectrum_path=new.original.relative(spectrum),spectrum_sha256=new.original.digest(spectrum)))
    two={'theory':{'camb':{'extra_args':{k:2 for k in old_design['numerical_controls']}}}}
    three=copy.deepcopy(two)
    for k in old_design['numerical_controls']:three['theory']['camb']['extra_args'][k]=3
    old_plan=directory/'old-plan.json';new.write_once(old_plan,dict(selected=selected,settings={'synthetic':True},design=old_design,doubled_native_configuration=two))
    screen=directory/'screen.json';screen.write_text(json.dumps(dict(plan_path=new.original.relative(old_plan),points=screen_rows)))
    review=directory/'review.json';review.write_text('{}')
    receipt=dict(reviewed_screen_status='numerical_sensitivity_requires_followup',reviewed_diagnostic_flags=['total_centered_RMS'],
                 reviewed_screen={'status':'numerical_sensitivity_requires_followup'},original_screen={'status':'incomplete_or_failed_numerical_screen'},
                 parent_target_identity='synthetic-qualified-parent')
    return screen,review,receipt,two,three,spectra


def cached_point(work,stage,index,plan,spectra,passed=True):
    path,log=new.file_paths(work,stage,index);log.write_text('Synthetic fixture only; no model calls.\n')
    for p in new.attempt_paths(work,stage,index):new.exclusive_payload(p,dict(plan_identity=plan['identity'],stage=stage,index=index))
    spec=path.with_name(path.stem+'-spectra.npz');np.savez(spec,**spectra)
    row=dict(stage=stage,audit_index=index,chain=index//8,point=plan['points'][index]['point'],plan_identity=plan['identity'],
             plan_sha256=new.original.digest(work/'plan.json'),native_logposterior_invocations=1,
             auxiliary_background_calls=1,accuracy2_reference_sha256=plan['points'][index]['accuracy2_record_sha256'],
             replay_gate_sha256=new.require_replay(work,plan) if stage=='refine' else None,
             spectrum_path=new.original.relative(spec),spectrum_sha256=new.original.digest(spec),
             status=('accuracy2_replay_passed' if passed else 'accuracy2_replay_failed') if stage=='replay' else 'finite_native_precision',
             failed_checks=[] if passed else ['replay'],comparison={'high_minus_declared_total_loglike':.01*(index%2),
             'high_minus_declared_component_loglikes':{'cmb':.01*(index%2)}})
    new.write_once(path,row)


def cache_contracts():
    with tempfile.TemporaryDirectory(dir=ROOT/'.work') as folder:
        d=Path(folder);screen,review,receipt,two,three,spectra=make_fixture(d)
        with patch.object(new.consumer,'verify',return_value=receipt),patch.object(new,'validation_guard',return_value='synthetic-validation'),\
             patch.object(new,'configurations',return_value=(two,three)),\
             patch.object(new.subprocess,'run',side_effect=AssertionError('No subprocess/native work permitted')):
            work=d/'new';plan=new.prepare(screen,review,work)
            assert new.verify_plan(work/'plan.json')==plan
            assert [p['point'] for p in plan['points']]==[p['point'] for p in new.original.sealed_read(d/'old-plan.json')['selected']]
            refuses(lambda:new.run_stage(work,'refine',True,1))
            assert not list(work.glob('refine*'))
            cached_point(work,'replay',0,plan,spectra)
            replay=new.run_stage(work,'replay',True,1)
            assert replay['status']=='accuracy2_replay_passed'
            new.require_replay(work,plan)
            for i in range(32):cached_point(work,'refine',i,plan,spectra)
            result=new.run_stage(work,'refine',True,2)
            assert result['status']=='no_large_variation_detected_on_fixed32'
            assert result['known_native_logposterior_invocations']==32 and result['cumulative_known_native_logposterior_invocations']==33
            assert result['unknown_count_attempts']==0 and result['posterior_qualification'] is False
            assert new.run_stage(work,'refine',False,1)==result
            # Missing/failed replay blocks the next stage before any attempt.
            failed=d/'failed';p2=new.prepare(screen,review,failed)
            cached_point(failed,'replay',0,p2,spectra,False)
            new.run_stage(failed,'replay',True,1)
            refuses(lambda:new.run_stage(failed,'refine',True,1))
            assert not list(failed.glob('refine*'))
            # Replay ledger cannot be omitted or silently replaced.
            ledger=work/'replay-record-hashes.json';backup=ledger.read_bytes();ledger.unlink()
            refuses(lambda:new.require_replay(work,plan));ledger.write_bytes(backup)
            # A killed attempt consumes the slot; direct worker retry is refused.
            path,log=new.file_paths(work,'replay',0);record=path.read_bytes();path.unlink()
            stopped=new.collect_or_execute(work,'replay',0,plan,True)
            assert stopped['unknown_native_invocation_count'] and stopped['status']=='incomplete_attempt_no_retry'
            refuses(lambda:new.worker(work/'plan.json','replay',0));path.write_bytes(record)
            # Source-map deletion, altered parent and cached numeric tampering.
            planpath=work/'plan.json';backup=planpath.read_bytes();bad=json.loads(backup)
            bad['source_sha256'].pop(next(iter(bad['source_sha256'])))
            bad.pop('payload_sha256');bad['identity']=new.original.identity({k:v for k,v in bad.items() if k!='identity'})
            planpath.unlink();new.write_once(planpath,bad);refuses(lambda:new.verify_plan(planpath));planpath.write_bytes(backup)
            original=screen.read_bytes();screen.write_bytes(original+b' ');refuses(lambda:new.verify_plan(planpath));screen.write_bytes(original)
            point=work/'refine-00.json';backup=point.read_bytes();bad=json.loads(backup);bad['comparison']['high_minus_declared_total_loglike']+=1
            point.write_text(json.dumps(bad));refuses(lambda:new.run_stage(work,'refine',False,1));point.write_bytes(backup)
            for status in ('no_large_variation_detected_on_fixed32','incomplete_or_failed_numerical_screen'):
                wrong=dict(receipt,reviewed_screen_status=status)
                with patch.object(new.consumer,'verify',return_value=wrong):
                    refuses(lambda:new.prepare(screen,review,d/'refused'))
                    assert not (d/'refused').exists()
    return {'full_synthetic_plan_and_33_record_cache_replayed':True,'original32_IDs_and_coordinates_preserved':True,
            'no_qualifying_posterior_claim':True,'missing_or_failed_replay_blocks_refinement':True,
            'interrupted_attempt_has_unknown_count_and_no_retry':True,'duplicate_worker_claim_refused':True,
            'missing_ledger_source_map_parent_numeric_tampering_refused':True,
            'parent_receipt_and_configurations_explicitly_mocked':True}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=new.VALIDATION)
    a=p.parse_args();assert not a.output.exists()
    import camb,cobaya.model
    def forbidden(*args,**kwargs):raise AssertionError('Physical calculations forbidden during refinement validation')
    with ExitStack() as guard:
        for owner,name in [(camb,'get_results'),(camb,'get_background'),(camb,'get_transfer_functions'),(cobaya.model,'get_model'),(cobaya.model.Model,'__init__')]:
            guard.enter_context(patch.object(owner,name,forbidden))
        start=new.dependencies()
        config_checks=[]
        for model in ('lcdm','cpl'):
            settings=dict(model=model,evolution='none',sample='dovekie',calibration='official_planck',fast_lensing=True,gpu=False)
            two,three=new.configurations(settings,json.loads(new.original.DESIGN.read_text()),json.loads(new.DESIGN.read_text()))
            assert set(two['theory'])==set(three['theory'])=={'camb'}
            config_checks.append(dict(model=model,only_three_boost_controls_changed=True,no_model_construction=True))
        result=dict(status='passed_synthetic_native_refinement_validation',physical_calls=0,
                    arithmetic=arithmetic(),cache=cache_contracts(),configuration_checks=config_checks,
                    source_sha256=start,scope='Synthetic payloads and explicit parent/configuration mocks for cache tests; actual factory dictionaries for configuration invariance. No physical model/background/spectrum evaluation.')
        assert new.dependencies()==start
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'physical_calls':0}))


if __name__=='__main__':main()
