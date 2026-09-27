"""Independent finite-control and failure-integrity checks; backgrounds only."""
import argparse
from contextlib import ExitStack
import copy
import importlib.metadata
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

import numpy as np

import native_precision_thermal_review as review

ROOT=review.ROOT
HERE=review.HERE
OUT=ROOT/'studies/unified_cosmology/results/inference/native-precision-thermal-review-validation.json'


def reject(call, contains=None):
    try:
        call()
    except (AssertionError,KeyError,ValueError) as error:
        if contains:assert contains in str(error),str(error)
        return
    raise AssertionError('Invalid thermal-review fixture accepted.')


def fixture(design, index=0):
    point={'H0':68.+index/100.}
    expansion={'q0':-.55,'q05':-.1,'q1':.15,'j0':1.,'j05':1.,'j1':1.}
    low={'rdrag':147.,'omegam':.3,'expansion':expansion,'H_km_s_Mpc':[68.,90.,120.,235.,1.58e6],
         'input_WantTransfer':False,'background_WantTransfer':True}
    values={'nominal':low,'high':copy.deepcopy(low)}
    stored={'point':point,'derived':dict(expansion,rdrag=147.,omegam=.3),
            'exact_loglikes':{'cmb':-100.,'bao':-7.,'sn':230.},'exact_logpost':126.}
    row={'audit_index':index,'chain':index//8,'point':point,'native_point_evaluations':1,
         'status':'failed_native_precision_checks','failed_checks':['nominal_background_rdrag_relative'],
         'high_loglikes':dict(stored['exact_loglikes']),'high_logpriors':[3.],'high_logpost':126.,
         'high_derived':dict(stored['derived']),'background':{'high_H_km_s_Mpc':low['H_km_s_Mpc']},
         'nominal_background_stored_native_closure':{'rdrag_relative':2e-6,'omegam_absolute':0.,'q_scaled':0.,'j_scaled':0.},
         'finalized_theory_extra_args':{},'comparison':{}}
    row['comparison']=review.original.density_comparison(stored,row['high_loglikes'],row['high_logpriors'],row['high_logpost'],design['diagnostic_thresholds'])
    return row,stored,values


def synthetic_checks(design, review_design):
    row,stored,values=fixture(design)
    initial=copy.deepcopy((row,stored,values))
    fixed,report=review.reviewed_row(row,stored,values,design,review_design)
    assert fixed['status']=='finite_native_precision' and fixed['failed_checks']==[]
    assert (row,stored,values)==initial
    assert fixed['comparison']==row['comparison']
    preserved=[]
    for label,change_prior,change_post in [('prior_sum_changed',.001,.001),('logdensity_accounting',0.,.001)]:
        bad=copy.deepcopy(row);bad['high_logpriors'][0]+=change_prior;bad['high_logpost']+=change_post
        bad['comparison']=review.original.density_comparison(stored,bad['high_loglikes'],bad['high_logpriors'],bad['high_logpost'],design['diagnostic_thresholds'])
        # Even a missing old flag is restored from unchanged density arithmetic.
        output,_=review.reviewed_row(bad,stored,values,design,review_design)
        assert label in output['failed_checks'] and output['status']=='failed_native_precision_checks'
        preserved.append(label)
    bad=copy.deepcopy(row);bad['failed_checks'].append('unrecognized_existing_failure')
    output,_=review.reviewed_row(bad,stored,values,design,review_design)
    assert output['failed_checks']==['unrecognized_existing_failure']
    for status in ['exception','nonfinite','incomplete_attempt_no_silent_retry','worker_failed_before_completed_record']:
        bad=copy.deepcopy(row);bad['status']=status
        reject(lambda:review.reviewed_row(bad,stored,values,design,review_design),'Incomplete')
    for which,key in [('nominal','q05'),('nominal','j1'),('high','q1'),('high','j05')]:
        bad=copy.deepcopy(values);bad[which]['expansion'][key]=float('nan')
        reject(lambda:review.reviewed_row(row,stored,bad,design,review_design),'Nonfinite')
    bad=copy.deepcopy(stored);bad['derived']['q05']=float('nan')
    reject(lambda:review.reviewed_row(row,bad,values,design,review_design),'Nonfinite')
    bad=copy.deepcopy(row);bad['high_loglikes']['cmb']=float('nan')
    reject(lambda:review.reviewed_row(bad,stored,values,design,review_design))
    for kind,threshold in [('H_relative','high_background_H_relative_max'),('rdrag_relative','high_background_rdrag_relative_max')]:
        bad=copy.deepcopy(values)
        if kind=='H_relative':bad['high']['H_km_s_Mpc'][2]*=1+review_design[threshold]*10
        else:bad['high']['rdrag']*=1+review_design[threshold]*10
        output,_=review.reviewed_row(row,stored,bad,design,review_design)
        assert 'high_background_'+kind in output['failed_checks']
    bad=copy.deepcopy(row);bad['comparison']['high_minus_declared_total_loglike']+=.01
    reject(lambda:review.reviewed_row(bad,stored,values,design,review_design),'arithmetic')
    # A successful thermal-path repair must not erase genuine precision spread.
    rows=[]
    delta=np.linspace(-.4,.4,32)
    for i,d in enumerate(delta):
        r,s,v=fixture(design,i);r['high_loglikes']['cmb']+=float(d);r['high_logpost']+=float(d)
        r['comparison']=review.original.density_comparison(s,r['high_loglikes'],r['high_logpriors'],r['high_logpost'],design['diagnostic_thresholds'])
        output,_=review.reviewed_row(r,s,v,design,review_design);rows.append(output)
    summary=review.original.summarize(rows,design)
    assert summary['status']=='numerical_sensitivity_requires_followup'
    assert {'total_centered_RMS','total_centered_max_absolute','component_variation:cmb'}<=set(summary['diagnostic_flags'])
    assert abs(summary['total_loglike_difference']['centered_RMS']-np.std(delta))<1e-13
    assert summary['posterior_qualification'] is False and summary['posterior_reweighting_performed'] is False
    return {'only_background_flags_repaired':True,'prior_and_accounting_failures_preserved':preserved,
            'unknown_failure_preserved':True,'incomplete_and_nonfinite_rejected':True,
            'high_H_and_rdrag_mismatches_rejected':True,'saved_density_tampering_rejected':True,
            'genuine_centered_density_variation_preserved':summary['diagnostic_flags']}


def consumer_checks(design):
    with tempfile.TemporaryDirectory(dir=ROOT/'.work') as directory:
        directory=Path(directory);folder=directory/'parent';exact=folder/'exact';exact.mkdir(parents=True)
        cache=directory/'screen';cache.mkdir()
        selection=exact/'selection.json';selection.write_text('{}')
        correction=directory/'correction.json';correction.write_text(json.dumps({'selection_path':review.original.relative(selection)}))
        parent={'input_sha256':{review.original.relative(p):review.original.digest(p) for p in [selection,correction]}}
        selected=[];source_rows=[];fake_backgrounds={}
        for i in range(32):
            row,stored,values=fixture(design,i)
            p=exact/f'{i:05d}.json';p.write_text(json.dumps(stored))
            selected.append({'point':row['point'],'native_record_path':review.original.relative(p),'native_record_sha256':review.original.digest(p)})
            row['native_record_sha256']=review.original.digest(p)
            source_rows.append(row);fake_backgrounds[row['point']['H0']]=values
        plan={'design':design,'selected':selected,'correction_summary_path':review.original.relative(correction),
              'correction_summary_sha256':review.original.digest(correction),
              'source_sha256':{review.original.relative(Path(review.__file__)):review.original.digest(review.__file__)},
              'qualified_parent_inputs':parent['input_sha256'],
              'frozen_target':{'source_sha256':{},'versions':{'numpy':importlib.metadata.version('numpy')}}}
        planpath=cache/'selection.json';review.original.sealed_write(planpath,plan)
        manifest={};augmented=[]
        for i,row in enumerate(source_rows):
            row['plan_sha256']=review.original.digest(planpath)
            p=cache/f'{i:02d}.json';review.original.sealed_write(p,row)
            record=review.original.sealed_read(p)
            record.update(record_path=review.original.relative(p),record_sha256=review.original.digest(p))
            augmented.append(record);manifest[review.original.relative(p)]=review.original.digest(p)
        ledger=cache/'record-hashes.json';ledger.write_text(json.dumps(manifest))
        screen=review.original.summarize(augmented,design)
        screen.update(points=augmented,plan_path=review.original.relative(planpath),plan_sha256=review.original.digest(planpath),
                      record_manifest_path=review.original.relative(ledger),record_manifest_sha256=review.original.digest(ledger))
        screenpath=directory/'screen.json';screenpath.write_text(json.dumps(screen))
        output=directory/'review.json'
        with patch.object(review,'summarize_run',side_effect=ValueError('unqualified-parent')), \
             patch.object(review,'backgrounds',side_effect=AssertionError('called-before-parent-qualification')):
            reject(lambda:review.actual(screenpath,output),'unqualified-parent')
        assert not output.exists()
        calls=[]
        def background(point,*args):calls.append(point);return fake_backgrounds[point['H0']]
        with patch.object(review,'summarize_run',return_value=parent),patch.object(review,'backgrounds',side_effect=background):
            result=review.actual(screenpath,output)
            assert len(calls)==32 and result['background_calls']==64 and result['CMB_spectrum_calls']==0
            assert result['reviewed_screen']['status']=='no_large_variation_detected_on_fixed32'
            assert result['original_screen_status']=='incomplete_or_failed_numerical_screen'
            assert result['posterior_qualification'] is False
            reject(lambda:review.actual(screenpath,output),'Preserve previous')
            assert len(calls)==32
            before=screenpath.read_bytes();bad=copy.deepcopy(screen);bad['points'][-1]['status']='exception'
            screenpath.write_text(json.dumps(bad))
            reject(lambda:review.actual(screenpath,directory/'incomplete.json'),'Incomplete')
            assert len(calls)==32;screenpath.write_bytes(before)
            before=ledger.read_bytes();ledger.write_bytes(before+b' ')
            reject(lambda:review.actual(screenpath,directory/'ledger.json'));assert len(calls)==32;ledger.write_bytes(before)
            p=cache/'31.json';before=p.read_bytes();p.write_bytes(before+b' ')
            reject(lambda:review.actual(screenpath,directory/'changed-record.json'));assert len(calls)==32;p.write_bytes(before)
            before=correction.read_bytes();correction.write_bytes(before+b' ')
            reject(lambda:review.actual(screenpath,directory/'changed-parent.json'));assert len(calls)==32;correction.write_bytes(before)
            assert screenpath.read_bytes()==json.dumps(screen).encode()
        return {'fake_consumer_records':32,'qualification_and_background_only_explicitly_mocked':True,
                'real_file_hashes_and_seals_checked':True,'unqualified_parent_refused_before_backgrounds':True,
                'late_incomplete_record_rejected_before_backgrounds':True,'changed_records_manifest_parent_rejected':True,
                'existing_review_not_overwritten':True,'original_screen_preserved':True}


def physical_checks(design):
    import camb
    from modern_fast import configuration
    from exact_correction import verify_record
    input_hashes={};checks=[]
    pinpath=ROOT/'studies/unified_cosmology/results/inference/expansion-history-validation.json'
    pins=json.loads(pinpath.read_text())['full_native_thermal_reference_checks']
    input_hashes[review.original.relative(pinpath)]=review.original.digest(pinpath)
    extra=configuration(model='lcdm')['theory']['camb']['extra_args']
    for pin in pins:
        index=pin['index'];path=ROOT/f'.work/unified-cosmology/inference/spectral-capture-validation-v2/{index:05d}.json'
        assert review.original.digest(path)==pin['native_record_sha256']
        record=json.loads(path.read_text());verify_record(record)
        input_hashes[review.original.relative(path)]=review.original.digest(path)
        point=record['point'];physics={k:point[k] for k in ['H0','ombh2','omch2','ns','tau']}
        physics.update(As=1e-10*np.exp(point['logA']),w=point.get('w',-1.),wa=point.get('wa',0.))
        plain=camb.get_background(camb.set_params(**physics,**extra))
        plain_error=abs(plain.get_derived_params()['rdrag']/record['derived']['rdrag']-1)
        tested=review.backgrounds(point,extra,design)
        low=tested['nominal'];reference=record['derived']
        errors={'rdrag_relative':abs(low['rdrag']/reference['rdrag']-1),
                'omegam_absolute':abs(low['omegam']-reference['omegam']),
                'q_scaled':float(np.max([abs(low['expansion'][k]-reference[k])/(1+abs(reference[k])) for k in ['q0','q05','q1']])),
                'j_scaled':float(np.max([abs(low['expansion'][k]-reference[k])/(1+abs(reference[k])) for k in ['j0','j05','j1']]))}
        assert all(error<=design['diagnostic_thresholds']['nominal_background_'+key+'_closure_max'] for key,error in errors.items())
        checks.append({'index':index,'plain_background_rdrag_relative_error':float(plain_error),
                       'plain_background_fails_original_rdrag_gate':plain_error>design['diagnostic_thresholds']['nominal_background_rdrag_relative_closure_max'],
                       'thermal_background_errors':errors,'source_sha256':review.original.digest(path)})
    assert len(checks)==4 and any(r['plain_background_fails_original_rdrag_gate'] for r in checks)

    # Independent earlier fixed-reference/holdout native pairs supply actual
    # high-accuracy H(z) and drag-scale references, not reconstructed backgrounds.
    auditpath=ROOT/'studies/unified_cosmology/results/inference/native-precision-audit.json'
    audit=json.loads(auditpath.read_text());input_hashes[review.original.relative(auditpath)]=review.original.digest(auditpath)
    records=[]
    for pin in audit['points']:
        path=ROOT/f".work/unified-cosmology/inference/native-precision/{pin['index']:02d}.json"
        assert review.original.digest(path)==pin['record_sha256']
        r=json.loads(path.read_text());assert r['status']=='finite_native'
        assert review.original.digest(ROOT/r['spectrum_file'])==r['spectrum_sha256']
        input_hashes[review.original.relative(path)]=review.original.digest(path)
        input_hashes[r['spectrum_file']]=r['spectrum_sha256'];records.append(r)
    paired=[]
    for nominal,high in zip(records[::2],records[1::2]):
        assert nominal['point']==high['point'] and nominal['label']==high['label']
        values=review.backgrounds(nominal['point'],high['finalized_theory_extra_args'],design)
        errors={}
        for name,reference in [('nominal',nominal),('high',high)]:
            errors[name]={'H_relative':float(np.max(abs(np.asarray(values[name]['H_km_s_Mpc'])/reference['H_km_s_Mpc']-1))),
                          'rdrag_relative':abs(values[name]['rdrag']/reference['derived']['rdrag']-1)}
            assert errors[name]['H_relative']<=1e-10 and errors[name]['rdrag_relative']<=1e-7
        paired.append({'label':nominal['label'],'errors':errors,
                       'saved_native_total_loglike_difference':float(sum(high['loglikes'].values())-sum(nominal['loglikes'].values()))})
    assert len(paired)==4
    return checks,paired,input_hashes


def validate():
    import camb
    import cobaya.model
    design=json.loads(review.original.DESIGN.read_text());review_design=json.loads(review.DESIGN.read_text())
    source_paths=[Path(__file__),Path(review.__file__),review.DESIGN,review.original.DESIGN,
                  HERE/'native_posterior_precision.py',HERE/'expansion_history.py',HERE/'likelihood.py',
                  HERE/'measurement_summary.py',HERE/'exact_correction.py',HERE/'modern_fast.py',
                  HERE/'modern_run.py',HERE.parent/'external_probes/modern_adapter.py']
    sources={review.original.relative(p):review.original.digest(p) for p in source_paths}
    forbidden_calls=[];background_calls=[];original_background=camb.get_background
    def forbidden(*args,**kwargs):forbidden_calls.append(True);raise AssertionError('Native spectra/models forbidden in validation.')
    def background(*args,**kwargs):background_calls.append(True);return original_background(*args,**kwargs)
    with ExitStack() as stack:
        for name in ['get_results','get_transfer_functions']:stack.enter_context(patch.object(camb,name,forbidden))
        stack.enter_context(patch.object(cobaya.model,'get_model',forbidden))
        stack.enter_context(patch.object(cobaya.model.Model,'__init__',forbidden))
        stack.enter_context(patch.object(camb,'get_background',background))
        synthetic=synthetic_checks(design,review_design)
        consumer=consumer_checks(design)
        controls,paired,inputs=physical_checks(design)
    assert not forbidden_calls and len(background_calls)==20
    assert all(review.original.digest(ROOT/p)==h for p,h in dict(sources,**inputs).items())
    return {'status':'passed_background_only_thermal_review_validation','source_sha256':sources,'input_sha256':inputs,
            'four_pinned_full_native_controls':controls,'four_saved_nominal_high_native_pairs':paired,
            'synthetic_failure_tests':synthetic,'synthetic_consumer_contract':consumer,
            'CMB_spectrum_calls':0,'actual_background_calls':len(background_calls),'Cobaya_model_constructions':0,
            'posterior_points_evaluated':0,'posterior_qualification':False,
            'versions':{p:importlib.metadata.version(p) for p in ['camb','numpy','scipy','cobaya']},
            'scope':'Fixed previously evaluated physical controls and synthetic failure/consumer checks only. No32-point posterior screen, likelihood reevaluation, numerical reweighting or posterior measurement was performed.'}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,default=OUT)
    args=p.parse_args();result=validate();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'background_calls':result['actual_background_calls'],'CMB_spectrum_calls':0}))


if __name__=='__main__':main()
