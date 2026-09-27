"""Synthetic lineage/status tests of the supplemental precision consumer."""
import argparse
from contextlib import ExitStack
import copy
import importlib.metadata
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp

import native_precision_review_consumer as consumer
from exact_correction import record_digest
from native_precision_thermal_review_validate import fixture as numerical_fixture

ROOT=consumer.ROOT; HERE=consumer.HERE
original=consumer.original; thermal=consumer.thermal


def write(path, value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def fixture(directory, mode='auxiliary_only'):
    """Real files/seals with synthetic densities and an explicit fake qualifier."""
    parentdir=directory/'parent'; exact=parentdir/'exact'; exact.mkdir(parents=True)
    cache=directory/'screen'; cache.mkdir()
    design=consumer.read(original.DESIGN)
    frozen={'identity':'explicit-synthetic-target','source_sha256':{},'versions':{'numpy':importlib.metadata.version('numpy')}}
    settings={'synthetic_fixture':True}
    run=parentdir/'run-0.json'; write(run,{'target_identity':frozen,'arguments':settings})
    points=[]; groups=[]; locations=[]; stored=[]
    for i in range(2000):
        _, row, _=numerical_fixture(design,i)
        row.update(status='finite',log_weight=0.)
        row['payload_sha256']=record_digest(row)
        path=exact/f'{i:05d}.json'; write(path,row)
        points.append(row['point']); groups.append(i//500)
        locations.append({'expanded_index':1000+i%500,'row':i%500,'chain':f'chain.{i//500+1}.txt'})
        stored.append(row)
    selection=exact/'selection.json';write(selection,{'points':points,'groups':groups,'locations':locations,'settings':settings})
    correction=directory/'correction.json';write(correction,{'selection_path':original.relative(selection)})
    paths=[run,selection,correction]+[exact/f'{i:05d}.json'for i in range(2000)]
    parent={'qualified_under_declared_numerical_gates':True,'target_identity':frozen['identity'],
            'input_sha256':{original.relative(p):original.digest(p)for p in paths}}
    selected=original.select_indices(groups,locations,8)
    weights=np.exp(np.zeros(2000)-logsumexp(np.zeros(2000)))
    for item in selected:
        i=item['parent_index'];p=exact/f'{i:05d}.json'
        item.update(point=points[i],native_record_path=original.relative(p),native_record_sha256=original.digest(p),
                    original_location=locations[i],original_native_logweight=0.,weight_in_full_qualified_parent=float(weights[i]))
    source_paths=[Path(original.__file__),original.DESIGN,HERE/'measurement_summary.py',HERE/'exact_correction.py',
                  HERE/'native_precision_audit.py',HERE/'likelihood.py',HERE/'target_identity.py']
    plan={'design':design,'settings':settings,'selected':selected,'frozen_target':frozen,
          'qualified_parent_inputs':parent['input_sha256'],
          'source_sha256':{original.relative(p):original.digest(p)for p in source_paths},
          'nominal_native_configuration':{'synthetic':'nominal'},'doubled_native_configuration':{'synthetic':'high'},
          'weight_context_only':{'selected_weight_mass_in_full_parent':float(sum(x['weight_in_full_qualified_parent']for x in selected)),
              'full_parent_chain_weight_masses':{str(c):float(weights[np.asarray(groups)==c].sum())for c in range(4)},
              'use':'Context only. Selection and numerical-screen centering are unweighted; no posterior reweighting.'},
          'correction_summary_path':original.relative(correction),'correction_summary_sha256':original.digest(correction)}
    plan['identity']=original.identity(plan)
    planpath=cache/'selection.json';original.sealed_write(planpath,plan)
    rows=[]; manifests={}; values_by_point={}
    for i,item in enumerate(selected):
        j=item['parent_index'];row,_,values=numerical_fixture(design,j)
        row.update(audit_index=i,parent_index=j,chain=i//8,identity=plan['identity'],
                   point=points[j],plan_sha256=original.digest(planpath),native_record_sha256=item['native_record_sha256'],
                   original_native_logweight=0.,weight_in_full_qualified_parent=item['weight_in_full_qualified_parent'],
                   stored_declared_loglikes=stored[j]['exact_loglikes'],stored_declared_logpost=stored[j]['exact_logpost'])
        if mode=='real_density_variation':
            delta=float(np.linspace(-.4,.4,32)[i]);row['high_loglikes']['cmb']+=delta;row['high_logpost']+=delta
        if mode=='prior_failure' and i==31:
            row['high_logpriors'][0]+=.001;row['high_logpost']+=.001;row['failed_checks'].append('prior_sum_changed')
        if mode=='unknown_failure' and i==31:row['failed_checks'].append('unrecognized_existing_failure')
        if mode=='high_background_failure' and i==31:values['high']['H_km_s_Mpc'][2]*=1.00001
        row['comparison']=original.density_comparison(stored[j],row['high_loglikes'],row['high_logpriors'],row['high_logpost'],design['diagnostic_thresholds'])
        path=cache/f'{i:02d}.json';spectrum=cache/f'{i:02d}-spectra.npz';log=cache/f'{i:02d}.log'
        np.savez(spectrum,synthetic=np.arange(4));log.write_text('Synthetic fixture; no native call.\n')
        row.update(spectrum_path=original.relative(spectrum),spectrum_sha256=original.digest(spectrum))
        original.sealed_write(path,row);row=original.sealed_read(path)
        row.update(record_path=original.relative(path),record_sha256=original.digest(path),
                   log_path=original.relative(log),log_sha256=original.digest(log),integrated_time_mismatch_warning_count=0)
        for kind in ['record','log','spectrum']:manifests[row[kind+'_path']]=row[kind+'_sha256']
        rows.append(row);values_by_point[row['point']['H0']]=values
    manifest=cache/'record-hashes.json';write(manifest,manifests)
    screen=original.summarize(rows,design)
    screen.update(points=rows,plan_path=original.relative(planpath),plan_sha256=original.digest(planpath),
        record_manifest_path=original.relative(manifest),record_manifest_sha256=original.digest(manifest),
        source_sha256=plan['source_sha256'],frozen_target_identity=frozen['identity'],declared_target_changed=False,
        known_native_point_evaluations=32,incomplete_attempts_with_unknown_native_call_count=0,
        weight_context_only=plan['weight_context_only'],numerical_controls=design['numerical_controls'])
    screenpath=directory/'screen.json';write(screenpath,screen)
    reviewpath=directory/'review.json'
    with patch.object(thermal,'summarize_run',return_value=parent),patch.object(thermal,'backgrounds',side_effect=lambda point,*args:values_by_point[point['H0']]):
        thermal.actual(screenpath,reviewpath)
    return screenpath,reviewpath,parent


def validate():
    import camb
    import cobaya.model
    forbidden_calls=[]
    def forbidden(*args,**kwargs):forbidden_calls.append(True);raise AssertionError('No physical call.')
    tests=[];mutations=[];requalification_calls=[]
    with ExitStack() as stack:
        for owner,name in [(camb,'get_background'),(camb,'get_results'),(camb,'get_transfer_functions'),
                            (cobaya.model,'get_model'),(cobaya.model.Model,'__init__')]:
            stack.enter_context(patch.object(owner,name,forbidden))
        stack.enter_context(patch.object(original,'native_configurations',return_value=({'synthetic':'nominal'},{'synthetic':'high'})))
        for mode in ['auxiliary_only','real_density_variation','prior_failure','unknown_failure','high_background_failure']:
            with tempfile.TemporaryDirectory(dir=ROOT/'.work') as temporary:
                folder=Path(temporary);screenpath,reviewpath,parent=fixture(folder,mode)
                def qualifier(*args):requalification_calls.append(args);return parent
                with patch.object(consumer,'summarize_run',side_effect=qualifier):
                    result=consumer.verify(screenpath,reviewpath)
                    assert result['original_screen_status']=='incomplete_or_failed_numerical_screen'
                    expected={'auxiliary_only':'no_large_variation_detected_on_fixed32','real_density_variation':'numerical_sensitivity_requires_followup'}.get(mode,'incomplete_or_failed_numerical_screen')
                    assert result['reviewed_screen_status']==expected
                    assert result['reviewed_precision_screen_supported']==(mode=='auxiliary_only')
                    assert result['posterior_qualification'] is False and result['posterior_reweighting_performed'] is False
                    tests.append({'mode':mode,'original_status':result['original_screen_status'],'reviewed_status':expected,
                                  'supported':result['reviewed_precision_screen_supported'],'flags':result['reviewed_diagnostic_flags'],
                                  'failures':result['reviewed_failures']})
                    if mode!='auxiliary_only':continue
                    originals={p:p.read_bytes()for p in folder.rglob('*')if p.is_file()}
                    def restore():
                        for path,body in originals.items():path.write_bytes(body)
                    def reject(label,change):
                        restore();change()
                        try:consumer.verify(screenpath,reviewpath)
                        except (AssertionError,KeyError,ValueError):mutations.append(label)
                        else:raise AssertionError('Accepted mutation '+label)
                    def alter_review(change,allow_nan=False):
                        obj=json.loads(originals[reviewpath]);change(obj)
                        reviewpath.write_text(json.dumps(obj,allow_nan=allow_nan))
                    reject('changed_review_status',lambda:alter_review(lambda x:x['reviewed_screen'].update(status='numerical_sensitivity_requires_followup')))
                    reject('erased_original_failure',lambda:alter_review(lambda x:x['point_reviews'][0].update(original_failed_checks=[])))
                    reject('changed_saved_background',lambda:alter_review(lambda x:x['point_reviews'][31]['backgrounds']['nominal'].update(rdrag=146.)))
                    reject('NaN_nonfirst_background',lambda:alter_review(lambda x:x['point_reviews'][31]['backgrounds']['nominal']['expansion'].update(q05=float('nan')),True))
                    reject('missing_source_binding',lambda:alter_review(lambda x:x['input_source_sha256'].pop(next(iter(x['input_source_sha256'])))))
                    reject('changed_parent_bytes',lambda:(folder/'parent/run-0.json').write_bytes(originals[folder/'parent/run-0.json']+b' '))
                    reject('changed_native_point_bytes',lambda:(folder/'screen/31.json').write_bytes(originals[folder/'screen/31.json']+b' '))
                    reject('changed_log_bytes',lambda:(folder/'screen/31.log').write_text('Changed log'))
                    reject('changed_spectrum_bytes',lambda:(folder/'screen/31-spectra.npz').write_bytes(b'changed'))
                    def altered_screen(kind):
                        obj=json.loads(originals[screenpath])
                        if kind=='point':obj['points'][31]['point']['H0']+=1
                        else:obj['points'][31]['status']='exception'
                        write(screenpath,obj)
                        r=json.loads(originals[reviewpath]);r['original_screen_sha256']=original.digest(screenpath)
                        r['input_source_sha256'][original.relative(screenpath)]=original.digest(screenpath);write(reviewpath,r)
                    reject('changed_point_even_with_screen_hash_refreshed',lambda:altered_screen('point'))
                    reject('incomplete_native_point_not_repaired',lambda:altered_screen('exception'))
                    restore()
                    for label,changed in [('unqualified_parent',dict(parent,qualified_under_declared_numerical_gates=False)),
                                          ('changed_parent_target',dict(parent,target_identity='different')),
                                          ('changed_parent_input_set',dict(parent,input_sha256={}))]:
                        with patch.object(consumer,'summarize_run',return_value=changed):
                            try:consumer.verify(screenpath,reviewpath)
                            except AssertionError:mutations.append(label)
                            else:raise AssertionError('Accepted '+label)
                    restore()
                    assert consumer.verify(screenpath,reviewpath)==result
                    assert all(path.read_bytes()==body for path,body in originals.items()),'Consumer wrote input evidence.'
    assert not forbidden_calls
    paths=[Path(__file__),Path(consumer.__file__),consumer.DESIGN,Path(thermal.__file__),thermal.DESIGN,
           Path(original.__file__),original.DESIGN,HERE/'native_precision_thermal_review_validate.py',
           HERE/'measurement_summary.py',HERE/'exact_correction.py',HERE/'target_identity.py']
    validation=consumer.read(consumer.THERMAL_VALIDATION)
    consumer.hashes(validation['source_sha256']);consumer.hashes(validation['input_sha256'])
    return {'status':'passed_synthetic_supplemental_precision_consumer_validation','genuine_status_cases':tests,
            'refused_mutations':mutations,'fresh_mock_qualifier_calls':len(requalification_calls),
            'real_background_native_or_model_calls':0,'input_bytes_preserved':True,
            'source_sha256':{original.relative(p):original.digest(p)for p in paths},
            'existing_thermal_validation_sha256':{original.relative(consumer.THERMAL_VALIDATION):original.digest(consumer.THERMAL_VALIDATION)},
            'posterior_qualification':False,
            'mock_scope':'Fresh parent qualifier and native configuration construction are explicitly synthetic dependencies;32saved backgrounds are synthetic. Actual files, seals, source hashes, selection, precision arithmetic and status consumers execute normally. No actual posterior screen was evaluated.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'studies/unified_cosmology/results/inference/native-precision-review-consumer-validation.json')
    args=p.parse_args();value=validate();write(args.output,value)
    print(json.dumps({k:value[k]for k in ['status','refused_mutations','real_background_native_or_model_calls']}))
