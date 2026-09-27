"""Synthetic receipt and lineage controls for the supplementary auditor."""
import argparse
from contextlib import ExitStack
from copy import deepcopy
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

import final_audit as base
from final_audit_validate import fixture as original_fixture, put
import supplementary_audit as auditor

OUTPUT = base.ROOT/'studies/unified_cosmology/results/inference/supplementary-audit-validation.json'


def fixture(root):
    spec, parent, records, progress, status_path, rewrite = original_fixture(root)
    progress['completed_results'] = progress['completed_results'][:3]
    progress.update(stage='stopped_no_retry', error="CalledProcessError: quantile_precision.py rejected chronological ties")
    put(status_path, progress)
    producer = root/'studies/unified_cosmology/code/inference'
    continuation = producer/'postprocessing_continuation.py'
    continuation.write_text('# Synthetic continuation writer, never executed.\n')
    ordered = producer/'quantile_precision_ordered.py'
    ordered.write_text('# Synthetic ordered quantile producer, never executed.\n')
    entry = {'original_spec': spec, 'original_status_path': str(status_path.relative_to(root)),
             'original_status_sha256': base.digest(status_path), 'receipts': [], 'precision_review': None}
    input_binding = {spec['plan_path']: spec['plan_sha256'], entry['original_status_path']: entry['original_status_sha256']}
    planned = base.read(root/spec['plan_path'])
    for suffix, names in [('derived', ['quantile_precision.py', 'expansion_history.py', 'luminosity_history.py', 'probe_omission.py:bao']),
                          ('precision', ['native_posterior_precision.py'])]:
        directory = root/'.work/supplemental'/suffix; directory.mkdir(parents=True)
        stages = []
        for requested in names:
            name, _, omission = requested.partition(':')
            row = deepcopy(next(s for s in planned['stages'] if s['script'] == name and (not omission or s['argv'][-1] == omission)))
            if name == 'quantile_precision.py':
                row.update(preserved_original_script=name, preserved_original_source_sha256=row['source_sha256'],
                           script=ordered.name, source_sha256=base.digest(ordered))
                row['argv'] = [str(ordered) if x == str(producer/name) else x for x in row['argv']]
                value = base.read(Path(row['output']))
                value['source_sha256'][str(ordered.relative_to(root))] = base.digest(ordered)
                put(Path(row['output']), value)
            stages.append(row)
        plan = {'schema': 'explicit-postprocessing-continuation-v1',
            'original_failed_plan_and_status_sha256': input_binding,
            'source_sha256': base.digest(continuation), 'folder': planned['folder'],
            'target_identity': spec['target_identity'], 'stages': stages}
        pp = directory/'plan.json'; put(pp, plan)
        entry['receipts'].append({'plan_path': str(pp.relative_to(root)), 'plan_sha256': base.digest(pp)})
        results = []
        for index, stage in enumerate(stages):
            put(directory/f'{index:02d}-attempt.json', {'script': stage['script'], 'argv': stage['argv'], 'started_utc': '2026-09-27T20:00:00Z'})
            output = Path(stage['output']); value = base.read(output)
            result = {'script': stage['script'], 'path': str(output.relative_to(root)), 'sha256': base.digest(output), 'status': value['status']}
            put(directory/f'{index:02d}-result.json', result); results.append(result)
            (directory/f'{index:02d}.log').write_text('Synthetic invocation completed, no model calls.\n')
        put(directory/'completed.json', {'results': results, 'status': 'execution_complete_scientific_statuses_retained', 'plan_sha256': base.digest(pp)})
    return entry, parent, records


def mutate_child(root, entry, script, callback):
    for spec in entry['receipts']:
        pp = root/spec['plan_path']; plan = base.read(pp)
        for index, row in enumerate(plan['stages']):
            if row['script'] == script:
                path = Path(row['output']); value = base.read(path); callback(value); put(path, value)
                result_path = pp.parent/f'{index:02d}-result.json'
                result = base.read(result_path); result.update(sha256=base.digest(path), status=value['status']); put(result_path, result)
                complete = base.read(pp.parent/'completed.json'); complete['results'][index] = result; put(pp.parent/'completed.json', complete)
                return
    raise AssertionError(script)


def alter_plan(root, entry, index, callback):
    p = root/entry['receipts'][index]['plan_path']; value = base.read(p); callback(value); put(p, value)
    entry['receipts'][index]['plan_sha256'] = base.digest(p)
    complete = base.read(p.parent/'completed.json'); complete['plan_sha256'] = base.digest(p); put(p.parent/'completed.json', complete)


def precision_fixture(root, entry, parent, mode):
    def modify(row):
        row['status'] = 'incomplete_or_failed_numerical_screen'
        row['points'][0].update(status='failed_native_precision_checks', failed_checks=['nominal_background_rdrag_relative'])
        row['failures'] = [{'audit_index': 0, 'status': 'failed_native_precision_checks'}]
    mutate_child(root, entry, 'native_posterior_precision.py', modify)
    p = root/entry['receipts'][1]['plan_path']
    screen = Path(base.read(p)['stages'][0]['output'])
    review = root/'.work/thermal-review.json'; put(review, {'status': 'explicit_synthetic_review', 'input_sha256': parent['input_sha256']})
    entry['precision_review'] = {'review_path': str(review.relative_to(root)), 'review_sha256': base.digest(review)}
    success = mode in ['supported', 'wrong_parent', 'contradictory_flags', 'changed_saved_receipt']
    receipt = {'schema': 'verified-supplemental-native-precision-v1', 'status': 'verified_saved_precision_review',
        'parent_freshly_qualified': True, 'parent_target_identity': parent['target_identity'],
        'original_screen_status': 'incomplete_or_failed_numerical_screen',
        'reviewed_screen_status': 'no_large_variation_detected_on_fixed32' if success else 'numerical_sensitivity_requires_followup',
        'original_diagnostic_flags': [], 'original_failures': [{'status': 'failed_native_precision_checks'}],
        'reviewed_diagnostic_flags': [] if success else ['component_variation:synthetic'],
        'reviewed_failures': [], 'reviewed_precision_screen_supported': success,
        'original_screen_ref': {'path': str(screen.relative_to(root)), 'sha256': base.digest(screen)},
        'supplemental_review_ref': {'path': str(review.relative_to(root)), 'sha256': base.digest(review)},
        'input_sha256': parent['input_sha256'], 'posterior_qualification': False,
        'posterior_reweighting_performed': False, 'background_or_native_or_model_calls': 0}
    if mode == 'wrong_parent': receipt['parent_target_identity'] = 'b'*64
    if mode == 'contradictory_flags': receipt['reviewed_diagnostic_flags'] = ['total_centered_RMS']
    if mode == 'real_failure':
        receipt.update(reviewed_screen_status='incomplete_or_failed_numerical_screen', reviewed_diagnostic_flags=[],
                       reviewed_failures=[{'status': 'high_background_failure'}], reviewed_precision_screen_supported=False)
    if mode == 'changed_saved_receipt':
        saved = root/'.work/consumer-receipt.json'; changed = deepcopy(receipt); changed['original_failures'] = []; put(saved, changed)
        entry['precision_review'].update(verified_receipt_path=str(saved.relative_to(root)), verified_receipt_sha256=base.digest(saved))
    return lambda *args: deepcopy(receipt)


def run():
    cases = []
    def case(name, modify=lambda *x: None, mode=None, predicate=None, rejected_parent=False):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); entry, parent, records = fixture(root)
            reviewer = precision_fixture(root, entry, parent, mode) if mode else None
            modify(root, entry, parent, records)
            original_status = (root/entry['original_status_path']).read_bytes()
            calls = []
            def qualifier(folder, correction):
                calls.append(str(correction))
                if rejected_parent: raise AssertionError('Synthetic parent no longer qualified')
                return deepcopy(parent)
            output = auditor.audit_cohort(root, entry, _qualifier=qualifier, _reviewer=reviewer)
            assert (root/entry['original_status_path']).read_bytes() == original_status
            if predicate: assert predicate(output), (name, output['integrity_errors'], output['effective_stages'])
            minimal = {'schema':'supplementary-scientific-audit-v1','status':auditor.aggregate_status([output]),'cohorts':[output],
                       'physical_model_or_background_calls':0,'execution_completion_is_scientific_qualification':False,'original_failure_records_unchanged':True}
            auditor.validate_report(minimal)
            if name == 'completed_omitted_BAO_still_scientifically_failed': assert minimal['status']=='complete_with_failed_gates'
            if name == 'genuine_density_variation_not_repaired': assert minimal['status']=='complete_with_failed_gates'
            cases.append({'case': name, 'passed': True, 'mock_parent_qualifications': len(calls),
                'native_accuracy_diagnostic_supported': output['native_accuracy_diagnostic_supported'],
                'typed_plot_histories': sorted(output['verified_children']), 'integrity_errors_expected': bool(output['integrity_errors'])})
    invalid = lambda r: bool(r['integrity_errors']) and not r['verified_children'] and r['posterior_at_declared_accuracy'] is None
    good = lambda r: not r['integrity_errors'] and r['native_accuracy_diagnostic_supported'] and set(r['verified_children']) == {'expansion_history', 'luminosity_history'} and r['original']['execution_snapshot']['stage'] == 'stopped_no_retry'
    case('completed_receipts_with_qualified_scientific_outputs', predicate=good)
    case('missing_completed_receipt_remains_pending',
         lambda root,e,*_: (root/e['receipts'][1]['plan_path']).with_name('completed.json').unlink(),
         predicate=lambda r: not r['integrity_errors'] and not r['native_accuracy_diagnostic_supported'] and not r['verified_children'])
    def failed_terminal(root, e, *_):
        pp = root/e['receipts'][1]['plan_path']; (pp.parent/'completed.json').unlink()
        put(pp.parent/'failure.json', {'error':'Synthetic failed process','completed':[], 'plan_sha256':base.digest(pp),'retry_authorized':False})
    case('failed_execution_not_promoted', failed_terminal,
         predicate=lambda r: not r['integrity_errors'] and r['supplemental_receipts'][1]['execution']=='failed' and not r['verified_children'])
    case('duplicate_receipt_rejected', lambda root,e,*_: e['receipts'].append(deepcopy(e['receipts'][0])), predicate=invalid)
    case('missing_log_rejected', lambda root,e,*_: (root/e['receipts'][0]['plan_path']).with_name('00.log').unlink(), predicate=invalid)
    case('missing_stage_result_rejected', lambda root,e,*_: (root/e['receipts'][0]['plan_path']).with_name('00-result.json').unlink(), predicate=invalid)
    case('changed_output_bytes_rejected', lambda root,e,p,r: Path(r['expansion_history.py'][1]).write_text('{}\n'), predicate=invalid)
    case('changed_source_rejected', lambda root,e,p,r: Path(r['expansion_history.py'][0]['argv'][1]).write_text('# Changed\n'), predicate=invalid)
    case('changed_attempt_argv_rejected', lambda root,e,*_: put((root/e['receipts'][0]['plan_path']).with_name('00-attempt.json'),
         {'script':'quantile_precision_ordered.py','argv':['not-the-plan'],'started_utc':'2026-09-27T20:00:00Z'}), predicate=invalid)
    case('unauthorized_plan_argv_rejected', lambda root,e,*_: alter_plan(root,e,0,lambda p:p['stages'][1]['argv'].append('--changed')), predicate=invalid)
    case('wrong_original_failure_hash_rejected', lambda root,e,*_: alter_plan(root,e,0,lambda p:p['original_failed_plan_and_status_sha256'].update({e['original_status_path']:'0'*64})), predicate=invalid)
    case('wrong_parent_target_rejected', lambda root,e,*_: mutate_child(root,e,'expansion_history.py',lambda r:r.update(target_identity='b'*64)), predicate=invalid)
    case('missing_parent_binding_rejected', lambda root,e,*_: mutate_child(root,e,'luminosity_history.py',lambda r:r.pop('input_sha256')), predicate=invalid)
    case('false_success_with_failed_gates_rejected', lambda root,e,*_: mutate_child(root,e,'expansion_history.py',lambda r:r.update(failed_gates=['stability'])), predicate=invalid)
    case('completed_omitted_BAO_still_scientifically_failed',
         lambda root,e,*_: mutate_child(root,e,'probe_omission.py',lambda r:r.update(status='insufficient_probe_omission_overlap_or_stability',qualified_under_declared_numerical_gates=False,posterior=None,failed_gates=['raw_weight_ESS','Pareto_k','weighted_chain_stability'])),
         predicate=lambda r:not r['integrity_errors'] and r['effective_stages']['probe_omission:bao']['state']=='failed' and all(x['execution']=='complete'for x in r['supplemental_receipts']))
    case('unqualified_parent_withheld', rejected_parent=True, predicate=invalid)
    case('auxiliary_thermal_repair_preserves_original_failure', mode='supported', predicate=lambda r: good(r) and r['effective_stages']['native_posterior_precision']['state']=='failed' and r['precision_review']['original_failures'])
    case('genuine_density_variation_not_repaired', mode='variation', predicate=lambda r:not r['integrity_errors'] and not r['native_accuracy_diagnostic_supported'] and not r['verified_children'] and bool(r['verified_conditional_histories']))
    case('genuine_high_background_failure_not_repaired', mode='real_failure', predicate=lambda r:not r['integrity_errors'] and not r['verified_children'] and r['effective_precision_diagnostic']['state']=='failed')
    case('wrong_thermal_review_parent_rejected', mode='wrong_parent', predicate=invalid)
    case('contradictory_supported_review_flags_rejected', mode='contradictory_flags', predicate=invalid)
    case('changed_saved_pure_consumer_receipt_rejected', mode='changed_saved_receipt', predicate=invalid)
    sources = [Path(__file__), Path(auditor.__file__), auditor.DESIGN, Path(base.__file__),
               base.CODE/'final_audit_validate.py',base.CODE/'final-audit-schema.json',
               base.CODE/'inference/postprocessing_continuation.py',base.CODE/'inference/measurement_summary.py',
               base.CODE/'inference/native_precision_review_consumer.py']
    return {'schema':'supplementary-audit-synthetic-validation-v1','status':'passed', 'cases':cases,
        'physical_model_calls':0,'background_calls':0,'observational_audit_runs':0,
        'scope':'Synthetic receipt trees, source/input hashes and explicit mocked parent/thermal-review consumers; real qualification/review remain delegated unchanged to frozen consumers.',
        'source_sha256':{base.relative(base.ROOT,p):base.digest(p) for p in sources}}


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,default=OUTPUT);args=parser.parse_args()
    sys_path = str(base.CODE/'inference')
    import sys
    sys.path.insert(0, sys_path)
    import camb
    import cobaya.model
    def forbidden(*args,**kwargs):raise AssertionError('Physical call forbidden in synthetic validation')
    with ExitStack() as stack:
        for owner,name in [(camb,'get_results'),(camb,'get_background'),(camb,'get_transfer_functions'),(cobaya.model,'get_model'),(cobaya.model.Model,'__init__')]:
            stack.enter_context(patch.object(owner,name,forbidden))
        result=run()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'cases':len(result['cases']),'sha256':base.digest(args.output)}))


if __name__=='__main__':main()
