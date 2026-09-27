"""Read-only audit of explicit postprocessing continuations and thermal reviews.

Original finite-pipeline failures stay immutable. Execution receipts certify a
completed invocation, never its scientific success. No cosmology is evaluated.
"""
import argparse
import copy
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

import final_audit as base

ROOT = base.ROOT
CODE = Path(__file__).resolve().parent
DESIGN = CODE/'supplementary-audit-design.json'
ALLOWED = {'quantile_precision.py', 'expansion_history.py', 'luminosity_history.py',
           'sn_predictive_check.py', 'native_posterior_precision.py',
           'probe_omission.py', 'luminosity_bridge.py'}


def state(script, record):
    if script == 'quantile_precision_ordered.py':
        return base.scientific_state('quantile_precision.py', record)
    return base.scientific_state(script, record)


def review_precision(screen, review):
    sys.path.insert(0, str(CODE/'inference'))
    from native_precision_review_consumer import verify
    return verify(screen, review)


def stage_match(root, row, originals):
    original_name = row.get('preserved_original_script', row['script'])
    assert original_name in ALLOWED, 'Unsupported continuation producer.'
    matches = [s for s in originals if s['script'] == original_name and s['output'] == row['output']]
    assert len(matches) == 1, 'Stage does not identify one original planned output.'
    expected = copy.deepcopy(matches[0])
    assert expected['condition'] is None, 'Conditional stages require a distinct reviewed launcher.'
    if original_name == 'quantile_precision.py':
        source = root/'studies/unified_cosmology/code/inference/quantile_precision_ordered.py'
        old_source = root/'studies/unified_cosmology/code/inference/quantile_precision.py'
        expected.update(preserved_original_script=original_name,
                        preserved_original_source_sha256=expected['source_sha256'],
                        script=source.name, source_sha256=base.digest(source))
        expected['argv'] = [str(source) if x == str(old_source) else x for x in expected['argv']]
    assert row == expected, 'Supplemental stage changed original argv/output/conditions or an unauthorized producer.'
    return matches[0]


def verify_child(root, stage, receipt, parent, spec, evidence):
    output = Path(stage['output'])
    assert receipt['script'] == stage['script']
    assert receipt['path'] == base.relative(root, output)
    child = base.Evidence(root)
    record = child.pin(output, receipt['sha256'], inspect=True)
    assert record['status'] == receipt['status']
    result = state(stage['script'], record)
    for path, expected in parent['input_sha256'].items():
        assert child.bindings.get(path) == expected, 'Child lacks complete qualified-parent binding: '+path
    for key in ['target_identity', 'parent_target_identity', 'parent_proposal_target_identity', 'frozen_target_identity']:
        if key in record:
            assert record[key] == spec['target_identity'], 'Supplemental output belongs to another target.'
    for key in ['settings', 'parent_settings', 'source_settings']:
        if key in record:
            assert record[key] == spec['settings'], 'Supplemental output settings differ from parent.'
    source = root/'studies/unified_cosmology/code/inference'/stage['script']
    assert child.bindings.get(base.relative(root, source)) == stage['source_sha256'], 'Child does not bind its actual producer.'
    if stage['script'] == 'native_posterior_precision.py':
        plan = base.read(root/record['plan_path'])
        base.check_precision(record, plan['design'])
    if stage['script'] == 'probe_omission.py':
        omission = stage['argv'][stage['argv'].index('--omit')+1]
        expected = {'sn': ['released_sn'], 'bao': ['bao.desi_dr2'],
                    'lensing': ['act_dr6_lenslike.ACTDR6LensLike', 'SPT2023_lensing']}[omission]
        assert record['omission'] == omission and record['removed_components'] == expected
        if omission == 'sn' and spec['settings']['evolution'] == 'linear':
            assert 'epsilon' in record['unused_auxiliary_coordinates']
            assert record.get('posterior') is None or 'epsilon' not in record['posterior']
    evidence.bindings.update(child.bindings)
    return record, dict(result, output_path=receipt['path'], output_sha256=receipt['sha256'],
                        script=stage['script'], execution='completed_invocation')


def audit_cohort(root, entry, *, _qualifier=None, _reviewer=None):
    root = Path(root).resolve()
    spec = entry['original_spec']
    evidence = base.Evidence(root)
    out = {'label': spec['label'], 'target_identity': spec['target_identity'], 'settings': spec['settings'],
           'original': None, 'supplemental_receipts': [], 'effective_stages': {},
           'parent_qualification': 'pending', 'posterior_at_declared_accuracy': None,
           'precision_review': None, 'native_accuracy_diagnostic_supported': False,
           'native_accuracy_supported_downstream_allowed': False,
           'verified_conditional_histories': {}, 'verified_children': {}, 'integrity_errors': []}
    parents = []
    def qualifier(folder, correction):
        parent = (_qualifier or base.default_qualifier)(folder, correction)
        parents.append(parent)
        return parent
    try:
        original = base.audit_cohort(root, spec, _qualifier=qualifier)
        out['original'] = original
        assert not original['integrity_errors'], 'Original cohort fails its unchanged integrity audit.'
        assert len(parents) == 1 and original['parent_qualification'] == 'qualified_at_declared_native_accuracy'
        parent = parents[0]
        assert parent['qualified_under_declared_numerical_gates'] is True
        assert parent['target_identity'] == spec['target_identity'] and parent['settings'] == spec['settings']
        out['parent_qualification'] = original['parent_qualification']
        out['posterior_at_declared_accuracy'] = original['posterior']
        out['effective_stages'] = copy.deepcopy(original['stages'])
        evidence.bindings.update(original['evidence_sha256'])
        plan_path = root/spec['plan_path']; plan = base.read(plan_path); folder = plan_path.parent
        status_path = root/entry['original_status_path']
        assert status_path == folder/'finite-pipeline-status.json'
        evidence.pin(status_path, entry['original_status_sha256'])
        status = base.read(status_path)
        assert status['stage'] == 'stopped_no_retry' and status['plan_sha256'] == spec['plan_sha256']
        assert 'quantile_precision.py' in status['error']
        assert [r['script'] for r in status['completed_results']] == ['diagnostics.py', 'spectral_correction.py', 'measurement_summary.py']
        assert status['completed_results'][-1]['status'] == 'qualified_conditional_measurements'
        original_bindings = {spec['plan_path']: spec['plan_sha256'], entry['original_status_path']: entry['original_status_sha256']}
        original_outputs = {r['path'] for r in status['completed_results']}
        receipt_specs = entry['receipts']
        receipt_paths = [str((root/s['plan_path']).resolve()) for s in receipt_specs]
        assert receipt_specs and len(set(receipt_paths)) == len(receipt_paths), 'Duplicate continuation receipt.'
        seen_outputs = set(original_outputs)
        screens = {}
        for requested in receipt_specs:
            receipt_plan_path = root/requested['plan_path']
            evidence.pin(receipt_plan_path, requested['plan_sha256'])
            receipt_plan = base.read(receipt_plan_path); receipt_folder = receipt_plan_path.parent
            assert receipt_plan_path.name == 'plan.json'
            assert receipt_plan['schema'] == 'explicit-postprocessing-continuation-v1'
            assert receipt_plan['original_failed_plan_and_status_sha256'] == original_bindings
            assert Path(receipt_plan['folder']).resolve() == folder
            assert receipt_plan['target_identity'] == spec['target_identity']
            continuation = root/'studies/unified_cosmology/code/inference/postprocessing_continuation.py'
            evidence.pin(continuation, receipt_plan['source_sha256'])
            rows = receipt_plan['stages']; assert rows
            for row in rows:
                previous = stage_match(root, row, plan['stages'])
                key = base.stage_key(previous)
                output_relative = base.relative(root, row['output'])
                assert output_relative not in seen_outputs, 'Duplicate or previously completed stage output.'
                seen_outputs.add(output_relative)
                evidence.pin(root/'studies/unified_cosmology/code/inference'/row['script'], row['source_sha256'])
                assert out['effective_stages'][key]['state'] == 'pending'
            completed_path, failure_path = receipt_folder/'completed.json', receipt_folder/'failure.json'
            assert not (completed_path.exists() and failure_path.exists()), 'Contradictory terminal continuation receipts.'
            receipt_state = {'plan_path': requested['plan_path'], 'plan_sha256': requested['plan_sha256'],
                             'execution': 'pending', 'outputs': []}
            out['supplemental_receipts'].append(receipt_state)
            # No partial invocation is promoted from the presence of an output.
            if not completed_path.exists():
                if failure_path.exists():
                    evidence.pin(failure_path, base.digest(failure_path), inspect=True)
                    failure = base.read(failure_path)
                    assert failure['plan_sha256'] == requested['plan_sha256'] and failure['retry_authorized'] is False
                    receipt_state.update(execution='failed', failure=failure)
                for row in rows:
                    key = base.stage_key(stage_match(root, row, plan['stages']))
                    out['effective_stages'][key].update(pending_reason='complete_supplemental_receipt_missing')
                continue
            evidence.pin(completed_path, base.digest(completed_path), inspect=True)
            completed = base.read(completed_path)
            assert completed['plan_sha256'] == requested['plan_sha256']
            assert completed['status'] == 'execution_complete_scientific_statuses_retained'
            assert len(completed['results']) == len(rows)
            expected_files = {'plan.json', 'completed.json'}
            for index, (row, final_receipt) in enumerate(zip(rows, completed['results'])):
                paths = [receipt_folder/f'{index:02d}-attempt.json', receipt_folder/f'{index:02d}-result.json', receipt_folder/f'{index:02d}.log']
                for p in paths:
                    evidence.pin(p, base.digest(p), inspect=p.suffix == '.json')
                    expected_files.add(p.name)
                attempt, result = [base.read(p) for p in paths[:2]]
                assert set(attempt) == {'script', 'argv', 'started_utc'}
                assert attempt['script'] == row['script'] and attempt['argv'] == row['argv']
                datetime.strptime(attempt['started_utc'], '%Y-%m-%dT%H:%M:%SZ')
                assert result == final_receipt
                record, child_state = verify_child(root, row, result, parent, spec, evidence)
                key = base.stage_key(stage_match(root, row, plan['stages']))
                out['effective_stages'][key] = child_state
                receipt_state['outputs'].append(dict(path=result['path'], sha256=result['sha256'], scientific_state=child_state['state']))
                if key == 'native_posterior_precision': screens[key] = (Path(row['output']), record)
            assert {p.name for p in receipt_folder.iterdir()} == expected_files, 'Unexpected or duplicate continuation artifact.'
            receipt_state['execution'] = 'complete'
        precision = out['effective_stages'].get('native_posterior_precision', {})
        supported = precision.get('state') == 'diagnostic' and precision.get('status') == base.SUCCESS['native_posterior_precision.py']
        review_spec = entry.get('precision_review')
        if review_spec is not None:
            assert 'native_posterior_precision' in screens, 'A reviewed screen needs a completed, verified supplemental screen receipt.'
            screen_path, screen = screens['native_posterior_precision']
            review_path = root/review_spec['review_path']
            evidence.pin(review_path, review_spec['review_sha256'], inspect=True)
            receipt = (_reviewer or review_precision)(screen_path, review_path)
            assert receipt['status'] == 'verified_saved_precision_review'
            assert receipt['parent_freshly_qualified'] is True and receipt['parent_target_identity'] == spec['target_identity']
            assert receipt['original_screen_ref'] == {'path': base.relative(root, screen_path), 'sha256': base.digest(screen_path)}
            assert receipt['supplemental_review_ref'] == {'path': review_spec['review_path'], 'sha256': review_spec['review_sha256']}
            assert receipt['original_screen_status'] == screen['status']
            assert receipt['posterior_qualification'] is False and receipt['posterior_reweighting_performed'] is False
            assert receipt['background_or_native_or_model_calls'] == 0
            supported = receipt['reviewed_precision_screen_supported']
            assert supported == (receipt['reviewed_screen_status'] == base.SUCCESS['native_posterior_precision.py'])
            if supported:
                assert not receipt['reviewed_diagnostic_flags'] and not receipt['reviewed_failures']
            review_evidence = base.Evidence(root)
            review_evidence.walk(receipt)
            for path, expected in parent['input_sha256'].items():
                assert review_evidence.bindings.get(path) == expected, 'Thermal consumer receipt lacks qualified-parent binding.'
            evidence.bindings.update(review_evidence.bindings)
            if 'verified_receipt_path' in review_spec:
                p = root/review_spec['verified_receipt_path']
                evidence.pin(p, review_spec['verified_receipt_sha256'], inspect=True)
                assert base.read(p) == receipt, 'Saved pure-consumer receipt differs from fresh verification.'
            out['precision_review'] = receipt
            # Keep the original failed screen intact in original/effective stage evidence.
            out['effective_precision_diagnostic'] = {'state': 'diagnostic' if supported else 'failed',
                'status': receipt['reviewed_screen_status'], 'posterior_qualification': False,
                'diagnostic_flags': receipt['reviewed_diagnostic_flags'], 'failures': receipt['reviewed_failures']}
        else:
            out['effective_precision_diagnostic'] = dict(precision, posterior_qualification=False)
        out['native_accuracy_diagnostic_supported'] = bool(supported)
        out['native_accuracy_supported_downstream_allowed'] = bool(supported)
        correction = next(s for s in plan['stages'] if s['script'] == 'spectral_correction.py')['output']
        for key in ['expansion_history', 'luminosity_history']:
            row = out['effective_stages'].get(key, {})
            if row.get('state') == 'qualified':
                child = {'kind': key, 'path': row['output_path'], 'sha256': row['output_sha256'],
                    'qualified': True, 'target_identity': spec['target_identity'],
                    'parent_correction_path': base.relative(root, correction), 'parent_correction_sha256': base.digest(correction),
                    'native_accuracy_diagnostic_supported': bool(supported)}
                out['verified_conditional_histories'][key] = child
                if supported: out['verified_children'][key] = child
        for path, expected in list(evidence.bindings.items()): evidence.pin(path, expected)
    except Exception as error:
        out['integrity_errors'].append(str(error) or type(error).__name__)
        out['posterior_at_declared_accuracy'] = None
        out['verified_children'] = {}; out['verified_conditional_histories'] = {}
        out['native_accuracy_diagnostic_supported'] = False
        out['native_accuracy_supported_downstream_allowed'] = False
    out['evidence_sha256'] = evidence.bindings
    return out


def audit_all(root=ROOT, design_path=DESIGN, *, _qualifier=None, _reviewer=None):
    root, design_path = Path(root).resolve(), Path(design_path).resolve()
    design = base.read(design_path)
    assert design['schema'] == 'supplementary-scientific-audit-design-v1'
    entries = design['cohorts']
    assert entries and len({e['original_spec']['label'] for e in entries}) == len(entries)
    cohorts = [audit_cohort(root, entry, _qualifier=_qualifier, _reviewer=_reviewer) for entry in entries]
    result = {'schema': 'supplementary-scientific-audit-v1', 'status': aggregate_status(cohorts),
        'audit_time_utc': datetime.now(timezone.utc).isoformat(), 'cohorts': cohorts,
        'original_failure_records_unchanged': True, 'execution_completion_is_scientific_qualification': False,
        'physical_model_or_background_calls': 0, 'limits': design['limits'],
        'source_sha256': {base.relative(ROOT, p): base.digest(p) for p in
            [Path(__file__), Path(base.__file__), CODE/'final-audit-schema.json', design_path,
             CODE/'inference/postprocessing_continuation.py', CODE/'inference/measurement_summary.py',
             CODE/'inference/native_precision_review_consumer.py']}}
    validate_report(result)
    return result


def aggregate_status(cohorts):
    invalid = any(c['integrity_errors'] for c in cohorts)
    pending = any(r['execution'] == 'pending' for c in cohorts for r in c['supplemental_receipts'])
    failed = any(r['execution'] == 'failed' for c in cohorts for r in c['supplemental_receipts'])
    scientific_failed = any(s['state'] in {'failed', 'invalid'} for c in cohorts for key,s in c['effective_stages'].items()
                            if key != 'native_posterior_precision')
    precision_failed = any(c.get('effective_precision_diagnostic', {}).get('state') in {'failed', 'invalid'} for c in cohorts)
    unsupported = any(not c['native_accuracy_diagnostic_supported'] for c in cohorts)
    failure = failed or scientific_failed or precision_failed
    return ('failed_integrity_audit' if invalid else 'incomplete_with_failed_gates' if pending and failure
            else 'pending_supplemental_completion' if pending else 'complete_with_failed_gates' if failure
            else 'pending_native_precision_support' if unsupported else 'completed_supplemental_results_audited')


def validate_report(report):
    """Typed plot-consumer invariants, independent of receipt completion labels."""
    assert report['schema'] == 'supplementary-scientific-audit-v1'
    assert report['physical_model_or_background_calls'] == 0
    assert report['execution_completion_is_scientific_qualification'] is False
    assert report['original_failure_records_unchanged'] is True
    cohorts = report['cohorts']
    assert cohorts and len({c['label'] for c in cohorts}) == len(cohorts)
    assert report['status'] == aggregate_status(cohorts)
    for cohort in cohorts:
        assert cohort['native_accuracy_supported_downstream_allowed'] == cohort['native_accuracy_diagnostic_supported']
        if cohort['integrity_errors']:
            assert cohort['posterior_at_declared_accuracy'] is None
            assert not cohort['verified_children'] and not cohort['verified_conditional_histories']
            assert not cohort['native_accuracy_diagnostic_supported']
        if not cohort['native_accuracy_diagnostic_supported']:
            assert not cohort['verified_children']
        for group in ['verified_conditional_histories', 'verified_children']:
            for kind, child in cohort[group].items():
                assert kind in {'expansion_history', 'luminosity_history'} and child['kind'] == kind
                assert child['qualified'] is True and child['target_identity'] == cohort['target_identity']
                assert cohort['parent_qualification'] == 'qualified_at_declared_native_accuracy'
                row = cohort['effective_stages'][kind]
                assert row['state'] == 'qualified' and child['path'] == row['output_path'] and child['sha256'] == row['output_sha256']
                assert child['native_accuracy_diagnostic_supported'] == cohort['native_accuracy_diagnostic_supported']
                if group == 'verified_children': assert child == cohort['verified_conditional_histories'][kind]
        for row in cohort['effective_stages'].values():
            assert row['state'] in {'pending', 'not_run', 'failed', 'invalid', 'diagnostic', 'qualified'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--design', type=Path, default=DESIGN,
                        help='Explicit pinned cohort/receipt manifest; no directory autodiscovery.')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = audit_all(design_path=args.design)
    assert not args.output.exists(), 'Write a new audit snapshot; preserve earlier outcomes.'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'output': str(args.output)}))


if __name__ == '__main__': main()
