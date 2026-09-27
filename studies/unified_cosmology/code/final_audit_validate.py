"""Synthetic decision/identity tests; no observational posterior or native call."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import numpy as np
from final_audit import (ROOT, CODE, DESIGN, SCHEMA, SUCCESS, audit_cohort, audit_all,
                         digest, scientific_state, check_precision, validate_report)


def put(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2)+'\n')


def fixture(root, evolution='none'):
    folder = root/'.work/synthetic-cohort'
    producer = root/'studies/unified_cosmology/code/inference'
    folder.mkdir(parents=True)
    producer.mkdir(parents=True)
    input_path = root/'.work/synthetic-input.txt'
    input_path.write_text('SYNTHETIC data, not an observed cosmology\n')
    inputs = {str(input_path.relative_to(root)): digest(input_path)}
    settings = {'model': 'cpl', 'evolution': evolution, 'sample': 'synthetic', 'calibration': 'synthetic'}
    target = 'a'*64
    parent = {'settings': settings, 'target_identity': target, 'qualified_under_declared_numerical_gates': True,
              'input_sha256': inputs, 'posterior': {'synthetic_coordinate': {'mean': 12345.}}}
    for i in range(4):
        put(folder/f'run-{i}.json', {'arguments': settings, 'target_identity': {'identity': target}})
    scripts = list(SUCCESS)
    scripts.remove('luminosity_bridge.py')
    scripts.remove('joint_lensing_bridge.py')
    stages = []
    records = {}
    precision_design = {'points': 32, 'diagnostic_thresholds': {'total_centered_RMS_loglike_max': .05,
                         'total_centered_max_absolute_loglike_max': .1, 'component_centered_max_absolute_loglike_max': .2}}
    precision_path = folder/'precision-plan.json'
    put(precision_path, {'design': precision_design, 'qualified_parent_inputs': inputs})
    def add(script, suffix='', omission=None):
        source = producer/script
        source.write_text('# Synthetic producer fixture; never executed.\n')
        output = folder/(Path(script).stem+suffix+'.json')
        stage = {'script': script, 'source_sha256': digest(source), 'output': str(output),
                 'argv': ['python', str(source)]+(['--omit', omission] if omission else []),
                 'condition': 'native_precision_screen_passed' if script == 'joint_lensing_bridge.py' else None}
        row = {'status': SUCCESS[script], 'input_sha256': inputs, 'source_sha256': {str(source.relative_to(root)): digest(source)},
               'target_identity': target, 'settings': settings, 'failed_gates': []}
        if script == 'measurement_summary.py':
            row = {'status': SUCCESS[script], 'source_sha256': digest(source), 'measurements': [deepcopy(parent)]}
        if script == 'expansion_history.py':
            row['history_weighted_stability'] = {'qualified_overlap': True, 'failed_gates': []}
        if script == 'sn_predictive_check.py':
            row.update(observations=10, replicate_degrees_of_freedom=9, subtracted_unpenalized_intercepts=1,
                       subtracted_cosmological_or_proper_Gaussian_mode_counts=0, posterior_tail_diagnostics={})
        if script == 'native_posterior_precision.py':
            row.update(posterior_qualification=False, points=[{'status': 'finite_native_precision',
                'comparison': {'high_minus_declared_total_loglike': 0., 'high_minus_declared_component_loglikes': {'A': 0.}}} for _ in range(32)],
                plan_path=str(precision_path.relative_to(root)), plan_sha256=digest(precision_path), diagnostic_flags=[])
        if script in {'probe_omission.py', 'luminosity_bridge.py'}:
            row.update(qualified_under_declared_numerical_gates=True, posterior={'synthetic_coordinate': {}}, points=2000)
        if omission:
            row.update(omission=omission, removed_components={'sn': ['released_sn'], 'bao': ['bao.desi_dr2'],
                       'lensing': ['act_dr6_lenslike.ACTDR6LensLike', 'SPT2023_lensing']}[omission],
                       unused_auxiliary_coordinates=['epsilon'] if omission == 'sn' and evolution == 'linear' else [])
        if script == 'joint_lensing_bridge.py':
            row['variants'] = {name: {'status': 'qualified_conditional_bridge', 'qualified_under_declared_numerical_gates': True,
                        'points': 2000, 'posterior': {'synthetic_coordinate': {}}, 'failed_gates': [], 'weighted_covariance': {'qualified': True}}
                        for name in ['baseline', 'extended']}
        put(output, row)
        records[script+suffix] = (stage, output)
        stages.append(stage)
    for script in scripts:
        if script == 'probe_omission.py':
            for omission in ['sn', 'bao', 'lensing']:
                add(script, '-'+omission, omission)
        else:
            add(script)
    if evolution == 'smooth01':
        add('luminosity_bridge.py')
    add('joint_lensing_bridge.py')
    plan = folder/'finite-pipeline-plan.json'
    put(plan, {'label': 'synthetic', 'folder': str(folder), 'stages': stages})
    spec = {'label': 'synthetic', 'plan_path': str(plan.relative_to(root)), 'plan_sha256': digest(plan),
            'settings': settings, 'target_identity': target}
    status_path = folder/'finite-pipeline-status.json'
    progress = {'stage': 'finite_execution_completed_result_review_required', 'plan_sha256': digest(plan),
                'completed_results': [{'script': s['script'], 'path': str(p.relative_to(root)), 'sha256': digest(p),
                                      'status': json.loads(p.read_text())['status']} for s, p in records.values()]}
    put(status_path, progress)
    def rewrite(name, function):
        stage, path = records[name]
        value = json.loads(path.read_text()); function(value); put(path, value)
        for receipt in progress['completed_results']:
            if receipt['path'] == str(path.relative_to(root)):
                receipt.update(sha256=digest(path), status=value['status'])
        put(status_path, progress)
    return spec, parent, records, progress, status_path, rewrite


def run_tests():
    results = []
    def case(name, mutation, predicate, evolution='none', reject_parent=False):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            spec, parent, records, progress, status_path, rewrite = fixture(root, evolution)
            mutation(spec, parent, records, progress, status_path, rewrite)
            calls = []
            def qualify(folder, correction):
                calls.append(str(correction))
                if reject_parent:
                    raise ValueError('SYNTHETIC current scientific identity rejected')
                return deepcopy(parent)
            actual = audit_cohort(root, spec, _qualifier=qualify)
            assert predicate(actual), (name, actual['integrity_errors'])
            results.append({'case': name, 'passed': True, 'mock_parent_calls': len(calls),
                            'states': {k: v['state'] for k, v in actual['stages'].items()},
                            'posterior_exposed': actual['posterior'] is not None})
    noop = lambda *args: None
    good = lambda x: not x['integrity_errors'] and x['posterior'] is not None and set(x['verified_children']) == {'expansion_history', 'luminosity_history'}
    invalid = lambda key: lambda x: x['stages'][key]['state'] == 'invalid'
    case('fully_qualified_synthetic_pipeline', noop, good)
    case('smooth_extra_stage', noop, good, 'smooth01')
    case('linear_auxiliary_not_reported', noop, good, 'linear')
    def pending(spec, parent, records, progress, path, rewrite):
        progress['completed_results'] = []; progress['stage'] = 'waiting_for_one_sampler_lifetime'; put(path, progress)
    case('present_files_without_receipts_are_pending', pending, lambda x: x['posterior'] is None and all(v['state'] == 'pending' for v in x['stages'].values()))
    def missing(spec, parent, records, progress, path, rewrite):
        pending(spec, parent, records, progress, path, rewrite)
        for _, output in records.values():
            output.unlink()
    case('missing_outputs_are_pending_not_passed', missing, lambda x: x['posterior'] is None and not x['integrity_errors'])
    case('current_parent_rejection_withholds_all_numbers', noop, lambda x: x['posterior'] is None and not x['verified_children'], reject_parent=True)
    case('wrong_child_target', lambda *a: a[-1]('expansion_history.py', lambda r: r.update(target_identity='b'*64)), invalid('expansion_history'))
    case('wrong_child_settings', lambda *a: a[-1]('luminosity_history.py', lambda r: r.update(settings={'model': 'wrong'})), invalid('luminosity_history'))
    case('child_missing_parent_binding_not_borrowed_from_sibling', lambda *a: a[-1]('luminosity_history.py', lambda r: r.pop('input_sha256')), invalid('luminosity_history'))
    case('successful_label_with_failed_gates', lambda *a: a[-1]('expansion_history.py', lambda r: r.update(failed_gates=['stability'])), invalid('expansion_history'))
    case('false_stability_boolean', lambda *a: a[-1]('expansion_history.py', lambda r: r['history_weighted_stability'].update(qualified_overlap=False)), invalid('expansion_history'))
    case('unknown_status_is_not_success', lambda *a: a[-1]('luminosity_history.py', lambda r: r.update(status='finished')), invalid('luminosity_history'))
    case('wrong_saved_measurement', lambda *a: a[-1]('measurement_summary.py', lambda r: r['measurements'][0]['posterior']['synthetic_coordinate'].update(mean=-999)),
         lambda x: x['posterior'] is None and x['stages']['measurement_summary']['state'] == 'invalid')
    def changed_bytes(spec, parent, records, progress, path, rewrite):
        output = records['luminosity_history.py'][1]
        output.write_text(output.read_text()+' ')
    case('changed_receipt_bound_output', changed_bytes, invalid('luminosity_history'))
    def changed_source(spec, parent, records, progress, path, rewrite):
        output = Path(records['luminosity_history.py'][0]['argv'][1]); output.write_text('# changed\n')
    case('changed_planned_producer', changed_source, invalid('luminosity_history'))
    case('bad_payload_seal', lambda *a: a[-1]('luminosity_history.py', lambda r: r.update(payload_sha256='0'*64)), invalid('luminosity_history'))
    case('wrong_SN_df', lambda *a: a[-1]('sn_predictive_check.py', lambda r: r.update(replicate_degrees_of_freedom=8)), invalid('sn_predictive_check'))
    case('wrong_omitted_component', lambda *a: a[-1]('probe_omission.py-bao', lambda r: r.update(removed_components=['released_sn'])), invalid('probe_omission:bao'))
    case('unused_epsilon_cannot_be_published', lambda *a: a[-1]('probe_omission.py-sn', lambda r: r['posterior'].update(epsilon={})), invalid('probe_omission:sn'), 'linear')
    case('missing_lensing_variant', lambda *a: a[-1]('joint_lensing_bridge.py', lambda r: r['variants'].pop('extended')), invalid('joint_lensing_bridge'))
    def failed_variant(*args):
        args[-1]('joint_lensing_bridge.py', lambda r: r['variants']['extended'].update(status='insufficient_bridge_overlap_or_stability',
                    qualified_under_declared_numerical_gates=False, posterior=None, failed_gates=['ESS']))
    case('one_lensing_variant_failed_not_hidden_by_top_status', failed_variant,
         lambda x: x['stages']['joint_lensing_bridge']['state'] == 'failed' and not x['integrity_errors'])
    case('false_variant_qualification', lambda *a: a[-1]('joint_lensing_bridge.py', lambda r: r['variants']['baseline'].update(qualified_under_declared_numerical_gates=False)), invalid('joint_lensing_bridge'))
    case('fake_precision_pass_recomputed_from_points', lambda *a: a[-1]('native_posterior_precision.py', lambda r: r['points'][0]['comparison'].update(high_minus_declared_total_loglike=10.)), invalid('native_posterior_precision'))
    def fail_precision(spec, parent, records, progress, path, rewrite):
        def alter(r):
            r.update(status='numerical_sensitivity_requires_followup', diagnostic_flags=['total_centered_RMS', 'total_centered_max_absolute'])
            r['points'][0]['comparison']['high_minus_declared_total_loglike'] = 10.
        rewrite('native_posterior_precision.py', alter)
        stage, output = records['joint_lensing_bridge.py']
        progress['completed_results'] = [v for v in progress['completed_results'] if v['script'] != 'joint_lensing_bridge.py']
        progress['completed_results'].append({'script': 'joint_lensing_bridge.py', 'status': 'not_run_native_precision_support_gate'})
        output.unlink(); put(path, progress)
    case('precision_failure_and_deliberate_not_run_are_distinct', fail_precision,
         lambda x: x['stages']['native_posterior_precision']['state'] == 'failed' and x['stages']['joint_lensing_bridge']['state'] == 'not_run' and not x['integrity_errors'])
    # Independent numerical oracle: numpy centered spread versus scalar audit.
    rng = np.random.default_rng(273331)
    design = {'points': 32, 'diagnostic_thresholds': {'total_centered_RMS_loglike_max': .05,
              'total_centered_max_absolute_loglike_max': .1, 'component_centered_max_absolute_loglike_max': .2}}
    for scale in [0., .01, .1, 1.]:
        values = rng.normal(size=(32, 2))*scale
        total = values.sum(axis=1); centered = total-total.mean()
        flags = []
        if np.sqrt(np.mean(centered**2)) > .05: flags.append('total_centered_RMS')
        if np.max(abs(centered)) > .1: flags.append('total_centered_max_absolute')
        flags += ['component_variation:'+str(j) for j in range(2) if np.max(abs(values[:, j]-values[:, j].mean())) > .2]
        row = {'status': 'numerical_sensitivity_requires_followup' if flags else 'no_large_variation_detected_on_fixed32',
               'diagnostic_flags': flags, 'points': [{'status': 'finite_native_precision', 'comparison': {
               'high_minus_declared_total_loglike': float(total[i]), 'high_minus_declared_component_loglikes': {str(j): float(values[i, j]) for j in range(2)}}} for i in range(32)]}
        check_precision(row, design)
    # A manually specified four-row pending report tests result-level withholding,
    # independently of audit_cohort construction and without a mock measurement.
    manual = {'schema': 'final-scientific-audit-v1', 'status': 'pending_scientific_completion',
              'new_model_or_background_evaluations': 0, 'execution_completion_is_scientific_qualification': False,
              'cohorts': [{'label': str(i), 'settings': {'synthetic_case': i}, 'target_identity': 'a'*64,
                'parent_qualification': 'pending', 'posterior': None, 'verified_children': {}, 'integrity_errors': [],
                'stages': {'measurement_summary': {'state': 'pending', 'failed_gates': []}}} for i in range(4)]}
    validate_report(manual)
    for mutation in [lambda r: r.update(status='completed_conditional_results_audited'),
                     lambda r: r['cohorts'][0].update(posterior={'forbidden_coordinate': 1.}),
                     lambda r: r.update(new_model_or_background_evaluations=1)]:
        altered = deepcopy(manual); mutation(altered)
        try:
            validate_report(altered)
        except AssertionError:
            pass
        else:
            raise AssertionError('Invalid output report was accepted')
    return results


def main():
    checks = run_tests()
    import ast
    for path in [CODE/'final_audit.py', Path(__file__)]:
        ast.parse(path.read_text())
    result = {'status': 'passed_synthetic_audit_validation', 'synthetic_cases': checks,
              'synthetic_case_count': len(checks), 'independent_numpy_precision_oracle_cases': 4,
              'scientific_model_or_background_calls': 0,
              'qualification_dependency': 'Mocked ONLY in temporary synthetic fixtures; real CLI has no qualification bypass.',
              'manual_result_schema_cases': 4,
              'source_sha256': {str(path.relative_to(ROOT)): digest(path) for path in [CODE/'final_audit.py', Path(__file__), DESIGN, SCHEMA]},
              'scope': 'Decision/identity/withholding validation, not a cosmological posterior or certification of pending chains.'}
    output = ROOT/'studies/unified_cosmology/results/final-audit-validation.json'
    put(output, result)
    print(json.dumps({'status': result['status'], 'synthetic_case_count': len(checks)}))


if __name__ == '__main__':
    main()
