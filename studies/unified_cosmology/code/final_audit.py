"""Read-only scientific audit of finite pipelines; never evaluates a cosmology."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
CODE = Path(__file__).resolve().parent
DESIGN = CODE/'final-audit-design.json'
SCHEMA = CODE/'final-audit-schema.json'
SUCCESS = {
    'diagnostics.py': 'passed', 'spectral_correction.py': 'passed_importance_weight_gates',
    'measurement_summary.py': 'qualified_conditional_measurements',
    'quantile_precision.py': 'conditional_quantile_precision_diagnostic',
    'expansion_history.py': 'qualified_pointwise_expansion_history',
    'luminosity_history.py': 'qualified_conditional_luminosity_history',
    'sn_predictive_check.py': 'computed_conditional_marginal_Gaussian_quadratic_check',
    'native_posterior_precision.py': 'no_large_variation_detected_on_fixed32',
    'probe_omission.py': 'qualified_conditional_probe_omission',
    'luminosity_bridge.py': 'qualified_conditional_bridge',
    'joint_lensing_bridge.py': 'separate_lensing_sensitivities',
}
DIAGNOSTIC = {'diagnostics.py', 'spectral_correction.py', 'quantile_precision.py',
              'sn_predictive_check.py', 'native_posterior_precision.py'}
UNITS = {'H0': 'km/s/Mpc', 'H_km_s_Mpc': 'km/s/Mpc', 'rdrag': 'Mpc',
         'cosmic_age': 'Gyr', 'B': 'mag', 'epsilon': 'mag', 'q': 'dimensionless',
         'j': 'dimensionless', 'w': 'dimensionless', 'wa': 'dimensionless'}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def relative(root, path):
    return str(Path(path).resolve().relative_to(root.resolve()))


def canonical_digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def gates(value, prefix=''):
    """Preserve every reported failed gate, including nested stability/variants."""
    failures = []
    if isinstance(value, dict):
        for name, item in value.items():
            path = prefix+'.'+name if prefix else name
            if name in {'failed_gates', 'failures', 'failed_checks', 'failed_numerical_gates', 'diagnostic_flags'} and item:
                failures.append({'field': path, 'value': item})
            elif isinstance(item, (dict, list)):
                failures.extend(gates(item, path))
    elif isinstance(value, list):
        for i, item in enumerate(value):
            if isinstance(item, (dict, list)):
                failures.extend(gates(item, f'{prefix}[{i}]'))
    return failures


def scientific_state(script, record):
    """Pure classification; successful execution is deliberately irrelevant."""
    status = record.get('status')
    failed = gates(record)
    if script == 'joint_lensing_bridge.py' and status == SUCCESS[script]:
        variants = record.get('variants', {})
        assert set(variants) == {'baseline', 'extended'}, 'Missing or unexpected joint-lensing variant.'
        states = {name: scientific_state('luminosity_bridge.py', item) for name, item in variants.items()}
        for name, item in variants.items():
            if states[name]['state'] == 'qualified':
                assert item['weighted_covariance']['qualified'] is True
        return {'state': 'qualified' if all(v['state'] == 'qualified' for v in states.values()) else 'failed',
                'status': status, 'failed_gates': failed, 'variants': states}
    if status == SUCCESS[script]:
        assert not failed, 'Successful status contradicts reported failed gates.'
        if script in {'probe_omission.py', 'luminosity_bridge.py'}:
            assert record.get('qualified_under_declared_numerical_gates') is True
            assert record.get('posterior') is not None
            assert record.get('points', 0) >= 2000
        if script == 'native_posterior_precision.py':
            assert record.get('posterior_qualification') is False
            assert len(record['points']) == 32
            assert all(row['status'] == 'finite_native_precision' for row in record['points'])
        if script == 'sn_predictive_check.py':
            assert record['replicate_degrees_of_freedom'] == record['observations']-1
            assert record['subtracted_unpenalized_intercepts'] == 1
            assert record['subtracted_cosmological_or_proper_Gaussian_mode_counts'] == 0
            assert record.get('posterior_tail_diagnostics') is not None
        if script == 'expansion_history.py':
            assert record['history_weighted_stability']['qualified_overlap'] is True
        return {'state': 'diagnostic' if script in DIAGNOSTIC else 'qualified', 'status': status, 'failed_gates': []}
    if isinstance(status, str) and (status.startswith(('failed_', 'insufficient_', 'incomplete_or_failed_')) or status in {'not_converged', 'numerical_sensitivity_requires_followup'}):
        return {'state': 'failed', 'status': status, 'failed_gates': failed or [{'field': 'status', 'value': status}]}
    raise ValueError('Unknown scientific status: '+repr(status))


def check_precision(record, design):
    """Independently reconstruct the finite-screen decision; no theory imports."""
    points = record['points']
    if len(points) != design['points'] or any(p['status'] != 'finite_native_precision' for p in points):
        assert record['status'] == 'incomplete_or_failed_numerical_screen'
        return
    def spread(values):
        assert all(math.isfinite(v) for v in values)
        mean = sum(values)/len(values)
        return math.sqrt(sum((v-mean)**2 for v in values)/len(values)), max(abs(v-mean) for v in values)
    rms, maximum = spread([p['comparison']['high_minus_declared_total_loglike'] for p in points])
    limits = design['diagnostic_thresholds']
    flags = []
    if rms > limits['total_centered_RMS_loglike_max']:
        flags.append('total_centered_RMS')
    if maximum > limits['total_centered_max_absolute_loglike_max']:
        flags.append('total_centered_max_absolute')
    components = points[0]['comparison']['high_minus_declared_component_loglikes']
    assert all(set(p['comparison']['high_minus_declared_component_loglikes']) == set(components) for p in points)
    for key in components:
        _, maximum = spread([p['comparison']['high_minus_declared_component_loglikes'][key] for p in points])
        if maximum > limits['component_centered_max_absolute_loglike_max']:
            flags.append('component_variation:'+key)
    assert record['diagnostic_flags'] == flags, 'Native precision flags do not recompute.'
    expected = 'numerical_sensitivity_requires_followup' if flags else 'no_large_variation_detected_on_fixed32'
    assert record['status'] == expected, 'Native precision status does not recompute.'


class Evidence:
    """Verify hashes recursively, including sealed cache records and manifests."""
    def __init__(self, root):
        self.root = root
        self.bindings = {}
        self.visited = set()

    def pin(self, path, expected, inspect=False):
        p = Path(path)
        if not p.is_absolute():
            p = self.root/p
        key = relative(self.root, p)
        assert p.is_file(), 'Missing evidence: '+key
        assert digest(p) == expected, 'Changed evidence: '+key
        if key in self.bindings:
            assert self.bindings[key] == expected, 'Conflicting hash binding: '+key
        self.bindings[key] = expected
        if inspect and p.suffix == '.json' and key not in self.visited:
            self.visited.add(key)
            value = read(p)
            self.walk(value)
            assert digest(p) == expected, 'Evidence changed during audit: '+key
            return value
        return None

    def walk(self, value):
        if isinstance(value, list):
            for item in value:
                self.walk(item)
            return
        if not isinstance(value, dict):
            return
        if 'payload_sha256' in value:
            assert canonical_digest({k: v for k, v in value.items() if k != 'payload_sha256'}) == value['payload_sha256'], 'Invalid sealed payload.'
        # Hash manifests can be named source_sha256, parent_inputs, or directly
        # contain path->hash pairs. Only actual path-like keys are interpreted.
        for key, item in value.items():
            if isinstance(item, str) and len(item) == 64 and all(c in '0123456789abcdef' for c in item) and '/' in key:
                self.pin(key, item, inspect=True)
            if key.endswith('_sha256') and isinstance(item, str):
                base = key[:-7]
                path = value.get(base+'_path', value.get(base))
                if isinstance(path, str) and ('/' in path or path.endswith('.json')):
                    self.pin(path, item, inspect=True)
            if isinstance(item, (dict, list)):
                self.walk(item)


def default_qualifier(folder, correction):
    sys.path.insert(0, str(CODE/'inference'))
    from measurement_summary import summarize_run
    return summarize_run(folder, correction)


def stage_key(stage):
    key = Path(stage['script']).stem
    if '--omit' in stage['argv']:
        key += ':'+stage['argv'][stage['argv'].index('--omit')+1]
    return key


def audit_cohort(root, spec, *, _qualifier=None):
    """Return plotting-safe verified_children; test injection is not a CLI option."""
    root = Path(root).resolve()
    evidence = Evidence(root)
    plan_path = root/spec['plan_path']
    result = {'label': spec['label'], 'settings': spec['settings'], 'target_identity': spec['target_identity'],
              'parent_qualification': 'pending', 'posterior': None, 'stages': {}, 'verified_children': {},
              'integrity_errors': [], 'units': UNITS}
    try:
        evidence.pin(plan_path, spec['plan_sha256'])
        plan = read(plan_path)
        folder = Path(plan['folder'])
        assert plan['label'] == spec['label'] and folder.resolve() == plan_path.parent.resolve()
        manifests = [read(folder/f'run-{i}.json') for i in range(4)]
        manifest = manifests[0]
        assert all(x['target_identity'] == manifest['target_identity'] for x in manifests)
        assert manifest['target_identity']['identity'] == spec['target_identity']
        assert {k: manifest['arguments'][k] for k in spec['settings']} == spec['settings']
        for i in range(4):
            evidence.pin(folder/f'run-{i}.json', digest(folder/f'run-{i}.json'))
        progress_path = folder/'finite-pipeline-status.json'
        progress = read(progress_path) if progress_path.exists() else {'stage': 'status_not_written', 'completed_results': []}
        if progress_path.exists():
            assert progress['plan_sha256'] == spec['plan_sha256']
        result['execution_snapshot'] = {'stage': progress['stage'], 'sha256': digest(progress_path) if progress_path.exists() else None,
                                        'time_utc': progress.get('time_utc'), 'error': progress.get('error')}
        stages = plan['stages']
        assert len(stages) == (13 if spec['settings']['evolution'] == 'smooth01' else 12)
        correction_stage = next(s for s in stages if s['script'] == 'spectral_correction.py')
        correction = Path(correction_stage['output'])
        receipts = {str((root/r['path']).resolve()): r for r in progress['completed_results'] if 'path' in r}
        skipped = {r['script']: r for r in progress['completed_results'] if r['status'].startswith('not_run_')}
        qualified = None
        if str(correction.resolve()) in receipts:
            receipt = receipts[str(correction.resolve())]
            evidence.pin(correction, receipt['sha256'])
            if read(correction).get('status') == SUCCESS['spectral_correction.py']:
                try:
                    qualified = (_qualifier or default_qualifier)(folder, correction)
                    assert qualified['qualified_under_declared_numerical_gates'] is True
                    assert qualified['target_identity'] == spec['target_identity'] and qualified['settings'] == spec['settings']
                    result['parent_qualification'] = 'qualified_at_declared_native_accuracy'
                except Exception as error:
                    result['parent_qualification'] = 'failed_current_parent_requalification'
                    result['integrity_errors'].append({'stage': 'parent', 'error': str(error) or type(error).__name__})
            else:
                result['parent_qualification'] = 'failed_native_correction_gates'
        for stage in stages:
            key = stage_key(stage)
            output = Path(stage['output'])
            row = {'script': stage['script'], 'output_path': relative(root, output), 'state': 'pending', 'failed_gates': []}
            result['stages'][key] = row
            try:
                assert stage['script'] in SUCCESS
                evidence.pin(root/'studies/unified_cosmology/code/inference'/stage['script'], stage['source_sha256'])
                receipt = receipts.get(str(output.resolve()))
                if receipt is None:
                    if stage['script'] in skipped:
                        row.update(state='not_run', status=skipped[stage['script']]['status'], gate=skipped[stage['script']])
                    else:
                        row['pending_reason'] = 'output_present_without_pipeline_receipt' if output.exists() else 'planned_output_missing'
                    continue
                assert receipt['script'] == stage['script']
                child_evidence = Evidence(root)
                record = child_evidence.pin(output, receipt['sha256'], inspect=True)
                assert record['status'] == receipt['status']
                row.update(scientific_state(stage['script'], record), output_sha256=receipt['sha256'])
                downstream = stage['script'] not in {'diagnostics.py', 'spectral_correction.py'}
                if downstream:
                    assert qualified is not None, 'Child output has no currently qualified parent.'
                    for path, digest_value in qualified['input_sha256'].items():
                        assert child_evidence.bindings.get(path) == digest_value, 'Missing/changed qualified-parent binding: '+path
                    for target_key in ['target_identity', 'parent_target_identity', 'parent_proposal_target_identity', 'frozen_target_identity']:
                        if target_key in record:
                            assert record[target_key] == spec['target_identity'], 'Child belongs to another target.'
                    for settings_key in ['settings', 'parent_settings', 'source_settings']:
                        if settings_key in record:
                            assert record[settings_key] == spec['settings'], 'Child settings differ from parent.'
                if stage['script'] == 'measurement_summary.py':
                    assert record['source_sha256'] == stage['source_sha256']
                    assert record['measurements'] == [qualified], 'Saved measurement differs from fresh qualification.'
                    result['posterior'] = qualified['posterior']
                if stage['script'] == 'probe_omission.py':
                    omission = stage['argv'][stage['argv'].index('--omit')+1]
                    expected = {'sn': ['released_sn'], 'bao': ['bao.desi_dr2'],
                                'lensing': ['act_dr6_lenslike.ACTDR6LensLike', 'SPT2023_lensing']}[omission]
                    assert record['omission'] == omission and record['removed_components'] == expected
                    if omission == 'sn' and spec['settings']['evolution'] == 'linear':
                        assert 'epsilon' in record['unused_auxiliary_coordinates']
                        assert record.get('posterior') is None or 'epsilon' not in record['posterior']
                if stage['script'] == 'native_posterior_precision.py':
                    p = read(root/record['plan_path'])
                    check_precision(record, p['design'])
                if stage['condition']:
                    assert result['stages']['native_posterior_precision'].get('status') == SUCCESS['native_posterior_precision.py'], 'Conditional stage ran without the declared precision support.'
                if key in {'expansion_history', 'luminosity_history'} and row['state'] == 'qualified':
                    result['verified_children'][key] = {'path': relative(root, output), 'sha256': receipt['sha256'],
                        'qualified': True, 'target_identity': spec['target_identity'],
                        'parent_correction_path': relative(root, correction), 'parent_correction_sha256': digest(correction)}
                evidence.bindings.update(child_evidence.bindings)
            except Exception as error:
                row.update(state='invalid', error=str(error) or type(error).__name__)
                result['integrity_errors'].append({'stage': key, 'error': row['error']})
        # Do not expose a parent table whose saved measurement did not verify.
        if result['stages'].get('measurement_summary', {}).get('state') != 'qualified':
            result['posterior'] = None
        result['native_precision_screen'] = result['stages']['native_posterior_precision']['state']
        result['luminosity_meaning'] = {'none': 'B(z)=0 fixed by assumption; zero width is not an empirical measurement.',
            'linear': 'Sampled epsilon*log(1+z)/log(2), in magnitudes; conditional sensitivity, not a host-age measurement.',
            'smooth01': 'Analytically marginalized Gaussian spline prior with 0.1-mag coefficient scale; retain conditional coefficient variance.'}[spec['settings']['evolution']]
        # Recheck all immutable bytes read during this pass, but not the explicitly
        # timestamped mutable execution-status snapshot.
        for path, expected in list(evidence.bindings.items()):
            evidence.pin(path, expected)
    except Exception as error:
        result['integrity_errors'].append({'stage': 'cohort', 'error': str(error) or type(error).__name__})
        result['posterior'] = None
        result['verified_children'] = {}
    result['evidence_sha256'] = evidence.bindings
    return result


def audit_all(root=ROOT, design_path=DESIGN, *, _qualifier=None):
    root, design_path = Path(root).resolve(), Path(design_path).resolve()
    design = read(design_path)
    cohorts = [audit_cohort(root, spec, _qualifier=_qualifier) for spec in design['cohorts']]
    states = [stage['state'] for cohort in cohorts for stage in cohort['stages'].values()]
    errors = any(c['integrity_errors'] for c in cohorts)
    pending = 'pending' in states
    failed = any(s in {'failed', 'not_run', 'invalid'} for s in states)
    status = ('failed_integrity_audit' if errors else 'incomplete_with_failed_gates' if pending and failed
              else 'pending_scientific_completion' if pending else 'complete_with_failed_gates' if failed
              else 'completed_conditional_results_audited')
    report = {'schema': 'final-scientific-audit-v1', 'status': status, 'audit_time_utc': datetime.now(timezone.utc).isoformat(),
            'cohorts': cohorts, 'limits': design['limits'], 'new_model_or_background_evaluations': 0,
            'fresh_parent_qualification': 'Current measurement_summary.summarize_run; numerical records and target/assets reverified, no density evaluation.',
            'source_sha256': {relative(root, p): digest(p) for p in [Path(__file__), design_path, SCHEMA]},
            'execution_completion_is_scientific_qualification': False}
    validate_report(report)
    return report


def validate_report(report):
    """Runtime cross-field schema invariants; no optional package dependency."""
    assert report['schema'] == 'final-scientific-audit-v1'
    assert report['new_model_or_background_evaluations'] == 0
    assert report['execution_completion_is_scientific_qualification'] is False
    cohorts = report['cohorts']
    assert len(cohorts) == 4 and len({c['label'] for c in cohorts}) == 4
    assert len({tuple(sorted(c['settings'].items())) for c in cohorts}) == 4
    states = []
    for cohort in cohorts:
        assert isinstance(cohort['target_identity'], str) and len(cohort['target_identity']) == 64
        assert isinstance(cohort['integrity_errors'], list)
        if cohort['parent_qualification'] != 'qualified_at_declared_native_accuracy':
            assert cohort['posterior'] is None and not cohort['verified_children']
        for name, child in cohort['verified_children'].items():
            assert name in {'expansion_history', 'luminosity_history'}
            assert child['qualified'] is True and child['target_identity'] == cohort['target_identity']
            assert cohort['stages'][name]['state'] == 'qualified'
            assert child['path'] == cohort['stages'][name]['output_path']
            assert child['sha256'] == cohort['stages'][name]['output_sha256']
        if cohort['posterior'] is not None:
            assert cohort['stages']['measurement_summary']['state'] == 'qualified'
        for stage in cohort['stages'].values():
            assert stage['state'] in {'pending', 'not_run', 'failed', 'invalid', 'diagnostic', 'qualified'}
            assert isinstance(stage['failed_gates'], list)
            states.append(stage['state'])
    errors = any(c['integrity_errors'] for c in cohorts)
    pending = 'pending' in states
    failed = any(s in {'failed', 'not_run', 'invalid'} for s in states)
    expected = ('failed_integrity_audit' if errors else 'incomplete_with_failed_gates' if pending and failed
                else 'pending_scientific_completion' if pending else 'complete_with_failed_gates' if failed
                else 'completed_conditional_results_audited')
    assert report['status'] == expected, 'Aggregate status hides pending/failed evidence.'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    report = audit_all()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix('.part')
    temporary.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    temporary.replace(args.output)
    print(json.dumps({'status': report['status'], 'cohorts': [{'label': r['label'], 'parent': r['parent_qualification'],
        'integrity_errors': r['integrity_errors']} for r in report['cohorts']]}))


if __name__ == '__main__':
    main()
