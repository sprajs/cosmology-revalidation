"""Synthetic tests only: no CAMB background, spectra or likelihood evaluations."""
import hashlib
import json
from pathlib import Path
import tempfile

import numpy as np

import native_posterior_precision as audit
from target_identity import canonical


def main():
    import camb
    original_results = camb.get_results; original_background = camb.get_background
    original_transfers = camb.get_transfer_functions
    calls = {'forbidden_attempts': 0}
    def forbidden(*args, **kwargs):
        calls['forbidden_attempts'] += 1
        raise AssertionError('Synthetic validation cannot evaluate native CAMB.')
    camb.get_results = forbidden; camb.get_background = forbidden; camb.get_transfer_functions = forbidden
    design = json.loads(audit.DESIGN.read_text())
    try:
        groups = np.repeat(np.arange(4), 500)
        locations = [{'expanded_index': int(i % 500)*31+7} for i in range(2000)]
        selected = audit.select_indices(groups, locations)
        expected = [31, 93, 156, 218, 281, 343, 406, 468]
        assert len(selected) == 32
        assert all([r['within_chain_ordered_index'] for r in selected if r['chain'] == c] == expected for c in range(4))
        permutation = np.random.default_rng(273021).permutation(2000)
        shuffled = audit.select_indices(groups[permutation], [locations[i] for i in permutation])
        assert [(r['chain'], r['expanded_index']) for r in shuffled] == [(r['chain'], r['expanded_index']) for r in selected]
        configurations = 0
        for model in ['lcdm', 'cpl']:
            for evolution in ['none', 'linear', 'smooth01', 'smooth03']:
                baseline = None
                for gpu in [False, True]:
                    settings = {'model': model, 'evolution': evolution, 'sample': 'dovekie',
                                'calibration': 'official_planck', 'fast_lensing': True, 'gpu': gpu}
                    low, high = audit.native_configurations(settings, design)
                    assert set(low['theory']) == {'camb'} == set(high['theory'])
                    assert canonical(low['params']) == canonical(high['params'])
                    assert canonical(low['likelihood']) == canonical(high['likelihood'])
                    differences = {k for k in low['theory']['camb']['extra_args']
                                   if low['theory']['camb']['extra_args'][k] != high['theory']['camb']['extra_args'][k]}
                    assert differences == set(design['numerical_controls'])
                    if baseline is None:
                        baseline = canonical(low)
                    else:
                        assert baseline == canonical(low), 'GPU must not alter a native-CAMB configuration.'
                    configurations += 1
        # Check the concrete native parameter object supports an independent copy
        # and the three numerical controls without invoking native calculations.
        parameters = camb.CAMBparams()
        for key in design['numerical_controls']:
            setattr(parameters.Accuracy, key, 2.)
        copied = parameters.copy()
        for key in design['numerical_controls']:
            setattr(copied.Accuracy, key, 1.)
        assert all(getattr(parameters.Accuracy, k) == 2. and getattr(copied.Accuracy, k) == 1. for k in design['numerical_controls'])
        # Instantiate the actual native pipeline, without evaluating it. Cobaya
        # inserts its native transfer helper only at this stage.
        from cobaya.model import get_model
        info = audit.native_configurations({'model': 'cpl', 'evolution': 'linear',
            'sample': 'dovekie', 'calibration': 'official_planck', 'fast_lensing': True}, design)[1]
        info['debug'] = 40
        with get_model(info) as model:
            model.add_requirements({'CAMBdata': None})
            native_components = list(model.theory)
            assert set(native_components) == {'camb', 'camb.transfers'}
        base = {'exact_loglikes': {'CMB': -123.4, 'BAO': -7.3, 'SN': -44.1},
                'exact_logpost': -123.4-7.3-44.1-2.1-1.9}
        def rows_for(deltas, cancellation=False):
            rows = []
            for i, delta in enumerate(deltas):
                high = dict(base['exact_loglikes'])
                high['CMB'] += float(delta)
                if cancellation:
                    high['BAO'] -= float(delta)
                prior = [-2.1, -1.9]
                value = audit.density_comparison(base, high, prior, sum(high.values())+sum(prior), design['diagnostic_thresholds'])
                assert not value['failed_density_checks']
                rows.append({'audit_index': i, 'status': 'finite_native_precision', 'chain': i//8,
                             'point': {'synthetic_x': i}, 'comparison': value})
            return rows
        constant = audit.summarize(rows_for(np.full(32, 7.25)), design)
        assert constant['status'] == 'no_large_variation_detected_on_fixed32'
        assert abs(constant['total_loglike_difference']['common_mean_offset']-7.25) < 1e-12
        assert constant['total_loglike_difference']['centered_RMS'] < 1e-12
        small = audit.summarize(rows_for(7.25+np.tile([-.02, .02], 16)), design)
        assert small['status'] == 'no_large_variation_detected_on_fixed32'
        assert abs(small['total_loglike_difference']['centered_RMS']-.02) < 1e-12
        large = audit.summarize(rows_for(np.tile([-.2, .2], 16)), design)
        assert large['status'] == 'numerical_sensitivity_requires_followup'
        assert 'total_centered_RMS' in large['diagnostic_flags']
        cancel = audit.summarize(rows_for(np.tile([-.3, .3], 16), True), design)
        assert cancel['status'] == 'numerical_sensitivity_requires_followup'
        assert cancel['total_loglike_difference']['centered_RMS'] < 1e-12
        assert {'component_variation:CMB', 'component_variation:BAO'} <= set(cancel['diagnostic_flags'])
        wrong_prior = audit.density_comparison(base, base['exact_loglikes'], [-2.1, -1.8],
                                               sum(base['exact_loglikes'].values())-3.9, design['diagnostic_thresholds'])
        assert wrong_prior['failed_density_checks'] == ['prior_sum_changed']
        failed = rows_for(np.zeros(32)); failed[6] = {'audit_index': 6, 'status': 'exception', 'error': 'synthetic_failure'}
        assert audit.summarize(failed, design)['status'] == 'incomplete_or_failed_numerical_screen'
        with tempfile.TemporaryDirectory() as directory:
            cache = Path(directory)/'cache.json'
            audit.sealed_write(cache, {'identity': 'synthetic', 'numerical': [1., 2.]})
            assert audit.sealed_read(cache)['numerical'] == [1., 2.]
            tampered = json.loads(cache.read_text()); tampered['numerical'][0] = 3.
            cache.write_text(json.dumps(tampered))
            try:
                audit.sealed_read(cache)
            except AssertionError:
                pass
            else:
                raise AssertionError('Cache numerical alteration was accepted.')
            import measurement_summary
            previous_gate = measurement_summary.summarize_run
            previous_factory = audit.factory
            def rejected_parent(*args):
                raise ValueError('synthetic_parent_not_qualified')
            def forbidden_factory(*args):
                raise AssertionError('Factory reached before parent gate.')
            try:
                measurement_summary.summarize_run = rejected_parent
                audit.factory = forbidden_factory
                child = Path(directory)/'must-not-create'
                try:
                    audit.prepare(Path(directory), Path(directory)/'missing-summary.json', child)
                except ValueError as error:
                    assert str(error) == 'synthetic_parent_not_qualified'
                else:
                    raise AssertionError('Unqualified parent was accepted.')
                assert not child.exists()
            finally:
                measurement_summary.summarize_run = previous_gate
                audit.factory = previous_factory
        assert calls['forbidden_attempts'] == 0
    finally:
        camb.get_results = original_results; camb.get_background = original_background
        camb.get_transfer_functions = original_transfers
    sources = [Path(__file__), Path(audit.__file__), audit.DESIGN,
               audit.HERE/'native_precision_audit.py', audit.HERE/'measurement_summary.py',
               audit.HERE/'modern_gpu.py', audit.HERE/'modern_fast.py', audit.HERE/'modern_run.py']
    report = {'status': 'passed_synthetic_posterior_precision_validation',
              'observational_points_used': 0, 'native_spectrum_calls': 0, 'native_background_calls': 0,
              'configuration_comparisons': configurations, 'deterministic_strata': expected,
              'selected_points': len(selected), 'chronological_permutation_invariance': True,
              'common_offset_not_mistaken_for_variable_error': True,
              'centered_RMS_analytic_absolute_error': abs(small['total_loglike_difference']['centered_RMS']-.02),
              'component_cancellation_flagged': True, 'changed_prior_rejected': True,
              'failed_point_withholds_screen': True, 'changed_cached_payload_rejected': True,
              'unqualified_parent_rejected_before_cache_or_factory': True,
              'copied_CAMB_accuracy_controls_independent': True,
              'initialized_native_theory_components_without_evaluation': native_components,
              'source_sha256': {audit.relative(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
              'scope': 'Pure arithmetic, selection, configuration and fail-closed tests. The32-point native run has not been performed.'}
    target = audit.ROOT/'studies/unified_cosmology/results/inference/native-posterior-precision-validation.json'
    target.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': report['status'], 'native_calls': 0, 'configurations': configurations}))


if __name__ == '__main__':
    main()
