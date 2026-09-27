"""Recheck auxiliary backgrounds without repeating a frozen native screen."""
import argparse
import copy
import json
import os
from pathlib import Path
from unittest.mock import patch

import numpy as np

import native_posterior_precision as original
from expansion_history import native_thermal_parameters
from likelihood import expansion_diagnostics
from measurement_summary import summarize_run

ROOT = original.ROOT
HERE = Path(__file__).resolve().parent
DESIGN = HERE/'native-precision-thermal-review-design.json'


def backgrounds(point, extra, design):
    import camb
    physics = {k: point[k] for k in ['H0', 'ombh2', 'omch2', 'ns', 'tau']}
    physics.update(As=1e-10*np.exp(point['logA']), w=point.get('w', -1.), wa=point.get('wa', 0.))
    def forbidden(*args, **kwargs):
        raise RuntimeError('A background review must not calculate spectra or transfer functions.')
    result = {}
    with patch.object(camb, 'get_results', forbidden), patch.object(camb, 'get_transfer_functions', forbidden):
        for label, accuracy in [('nominal', design['declared_value']), ('high', design['comparison_value'])]:
            options = dict(extra)
            for key in design['numerical_controls']:
                options[key] = accuracy
            parameters = camb.set_params(**physics, **options)
            adjusted = native_thermal_parameters(parameters)
            background = camb.get_background(adjusted)
            z = np.concatenate([np.arange(5)*.001+v for v in [0., .5, 1.]])
            result[label] = {
                'rdrag': float(background.get_derived_params()['rdrag']),
                'omegam': float(adjusted.omegam),
                'expansion': expansion_diagnostics(z, background.hubble_parameter(z)),
                'H_km_s_Mpc': background.hubble_parameter([0., .5, 1., 2.33, 1100.]).tolist(),
                'input_WantTransfer': bool(parameters.WantTransfer),
                'background_WantTransfer': bool(adjusted.WantTransfer)}
    return result


def reviewed_row(row, stored, values, design, review_design):
    assert row['status'] in {'finite_native_precision', 'failed_native_precision_checks'}, 'Incomplete native attempts cannot be repaired.'
    assert row['native_point_evaluations'] == 1
    assert row['point'] == stored['point']
    low, high = values['nominal'], values['high']
    # Python max can hide a non-first NaN. Reject every referenced quantity
    # before reductions, rather than checking only the final maxima.
    checked = [*stored['derived'].values(), *row['high_derived'].values(),
               *row['background']['high_H_km_s_Mpc']]
    for value in [low, high]:
        checked += [value['rdrag'], value['omegam'], *value['expansion'].values(), *value['H_km_s_Mpc']]
    assert np.isfinite(np.asarray(checked, dtype=float)).all(), 'Nonfinite background/derived quantity.'
    closure = {
        'rdrag_relative': abs(low['rdrag']/stored['derived']['rdrag']-1),
        'omegam_absolute': abs(low['omegam']-stored['derived']['omegam']),
        'q_scaled': max(abs(low['expansion'][k]-stored['derived'][k])/(1+abs(stored['derived'][k])) for k in ['q0','q05','q1']),
        'j_scaled': max(abs(low['expansion'][k]-stored['derived'][k])/(1+abs(stored['derived'][k])) for k in ['j0','j05','j1'])}
    high_closure = {'H_relative': float(np.max(abs(np.asarray(high['H_km_s_Mpc'])/row['background']['high_H_km_s_Mpc']-1))),
                    'rdrag_relative': abs(high['rdrag']/row['high_derived']['rdrag']-1)}
    comparison = original.density_comparison(stored, row['high_loglikes'], row['high_logpriors'],
                                            row['high_logpost'], design['diagnostic_thresholds'])
    assert comparison == row['comparison'], 'Saved likelihood arithmetic does not replay.'
    permitted = {'nominal_background_'+k for k in closure}
    failures = [f for f in row['failed_checks'] if f not in permitted]
    for f in comparison['failed_density_checks']:
        if f not in failures:
            failures.append(f)
    for key, value in closure.items():
        if not np.isfinite(value) or value > design['diagnostic_thresholds']['nominal_background_'+key+'_closure_max']:
            failures.append('nominal_background_'+key)
    for key, value in high_closure.items():
        if not np.isfinite(value) or value > review_design['high_background_'+key+'_max']:
            failures.append('high_background_'+key)
    adjusted = copy.deepcopy(row)
    adjusted.update(status='finite_native_precision' if not failures else 'failed_native_precision_checks',
                    failed_checks=failures, nominal_background_stored_native_closure=closure)
    return adjusted, {'original_status': row['status'], 'original_failed_checks': row['failed_checks'],
                     'reviewed_status': adjusted['status'], 'reviewed_failed_checks': failures,
                     'original_background_closure': row['nominal_background_stored_native_closure'],
                     'reviewed_background_closure': closure, 'high_background_closure': high_closure,
                     'backgrounds': values, 'unchanged_native_comparison': comparison}


def actual(screen_path, output):
    assert not output.exists(), 'Preserve previous reviews; use a fresh output.'
    assert all(os.environ.get(k) == '1' for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','CLIPY_NOJAX'])
    screen = json.loads(screen_path.read_text())
    plan_path = ROOT/screen['plan_path']
    assert original.digest(plan_path) == screen['plan_sha256']
    plan = original.sealed_read(plan_path)
    correction = ROOT/plan['correction_summary_path']
    assert original.digest(correction) == plan['correction_summary_sha256']
    parent_selection = ROOT/json.loads(correction.read_text())['selection_path']
    parent = summarize_run(parent_selection.parent.parent, correction)
    original.verify_bindings(plan)
    manifest_path = ROOT/screen['record_manifest_path']
    assert original.digest(manifest_path) == screen['record_manifest_sha256']
    bindings = dict(parent['input_sha256'], **json.loads(manifest_path.read_text()))
    for p in [screen_path, plan_path, correction, manifest_path, Path(__file__), DESIGN, HERE/'expansion_history.py', HERE/'native_posterior_precision.py']:
        bindings[original.relative(p)] = original.digest(p)
    for path, checksum in bindings.items():
        assert original.digest(ROOT/path) == checksum
    assert len(screen['points']) == plan['design']['points'] == 32
    for index, row in enumerate(screen['points']):
        assert row['status'] in {'finite_native_precision', 'failed_native_precision_checks'}, 'Incomplete native attempts cannot be repaired.'
        assert row['audit_index'] == index and row['native_point_evaluations'] == 1
        assert all(k in row for k in ['point', 'high_loglikes', 'high_logpriors', 'high_logpost',
                   'high_derived', 'background', 'failed_checks', 'comparison',
                   'nominal_background_stored_native_closure', 'finalized_theory_extra_args'])
    rows, reviews = [], []
    review_design = json.loads(DESIGN.read_text())
    for index, row in enumerate(screen['points']):
        selected = plan['selected'][index]
        assert row['audit_index'] == index and row['point'] == selected['point']
        assert row['plan_sha256'] == original.digest(plan_path)
        assert row['native_record_sha256'] == selected['native_record_sha256']
        source = ROOT/row['record_path']
        assert original.digest(source) == row['record_sha256']
        sealed = original.sealed_read(source)
        assert all(row[k] == v for k, v in sealed.items())
        native_path = ROOT/selected['native_record_path']
        assert original.digest(native_path) == selected['native_record_sha256']
        stored = json.loads(native_path.read_text())
        values = backgrounds(row['point'], row['finalized_theory_extra_args'], plan['design'])
        adjusted, review = reviewed_row(row, stored, values, plan['design'], review_design)
        rows.append(adjusted)
        reviews.append(dict(review, audit_index=index, original_record_path=original.relative(source),
                            original_record_sha256=original.digest(source)))
    corrected = original.summarize(rows, plan['design'])
    for path, checksum in bindings.items():
        assert original.digest(ROOT/path) == checksum
    result = {'status': 'completed_auxiliary_thermal_path_review', 'original_screen_status': screen['status'],
              'reviewed_screen': corrected, 'point_reviews': reviews, 'input_source_sha256': bindings,
              'original_screen_path': original.relative(screen_path), 'original_screen_sha256': original.digest(screen_path),
              'CMB_spectrum_calls': 0, 'background_calls': 64, 'posterior_qualification': False,
              'interpretation': review_design['interpretation']}
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--screen', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    result = actual(args.screen.resolve(), args.output.resolve())
    print(json.dumps({'status': result['status'], 'reviewed_screen_status': result['reviewed_screen']['status']}))


if __name__ == '__main__':
    main()
