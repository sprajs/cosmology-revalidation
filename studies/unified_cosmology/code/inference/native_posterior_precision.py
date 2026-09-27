"""Native numerical precision on32 fixed, already qualified posterior points.

Preparation always invokes measurement_summary.summarize_run. No provisional
input bypass is provided. Workers evaluate native CAMB only, once per point.
"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import argparse
from concurrent.futures import ThreadPoolExecutor
import copy
import hashlib
import importlib.metadata
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
from scipy.special import logsumexp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'native-posterior-precision-design.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def identity(value):
    payload = {k: v for k, v in value.items() if k != 'payload_sha256'}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def sealed_write(path, value):
    value = dict(value, payload_sha256=identity(value))
    temporary = Path(path).with_suffix('.part')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def sealed_read(path):
    value = json.loads(Path(path).read_text())
    assert value['payload_sha256'] == identity(value), 'Cached numerical payload changed.'
    return value


def factory(settings):
    if settings.get('gpu'):
        from modern_gpu import configuration, identify
    elif settings.get('fast_lensing'):
        from modern_fast import configuration, identify
    else:
        from modern_run import configuration
        from target_identity import identify
    return configuration, identify


def native_configurations(settings, design):
    from target_identity import canonical
    configuration, _ = factory(settings)
    kwargs = {k: settings[k] for k in ['model', 'evolution', 'sample', 'calibration']}
    nominal = configuration(**kwargs)
    assert set(nominal['theory']) == {'camb'}, 'No spectral emulator may enter this check.'
    extra = nominal['theory']['camb']['extra_args']
    assert extra['lens_potential_accuracy'] == 4
    assert all(extra[k] == design['declared_value'] for k in design['numerical_controls'])
    doubled = copy.deepcopy(nominal)
    for key in design['numerical_controls']:
        doubled['theory']['camb']['extra_args'][key] = design['comparison_value']
    restored = copy.deepcopy(doubled)
    for key in design['numerical_controls']:
        restored['theory']['camb']['extra_args'][key] = extra[key]
    assert canonical(restored) == canonical(nominal)
    return nominal, doubled


def select_indices(groups, locations, points_per_chain=8):
    groups = np.asarray(groups)
    assert set(groups) == {0, 1, 2, 3}
    assert len(groups) == len(locations)
    output = []
    for chain in range(4):
        eligible = np.flatnonzero(groups == chain)
        assert len(eligible) >= points_per_chain
        ordered = sorted(eligible.tolist(), key=lambda i: locations[i]['expanded_index'])
        positions = np.floor((np.arange(points_per_chain)+.5)*len(ordered)/points_per_chain).astype(int)
        assert len(set(positions)) == points_per_chain
        output.extend({'parent_index': int(ordered[position]), 'chain': chain,
                       'stratum': j, 'within_chain_ordered_index': int(position),
                       'expanded_index': locations[ordered[position]]['expanded_index']}
                      for j, position in enumerate(positions))
    return output


def density_comparison(stored, high, prior_values, high_logpost, thresholds):
    low = stored['exact_loglikes']
    assert set(high) == set(low), 'Native likelihood component set changed.'
    delta = {k: float(high[k]-low[k]) for k in low}
    prior_inferred = float(stored['exact_logpost']-sum(low.values()))
    prior_difference = float(sum(prior_values)-prior_inferred)
    total = float(sum(delta.values()))
    closure = float(high_logpost-stored['exact_logpost']-total-prior_difference)
    assert all(np.isfinite(v) for v in [*delta.values(), prior_difference, closure])
    failed = []
    if abs(prior_difference) > thresholds['prior_sum_absolute_difference_max']:
        failed.append('prior_sum_changed')
    if abs(closure) > thresholds['logposterior_minus_loglike_difference_absolute_max']:
        failed.append('logdensity_accounting')
    return {'high_minus_declared_component_loglikes': delta,
            'high_minus_declared_total_loglike': total,
            'declared_prior_sum_inferred_from_stored_logpost': prior_inferred,
            'high_minus_declared_prior_sum': prior_difference,
            'logposterior_accounting_closure': closure, 'failed_density_checks': failed}


def spread(values):
    values = np.asarray(values, dtype=float)
    assert len(values) and np.isfinite(values).all()
    offset = float(np.mean(values)); centered = values-offset
    return {'common_mean_offset': offset, 'common_median_offset': float(np.median(values)),
            'centered_RMS': float(np.sqrt(np.mean(centered**2))),
            'centered_max_absolute': float(np.max(abs(centered))),
            'range': float(np.ptp(values)), 'minimum': float(values.min()), 'maximum': float(values.max())}


def summarize(rows, design):
    bad = [{'audit_index': r['audit_index'], 'status': r['status'], 'error': r.get('error'),
            'failed_checks': r.get('failed_checks', [])}
           for r in rows if r['status'] != 'finite_native_precision']
    if len(rows) != design['points'] or bad:
        return {'status': 'incomplete_or_failed_numerical_screen', 'failures': bad,
                'completed_records': len(rows), 'posterior_qualification': False}
    total = spread([r['comparison']['high_minus_declared_total_loglike'] for r in rows])
    keys = rows[0]['comparison']['high_minus_declared_component_loglikes']
    components = {k: spread([r['comparison']['high_minus_declared_component_loglikes'][k] for r in rows]) for k in keys}
    chains = {str(c): spread([r['comparison']['high_minus_declared_total_loglike'] for r in rows if r['chain'] == c]) for c in range(4)}
    limits = design['diagnostic_thresholds']; flags = []
    if total['centered_RMS'] > limits['total_centered_RMS_loglike_max']:
        flags.append('total_centered_RMS')
    if total['centered_max_absolute'] > limits['total_centered_max_absolute_loglike_max']:
        flags.append('total_centered_max_absolute')
    flags.extend('component_variation:'+k for k, s in components.items()
                 if s['centered_max_absolute'] > limits['component_centered_max_absolute_loglike_max'])
    return {'status': 'numerical_sensitivity_requires_followup' if flags else 'no_large_variation_detected_on_fixed32',
            'diagnostic_flags': flags, 'total_loglike_difference': total,
            'component_loglike_differences': components, 'per_chain_loglike_differences': chains,
            'posterior_qualification': False, 'posterior_reweighting_performed': False,
            'unique_physical_nuisance_points': len({identity(r['point']) for r in rows}),
            'interpretation': design['threshold_interpretation'], 'limitations': design['limitations']}


def verify_bindings(plan):
    for mapping in [plan['source_sha256'], plan['qualified_parent_inputs'], plan['frozen_target']['source_sha256']]:
        assert all(digest(ROOT/path) == value for path, value in mapping.items()), 'Frozen source/input changed.'
    for package, version in plan['frozen_target']['versions'].items():
        assert importlib.metadata.version(package) == version, 'Frozen numerical environment changed.'


def prepare(folder, summary_path, cache):
    from measurement_summary import summarize_run
    from late_geometry import sample_path
    from target_identity import canonical
    # This must precede cache creation or any native model construction.
    parent = summarize_run(folder, summary_path)
    summary = json.loads(summary_path.read_text())
    selection_path = ROOT/summary['selection_path']
    selection = json.loads(selection_path.read_text())
    settings = selection['settings']
    frozen = json.loads((folder/'run-0.json').read_text())['target_identity']
    design = json.loads(DESIGN.read_text())
    configuration, identify_target = factory(settings)
    kwargs = {k: settings[k] for k in ['model', 'evolution', 'sample', 'calibration']}
    proposal = configuration(**kwargs, surrogate=Path(settings['surrogate']))
    actual = identify_target(proposal, sample_path(settings['sample']), Path(settings['surrogate']))
    assert actual == frozen, 'Qualified parent scientific target identity changed.'
    nominal, doubled = native_configurations(settings, design)
    selected = select_indices(selection['groups'], selection['locations'], design['points_per_chain'])
    assert len(selected) == design['points']
    all_logweights = np.array([json.loads((selection_path.parent/f'{i:05d}.json').read_text())['log_weight']
                              for i in range(len(selection['points']))])
    full_weights = np.exp(all_logweights-logsumexp(all_logweights))
    for item in selected:
        index = item['parent_index']; path = selection_path.parent/f'{index:05d}.json'
        item.update(point=selection['points'][index], native_record_path=relative(path),
                    native_record_sha256=digest(path), original_location=selection['locations'][index],
                    original_native_logweight=float(all_logweights[index]),
                    weight_in_full_qualified_parent=float(full_weights[index]))
    sources = [Path(__file__), DESIGN, HERE/'measurement_summary.py', HERE/'exact_correction.py',
               HERE/'native_precision_audit.py', HERE/'likelihood.py', HERE/'target_identity.py']
    plan = {'design': design, 'settings': settings, 'selected': selected,
            'frozen_target': frozen, 'qualified_parent_inputs': parent['input_sha256'],
            'source_sha256': {relative(p): digest(p) for p in sources},
            'nominal_native_configuration': canonical(nominal),
            'doubled_native_configuration': canonical(doubled),
            'weight_context_only': {
                'selected_weight_mass_in_full_parent': float(sum(x['weight_in_full_qualified_parent'] for x in selected)),
                'full_parent_chain_weight_masses': {str(c): float(full_weights[np.asarray(selection['groups']) == c].sum()) for c in range(4)},
                'use': 'Context only. Selection and numerical-screen centering are unweighted; no posterior reweighting.'},
            'correction_summary_path': relative(summary_path), 'correction_summary_sha256': digest(summary_path)}
    plan['identity'] = identity(plan)
    verify_bindings(plan)
    cache.mkdir(parents=True, exist_ok=True)
    plan_path = cache/'selection.json'
    if plan_path.exists():
        old = sealed_read(plan_path)
        assert {k: v for k, v in old.items() if k != 'payload_sha256'} == plan
    else:
        sealed_write(plan_path, plan)
    return sealed_read(plan_path)


def worker(plan_path, index):
    from exact_correction import verify_record
    from target_identity import canonical
    from likelihood import expansion_diagnostics
    from cobaya.model import get_model
    import camb
    plan_path = Path(plan_path).resolve(); plan = sealed_read(plan_path)
    verify_bindings(plan)
    assert 0 <= index < plan['design']['points']
    item = plan['selected'][index]
    path = plan_path.parent/f'{index:02d}.json'
    assert not path.exists(), 'No silent repeat of a native precision point.'
    stored_path = ROOT/item['native_record_path']
    assert digest(stored_path) == item['native_record_sha256']
    stored = json.loads(stored_path.read_text()); verify_record(stored)
    assert stored['point'] == item['point'] and stored['status'] == 'finite'
    record = {'audit_index': index, 'parent_index': item['parent_index'], 'chain': item['chain'],
              'point': item['point'], 'identity': plan['identity'], 'plan_sha256': digest(plan_path),
              'native_record_sha256': item['native_record_sha256'], 'native_point_evaluations': 0,
              'original_native_logweight': item['original_native_logweight'],
              'weight_in_full_qualified_parent': item['weight_in_full_qualified_parent']}
    started = time.monotonic()
    try:
        nominal, doubled = native_configurations(plan['settings'], plan['design'])
        assert canonical(nominal) == plan['nominal_native_configuration']
        assert canonical(doubled) == plan['doubled_native_configuration']
        with get_model(doubled) as model:
            model.add_requirements({'CAMBdata': None})
            assert 'camb' in model.theory and set(model.theory) <= {'camb', 'camb.transfers'}, 'Only native CAMB and its transfer helper are permitted.'
            record['native_point_evaluations'] = 1
            result = model.logposterior(item['point'])
            assert np.isfinite(result.logpost), 'Nonfinite doubled-accuracy native density.'
            high = dict(zip(model.likelihood, map(float, result.loglikes)))
            comparison = density_comparison(stored, high, list(map(float, result.logpriors)),
                                            float(result.logpost), plan['design']['diagnostic_thresholds'])
            background = model.provider.get_CAMBdata()
            raw = model.provider.get_Cl(ell_factor=True)
            spectra_path = path.with_name(f'{index:02d}-spectra.npz')
            np.savez_compressed(spectra_path, ell=raw['ell'], **{k: raw[k] for k in ['tt', 'ee', 'bb', 'te', 'pp']})
            # Same finalized physical parameters, declared numerical background.
            baseline_parameters = background.Params.copy()
            for key in plan['design']['numerical_controls']:
                setattr(baseline_parameters.Accuracy, key, plan['design']['declared_value'])
            nominal_background = camb.get_background(baseline_parameters)
            grid = np.array([0., .5, 1., 2.33, 1100.])
            high_h = background.hubble_parameter(grid)
            low_h = nominal_background.hubble_parameter(grid)
            high_derived = background.get_derived_params(); low_derived = nominal_background.get_derived_params()
            check_grid = np.concatenate([np.arange(5)*.001+z for z in [0., .5, 1.]])
            low_expansion = expansion_diagnostics(check_grid, nominal_background.hubble_parameter(check_grid))
            gates = plan['design']['diagnostic_thresholds']
            closure = {'rdrag_relative': abs(low_derived['rdrag']/stored['derived']['rdrag']-1),
                       'omegam_absolute': abs(float(baseline_parameters.omegam)-stored['derived']['omegam']),
                       'q_scaled': max(abs(low_expansion[k]-stored['derived'][k])/(1+abs(stored['derived'][k])) for k in ['q0', 'q05', 'q1']),
                       'j_scaled': max(abs(low_expansion[k]-stored['derived'][k])/(1+abs(stored['derived'][k])) for k in ['j0', 'j05', 'j1'])}
            failed = list(comparison['failed_density_checks'])
            for name in closure:
                key = {'rdrag_relative': 'rdrag_relative', 'omegam_absolute': 'omegam_absolute', 'q_scaled': 'q_scaled', 'j_scaled': 'j_scaled'}[name]
                if not np.isfinite(closure[name]) or closure[name] > gates['nominal_background_'+key+'_closure_max']:
                    failed.append('nominal_background_'+name)
            record.update(status='finite_native_precision' if not failed else 'failed_native_precision_checks',
                failed_checks=failed, comparison=comparison, high_loglikes=high,
                high_logpriors=list(map(float, result.logpriors)), high_logpost=float(result.logpost),
                high_derived=dict(zip(model.parameterization.derived_params(), map(float, result.derived))),
                nominal_background_stored_native_closure=closure,
                background={'redshift': grid.tolist(), 'declared_H_km_s_Mpc': low_h.tolist(),
                    'high_H_km_s_Mpc': high_h.tolist(), 'H_relative_difference': (high_h/low_h-1).tolist(),
                    'thetaMC_times100_difference': float(100*(background.cosmomc_theta()-nominal_background.cosmomc_theta())),
                    'recombination_differences': {k: float(high_derived[k]-low_derived[k]) for k in ['age', 'zstar', 'thetastar', 'rstar', 'zdrag', 'rdrag']}},
                finalized_theory_extra_args=dict(model.theory['camb'].extra_args),
                CAMB_Params_max_l=int(background.Params.max_l), provider_Dl_length=len(raw['ell']),
                spectrum_path=relative(spectra_path), spectrum_sha256=digest(spectra_path),
                stored_declared_loglikes=stored['exact_loglikes'], stored_declared_logpost=stored['exact_logpost'])
            json.dumps(record, allow_nan=False)
    except Exception as error:
        # Keep the failed point without inventing a finite numerical comparison.
        record = {k: record[k] for k in ['audit_index', 'parent_index', 'chain', 'point', 'identity',
                   'plan_sha256', 'native_record_sha256', 'native_point_evaluations']}
        record.update(status='exception', error=repr(error))
    record['seconds'] = time.monotonic()-started
    verify_bindings(plan)
    sealed_write(path, record)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--chain-folder', type=Path); p.add_argument('--correction-summary', type=Path)
    p.add_argument('--cache', type=Path); p.add_argument('--output', type=Path)
    p.add_argument('--workers', type=int, default=1); p.add_argument('--plan-only', action='store_true')
    p.add_argument('--worker-plan', type=Path); p.add_argument('--worker-index', type=int)
    args = p.parse_args()
    if args.worker_plan:
        worker(args.worker_plan, args.worker_index); return
    if not all([args.chain_folder, args.correction_summary, args.cache, args.output]):
        p.error('Require qualified chain folder, correction summary, cache and output.')
    assert 1 <= args.workers <= 4
    cache = args.cache.resolve(); assert cache.is_relative_to(ROOT/'.work')
    plan = prepare(args.chain_folder.resolve(), args.correction_summary.resolve(), cache)
    plan_path = cache/'selection.json'
    if args.plan_only:
        print(json.dumps({'status': 'qualified_parent_plan_prepared', 'points': len(plan['selected']),
                          'native_point_evaluations': 0, 'plan_sha256': digest(plan_path)})); return
    manifest = cache/'record-hashes.json'
    if manifest.exists():
        for path, expected in json.loads(manifest.read_text()).items():
            assert digest(ROOT/path) == expected, 'Previously summarized precision cache changed.'
    def execute(index):
        path = cache/f'{index:02d}.json'; log = cache/f'{index:02d}.log'
        if not path.exists():
            if log.exists():
                return {'audit_index': index, 'status': 'incomplete_attempt_no_silent_retry', 'log_path': relative(log)}
            with log.open('w') as output:
                process = subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker-plan', str(plan_path),
                                          '--worker-index', str(index)], stdout=output, stderr=subprocess.STDOUT)
            if not path.exists():
                return {'audit_index': index, 'status': 'worker_failed_before_completed_record',
                        'returncode': process.returncode, 'log_path': relative(log)}
        row = sealed_read(path)
        assert row['identity'] == plan['identity'] and row['audit_index'] == index
        assert row['plan_sha256'] == digest(plan_path)
        assert row['native_record_sha256'] == plan['selected'][index]['native_record_sha256']
        assert row['point'] == plan['selected'][index]['point']
        if 'spectrum_path' in row:
            assert digest(ROOT/row['spectrum_path']) == row['spectrum_sha256']
        row['record_path'] = relative(path); row['record_sha256'] = digest(path)
        row['log_path'] = relative(log); row['log_sha256'] = digest(log)
        row['integrated_time_mismatch_warning_count'] = log.read_text().count('mismatch in integrated times')
        return row
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(execute, range(plan['design']['points'])))
    hashes = {}
    for row in rows:
        for key in ['record', 'log', 'spectrum']:
            path = row.get(key+'_path')
            if path:
                hashes[path] = digest(ROOT/path)
    manifest.write_text(json.dumps(hashes, indent=2)+'\n')
    result = summarize(rows, plan['design'])
    result.update(points=rows, plan_path=relative(plan_path), plan_sha256=digest(plan_path),
                  record_manifest_path=relative(manifest), record_manifest_sha256=digest(manifest),
                  source_sha256=plan['source_sha256'], frozen_target_identity=plan['frozen_target']['identity'],
                  declared_target_changed=False,
                  known_native_point_evaluations=sum(r.get('native_point_evaluations', 0) for r in rows),
                  incomplete_attempts_with_unknown_native_call_count=sum('native_point_evaluations' not in r for r in rows),
                  weight_context_only=plan['weight_context_only'],
                  numerical_controls=plan['design']['numerical_controls'])
    verify_bindings(plan)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'points': len(rows), 'posterior_qualification': False}))


if __name__ == '__main__':
    main()
