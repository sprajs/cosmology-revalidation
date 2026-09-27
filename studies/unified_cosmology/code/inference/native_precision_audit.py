"""Eight fixed-point native precision comparisons; never alter the active target."""
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
EXTERNAL = HERE.parent/'external_probes'
sys.path.insert(0, str(EXTERNAL))
WORK = ROOT/'.work/unified-cosmology/inference/native-precision'
DESIGN = HERE/'native-precision-design.json'
RESULT = ROOT/'studies/unified_cosmology/results/inference/native-precision-audit.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def freeze():
    reference_path = ROOT/'studies/unified_cosmology/results/external_probes/modern-validation.json'
    reference = json.loads(reference_path.read_text())['reference_point']
    original_design = EXTERNAL/'spectral-training-design.json'
    theta = json.loads(original_design.read_text())['centre'][0]
    points = [{'label': 'fixed_reference', 'point': reference, 'theta_target': theta,
               'reason': 'Existing fixed modern validation reference, not a best fit.',
               'input_sha256': {relative(reference_path): digest(reference_path), relative(original_design): digest(original_design)}}]
    for number in [34, 37, 90]:
        path = ROOT/f'.work/unified-cosmology/inference/broad-spectral/holdout/{number:04d}.json'
        record = json.loads(path.read_text()); array = ROOT/record['file']
        assert digest(array) == record['sha256']
        points.append({'label': f'broad_holdout_{number:04d}', 'point': dict(reference, **record['physical_point']),
                       'theta_target': record['requested_coordinates'][0],
                       'reason': 'Large original cubic holdout error near early-w boundary.' if number != 90 else 'Near-baseline support control.',
                       'original_native_loglikes': dict(zip(record['likelihood_names'], record['loglikes'])),
                       'input_sha256': {relative(path): digest(path), relative(array): digest(array)}})
    sources = [Path(__file__), EXTERNAL/'modern_adapter.py', EXTERNAL/'adapter.py', EXTERNAL/'fast_lensing.py']
    value = {'declared_before_new_native_calls': True, 'maximum_native_point_evaluations': 8,
             'workers': 1, 'OMP_threads': 1, 'points': points,
             'accuracy_levels': [{'label': 'declared', 'AccuracyBoost': 1, 'lAccuracyBoost': 1, 'lSampleBoost': 1},
                                 {'label': 'doubled', 'AccuracyBoost': 2, 'lAccuracyBoost': 2, 'lSampleBoost': 2}],
             'comparison': 'Identical physical parameters, priors and data. Only all three numerical accuracy controls double. Exact validated fast-lensing reassociation in both arms. No spectral interpolation used.',
             'warning_source': {'url': 'https://raw.githubusercontent.com/cmbant/CAMB/1.6.6/fortran/cmbmain.f90',
                                'lines': [1024, 1033], 'threshold': 'abs(tau-Transfer_Times) >5e-5 in conformal Mpc at final source step; warning may repeat for multiple wavenumbers.'},
             'source_sha256': {relative(p): digest(p) for p in sources},
             'scope': 'Bounded numerical sensitivity after seeing broad holdout errors. No change/tuning of heldout model, scientific prior, active target or sampling weights. Not a posterior convergence study.'}
    if DESIGN.exists():
        assert json.loads(DESIGN.read_text()) == value
    else:
        DESIGN.write_text(json.dumps(value, indent=2)+'\n')
    WORK.mkdir(parents=True, exist_ok=True)
    return value


def worker(index):
    import camb
    from cobaya.model import get_model
    from modern_adapter import modern_info
    from fast_lensing import use_fast_lensing
    design = json.loads(DESIGN.read_text())
    for path, expected in design['source_sha256'].items():
        assert digest(ROOT/path) == expected
    item = design['points'][index//2]; accuracy = design['accuracy_levels'][index % 2]
    for path, expected in item['input_sha256'].items():
        assert digest(ROOT/path) == expected
    destination = WORK/f'{index:02d}.json'
    assert not destination.exists(), 'Do not silently repeat a native point.'
    info = use_fast_lensing(modern_info())
    info['theory']['camb']['extra_args'].update({k: accuracy[k] for k in ['AccuracyBoost', 'lAccuracyBoost', 'lSampleBoost']})
    started = time.monotonic()
    with get_model(info) as model:
        model.add_requirements({'CAMBdata': None})
        t = time.monotonic(); result = model.logposterior(item['point']); native_seconds = time.monotonic()-t
        assert np.isfinite(result.logpost)
        background = model.provider.get_CAMBdata()
        dls = model.provider.get_Cl(ell_factor=True)
        spectrum = WORK/f'{index:02d}-spectra.npz'
        np.savez_compressed(spectrum, ell=dls['ell'], **{k: dls[k] for k in ['tt', 'ee', 'bb', 'te', 'pp']})
        extra = dict(model.theory['camb'].extra_args)
        extra['halofit_version'] = 'mead2016'
        other = {k: item['point'][k] for k in ['ombh2', 'omch2', 'ns', 'tau', 'w', 'wa']}
        # Pure background inversion at the same requested acoustic coordinate.
        inverted = camb.set_params(cosmomc_theta=item['theta_target']/100,
                                   As=1e-10*np.exp(item['point']['logA']), **other, **extra)
        inverse_background = camb.get_background(inverted)
        z = np.array([0., .5, 1., 2.33, 1100.])
        derived = background.get_derived_params()
        record = {'status': 'finite_native', 'index': index, 'label': item['label'], 'accuracy': accuracy,
                  'point': item['point'], 'native_calls': 1, 'native_seconds': native_seconds,
                  'total_seconds': time.monotonic()-started,
                  'logposterior': float(result.logpost), 'logpriors': list(map(float, result.logpriors)),
                  'loglikes': dict(zip(model.likelihood, map(float, result.loglikes))),
                  'H_redshift': z.tolist(), 'H_km_s_Mpc': background.hubble_parameter(z).tolist(),
                  'thetaMC_times100': float(100*background.cosmomc_theta()),
                  'derived': {key: float(derived[key]) for key in ['age', 'zstar', 'thetastar', 'rstar', 'zdrag', 'rdrag']},
                  'inversion_at_fixed_requested_theta': {'target': item['theta_target'], 'H0': float(inverted.H0),
                      'actual_thetaMC_times100': float(100*inverse_background.cosmomc_theta())},
                  'finalized_theory_extra_args': dict(model.theory['camb'].extra_args),
                  'CAMB_Params_max_l': int(background.Params.max_l), 'provider_Dl_length': len(dls['ell']),
                  'spectrum_file': relative(spectrum), 'spectrum_sha256': digest(spectrum),
                  'design_sha256': digest(DESIGN), 'source_sha256': design['source_sha256'],
                  'camb_version': camb.__version__}
        destination.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('--freeze', action='store_true')
    parser.add_argument('--worker', type=int); args = parser.parse_args()
    if args.worker is not None:
        assert 0 <= args.worker < 8
        worker(args.worker); return
    design = freeze()
    if args.freeze:
        print(json.dumps({'design_sha256': digest(DESIGN), 'planned_native_calls': 8})); return
    rows = []
    for index in range(8):
        path = WORK/f'{index:02d}.json'; log = WORK/f'{index:02d}.log'
        if not path.exists():
            # A separate process flushes native Fortran output, attributing each
            # warning to one fixed point instead of an interleaved worker log.
            assert not log.exists(), 'A previous native attempt lacks a completed record; inspect it, do not silently retry.'
            with log.open('w') as output:
                subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker', str(index)],
                               stdout=output, stderr=subprocess.STDOUT, check=True)
        row = json.loads(path.read_text())
        assert row['design_sha256'] == digest(DESIGN)
        assert digest(ROOT/row['spectrum_file']) == row['spectrum_sha256']
        row['mismatch_warning_count'] = log.read_text().count('mismatch in integrated times')
        row['record_sha256'] = digest(path); row['log_sha256'] = digest(log)
        rows.append(row)
        print(json.dumps({'completed_native_points': index+1, 'label': row['label'],
                          'accuracy': row['accuracy']['label'], 'warnings': row['mismatch_warning_count']}), flush=True)
    comparisons = []
    for pair in range(4):
        low, high = rows[2*pair:2*pair+2]
        assert low['point'] == high['point'] and low['logpriors'] == high['logpriors']
        with np.load(ROOT/low['spectrum_file']) as first, np.load(ROOT/high['spectrum_file']) as second:
            spectra = {}
            for name in ['tt', 'ee', 'bb', 'te', 'pp']:
                a, b = first[name][2:], second[name][2:]
                scale = np.sqrt(abs(first['tt'][2:]*first['ee'][2:])) if name == 'te' else abs(a)
                scale = np.maximum(scale, 1e-6*np.max(scale))
                fractional = abs(b-a)/scale
                spectra[name] = {'absolute_rms': float(np.sqrt(np.mean((b-a)**2))),
                                 'scaled_absolute_error_quantiles': np.quantile(fractional, [.5, .95, 1.]).tolist(),
                                 'scale': 'sqrt(abs(TT*EE)), floored at1e-6 of its peak' if name == 'te' else 'abs(nominal spectrum), floored at1e-6 of its peak'}
        likelihood = {key: high['loglikes'][key]-value for key, value in low['loglikes'].items()}
        prior_native = design['points'][pair].get('original_native_loglikes')
        comparisons.append({'label': low['label'], 'w_plus_wa': low['point']['w']+low['point']['wa'],
            'high_minus_declared_loglike_by_component': likelihood,
            'high_minus_declared_total_loglike': sum(likelihood.values()),
            'original_nominal_reproduction_max_loglike_error': None if prior_native is None else max(abs(low['loglikes'][k]-v) for k, v in prior_native.items()),
            'H_relative_difference': (np.array(high['H_km_s_Mpc'])/low['H_km_s_Mpc']-1).tolist(),
            'thetaMC_times100_difference': high['thetaMC_times100']-low['thetaMC_times100'],
            'fixed_theta_H0_difference': high['inversion_at_fixed_requested_theta']['H0']-low['inversion_at_fixed_requested_theta']['H0'],
            'derived_differences': {k: high['derived'][k]-value for k, value in low['derived'].items()},
            'warning_counts': [low['mismatch_warning_count'], high['mismatch_warning_count']], 'spectra': spectra})
    report = {'status': 'completed_bounded_numerical_sensitivity', 'native_point_evaluations': 8,
              'active_target_changed': False, 'design_sha256': digest(DESIGN), 'source_sha256': design['source_sha256'],
              'points': rows, 'comparisons': comparisons,
              'scope': design['scope'],
              'warning_interpretation': 'Native final perturbation-source time differs from the separately tabulated transfer time by more than5e-5 conformal Mpc. Counts can repeat across wavenumbers and are not counts of bad cosmologies. Printed warnings provide no mismatch amplitude or automatic likelihood-validity verdict.'}
    RESULT.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
