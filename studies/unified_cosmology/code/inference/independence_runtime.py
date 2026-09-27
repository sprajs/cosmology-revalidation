"""Explicit one-time native CAMB initialization, preserving the scientific target."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
from independence_proposal import ROOT, digest, relative

HERE = Path(__file__).resolve().parent
# A fixed numerical initialization point, never a prior or an observational fit.
WARMUP = {'H0': 67.36, 'ombh2': .02237, 'omch2': .12, 'As': 2.1e-9,
          'ns': .9649, 'tau': .0544, 'w': -1., 'wa': 0.}


def numerical_runtime():
    jax = sys.modules.get('jax')
    return {'coordinate_numpy_dtype': str(np.dtype(float)),
            'jax_loaded': jax is not None,
            'jax_enable_x64': bool(jax.config.read('jax_enable_x64')) if jax is not None else None,
            'third_party_internal_dtypes': 'Not inferred from Python/NumPy return types; no dtype setting is changed.',
            'environment': {key: os.environ.get(key) for key in
                            ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'CLIPY_NOJAX', 'JAX_ENABLE_X64']}}


def native_initialize(info):
    """Initialize CAMB once before any sampled density, consuming no RNG draws."""
    assert os.environ.get('CLIPY_NOJAX') == '1'
    assert 'spectral_surrogate' in info['theory']
    extra = dict(info['theory']['spectral_surrogate']['extra_args'])
    # Matches the actual modern target's native spectrum request. This is only
    # an initialization calculation; the frozen theory class is not modified.
    extra['lmax'] = max(extra.get('lmax', 0), 9001)
    import camb
    started = time.perf_counter()
    pars = camb.set_params(**WARMUP, **extra)
    full = camb.get_results(pars)
    result = {'method': 'one explicit camb.get_results before sampled density evaluation',
              'fixed_physical_parameters': WARMUP, 'native_extra_args': extra,
              'rdrag_at_initialization_point_Mpc': float(full.get_derived_params()['rdrag']),
              'seconds': time.perf_counter()-started, 'source_sha256': digest(__file__),
              'runtime_after_initialization': numerical_runtime(),
              'interpretation': 'Numerical library initialization only; no prior, target factor or sampled point is changed.'}
    return result


def validate(pilot_folder, output):
    from modern_fast import configuration, identify
    from measurement_summary import verify_current_target
    from late_geometry import sample_path
    from cobaya.model import get_model
    pilot_folder, output = Path(pilot_folder).resolve(), Path(output).resolve()
    assert not output.exists(); output.mkdir(parents=True)
    manifest = json.loads((pilot_folder/'run-0.json').read_text())
    verify_current_target(manifest)
    records = [json.loads(line) for line in (pilot_folder/'pilot-records.jsonl').read_text().splitlines()]
    # The selection is identical to the previously failed cold closure; no tuning.
    chosen = []
    for label in ['Gaussian', 't5']:
        candidates = [r for r in records if r['finite_target'] and r['component'] == label]
        chosen += sorted(candidates, key=lambda r: hashlib.sha256(str(r['index']).encode()).hexdigest())[:4]
    design = {'selection': 'same four finite/component SHA256(index) points as failed cold closure',
              'indices': [r['index'] for r in chosen], 'orders': ['forward', 'reverse'],
              'fixed_native_initialization': WARMUP, 'absolute_density_tolerance': 1e-8,
              'source_sha256': digest(__file__), 'pilot_records_sha256': digest(pilot_folder/'pilot-records.jsonl')}
    design_path = output/'design.json'; design_path.write_text(json.dumps(design, indent=2)+'\n')
    settings = manifest['arguments']
    info = configuration(settings['model'], settings['evolution'], settings['sample'],
                         settings['calibration'], Path(settings['surrogate']))
    assert identify(info, sample_path(settings['sample']), Path(settings['surrogate'])) == manifest['target_identity']
    comparisons = []; initial_cold = None
    with get_model(info) as model:
        p = chosen[0]
        cold = model.logposterior(p['point'], cached=False)
        initial_cold = {'index': p['index'], 'rdrag': float(model.provider.get_param('rdrag')),
                        'target_logpost': float(cold.logpost),
                        'difference_from_warm_pilot': float(cold.logpost-p['target_logpost'])}
        initialization = native_initialize(info)
        for order, cohort in [('forward', chosen), ('reverse', chosen[::-1])]:
            for record in cohort:
                value = model.logposterior(record['point'], cached=False)
                comparisons.append({'order': order, 'index': record['index'],
                    'target_difference_from_warm_pilot': float(value.logpost-record['target_logpost']),
                    'component_differences': dict(zip(model.likelihood,
                       map(float, np.asarray(value.loglikes)-record['loglikes']))),
                    'rdrag': float(model.provider.get_param('rdrag')),
                    'spectrum_hashes': {key: hashlib.sha256(np.asarray(values).tobytes()).hexdigest()
                       for key, values in model.theory['spectral_surrogate'].get_Cl(ell_factor=True).items() if key in ['tt', 'te', 'ee', 'pp']}})
    maximum = max(abs(c['target_difference_from_warm_pilot']) for c in comparisons)
    component_max = max(abs(v) for c in comparisons for v in c['component_differences'].values())
    for a in comparisons[:8]:
        b = next(r for r in comparisons[8:] if r['index'] == a['index'])
        assert {k:v for k,v in a.items() if k != 'order'} == {k:v for k,v in b.items() if k != 'order'}
    result = {'status': 'passed_declared_initialized_runtime' if max(maximum, component_max) <= 1e-8 else 'failed_initialized_runtime',
              'initial_cold_difference_retained': initial_cold, 'native_initialization': initialization,
              'comparisons': comparisons, 'maximum_absolute_logpost_difference': maximum,
              'maximum_component_difference': component_max, 'reversed_order_exactly_identical': True,
              'full_native_initialization_calls': 1, 'surrogate_target_calls': 17,
              'design_sha256': digest(design_path), 'source_sha256': digest(__file__),
              'runtime_after_target_evaluations': numerical_runtime(),
              'scope': 'Numerical initial-state contract on eight fixed points, not a posterior-convergence or universal precision guarantee.'}
    (output/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    assert result['status'] == 'passed_declared_initialized_runtime', result
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('pilot_folder', type=Path); p.add_argument('output', type=Path)
    a = p.parse_args(); r = validate(a.pilot_folder, a.output)
    print(json.dumps({k:r[k] for k in ['status', 'maximum_absolute_logpost_difference', 'maximum_component_difference', 'runtime_after_target_evaluations']}, indent=2))
