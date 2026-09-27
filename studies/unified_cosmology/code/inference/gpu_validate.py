"""CPU/GPU proposal agreement, including every joint-likelihood contribution."""
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from cobaya.model import get_model
from modern_fast import configuration as cpu_configuration
from modern_gpu import configuration as gpu_configuration, GPUMatrix
from spectral_surrogate import polynomial
from target_identity import ROOT, canonical


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    folder = ROOT/'.work/unified-cosmology/inference/quartic-spectral'
    model_path = folder/'surrogate-quartic-v1.npz'
    model = np.load(model_path, allow_pickle=False)
    inputs = {str(model_path.relative_to(ROOT)): digest(model_path)}
    points = []
    # Fixed index ordering and a coordinate-only interior criterion avoid native
    # fallbacks in this bounded hardware comparison. It is not a data selection.
    for path in sorted((folder/'holdout_original').glob('*.json')):
        record = json.loads(path.read_text())
        if record['status'] != 'finite_exact':
            continue
        x = np.linalg.solve(model['coordinate_cholesky'],
                            np.array(record['actual_coordinates'])-model['centre'])
        if np.max(abs(x)) < 2.5:
            points.append(record['physical_point'])
            inputs[str(path.relative_to(ROOT))] = digest(path)
        if len(points) == 8:
            break
    assert len(points) == 8
    rng = np.random.default_rng(272899)
    features = polynomial(rng.normal(size=(100, 8)), model['exponents'])
    coefficients = model['coefficients']
    gpu = GPUMatrix(coefficients)
    _ = features[0] @ gpu
    started = time.perf_counter()
    cpu_values = np.array([x @ coefficients for x in features])
    cpu_seconds = time.perf_counter()-started
    started = time.perf_counter()
    gpu_values = np.array([x @ gpu for x in features])
    gpu_seconds = time.perf_counter()-started
    scaled_error = float(np.max(abs(cpu_values-gpu_values)/(1+abs(cpu_values))))
    assert scaled_error < 1e-11
    rows = []
    for evolution in ['none', 'linear', 'smooth01']:
        args = ('cpl', evolution, 'dovekie', 'official_planck', model_path)
        cpu_info = cpu_configuration(*args)
        gpu_info = gpu_configuration(*args)
        old_external = gpu_info['theory']['spectral_surrogate']['external']
        gpu_info['theory']['spectral_surrogate']['external'] = cpu_info['theory']['spectral_surrogate']['external']
        assert canonical(cpu_info) == canonical(gpu_info)
        gpu_info['theory']['spectral_surrogate']['external'] = old_external
        cpu = get_model(cpu_info)
        accelerated = get_model(gpu_info)
        names = cpu.parameterization.sampled_params()
        reference = {name: (definition['ref']['loc'] if isinstance(definition.get('ref'), dict)
                            else definition['ref']) for name, definition in cpu_info['params'].items()
                     if name in names}
        for index, physical in enumerate(points):
            point = dict(reference, **physical)
            if evolution == 'linear':
                point['epsilon'] = float(np.linspace(-.03, .03, len(points))[index])
            left = cpu.logposterior(point)
            right = accelerated.logposterior(point)
            assert np.isfinite(left.logpost) and np.isfinite(right.logpost)
            errors = {'logpost': abs(float(left.logpost-right.logpost)),
                      'loglikes': float(np.max(abs(np.array(left.loglikes)-right.loglikes))),
                      'logpriors': float(np.max(abs(np.array(left.logpriors)-right.logpriors))),
                      'derived': float(np.max(abs(np.array(left.derived)-right.derived)))}
            assert max(errors.values()) < 1e-6, errors
            assert cpu.theory['spectral_surrogate'].exact_calls == 0
            assert accelerated.theory['spectral_surrogate'].exact_calls == 0
            rows.append({'evolution': evolution, 'point': point, 'maximum_errors': errors})
        cpu.close(); accelerated.close()
    source = [Path(__file__), Path(__file__).with_name('modern_gpu.py'),
              Path(__file__).with_name('spectral_surrogate.py')]
    result = {'status': 'passed_numerical_equivalence', 'cases': rows,
              'component_benchmark': {'cases': 100, 'dtype': 'float64',
                  'cpu_seconds': cpu_seconds, 'gpu_seconds_including_transfers': gpu_seconds,
                  'speedup': cpu_seconds/gpu_seconds, 'relative_to_one_plus_abs_error': scaled_error},
              'input_sha256': inputs,
              'source_sha256': {str(path.relative_to(ROOT)): digest(path) for path in source},
              'scope': 'Three unchanged scientific targets at eight fixed interior coordinates. '
                       'No new native spectra, posterior qualification, full-chain speed claim '
                       'or validation of interpolation accuracy. Native correction remains mandatory.'}
    destination = ROOT/'studies/unified_cosmology/results/inference/gpu-validation.json'
    destination.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'cases': len(rows),
                      'component_benchmark': result['component_benchmark']}, indent=2))


if __name__ == '__main__':
    main()
