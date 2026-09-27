"""Past expansion from qualified native correction points, using background only."""
import argparse
from contextlib import contextmanager
import hashlib
import importlib.metadata
import json
from functools import lru_cache
from pathlib import Path
import time

import numpy as np
from scipy.special import logsumexp
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'expansion-history-design.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def payload_digest(record):
    return hashlib.sha256(json.dumps(record, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def read_cached_background(path, expected_identity, index, native_hash):
    record = json.loads(Path(path).read_text())
    checksum = record.pop('payload_sha256')
    assert payload_digest(record) == checksum, 'Background cache payload changed.'
    record['payload_sha256'] = checksum
    assert record['identity'] == expected_identity and record['index'] == index
    assert record['native_record_sha256'] == native_hash
    return record


def verify_cache_manifest(path):
    if Path(path).exists():
        for name, expected in json.loads(Path(path).read_text()).items():
            assert digest(ROOT/name) == expected, 'Previously recorded background cache bytes changed.'


@contextmanager
def no_spectra():
    import camb
    previous = camb.get_results
    def forbidden(*args, **kwargs):
        raise RuntimeError('Expansion history forbids CMB spectrum evaluation.')
    camb.get_results = forbidden
    try:
        yield
    finally:
        camb.get_results = previous


def differentiate_hubble(hubble, z, step):
    """q=-a a''/a'^2 and j=a^2 a'''/a'^3 from H(z), including radiation."""
    z = np.asarray(z, dtype=float)
    assert np.all(z >= 0.) and step > 0.
    nodes = np.tile(np.arange(-2., 3.), (len(z), 1))
    nodes[z < 2*step] = np.arange(5.)
    positions = z[:, None] + step*nodes
    assert np.min(positions) >= 0.
    h = np.asarray(hubble(positions.reshape(-1))).reshape(positions.shape)
    centre = np.asarray(hubble(z))
    assert np.isfinite(h).all() and np.all(h > 0.)
    # Solving for local polynomial coefficients independently for each stencil
    # avoids applying a one-sided derivative where a central stencil is valid.
    coefficients = np.array([np.linalg.solve(np.vander(t, 5, increasing=True),
                             values/base) for t, values, base in zip(nodes, h, centre)])
    hp = coefficients[:, 1]/step
    hpp = 2*coefficients[:, 2]/step**2
    q = (1+z)*hp-1
    j = 1-2*(1+z)*hp+(1+z)**2*(hpp+hp**2)
    return {'H_km_s_Mpc': centre, 'q': q, 'j': j}


@lru_cache(maxsize=None)
def age_nodes(count):
    nodes, weights = np.polynomial.legendre.leggauss(count)
    return (nodes+1)/2, weights/2


def integrate_age(background):
    """Independent quadrature rules on the identical physical H(a), in Gyr."""
    from camb import constants
    conversion = constants.Mpc/1000/constants.Gyr
    integral, error = quad(lambda a: 1/(a*background.hubble_parameter(1/a-1)),
                           0., 1., epsabs=1e-12, epsrel=1e-10, limit=200)
    t, weights = age_nodes(256)
    gaussian = weights@(2/(t*background.hubble_parameter(1/t**2-1)))
    return float(conversion*integral), float(conversion*max(error, abs(integral-gaussian)))


def background_history(background, z, step):
    import camb.constants as constants
    history = differentiate_hubble(background.hubble_parameter, z, step)
    half = differentiate_hubble(background.hubble_parameter, z, step/2)
    a = 1/(1+np.asarray(z))
    density = background.get_background_densities(a)
    friedmann = np.sqrt(density['tot']/(3*a**4))*constants.c/1000
    derived = background.get_derived_params()
    default_age = float(background.physical_time(0.))
    age, age_error = integrate_age(background)
    scalar = {'age_Gyr': age, 'rdrag_Mpc': float(derived['rdrag']),
              'H0_km_s_Mpc': float(background.Params.H0),
              'H0_rdrag_km_s': float(background.Params.H0*derived['rdrag']),
              'omegam': float(background.Params.omegam)}
    checks = {
        'q_step_scaled_error': float(np.max(abs(history['q']-half['q'])/(1+abs(half['q'])))),
        'j_step_scaled_error': float(np.max(abs(history['j']-half['j'])/(1+abs(half['j'])))),
        'friedmann_relative_error': float(np.max(abs(friedmann/history['H_km_s_Mpc']-1))),
        'age_internal_absolute_Gyr': abs(default_age-float(derived['age'])),
        'age_integral_absolute_Gyr_error': age_error,
        'age_CAMB_default_minus_integral_Gyr': default_age-age}
    return {'history': {k: v.tolist() for k, v in history.items()},
            'scalar': scalar, 'numerical_checks': checks}


def native_thermal_parameters(parameters):
    """Match the full CAMB nonlinear-lensing thermal path, without spectra.

    CAMB_GetResults temporarily enables WantTransfer before InitVars, then
    restores its input value. InitVars uses that flag in its starting-time/grid
    choice, which also affects the finite-tolerance drag-redshift calculation.
    get_background does not perform this temporary switch by itself.
    """
    from camb import model
    adjusted = parameters.copy()
    if (adjusted.WantCls and adjusted.WantScalars
            and adjusted.NonLinear in {model.NonLinear_lens, model.NonLinear_both}
            and (adjusted.DoLensing or len(adjusted.SourceWindows) > 0)):
        adjusted.WantTransfer = True
    return adjusted


def evaluate_background(point, extra, design, native):
    import camb
    cosmology = {k: point[k] for k in ['H0', 'ombh2', 'omch2', 'ns', 'tau']}
    cosmology.update(As=1e-10*np.exp(point['logA']), w=point.get('w', -1.), wa=point.get('wa', 0.))
    with no_spectra():
        parameters = camb.set_params(**cosmology, **extra)
        adjusted = native_thermal_parameters(parameters)
        background = camb.get_background(adjusted)
        result = background_history(background, design['redshift_grid'], design['finite_difference_step'])
    result['thermal_adapter'] = {
        'input_WantTransfer': bool(parameters.WantTransfer),
        'background_WantTransfer': bool(adjusted.WantTransfer),
        'CMB_spectrum_calls': 0,
        'meaning': 'Reproduce the full native nonlinear-lensing thermal initialization; no transfer functions or spectra are computed.'}
    if (not all(np.isfinite(values).all() for values in result['history'].values())
            or not all(np.isfinite(value) for value in result['scalar'].values())
            or not all(np.isfinite(value) for value in result['numerical_checks'].values())):
        raise ArithmeticError('Nonfinite background output; no history is qualified.')
    z = np.array(design['redshift_grid'])
    checks = result['numerical_checks']
    for name in ['q', 'j']:
        comparison = []
        for redshift, label in [(0., '0'), (.5, '05'), (1., '1')]:
            index = int(np.flatnonzero(z == redshift)[0])
            reference = native[name+label]
            comparison.append(abs(result['history'][name][index]-reference)/(1+abs(reference)))
        checks[name+'_native_scaled_error'] = float(max(comparison))
    checks['rdrag_relative_error'] = abs(result['scalar']['rdrag_Mpc']/native['rdrag']-1)
    checks['omegam_absolute_error'] = abs(result['scalar']['omegam']-native['omegam'])
    gates = design['numerical_gates']
    failed = []
    for field, threshold in [
        ('q_step_scaled_error', 'q_step_and_native_scaled_error_max'),
        ('q_native_scaled_error', 'q_step_and_native_scaled_error_max'),
        ('j_step_scaled_error', 'j_step_and_native_scaled_error_max'),
        ('j_native_scaled_error', 'j_step_and_native_scaled_error_max'),
        ('rdrag_relative_error', 'rdrag_relative_error_max'),
        ('omegam_absolute_error', 'omegam_absolute_error_max'),
        ('friedmann_relative_error', 'friedmann_relative_error_max'),
        ('age_internal_absolute_Gyr', 'age_internal_absolute_Gyr_max'),
        ('age_integral_absolute_Gyr_error', 'age_integral_absolute_Gyr_error_max')]:
        if not np.isfinite(checks[field]) or checks[field] > gates[threshold]:
            failed.append(field)
    result['failed_numerical_gates'] = failed
    return result


def pointwise_summary(array, weights, probabilities):
    array = np.asarray(array)
    if array.ndim == 1:
        array = array[:, None]
    quantiles = []
    for column in array.T:
        order = np.argsort(column)
        mass = weights[order]
        quantiles.append(np.interp(probabilities, np.cumsum(mass)-mass/2, column[order]).tolist())
    return {'quantile_probabilities': probabilities, 'quantiles_by_redshift_or_scalar': quantiles,
            'mean': (weights@array).tolist(),
            'sd': np.sqrt(weights@((array-weights@array)**2)).tolist()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--chain-folder', type=Path, required=True)
    parser.add_argument('--correction-summary', type=Path, required=True)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    from measurement_summary import summarize_run, weighted_fraction
    from luminosity_sensitivity import weight_diagnostics
    # This is intentionally the first operation. There is no provisional bypass.
    parent = summarize_run(args.chain_folder, args.correction_summary)
    summary = json.loads(args.correction_summary.read_text())
    selection_path = ROOT/summary['selection_path']
    assert selection_path.parent.parent == args.chain_folder.resolve()
    selection = json.loads(selection_path.read_text())
    manifest = json.loads((args.chain_folder/'run-0.json').read_text())
    target = manifest['target_identity']
    for path, expected in target['source_sha256'].items():
        assert digest(ROOT/path) == expected, 'Scientific target source changed.'
    for name, version in target['versions'].items():
        assert importlib.metadata.version(name) == version, 'Scientific environment changed.'
    theory = target['configuration']['theory']
    assert len(theory) == 1
    extra = next(iter(theory.values()))['extra_args']
    assert extra['dark_energy_model'] == 'ppf' and extra['omk'] == 0.
    design = json.loads(DESIGN.read_text())
    sources = [Path(__file__), DESIGN, HERE/'measurement_summary.py', HERE/'exact_correction.py',
               HERE/'luminosity_sensitivity.py', HERE/'luminosity-sensitivity-design.json']
    hashes = {relative(p): digest(p) for p in sources}
    lineage = {'parent': parent['input_sha256'], 'source_sha256': hashes,
               'target_identity': target['identity'], 'versions': target['versions']}
    identity = hashlib.sha256(json.dumps(lineage, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    cache = args.cache.resolve()
    assert cache.is_relative_to(ROOT/'.work'), 'Use an ignored .work cache.'
    cache.mkdir(parents=True, exist_ok=True)
    lineage_path = cache/'lineage.json'
    if lineage_path.exists():
        assert json.loads(lineage_path.read_text()) == lineage, 'Cache belongs to different inputs or implementation.'
    else:
        lineage_path.write_text(json.dumps(lineage, indent=2)+'\n')
    cache_manifest = cache/'background-records.json'
    verify_cache_manifest(cache_manifest)
    records = []; cache_hashes = {}; started = time.monotonic()
    for index, point in enumerate(selection['points']):
        native_path = selection_path.parent/f'{index:05d}.json'
        native = json.loads(native_path.read_text())
        path = cache/f'{index:05d}.json'
        if path.exists():
            record = read_cached_background(path, identity, index, digest(native_path))
        else:
            try:
                record = evaluate_background(point, extra, design, native['derived'])
                record['status'] = 'finite_background' if not record['failed_numerical_gates'] else 'failed_numerical_gates'
            except Exception as error:
                record = {'status': 'exception', 'error': repr(error)}
            record.update(identity=identity, index=index, native_record_sha256=digest(native_path))
            record['payload_sha256'] = payload_digest(record)
            temporary = path.with_suffix('.part')
            temporary.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
            temporary.replace(path)
        records.append(record); cache_hashes[relative(path)] = digest(path)
        if (index+1) % 100 == 0:
            print(json.dumps({'background_points': index+1, 'seconds': time.monotonic()-started}), flush=True)
    cache_manifest.write_text(json.dumps(cache_hashes, indent=2)+'\n')
    failures = [{'index': r['index'], 'status': r['status'], 'gates': r.get('failed_numerical_gates'),
                 'error': r.get('error')} for r in records if r['status'] != 'finite_background']
    result = {'status': 'failed_background_checks' if failures else 'qualified_pointwise_expansion_history',
              'failures': failures, 'points': len(records), 'target_identity': target['identity'],
              'settings': parent['settings'], 'redshift': design['redshift_grid'],
              'correction_summary': relative(args.correction_summary),
              'correction_summary_sha256': digest(args.correction_summary),
              'lineage_path': relative(lineage_path), 'lineage_sha256': digest(lineage_path),
              'cache_manifest_path': relative(cache_manifest), 'cache_manifest_sha256': digest(cache_manifest),
              'source_sha256': hashes, 'identity': identity, 'CMB_spectrum_calls': 0}
    if not failures:
        native_records = [json.loads((selection_path.parent/f'{i:05d}.json').read_text()) for i in range(len(records))]
        lw = np.array([r['log_weight'] for r in native_records])
        weights = np.exp(lw-logsumexp(lw)); groups = np.array(selection['groups'])
        arrays = {name: np.array([r['history'][name] for r in records]) for name in ['H_km_s_Mpc', 'q', 'j']}
        scalars = {name: np.array([r['scalar'][name] for r in records]) for name in records[0]['scalar']}
        values = dict(scalars)
        for name, array in arrays.items():
            values.update({f'{name}_grid{i}': column for i, column in enumerate(array.T)})
        stability = weight_diagnostics(lw, values, groups, summary['weighted_stability_gates'])
        if not stability['qualified_overlap']:
            result['status'] = 'failed_history_weighted_stability'
        result['history_weighted_stability'] = stability
        result['history'] = {name: pointwise_summary(x, weights, design['quantile_probabilities']) for name, x in arrays.items()}
        result['scalar'] = {name: pointwise_summary(x, weights, design['quantile_probabilities']) for name, x in scalars.items()}
        result['acceleration_fraction_by_redshift'] = [weighted_fraction(x < 0., weights, groups) for x in arrays['q'].T]
        result['any_acceleration_on_reported_past_grid'] = weighted_fraction(np.any(arrays['q'] < 0., axis=1), weights, groups)
        result['numerical_check_maxima'] = {key: max(r['numerical_checks'][key] for r in records)
            for key in records[0]['numerical_checks'] if key != 'age_CAMB_default_minus_integral_Gyr'}
        age_difference = [r['numerical_checks']['age_CAMB_default_minus_integral_Gyr'] for r in records]
        result['age_CAMB_default_minus_integral_Gyr_range'] = [min(age_difference), max(age_difference)]
    result['interpretation'] = [
        'H is in km/s/Mpc; q<0 means scale-factor acceleration; j is the dimensionless third time derivative of scale factor, not dq/dt.',
        'Quantiles are pointwise equal-tail intervals conditional on the named model, priors and probe factorization, not simultaneous confidence bands.',
        'Any acceleration means at least one listed past redshift; no future extrapolation, global mode guarantee or Gaussian-sigma conversion.',
        'Cosmic age is time since the hot big bang under the shared CMB cosmological model, not an independent stellar-population age.',
        'Weights are the original untrimmed exact/proposal ratios; all failures are retained and withhold qualification.']
    assert all(digest(ROOT/path) == expected for path, expected in hashes.items()), 'Summary source changed during computation.'
    assert all(digest(ROOT/path) == expected for path, expected in target['source_sha256'].items()), 'Scientific source changed during computation.'
    assert all(digest(ROOT/path) == expected for path, expected in lineage['parent'].items()), 'Qualified parent inputs changed during computation.'
    assert all(importlib.metadata.version(name) == version for name, version in target['versions'].items()), 'Scientific environment changed during computation.'
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'points': len(records), 'failures': len(failures)}))


if __name__ == '__main__':
    main()
