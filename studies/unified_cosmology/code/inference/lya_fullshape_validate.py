"""Independent array-only validation of the rounded Ly-alpha Gaussian replacement.

Uses released measurements with synthetic positive distance predictions. No CAMB
background, spectra, cosmological model, sampler, or posterior is evaluated.
"""
import argparse
from contextlib import ExitStack
import hashlib
import importlib.metadata
import inspect
import itertools
import json
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.stats import multivariate_normal

import lya_fullshape as kernel

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
SOURCES = ROOT/'studies/unified_cosmology/results/external_probes/lya-fullshape-sources.json'
OUTPUT = ROOT/'studies/unified_cosmology/results/inference/lya-fullshape-validation.json'
SEED = 273692


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def forbidden(*args, **kwargs):
    raise AssertionError('Physical model/background/spectrum calls forbidden in this validation')


def independent_variants():
    # Transcribed independently from the versioned Eq26 and its printed precision.
    center = np.array([39.32, 8.600, .33, .066, .225])
    half = np.array([.005, .0005, .005, .0005, .0005])
    output = {'nominal': center}
    for integer in range(32):
        bits = f'{integer:05b}'
        offset = np.array([1 if digit == '1' else -1 for digit in bits])*half
        output['rounding_'+bits] = center+offset
    return output


def dense(mean, covariance, prediction):
    residual = prediction-mean
    return float(residual @ np.linalg.solve(covariance, residual))


def independent_replacement(model, parameters):
    # Name lookup is independent of both the kernel's stored indices and file order.
    dm = next(i for i, pair in enumerate(zip(model.z, model.row_types))
              if tuple(pair) == (2.33, 'DM_over_rs'))
    dh = next(i for i, pair in enumerate(zip(model.z, model.row_types))
              if tuple(pair) == (2.33, 'DH_over_rs'))
    m = model.mean.copy()
    c = model.C.copy()
    m[dm], m[dh] = parameters[:2]
    sdm, sdh, rho = parameters[2:]
    c[dm, dm] = sdm**2
    c[dh, dh] = sdh**2
    c[dm, dh] = c[dh, dm] = rho*sdm*sdh
    return m, c


def expect_refusal(function):
    try:
        function()
    except (AssertionError, ValueError, np.linalg.LinAlgError):
        return True
    raise AssertionError('Required invalid-input refusal did not occur')


def source_check():
    report = json.loads(SOURCES.read_text())
    assert report['status'] == 'versioned_sources_acquired_printed_Gaussian_is_explicit_approximation'
    np.testing.assert_array_equal(report['mean_sigma_rho'], independent_variants()['nominal'])
    # In particular 8.600 has a half-unit rounding interval of .0005, not .05.
    equation = report['printed_equation_26']
    for token in ['39.32', '8.600', '0.33', '0.066', '0.225']:
        assert token in equation
    identities = {}
    for item in report['resources'].values():
        assert sha(ROOT/item['path']) == item['sha256']
        identities[item['path']] = item['sha256']
    for mapping in ['source_sha256', 'input_sha256']:
        for name, digest in report[mapping].items():
            assert sha(ROOT/name) == digest
            identities[name] = digest
    identities[relative(SOURCES)] = sha(SOURCES)
    design = json.loads(kernel.DESIGN.read_text())
    np.testing.assert_array_equal(design['printed_summary']['mean_sigma_rho'], independent_variants()['nominal'])
    np.testing.assert_array_equal(design['printed_summary']['rounding_half_steps'], [.005, .0005, .005, .0005, .0005])
    return identities


class FakeProvider:
    """Positive analytic test functions; not a cosmological background."""
    def __init__(self, coefficients):
        self.coefficients = coefficients
        self.calls = []

    def dm(self, z):
        a, b, _, _, _ = self.coefficients
        return a*np.log1p(np.asarray(z))*(1+b*np.asarray(z))

    def hubble(self, z):
        _, _, h, linear, _ = self.coefficients
        z = np.asarray(z)
        return h*(1+linear*z+.09*z*z)

    def get_angular_diameter_distance(self, z):
        self.calls.append('DA')
        return np.atleast_1d(self.dm(z)/(1+np.asarray(z)))

    def get_Hubble(self, z, units):
        self.calls.append(units)
        h = self.hubble(z)
        assert units in ['km/s/Mpc', '1/Mpc']
        return np.atleast_1d(h if units == 'km/s/Mpc' else h/299792.458)

    def get_param(self, name):
        assert name == 'rdrag'
        self.calls.append(name)
        return self.coefficients[-1]


def native_interface(model, rng):
    from cobaya.likelihoods.base_classes.bao import BAO
    native = BAO({'path': str(kernel.RELEASE), 'measurements_file': kernel.MEAN.name,
                  'cov_file': kernel.COVARIANCE.name},
                 packages_path=str(kernel.RELEASE.parents[2]))
    np.testing.assert_array_equal(native.data['observable'].to_numpy(), model.row_types)
    np.testing.assert_array_equal(native.data['z'].to_numpy(), model.z)
    np.testing.assert_array_equal(native.data['value'].to_numpy(), model.mean)
    np.testing.assert_array_equal(native.cov, model.C)
    differences, prediction_errors, dv_errors = [], [], []
    for _ in range(24):
        p = FakeProvider([rng.uniform(3000, 4000), rng.uniform(.01, .3),
                          rng.uniform(55, 85), rng.uniform(.5, 1.1), rng.uniform(130, 160)])
        native.provider = p
        dm, h, rd = p.dm(model.z), p.hubble(model.z), p.coefficients[-1]
        values = model.evaluate(dm, h, rd)
        actual = float(native.logp())  # Actual installed code, no get_model/theory object.
        differences.append(abs(actual-values['old_loglike']))
        native_prediction = np.array([native.theory_fun(z, typ)[0]
                                      for z, typ in zip(model.z, model.row_types)])
        prediction_errors.append(float(np.max(abs(native_prediction-values['prediction']))))
        i = int(np.flatnonzero(model.row_types == 'DV_over_rs')[0])
        # Independent scalar DV expression and explicit angular-distance conversion.
        da = float(p.dm(model.z[i])/(1+model.z[i]))
        dv = (((1+model.z[i])*da)**2*299792.458*model.z[i]/h[i])**(1/3)/rd
        dv_errors.append(abs(values['prediction'][i]-dv))
        assert {'DA', 'km/s/Mpc', '1/Mpc', 'rdrag'} <= set(p.calls)
    assert max(differences) < 1e-9
    assert max(prediction_errors) < 1e-12
    assert max(dv_errors) < 1e-12
    return {'cases': 24, 'native_old_loglike_max_abs_error': max(differences),
            'prediction_max_abs_error': max(prediction_errors),
            'independent_DV_max_abs_error': max(dv_errors),
            'bao_source': relative(inspect.getfile(BAO)),
            'bao_source_sha256': sha(inspect.getfile(BAO)),
            'interface': 'Installed BAO.initialize/theory_fun/logp with fake DA/H/rdrag provider; no cosmological model.'}


def run():
    source_inputs = source_check()
    rng = np.random.default_rng(SEED)
    model = kernel.GaussianBAOReplacement.from_release()
    variants = independent_variants()
    assert list(variants) == list(kernel.variants())
    for name, parameters in variants.items():
        np.testing.assert_array_equal(parameters, kernel.variants()[name])
    assert len(variants) == 33
    assert list(model.row_types[-2:]) == ['DH_over_rs', 'DM_over_rs']
    assert not np.count_nonzero(model.C[np.ix_(model.lya, ~model.lya)])
    predictions = model.mean+rng.normal(size=(64, 13))@np.linalg.cholesky(model.C).T*1.7
    assert np.all(predictions > 0)
    errors = {'dense_chi2': 0., 'block_closure': 0., 'delta_loglike': 0.,
              'normalized_density_constant': 0., 'row_permutation': 0.}
    constants = {}
    for prediction in predictions:
        result = model.evaluate_prediction(prediction)
        old = dense(model.mean, model.C, prediction)
        errors['dense_chi2'] = max(errors['dense_chi2'], abs(old-result['old_chi2']))
        errors['block_closure'] = max(errors['block_closure'], abs(result['old_block_closure']))
        normalized_old = multivariate_normal.logpdf(prediction, model.mean, model.C)
        for name, parameters in variants.items():
            m, c = independent_replacement(model, parameters)
            expected = dense(m, c, prediction)
            actual = result['variants'][name]
            errors['dense_chi2'] = max(errors['dense_chi2'], abs(expected-actual['new_chi2']))
            errors['block_closure'] = max(errors['block_closure'], abs(actual['new_block_closure']))
            errors['delta_loglike'] = max(errors['delta_loglike'], abs(-.5*(expected-old)-actual['delta_loglike']))
            # A normalized full 13-dimensional Gaussian differs by a fixed logdet ratio.
            constant = -.5*(np.linalg.slogdet(c)[1]-np.linalg.slogdet(model.C)[1])
            normalized_delta = multivariate_normal.logpdf(prediction, m, c)-normalized_old
            errors['normalized_density_constant'] = max(errors['normalized_density_constant'],
                                                       abs(normalized_delta-actual['delta_loglike']-constant))
            constants[name] = float(constant)
    for _ in range(16):
        order = rng.permutation(13)
        changed = kernel.GaussianBAOReplacement(model.z[order], model.mean[order],
                                               model.row_types[order], model.C[np.ix_(order, order)])
        for prediction in predictions[:4]:
            first, second = model.evaluate_prediction(prediction), changed.evaluate_prediction(prediction[order])
            errors['row_permutation'] = max(errors['row_permutation'], abs(first['old_chi2']-second['old_chi2']))
            for name in variants:
                errors['row_permutation'] = max(errors['row_permutation'],
                    abs(first['variants'][name]['new_chi2']-second['variants'][name]['new_chi2']))
    # Wrongly treating the file's last pair as DM,DH must be detectable numerically.
    wrong_mean = model.mean.copy(); wrong_mean[-2:] = variants['nominal'][:2]
    wrong_cov = model.C.copy()
    sdm, sdh, rho = variants['nominal'][2:]
    wrong_cov[-2:, -2:] = [[sdm**2, rho*sdm*sdh], [rho*sdm*sdh, sdh**2]]
    correct_mean, correct_cov = independent_replacement(model, variants['nominal'])
    correct = dense(correct_mean, correct_cov, model.mean)
    wrong = dense(wrong_mean, wrong_cov, model.mean)
    assert abs(wrong-correct) > 100
    # Correlations inside the new pair must not be discarded.
    diagonal_cov = correct_cov.copy()
    i, j = model.order; diagonal_cov[i, j] = diagonal_cov[j, i] = 0
    pair_correlations_effect = max(abs(dense(correct_mean, diagonal_cov, p)-dense(correct_mean, correct_cov, p)) for p in predictions)
    assert pair_correlations_effect > .1
    refusals = {}
    for name, value in [('nonzero_cross_covariance', 1e-4), ('tiny_nonzero_cross_covariance', 1e-300)]:
        covariance = model.C.copy(); covariance[0, i] = covariance[i, 0] = value
        refusals[name] = expect_refusal(lambda: kernel.GaussianBAOReplacement(model.z, model.mean, model.row_types, covariance))
    for name, mutation in [('negative_prediction', lambda p: p.__setitem__(0, -1)),
                           ('nonfinite_prediction', lambda p: p.__setitem__(0, np.nan))]:
        p = predictions[0].copy(); mutation(p)
        refusals[name] = expect_refusal(lambda: model.evaluate_prediction(p))
    for name, rd in [('zero_rdrag', 0.), ('nonfinite_rdrag', np.nan)]:
        refusals[name] = expect_refusal(lambda: model.evaluate(np.ones(13), np.ones(13), rd))
    duplicate_types = model.row_types.copy(); duplicate_types[i] = 'DH_over_rs'
    refusals['duplicate_named_Lya_observable'] = expect_refusal(lambda: kernel.GaussianBAOReplacement(model.z, model.mean, duplicate_types, model.C))
    for value in errors.values():
        assert value < 1e-9, errors
    native = native_interface(model, rng)
    source_paths = [Path(__file__), Path(kernel.__file__), kernel.DESIGN]
    inputs = {relative(kernel.MEAN): sha(kernel.MEAN), relative(kernel.COVARIANCE): sha(kernel.COVARIANCE), **source_inputs}
    return {'schema': 'lya-fullshape-independent-kernel-validation-v1', 'status': 'passed',
            'seed': SEED, 'scope': 'Synthetic positive predictions against released data; numerical/interface validation only, no cosmological inference or Gaussian-tail validation.',
            'physical_model_calls': 0, 'background_calls': 0, 'spectrum_calls': 0,
            'rows': 13, 'retained_rows': 11, 'variant_count': 33,
            'synthetic_prediction_vectors': 64, 'dense_variant_checks': 64*33,
            'row_permutations': 16, 'error_maxima': errors,
            'normalized_full_density_minus_unnormalized_ratio_constants': constants,
            'ordering_control': {'released_pair_order': model.row_types[-2:].tolist(),
                                 'correct_named_chi2': correct, 'wrong_positional_chi2': wrong},
            'dropping_pair_correlation_max_chi2_change': pair_correlations_effect,
            'invalid_inputs_refused': refusals, 'installed_cobaya_interface': native,
            'source_sha256': {relative(p): sha(p) for p in source_paths},
            'input_sha256': inputs,
            'software': {name: importlib.metadata.version(name) for name in ['numpy', 'scipy', 'cobaya', 'pandas']}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    args = parser.parse_args()
    # Construction and physical solvers are prohibited, even if a future import changes.
    import camb
    import cobaya.model
    import cobaya
    with ExitStack() as stack:
        for module, name in [(camb, 'get_background'), (camb, 'get_results'),
                             (camb, 'get_transfer_functions'), (cobaya, 'get_model'),
                             (cobaya.model, 'get_model')]:
            stack.enter_context(patch.object(module, name, forbidden))
        result = run()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'output': str(args.output),
                      'sha256': sha(args.output), 'error_maxima': result['error_maxima']}))


if __name__ == '__main__':
    main()
