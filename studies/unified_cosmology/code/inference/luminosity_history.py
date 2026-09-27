"""Recover conditional grey-luminosity histories from qualified joint fits.

Gaussian spline coefficients were integrated out of the sampling density. Their
conditional Gaussian distribution is recovered analytically, retaining its
variance as well as variation between native-corrected cosmological points.
Nothing here identifies stellar age as the physical cause of luminosity drift.
"""
import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.linalg import cho_solve
from scipy.optimize import brentq
from scipy.special import logsumexp, ndtr
from scipy.stats import norm

from measurement_summary import ROOT, summarize_run
from luminosity_sensitivity import IntegratedLuminosity, background_only
from exact_correction import weighted_summary
from late_geometry import sample_path

PROBABILITIES = [.025, .16, .5, .84, .975]
REDSHIFTS = np.array([0., .05, .1, .2, .4, .6, .8, 1.])


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def gaussian_mixture_summary(means, variance, weights):
    means = np.asarray(means)
    assert variance >= 0 and np.isfinite(means).all()
    if variance == 0:
        return weighted_summary(means, weights)
    sigma = np.sqrt(variance)
    average = float(weights @ means)
    total_variance = float(variance + weights @ ((means-average)**2))
    lower, upper = float(means.min()-12*sigma), float(means.max()+12*sigma)
    quantiles = [brentq(lambda x: float(weights @ ndtr((x-means)/sigma))-p,
                       lower, upper, xtol=1e-12) for p in PROBABILITIES]
    return {'mean': average, 'sd': float(np.sqrt(total_variance)),
            'quantiles_025_16_50_84_975': quantiles}


def validate():
    weights = np.array([.15, .35, .5])
    equal = gaussian_mixture_summary(np.full(3, .123), .027**2, weights)
    error = float(np.max(abs(np.array(equal['quantiles_025_16_50_84_975'])
                            - (.123 + .027*norm.ppf(PROBABILITIES)))))
    assert error < 1e-10 and abs(equal['sd']-.027) < 1e-14
    means = np.array([-.2, .05, .3])
    variance = .04**2
    mixture = gaussian_mixture_summary(means, variance, weights)
    # Independent inverse-CDF numerical integration through a fine density grid.
    grid = np.linspace(-.7, .8, 200001)
    density = sum(w*norm.pdf(grid, mu, np.sqrt(variance))
                  for w, mu in zip(weights, means))
    cumulative = np.cumsum((density[1:]+density[:-1])*.5*np.diff(grid))
    quantiles = np.interp(PROBABILITIES, cumulative, grid[1:])
    mixture_error = float(max(abs(quantiles-mixture['quantiles_025_16_50_84_975'])))
    assert mixture_error < 1e-7
    assert abs(mixture['mean'] - weights @ means) < 1e-14
    assert abs(mixture['sd']**2 - variance - weights @ ((means-weights@means)**2)) < 1e-14
    # Independent full GLS block inversion includes the unpenalized magnitude
    # intercept explicitly, instead of using the production whitened projection.
    rng = np.random.default_rng(272809)
    gls_rows = []
    for case in range(16):
        n = 24 + case
        z = np.sort(rng.uniform(.01, 1.15, n))
        factor = rng.normal(size=(n, n))
        covariance = (factor @ factor.T/n + np.eye(n)) * .1**2
        residual = rng.normal(size=n) * .2
        sigma = [.1, .3][case % 2]
        sn = IntegratedLuminosity(z, z, -residual, covariance)
        label = 'smooth01' if sigma == .1 else 'smooth03'
        conditional_cov = cho_solve(sn.gaussian[label][0], np.eye(4))
        conditional_mean = -conditional_cov @ sn.spline.T @ sn.project(residual)
        design = np.column_stack([np.ones(n), sn.basis])
        precision = np.linalg.inv(covariance)
        joint_precision = design.T @ precision @ design + np.diag([0.] + [1/sigma**2]*4)
        joint_cov = np.linalg.inv(joint_precision)
        joint_mean = -joint_cov @ design.T @ precision @ residual
        shifted = -conditional_cov @ sn.spline.T @ sn.project(residual + 53.)

        def direct_loglike(c):
            p = np.linalg.inv(c)
            u = p @ np.ones(n)
            norm = np.linalg.slogdet(c)[1] + np.log(u.sum()) + (n-1)*np.log(2*np.pi)
            return -.5 * (residual @ p @ residual - (u @ residual)**2/u.sum() + norm)

        direct_ratio = direct_loglike(covariance + sigma**2*sn.basis @ sn.basis.T) - direct_loglike(covariance)
        projected_ratio = sn.from_prediction(np.zeros(n))['log_likelihood_ratios'][label]
        errors = {
            'conditional_mean_error': float(np.max(abs(conditional_mean-joint_mean[1:]))),
            'conditional_covariance_error': float(np.max(abs(conditional_cov-joint_cov[1:, 1:]))),
            'global_offset_invariance_error': float(np.max(abs(shifted-conditional_mean))),
            'marginal_density_log_ratio_error': float(abs(direct_ratio-projected_ratio))}
        assert errors['conditional_mean_error'] < 1e-10
        assert errors['conditional_covariance_error'] < 1e-12
        assert errors['global_offset_invariance_error'] < 1e-11
        assert errors['marginal_density_log_ratio_error'] < 1e-9
        gls_rows.append({'case': case, 'sigma_mag': sigma, **errors})
    return {'status': 'passed', 'equal_component_normal_quantile_max_error': error,
            'independent_density_grid_quantile_max_error': mixture_error,
            'independent_full_intercept_GLS_cases': gls_rows,
            'independent_full_intercept_GLS_maximum_errors':
                {key: max(row[key] for row in gls_rows) for key in errors},
            'scope': 'Synthetic Gaussian-mixture recovery and independent full-intercept '
                     'GLS sign, covariance, normalized density and offset-invariance checks; no cosmological posterior.',
            'dependency_sha256': {str(path.relative_to(ROOT)): digest(path) for path in
                [Path(__file__).with_name('luminosity_sensitivity.py'),
                 Path(__file__).with_name('exact_correction.py')]},
            'source_sha256': digest(__file__)}


def actual(folder, summary_path):
    import camb
    qualified = summarize_run(folder, summary_path)
    manifest = json.loads((Path(folder)/'run-0.json').read_text())
    target = manifest['target_identity']
    for path, expected in target['source_sha256'].items():
        assert digest(ROOT/path) == expected, 'Scientific source changed before new background calculation.'
    for name, version in target['versions'].items():
        assert importlib.metadata.version(name) == version, 'Scientific environment changed.'
    theory = target['configuration']['theory']
    assert len(theory) == 1
    extra = next(iter(theory.values()))['extra_args']
    sources = [Path(__file__), Path(__file__).with_name('measurement_summary.py'),
               Path(__file__).with_name('luminosity_sensitivity.py'),
               Path(__file__).with_name('exact_correction.py')]
    hashes = {str(p.relative_to(ROOT)): digest(p) for p in sources}
    parent = json.loads(Path(summary_path).read_text())
    selection_path = ROOT / parent['selection_path']
    selection = json.loads(selection_path.read_text())
    settings = selection['settings']
    evolution = settings['evolution']
    records = [json.loads((selection_path.parent/f'{i:05d}.json').read_text())
               for i in range(len(selection['points']))]
    weights = np.exp(np.array([r['log_weight'] for r in records])
                     - logsumexp([r['log_weight'] for r in records]))
    data_path = sample_path(settings['sample'])
    assert digest(data_path) == target['sample_sha256']
    with np.load(data_path) as data:
        sn = IntegratedLuminosity(*(data[key] for key in ['zHD', 'zHEL', 'MU', 'covariance']))
    basis = CubicSpline(sn.knots, np.eye(5)[:, 1:], bc_type='natural')(REDSHIFTS)
    covariance = np.zeros((len(REDSHIFTS), len(REDSHIFTS)))
    coefficient_covariance = None
    if evolution.startswith('smooth'):
        factor = sn.gaussian[evolution][0]
        coefficient_covariance = cho_solve(factor, np.eye(4))
        covariance = basis @ coefficient_covariance @ basis.T
    means = []
    errors = []
    for index, record in enumerate(records):
        point = record['point']
        with background_only():
            pars = camb.set_params(**{key: point[key] for key in ['H0', 'ombh2', 'omch2', 'ns', 'tau']},
                                   As=1e-10*np.exp(point['logA']), w=point.get('w', -1.),
                                   wa=point.get('wa', 0.), **extra)
            background = camb.get_background(pars)
        z, inverse = np.unique(sn.z, return_inverse=True)
        da = background.angular_diameter_distance(z)[inverse]
        prediction = 5*np.log10(da*(1+sn.z)*(1+sn.zhel))+25
        base = sn.from_prediction(prediction)
        if evolution.startswith('smooth'):
            residual = sn.project(prediction-sn.observed)
            mean = -coefficient_covariance @ (sn.spline.T @ residual)
            means.append(basis @ mean)
            loglike = base['baseline_SN_loglike'] + base['log_likelihood_ratios'][evolution]
        else:
            epsilon = point.get('epsilon', 0.)
            means.append(epsilon*np.log1p(REDSHIFTS)/np.log(2))
            residual = sn.project(prediction+epsilon*np.log1p(sn.z)/np.log(2)-sn.observed)
            loglike = -.5*(float(residual@residual)+sn.lognorm)
        errors.append(abs(loglike-record['exact_loglikes']['released_sn']))
        if errors[-1] > 1e-5:
            raise ArithmeticError(f'Native SN density does not close at point {index}: {errors[-1]}')
        if (index+1) % 200 == 0:
            print(json.dumps({'background_points': index+1, 'total': len(records)}), flush=True)
    means = np.asarray(means)
    marginal = [gaussian_mixture_summary(means[:, j], float(covariance[j, j]), weights)
                for j in range(len(REDSHIFTS))]
    average = weights @ means
    deviations = means-average
    total_covariance = covariance+(deviations*weights[:, None]).T@deviations
    q = np.array([r['derived']['q0'] for r in records])
    qdev = q-weights@q
    cross = (deviations*weights[:, None]).T@qdev
    assert all(digest(ROOT/path) == expected for path, expected in hashes.items())
    assert all(digest(ROOT/path) == expected for path, expected in qualified['input_sha256'].items()), 'Qualified parent changed during background calculation.'
    assert all(digest(ROOT/path) == expected for path, expected in target['source_sha256'].items())
    assert all(importlib.metadata.version(name) == version for name, version in target['versions'].items())
    return {'status': 'qualified_conditional_luminosity_history',
            'settings': qualified['settings'], 'target_identity': qualified['target_identity'],
            'redshift': REDSHIFTS.tolist(), 'pointwise_magnitude_posterior': marginal,
            'covariance_mag2': total_covariance.tolist(),
            'conditional_coefficient_covariance_mag2': None if coefficient_covariance is None
                    else coefficient_covariance.tolist(),
            'covariance_with_q0_mag': cross.tolist(),
            'native_SN_loglike_max_absolute_difference': float(max(errors)),
            'data_redshift_range': [float(sn.z.min()), float(sn.z.max())],
            'outside_observed_redshift_range': ((REDSHIFTS < sn.z.min()) | (REDSHIFTS > sn.z.max())).tolist(),
            'interpretation': 'Positive B means dimmer standardized SNe at fixed distance. '
                'B(0)=0 is imposed. Baseline zero-width B is fixed, not independently measured. '
                'Smooth models retain the conditional coefficient variance analytically. '
                'These are pointwise conditional credible intervals, not simultaneous bands, '
                'a progenitor-age effect or a calibrated replacement survey correction.',
            'input_sha256': dict(qualified['input_sha256'], **{str(data_path.relative_to(ROOT)): digest(data_path)}),
            'source_sha256': hashes, 'scientific_versions': target['versions']}


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--validate', action='store_true')
    p.add_argument('--chain-folder', type=Path)
    p.add_argument('--correction-summary', type=Path)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if a.validate:
        if a.chain_folder or a.correction_summary:
            p.error('Validation uses synthetic inputs only.')
        result = validate()
    else:
        if not a.chain_folder or not a.correction_summary:
            p.error('A qualified chain folder and correction summary are required.')
        result = actual(a.chain_folder, a.correction_summary)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    a.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'output': str(a.output)}))


if __name__ == '__main__':
    main()
