"""Deterministic, anchored SN-only flat-LambdaCDM reconstruction.

H0 is integrated analytically conditional on Omega_m, with the Jacobian for a
flat H0 measure. This is a late-time check, not the shared CMB calculation.
"""
import argparse
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import time

import numpy as np
from scipy import linalg
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import log_ndtr, logsumexp
from numpy.polynomial.legendre import leggauss

from calibration_interface import ReleasedCalibration, ROOT, OUT as INTERFACE_RESULT

HERE = Path(__file__).resolve().parent
DESIGN = HERE/'calibration-lcdm-design.json'
PROBABILITIES = [.025, .16, .5, .84, .975]
K = np.log(10)/5
REFERENCE_H0 = 70.


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(path):
    return str(Path(path).resolve().relative_to(ROOT))


def normal_interval_logmass(lo, hi):
    """Stable log[Phi(hi)-Phi(lo)], including intervals in either tail."""
    lo, hi = np.broadcast_arrays(np.asarray(lo), np.asarray(hi))
    assert np.all(hi > lo)
    a = np.where(lo > 0, -lo, hi)
    b = np.where(lo > 0, -hi, lo)
    la, lb = log_ndtr(a), log_ndtr(b)
    return la + np.log(-np.expm1(lb-la))


@lru_cache(maxsize=None)
def nodes(count):
    x, w = leggauss(count)
    return (x+1)/2, w/2


class CalibrationLCDM:
    def __init__(self, released, distance_order=128, h_bounds=(50., 90.), omega_bounds=(.01, .99)):
        self.data = released
        self.distance_order = distance_order
        self.h_bounds = np.array(h_bounds)
        self.eta_bounds = np.log(self.h_bounds/REFERENCE_H0)/K
        self.omega_bounds = np.array(omega_bounds)
        self.direction = self.project((~released.calibrator).astype(float))
        self.a = float(self.direction @ self.direction)
        self.sigma = self.a**-.5
        assert self.a > 0

    def project(self, values):
        # Subtracting a common offset before whitening avoids large cancellation.
        values = values-np.mean(values)
        white = linalg.solve_triangular(self.data.chol, values, lower=True)
        return white-self.data.one_white*(self.data.one_white@white)/self.data.flat_m_precision

    def distances(self, omega, H0=REFERENCE_H0):
        z = self.data.z_hd_noncalibrator
        x, w = nodes(self.distance_order)
        inverse_e = 1/np.sqrt(omega*(1+z[:, None]*x)**3+1-omega)
        return 299792.458/H0*z*(inverse_e@w)/(1+z)

    @lru_cache(maxsize=20000)
    def coefficients(self, omega):
        r = self.data.theory_mu(self.distances(omega))-self.data.data.m_b_corr.to_numpy()
        white = self.project(r)
        b, c = float(self.direction@white), float(white@white)
        mu = b/self.a
        shifted = mu+K/self.a
        lo, hi = (self.eta_bounds-shifted)/self.sigma
        log_mass = float(normal_interval_logmass(lo, hi))
        log_integral = (-.5*(c-b*b/self.a+self.data.log_normalization)
                        + np.log(K*REFERENCE_H0/np.diff(self.h_bounds)[0])
                        + K*mu+.5*K*K/self.a + .5*np.log(2*np.pi/self.a)+log_mass)
        return {'b': b, 'c': c, 'mean_eta_flat_H0': shifted,
                'log_mass': log_mass, 'log_H0_integral': float(log_integral)}

    def H0_moment(self, omega, power):
        row = self.coefficients(float(omega))
        mu = row['mean_eta_flat_H0']
        moved = mu+power*K/self.a
        lo, hi = (self.eta_bounds-moved)/self.sigma
        return float(REFERENCE_H0**power*np.exp(power*K*mu+.5*(power*K)**2/self.a
                     + normal_interval_logmass(lo, hi)-row['log_mass']))

    def H0_cdf(self, omega, value):
        if value <= self.h_bounds[0]:
            return 0.
        if value >= self.h_bounds[1]:
            return 1.
        row = self.coefficients(float(omega))
        lo = (self.eta_bounds[0]-row['mean_eta_flat_H0'])/self.sigma
        hi = (np.log(value/REFERENCE_H0)/K-row['mean_eta_flat_H0'])/self.sigma
        return float(np.exp(normal_interval_logmass(lo, hi)-row['log_mass']))

    def quadrature(self, order, lower=None, upper=None):
        lower = self.omega_bounds[0] if lower is None else lower
        upper = self.omega_bounds[1] if upper is None else upper
        x, w = nodes(order)
        omega = lower+(upper-lower)*x
        log_terms = np.array([self.coefficients(float(o))['log_H0_integral'] for o in omega])
        log_terms += np.log(w*(upper-lower)/np.diff(self.omega_bounds)[0])
        log_total = float(logsumexp(log_terms))
        weights = np.exp(log_terms-log_total)
        return omega, weights, log_total

    def result(self, order):
        omega, weights, norm = self.quadrature(order)
        hmean = np.array([self.H0_moment(o, 1) for o in omega])
        hsecond = np.array([self.H0_moment(o, 2) for o in omega])
        mean_o, mean_h = float(weights@omega), float(weights@hmean)
        sd_o = float(np.sqrt(weights@((omega-mean_o)**2)))
        sd_h = float(np.sqrt(weights@hsecond-mean_h**2))
        # Adaptive CDF integration independently avoids assigning posterior mass
        # to quadrature nodes or interpolating a discretized empirical CDF.
        offset = self.coefficients(.33)['log_H0_integral']
        def omega_cdf(upper):
            integral = quad(lambda o: np.exp(self.coefficients(float(o))['log_H0_integral']-offset),
                            self.omega_bounds[0], upper, epsabs=1e-11, epsrel=1e-9,
                            points=[o for o in [.15, .25, .35, .5, .8] if self.omega_bounds[0] < o < upper],
                            limit=100)[0]
            return integral*np.exp(offset-norm)/np.diff(self.omega_bounds)[0]
        oq = [brentq(lambda v: omega_cdf(v)-p, self.omega_bounds[0]+1e-10,
                     self.omega_bounds[1], xtol=1e-11) for p in PROBABILITIES]
        hq = [brentq(lambda v: weights@np.array([self.H0_cdf(o, v) for o in omega])-p,
                     *self.h_bounds, xtol=1e-9) for p in PROBABILITIES]
        return {'Omega_m': {'mean': mean_o, 'sd': sd_o, 'quantiles': oq},
                'H0_km_s_Mpc': {'mean': mean_h, 'sd': sd_h, 'quantiles': hq},
                'q0': {'mean': 1.5*mean_o-1, 'sd': 1.5*sd_o,
                       'quantiles': (1.5*np.array(oq)-1).tolist()},
                'H0_Omega_m_covariance': float(weights@((omega-mean_o)*(hmean-mean_h))),
                'quantile_probabilities': PROBABILITIES, 'quadrature_order': order}


def validate_numeric(model):
    distance_checks = []
    chi2_errors = []
    integration_checks = []
    for omega in [.1, .3, .5, .8]:
        z = np.array([.01, .1, .5, 1., 2.26])
        adaptive = np.array([quad(lambda zz: 1/np.sqrt(omega*(1+zz)**3+1-omega),
                                 0, zz, epsabs=1e-12, epsrel=1e-12)[0] for zz in z])
        for count in [64, 128]:
            x, w = nodes(count)
            direct = z*((1/np.sqrt(omega*(1+z[:, None]*x)**3+1-omega))@w)
            distance_checks.append(float(np.max(abs(direct/adaptive-1))))
        row = model.coefficients(omega)
        for H0 in [50., 65., 75., 90.]:
            eta = np.log(H0/REFERENCE_H0)/K
            quadratic = row['c']-2*row['b']*eta+model.a*eta**2
            official = model.data.evaluate(model.distances(omega, H0))
            chi2_errors.append(abs(quadratic-official['chi2']))
        # Independent numerical integration in H0, not in transformed eta.
        x, w = nodes(256)
        h = model.h_bounds[0]+np.diff(model.h_bounds)[0]*x
        eta = np.log(h/REFERENCE_H0)/K
        q = row['c']-2*row['b']*eta+model.a*eta**2
        logs = -.5*(q+model.data.log_normalization)+np.log(w)
        numerical_log_integral = float(logsumexp(logs))
        normalized = np.exp(logs-numerical_log_integral)
        integration_checks.append({'Omega_m': omega,
            'log_integral_absolute_difference': abs(numerical_log_integral-row['log_H0_integral']),
            'first_moment_relative_difference': abs((normalized@h)/model.H0_moment(omega, 1)-1),
            'second_moment_relative_difference': abs((normalized@(h*h))/model.H0_moment(omega, 2)-1)})
    assert max(distance_checks) < 1e-10
    assert max(chi2_errors) < 1e-7
    assert max(v['log_integral_absolute_difference'] for v in integration_checks) < 1e-9
    assert max(v[k] for v in integration_checks for k in ['first_moment_relative_difference', 'second_moment_relative_difference']) < 1e-9
    return {'distance_max_relative_error': max(distance_checks),
            'direct_likelihood_chi2_max_absolute_difference': max(chi2_errors),
            'independent_H0_quadrature': integration_checks}


def compare(left, right):
    return {key: float(max(abs(np.array([left[key]['mean'], left[key]['sd']]+left[key]['quantiles'])
                             -np.array([right[key]['mean'], right[key]['sd']]+right[key]['quantiles']))))
            for key in ['H0_km_s_Mpc', 'Omega_m']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert not args.output.exists(), 'Preserve previous calculations; use a fresh output.'
    start = time.monotonic()
    interface = json.loads(INTERFACE_RESULT.read_text())
    assert interface['status'] == 'passed_released_Gaussian_interface_checks_no_cosmological_fit'
    assert interface['released_target_usable'] is True
    sources = {relative(p): sha(p) for p in [Path(__file__), DESIGN, INTERFACE_RESULT]}
    for path, expected in interface['source_sha256'].items():
        assert sha(ROOT/path) == expected
        sources[path] = expected
    released = ReleasedCalibration()
    inputs = {path: record['sha256'] for path, record in released.input_records.items()}
    model = CalibrationLCDM(released)
    checks = validate_numeric(model)
    baseline = model.result(512)
    lower = CalibrationLCDM(released, distance_order=64).result(256)
    refinement = compare(baseline, lower)
    assert refinement['H0_km_s_Mpc'] < 1e-5 and refinement['Omega_m'] < 1e-7
    omega, weight, norm = model.quadrature(512)
    offset = model.coefficients(.33)['log_H0_integral']
    def density(o):
        return np.exp(model.coefficients(float(o))['log_H0_integral']-offset)/np.diff(model.omega_bounds)[0]
    adaptive, adaptive_error = quad(density, *model.omega_bounds, epsabs=1e-11, epsrel=1e-9,
                                    points=[.15, .25, .35, .5, .8], limit=100)
    normalization_difference = abs(np.log(adaptive)+offset-norm)
    assert normalization_difference < 1e-8
    expanded = CalibrationLCDM(released, h_bounds=(30., 110.), omega_bounds=(.001, .999)).result(512)
    boundary_sensitivity = compare(baseline, expanded)
    for path, expected in dict(inputs, **sources).items():
        assert sha(ROOT/path) == expected, 'Scientific input changed during calculation.'
    result = {'status': 'qualified_conditional_SN_only_LCDM_quadrature',
        'posterior': baseline, 'model': json.loads(DESIGN.read_text())['physics'],
        'priors': json.loads(DESIGN.read_text())['priors'], 'selected_rows': len(released.data),
        'calibrator_rows': int(np.sum(released.calibrator)),
        'checks': checks, 'quadrature_refinement_max_absolute_difference': refinement,
        'independent_adaptive_Omega_log_normalization_difference': normalization_difference,
        'independent_adaptive_Omega_relative_error_estimate': adaptive_error/adaptive,
        'expanded_prior_result': expanded, 'prior_boundary_sensitivity_max_absolute_difference': boundary_sensitivity,
        'comparison': {'source': 'https://arxiv.org/html/2202.04077v2#S4.T3',
                       'published_H0_mean_sd': [73.6, 1.1], 'published_Omega_m_mean_sd': [.334, .018]},
        'source_sha256': sources, 'input_sha256': inputs, 'seconds': time.monotonic()-start,
        'CMB_spectrum_calls': 0, 'sampling_steps': 0,
        'limits': ['Conditional on the released corrected magnitudes, Cepheid means and total covariance.',
                   'The exact embedded Cepheid covariance construction is not independently recovered.',
                   'Matter plus Lambda at late times; no CMB, BAO or independent stellar-age likelihood.',
                   'One flat absolute-magnitude measure, no additional H0 prior; no absolute model evidence.',
                   'An independent calculation of shared measurements, not independent observations.']}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status': result['status'], 'posterior': baseline, 'seconds': result['seconds']}))


if __name__ == '__main__':
    main()
