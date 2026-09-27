"""Correlated Gaussian calibration prediction after excluding calibrator values."""
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cholesky, cho_solve, solve_triangular
from scipy.special import log_ndtr, logsumexp

HERE = Path(__file__).resolve().parent
DESIGN = HERE/'calibration-holdout-design.json'


class GaussianHoldout:
    """One flat common magnitude, full covariance, and fixed withheld row mask."""
    def __init__(self, covariance, calibrator):
        self.C = np.array(covariance, dtype=float, copy=True)
        self.mask = np.array(calibrator, copy=True)
        assert self.mask.dtype == bool and self.mask.ndim == 1
        assert self.C.shape == (len(self.mask), len(self.mask))
        assert np.isfinite(self.C).all() and np.array_equal(self.C, self.C.T)
        self.nc, self.nn = int(self.mask.sum()), int((~self.mask).sum())
        assert self.nc > 0 and self.nn >= 2
        self.L = cholesky(self.C, lower=True)
        self.Ln = cholesky(self.C[np.ix_(~self.mask, ~self.mask)], lower=True)
        self.on = solve_triangular(self.Ln, np.ones(self.nn), lower=True)
        self.an = float(self.on@self.on)
        self.ofull = solve_triangular(self.L, np.ones(len(self.mask)), lower=True)
        self.afull = float(self.ofull@self.ofull)
        self.lognorm_n = float(2*np.log(np.diag(self.Ln)).sum()+np.log(self.an)+(self.nn-1)*np.log(2*np.pi))
        self.lognorm_full = float(2*np.log(np.diag(self.L)).sum()+np.log(self.afull)+(len(self.mask)-1)*np.log(2*np.pi))
        Ccn = self.C[np.ix_(self.mask, ~self.mask)]
        self.A = cho_solve((self.Ln, True), Ccn.T).T
        self.b = np.ones(self.nc)-self.A@np.ones(self.nn)
        V0 = self.C[np.ix_(self.mask, self.mask)]-self.A@Ccn.T
        self.V0 = (V0+V0.T)/2
        cholesky(self.V0, lower=True)
        self.V = self.V0+np.outer(self.b, self.b)/self.an
        self.Lv = cholesky(self.V, lower=True)
        self.lognorm_c = float(2*np.log(np.diag(self.Lv)).sum()+self.nc*np.log(2*np.pi))
        # A common offset of the calibrator block, or a constant change in
        # the noncalibrator distance modulus, moves t along 1_C. The vector
        # b above propagates latent-M uncertainty; it is not this offset mode.
        direction = np.ones(self.nc)
        vb = cho_solve((self.Lv, True), direction)
        precision = float(direction@vb)
        assert np.isfinite(precision) and precision > 0, 'No identifiable common-calibration contrast.'
        self.v = vb/precision
        self.sigma = float(precision**-.5)

    def evaluate(self, residual):
        residual = np.asarray(residual, dtype=float)
        assert residual.shape == self.mask.shape and np.isfinite(residual).all()
        # A common subtraction leaves both flat-M factors and the conditional
        # predictive factor unchanged, while reducing cancellation.
        residual = residual-residual[np.flatnonzero(~self.mask)[0]]
        rn, rc = residual[~self.mask], residual[self.mask]
        wn = solve_triangular(self.Ln, rn, lower=True)
        mhat = float(self.on@wn/self.an)
        pn = wn-self.on*mhat
        chi_n = float(pn@pn)
        t = rc-self.A@rn-self.b*mhat
        wc = solve_triangular(self.Lv, t, lower=True)
        chi_c = float(wc@wc)
        wf = solve_triangular(self.L, residual, lower=True)
        pf = wf-self.ofull*(self.ofull@wf/self.afull)
        full = float(-.5*(pf@pf+self.lognorm_full))
        logn, logc = -.5*(chi_n+self.lognorm_n), -.5*(chi_c+self.lognorm_c)
        contrast = float(self.v@t)
        z = contrast/self.sigma
        return {'noncalibrator_loglike': float(logn), 'noncalibrator_chi2': chi_n,
                'conditional_calibrator_loglike': float(logc), 'conditional_calibrator_chi2': chi_c,
                'full_loglike': full, 'full_loglike_closure': float(full-logn-logc),
                'calibration_contrast_mag': contrast, 'contrast_sigma_mag': self.sigma, 'contrast_z': float(z),
                'conditional_logcdf': float(log_ndtr(z)), 'conditional_logsf': float(log_ndtr(-z)),
                'rows_noncalibrator': self.nn, 'rows_calibrator': self.nc}


class ReleasedHoldout:
    def __init__(self):
        from anchored_adapter import ReleasedCalibration
        self.release = ReleasedCalibration()
        self.kernel = GaussianHoldout(self.release.C, self.release.calibrator)
        assert (self.kernel.nn, self.kernel.nc) == (1580, 77)

    def evaluate(self, angular_diameter_distance_mpc):
        mu = self.release.theory_mu(angular_diameter_distance_mpc)
        result = self.kernel.evaluate(self.release.data.m_b_corr.to_numpy()-mu)
        full = self.release.evaluate(angular_diameter_distance_mpc)
        result['released_full_loglike_closure'] = float(result['full_loglike']-full['loglike'])
        return result


def contribution_precision(logweights, logintegrand, groups, limits):
    """Self-normalized importance influence, scaled to the predictive mean."""
    logweights, logintegrand, groups = map(np.asarray, [logweights, logintegrand, groups])
    assert logweights.shape == logintegrand.shape == groups.shape and logweights.ndim == 1
    assert len(groups) >= 2000 and set(groups) == {0, 1, 2, 3}
    assert np.isfinite(logweights).all() and np.isfinite(logintegrand).all()
    lw = logweights-logsumexp(logweights)
    mean = float(logsumexp(lw+logintegrand))
    weights = np.exp(lw)
    rho = np.exp(lw+logintegrand-mean)
    ess = float(1/(rho@rho))
    mass = {str(g): float(rho[groups == g].sum()) for g in range(4)}
    errors = {}
    for count in [10, 20, 40]:
        blocks = [idx for g in range(4) for idx in np.array_split(np.flatnonzero(groups == g), count)]
        assert min(map(len, blocks)) >= 10
        influence = [float((rho[idx]-weights[idx]).sum()) for idx in blocks]
        errors[str(count)] = float(np.sqrt(len(blocks)*np.var(influence, ddof=1)))
    chain_influence = [float((rho[groups == g]-weights[groups == g]).sum()) for g in range(4)]
    chain_error = float(np.sqrt(4*np.var(chain_influence, ddof=1)))
    failed = []
    if ess < limits['minimum_contribution_ESS']: failed.append('contribution_ESS')
    for key, value in mass.items():
        if not limits['chain_mass_min'] <= value <= limits['chain_mass_max']: failed.append('contribution_chain_mass:'+key)
    for key, value in errors.items():
        if value > limits['maximum_relative_batch_MCSE']: failed.append('relative_batch_MCSE:'+key)
    if chain_error > limits['maximum_relative_batch_MCSE']: failed.append('relative_independent_chain_MCSE')
    return {'qualified_predictive_precision': not failed, 'failed_gates': failed,
            'diagnostic_log_mean': mean, 'contribution_weight_ESS': ess,
            'contribution_chain_weight_fractions': mass,
            'relative_batch_MCSE': errors, 'relative_independent_chain_MCSE': chain_error,
            'largest_contribution_fraction': float(rho.max())}


def predictive_summary(rows, logweights, groups):
    """Call only after the HF-only posterior has passed all original gates."""
    design = json.loads(DESIGN.read_text())
    keys = {'CDF': 'conditional_logcdf', 'SF': 'conditional_logsf',
            'joint_density': 'conditional_calibrator_loglike'}
    checks = {name: contribution_precision(logweights, [r[key] for r in rows], groups,
                                         design['predictive_precision']) for name, key in keys.items()}
    lcdf, lsf = [checks[k]['diagnostic_log_mean'] for k in ['CDF', 'SF']]
    assert abs(float(np.logaddexp(lcdf, lsf))) < 1e-10
    side = 'CDF' if lcdf <= lsf else 'SF'
    logtail = min(0., float(np.log(2)+min(lcdf, lsf)))
    tail = None
    if checks[side]['qualified_predictive_precision']:
        tail = {'log_two_sided_equal_tail_probability': logtail,
                'two_sided_equal_tail_probability': float(np.exp(logtail)) if logtail > -700 else None,
                'probability_underflow_avoided': logtail <= -700,
                'smaller_tail': side, 'log_CDF': lcdf, 'log_SF': lsf,
                'scope': design['predictive_interpretation']}
    return {'status': 'qualified_predictive_precision' if tail is not None and checks['joint_density']['qualified_predictive_precision']
                       else 'predictive_precision_requires_followup',
            'calibration_contrast_predictive_tail': tail,
            'conditional_calibrator_log_predictive_density': checks['joint_density']['diagnostic_log_mean']
                if checks['joint_density']['qualified_predictive_precision'] else None,
            'predictive_precision_diagnostics': checks,
            'diagnostic_values_are_not_qualified_when_their_gate_fails': True,
            'held_out_rows': 77, 'model_comparison_or_frequentist_sigma': False,
            'scope': design['predictive_interpretation']}
