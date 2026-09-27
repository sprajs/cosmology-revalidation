"""Released SN distance likelihood sharing the external probes' CAMB background.

This does not replace the survey's light-curve, selection or CC inference.
The one unknown absolute magnitude is analytically integrated with a flat prior.
Optional smooth luminosity modes have an explicitly chosen Gaussian prior;
they are a sensitivity experiment, not observed host-age corrections.
"""
from pathlib import Path
import numpy as np
from scipy.linalg import cho_factor, cho_solve
from scipy.interpolate import CubicSpline
from cobaya.likelihood import Likelihood

ROOT = Path(__file__).resolve().parents[4]
DEFAULT_DATA = ROOT / '.work/unified-cosmology/survey-selection/normalized/dovekie-total.npz'


class ReleasedDistances(Likelihood):
    data_file: str = str(DEFAULT_DATA)
    smooth_sigma: float = 0.
    input_params = ['epsilon']
    type = 'SN'
    speed = 100.
    params = {'sn_chi2': {'derived': True}}

    def initialize(self):
        with np.load(self.data_file, allow_pickle=False) as f:
            self.z = np.asarray(f['zHD'])
            self.zhel = np.asarray(f['zHEL'])
            self.observed = np.asarray(f['MU'])
            covariance = np.asarray(f['covariance'])
        assert self.z.shape == self.observed.shape == self.zhel.shape
        assert covariance.shape == (len(self.z), len(self.z))
        assert np.all(self.z > 0) and np.isfinite(covariance).all()
        # A common zero point is unidentifiable and is integrated separately.
        # Four independent shape coefficients: zero at z=0, knots at .1,.4,.8,1.3.
        knots = np.array([0., .1, .4, .8, 1.3])
        self.basis = CubicSpline(knots, np.eye(5)[:, 1:], bc_type='natural')(self.z)
        if self.smooth_sigma:
            covariance = covariance + self.smooth_sigma**2 * self.basis @ self.basis.T
        self.factor = cho_factor(covariance, lower=True)
        self.precision = cho_solve(self.factor, np.eye(len(self.z)))
        self.u = self.precision @ np.ones(len(self.z))
        self.A = self.u.sum()
        self.logdet = 2*np.log(np.diag(self.factor[0])).sum()
        self.evolution = np.log1p(self.z) / np.log(2.)
        self.background_z, self.inverse_z = np.unique(self.z, return_inverse=True)

    def get_requirements(self):
        return {'angular_diameter_distance': {'z': self.background_z}}

    def logp(self, epsilon=0., _derived=None):
        da = self.provider.get_angular_diameter_distance(self.background_z)[self.inverse_z]
        dl = da * (1+self.z) * (1+self.zhel)
        prediction = 5*np.log10(dl) + 25 + epsilon*self.evolution
        residual = prediction - self.observed
        # Subtract a constant before projection to avoid cancellation for large M.
        residual = residual - residual.mean()
        score = float(residual @ self.precision @ residual - (residual @ self.u)**2/self.A)
        if _derived is not None:
            _derived['sn_chi2'] = score
        # Includes all fixed-covariance factors and the flat-M integral; no evidence
        # claim is made because the absolute normalisation of that prior is improper.
        return -.5*(score+self.logdet+np.log(self.A)+(len(self.z)-1)*np.log(2*np.pi))


def expansion_diagnostics(z, hubble):
    """Differentiate local H(z) in physical units; no matter-only q shortcut."""
    vals = {}
    for centre, label in [(0., '0'), (.5, '05'), (1., '1')]:
        pick = (z >= centre-1e-10) & (z <= centre+.00400001)
        coeff = np.polynomial.polynomial.polyfit(z[pick]-centre, hubble[pick]/hubble[pick][0], 4)
        hp, hpp = coeff[1], 2*coeff[2]
        q = (1+centre)*hp-1
        qp = hp+(1+centre)*(hpp-hp*hp)
        vals['q'+label] = float(q)
        vals['j'+label] = float(q*(2*q+1)+(1+centre)*qp)
    return vals


class ExpansionDiagnostics(Likelihood):
    input_params = []
    params = {x: {'derived': True} for x in ['q0', 'q05', 'q1', 'j0', 'j05', 'j1']}

    def initialize(self):
        self.z = np.concatenate([np.arange(5)*.001+c for c in [0., .5, 1.]])

    def get_requirements(self):
        return {'Hubble': {'z': self.z}}

    def logp(self, _derived=None):
        if _derived is not None:
            _derived.update(expansion_diagnostics(self.z, self.provider.get_Hubble(self.z)))
        return 0.
