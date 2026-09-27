"""Explicit Gaussian reconstruction of the rounded 2026 Ly-alpha distances.

This is not an author likelihood or a reconstruction of the forest spectra.
The overlapping old Ly-alpha pair is replaced, never counted a second time.
"""
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.linalg import cholesky, solve_triangular

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
DESIGN = HERE/'lya-fullshape-design.json'
RELEASE = ROOT/'.work/unified-cosmology/external-probes/packages/data/bao_data/desi_bao_dr2'
MEAN = RELEASE/'desi_gaussian_bao_ALL_GCcomb_mean.txt'
COVARIANCE = RELEASE/'desi_gaussian_bao_ALL_GCcomb_cov.txt'
C_KM_S = 299792.458


def variants():
    design = json.loads(DESIGN.read_text())
    nominal = np.asarray(design['printed_summary']['mean_sigma_rho'], dtype=float)
    half = np.asarray(design['printed_summary']['rounding_half_steps'], dtype=float)
    result = {'nominal': nominal}
    for bits in itertools.product([0, 1], repeat=5):
        result['rounding_'+''.join(map(str, bits))] = nominal+(2*np.asarray(bits)-1)*half
    return result


class GaussianBAOReplacement:
    def __init__(self, redshifts, mean, row_types, covariance):
        self.z, self.mean = [np.array(x, dtype=float, copy=True) for x in [redshifts, mean]]
        self.row_types = np.array(row_types, dtype=str, copy=True)
        self.C = np.array(covariance, dtype=float, copy=True)
        n = len(self.z)
        assert self.z.shape == self.mean.shape == self.row_types.shape == (n,)
        assert n >= 3 and self.C.shape == (n, n)
        assert np.isfinite(self.z).all() and (self.z > 0).all()
        assert np.isfinite(self.mean).all() and np.isfinite(self.C).all()
        assert np.array_equal(self.C, self.C.T)
        assert set(self.row_types) <= {'DV_over_rs', 'DM_over_rs', 'DH_over_rs'}
        self.lya = self.z == 2.33
        assert self.lya.sum() == 2
        self.dm = np.flatnonzero(self.lya & (self.row_types == 'DM_over_rs'))
        self.dh = np.flatnonzero(self.lya & (self.row_types == 'DH_over_rs'))
        assert len(self.dm) == len(self.dh) == 1
        self.order = np.array([int(self.dm[0]), int(self.dh[0])])
        assert not np.count_nonzero(self.C[np.ix_(self.lya, ~self.lya)]), 'Nonzero old cross-block covariance requires a different replacement model.'
        self.L = cholesky(self.C, lower=True)
        self.Lkeep = cholesky(self.C[np.ix_(~self.lya, ~self.lya)], lower=True)
        self.Loldlya = cholesky(self.C[np.ix_(self.order, self.order)], lower=True)
        self.new = {}
        for name, numbers in variants().items():
            dm, dh, sm, sh, rho = numbers
            assert sm > 0 and sh > 0 and abs(rho) < 1
            covariance = np.outer([sm, sh], [sm, sh])*np.array([[1., rho], [rho, 1.]])
            total = self.C.copy(); total[np.ix_(self.order, self.order)] = covariance
            m = self.mean.copy(); m[self.order] = [dm, dh]
            self.new[name] = (m, cholesky(total, lower=True), cholesky(covariance, lower=True))

    @classmethod
    def from_release(cls, mean_path=MEAN, covariance_path=COVARIANCE):
        rows = np.loadtxt(mean_path, dtype=str)
        result = cls(rows[:, 0].astype(float), rows[:, 1].astype(float), rows[:, 2], np.loadtxt(covariance_path))
        assert len(result.z) == 13 and np.sum(~result.lya) == 11
        return result

    def prediction(self, transverse_comoving_distance_mpc, hubble_km_s_Mpc, rdrag_mpc):
        dm, hubble = [np.asarray(x, dtype=float) for x in [transverse_comoving_distance_mpc, hubble_km_s_Mpc]]
        assert dm.shape == hubble.shape == self.z.shape
        assert np.isfinite(dm).all() and np.isfinite(hubble).all()
        assert (dm > 0).all() and (hubble > 0).all() and np.isfinite(rdrag_mpc) and rdrag_mpc > 0
        dh = C_KM_S/hubble
        lookup = {'DM_over_rs': dm/rdrag_mpc, 'DH_over_rs': dh/rdrag_mpc,
                  'DV_over_rs': np.cbrt(dm*dm*dh*self.z)/rdrag_mpc}
        return np.array([lookup[k][i] for i, k in enumerate(self.row_types)])

    @staticmethod
    def quadratic(L, residual):
        whitened = solve_triangular(L, residual, lower=True)
        return float(whitened@whitened)

    def evaluate_prediction(self, prediction):
        prediction = np.asarray(prediction, dtype=float)
        assert prediction.shape == self.mean.shape and np.isfinite(prediction).all() and (prediction > 0).all()
        old_residual = prediction-self.mean
        old = self.quadratic(self.L, old_residual)
        kept = self.quadratic(self.Lkeep, old_residual[~self.lya])
        old_lya = self.quadratic(self.Loldlya, old_residual[self.order])
        outcomes = {}
        for name, (mean, L, Lpair) in self.new.items():
            residual = prediction-mean
            total = self.quadratic(L, residual)
            pair = self.quadratic(Lpair, residual[self.order])
            outcomes[name] = {'new_loglike': -.5*total, 'new_chi2': total,
                              'new_lya_chi2': pair, 'new_block_closure': total-kept-pair,
                              'delta_loglike': -.5*(total-old)}
        return {'old_loglike': -.5*old, 'old_chi2': old, 'retained_chi2': kept,
                'old_lya_chi2': old_lya, 'old_block_closure': old-kept-old_lya,
                'prediction': prediction.tolist(), 'variants': outcomes}

    def evaluate(self, transverse_comoving_distance_mpc, hubble_km_s_Mpc, rdrag_mpc):
        return self.evaluate_prediction(self.prediction(transverse_comoving_distance_mpc, hubble_km_s_Mpc, rdrag_mpc))
