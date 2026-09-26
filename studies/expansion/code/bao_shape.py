"""Registered BAO shape-cone test; no supernova or dark-energy law enters.

See docs/research-2026-09-26/plan.md for mathematical and observational scope.
The 12D cone is represented independently by inequalities and 12 step rays.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.linalg import solve_triangular
from scipy.optimize import nnls, minimize, LinearConstraint
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / 'sources/repos/CobayaSampler__bao_data/desi_bao_dr2'
OUT = ROOT / 'runs/research_2026_09_26/bao_shape'
SOURCE_BYTES = Path(__file__).read_bytes()
PLAN = ROOT / 'docs/research-2026-09-26/plan.md'
PLAN_BYTES = PLAN.read_bytes()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def load_data():
    names = ['desi_gaussian_bao_ALL_GCcomb_mean.txt',
             'desi_gaussian_bao_ALL_GCcomb_cov.txt']
    paths = [INPUT / name for name in names]
    hashes = {str(p.relative_to(ROOT)): digest(p.read_bytes()) for p in paths}
    df = pd.read_csv(paths[0], sep=r'\s+', comment='#', names=['z', 'value', 'kind'])
    cov = np.loadtxt(paths[1])
    zs = np.sort(df.loc[df.kind == 'DH_over_rs', 'z'].unique())
    order = []
    for kind in ['DM_over_rs', 'DH_over_rs']:
        for z in zs:
            idx = np.flatnonzero((df.z == z) & (df.kind == kind))
            assert len(idx) == 1
            order.append(int(idx[0]))
    scale = np.r_[np.ones(len(zs)), 1 + zs]
    y = df.value.to_numpy()[order] * scale
    C = cov[np.ix_(order, order)] * np.outer(scale, scale)
    return zs, y, C, hashes, order, df, cov


def cone_matrices(t):
    n = len(t)
    rows, labels, rays = [], [], []
    for j in range(n):
        r = np.zeros(2*n); r[n+j] = 1
        rows.append(r); labels.append(f'positive_g_{j}')
        r = np.zeros(2*n); r[j] = 1; r[n+j] = -(t[j]-(t[j-1] if j else 0))
        if j:
            r[j-1] = -1
        rows.append(r); labels.append(f'lower_integral_{j}')
        if j:
            r = np.zeros(2*n); r[n+j-1] = 1; r[n+j] = -1
            rows.append(r); labels.append(f'monotone_g_{j}')
            r = np.zeros(2*n); r[j-1] = 1; r[j] = -1; r[n+j-1] = t[j]-t[j-1]
            rows.append(r); labels.append(f'upper_integral_{j}')
        for side in ['before', 'after']:
            g = (np.arange(n) < j if side == 'before' else np.arange(n) <= j).astype(float)
            rays.append(np.r_[np.minimum(t, t[j]), g])
    return np.array(rows), np.array(rays).T, labels


class Cone:
    def __init__(self, cov, rays, inequalities):
        self.L = np.linalg.cholesky(cov)
        self.rays = rays
        self.A = inequalities
        W = solve_triangular(self.L, rays, lower=True)
        self.norm = np.linalg.norm(W, axis=0)
        self.W = W/self.norm
        self.B = inequalities @ self.L
        self.B /= np.linalg.norm(self.B, axis=1)[:, None]

    def fit(self, y):
        v = solve_triangular(self.L, y, lower=True)
        coef, norm = nnls(self.W, v, maxiter=1000)
        fitted = self.rays @ (coef/self.norm)
        assert np.min(self.A @ fitted) > -1e-8
        resid = self.W @ coef - v
        grad = self.W.T @ resid
        kkt = max(float(np.max(np.abs(grad[coef > 1e-8]), initial=0)),
                  float(np.max(-grad[coef <= 1e-8], initial=0)))
        return norm*norm, fitted, kkt

    def independent(self, y):
        v = solve_triangular(self.L, y, lower=True)
        fit = minimize(lambda x: .5*np.sum((x-v)**2), np.zeros(len(y)),
                       jac=lambda x: x-v, method='SLSQP',
                       constraints=[LinearConstraint(self.B, 0, np.inf)],
                       options={'ftol': 1e-11, 'maxiter': 2000})
        feasibility = float(np.min(self.B @ fit.x))
        if not fit.success and feasibility < -1e-7:
            raise RuntimeError(str(fit.message))
        return 2*float(fit.fun), feasibility, bool(fit.success)

    def statistics(self, means, standard_noise):
        return np.array([self.fit(mean + self.L @ eps)[0]
                         for mean, eps in zip(np.broadcast_to(means, standard_noise.shape), standard_noise)])


def tail(values, threshold):
    n = len(values); k = int(np.sum(values >= threshold-1e-10))
    lo = 0 if k == 0 else float(beta.ppf(.025, k, n-k+1))
    hi = 1 if k == n else float(beta.ppf(.975, k+1, n-k))
    return {'exceedances': k, 'draws': n, 'fraction': k/n,
            'binomial_95_interval': [lo, hi], 'plus_one_tail': (k+1)/(n+1)}


def lcdm_vector(z, om, scale):
    e = lambda x: np.sqrt(om*(1+x)**3+1-om)
    d = np.array([scale*quad(lambda x: 1/e(x), 0, zz, epsabs=1e-12, epsrel=1e-12)[0] for zz in z])
    g = scale*(1+z)/e(z)
    return np.r_[d, g]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--draws', type=int, default=1000)
    ap.add_argument('--seed', type=int, default=260926)
    ap.add_argument('--name', default='primary')
    args = ap.parse_args()
    out = OUT/args.name
    if out.exists():
        raise FileExistsError(f'Preserve existing run: {out}')
    out.mkdir(parents=True)
    (out/'executed_source.py').write_bytes(SOURCE_BYTES)
    (out/'registered_plan.md').write_bytes(PLAN_BYTES)
    z, y, C, hashes, order, df, fullcov = load_data()
    n = len(z); t = np.log1p(z)
    A, rays, labels = cone_matrices(t)
    cone = Cone(C, rays, A)
    observed, fitted, kkt = cone.fit(y)
    independent = cone.independent(y)
    rng = np.random.default_rng(args.seed)
    checks = {'ray_feasibility_min': float(np.min(A @ rays)),
              'independent_observed_statistic': independent[0],
              'independent_observed_success': independent[2],
              'observed_kkt_error': kkt,
              'marginal_row_order': order}
    deltas, maxkkt, maxdom = [], kkt, -np.inf
    # Both solvers see realizations near the apex and near physical distances.
    for i in range(40):
        point = rng.normal(size=2*n)*5 if i < 20 else y + cone.L @ rng.normal(size=2*n)
        s, _, kk = cone.fit(point)
        ss, _, _ = cone.independent(point)
        deltas.append(abs(s-ss)); maxkkt = max(maxkkt, kk)
        eps = cone.L @ rng.normal(size=2*n)
        nullmean = rays @ rng.exponential(4, size=2*n)
        maxdom = max(maxdom, cone.fit(nullmean+eps)[0]-cone.fit(eps)[0])
        assert cone.fit(nullmean)[0] < 1e-16
    checks['max_independent_statistic_difference'] = max(deltas)
    checks['max_kkt_error'] = maxkkt
    checks['max_pointwise_null_domination_violation'] = maxdom
    scaled = Cone(C*2.7**2, rays, A).fit(y*2.7)[0]
    checks['common_scale_statistic_difference'] = abs(scaled-observed)
    assert max(deltas) < 1e-6 and maxkkt < 1e-6 and maxdom < 1e-7
    assert abs(scaled-observed) < 1e-8
    # Exact constant-q endpoints: positive q must be feasible; negative q violates.
    cqchecks = {}
    for q in [-.5, 0, .5, 1.0]:
        g = 30*np.exp(-q*t)
        d = 30*t if q == 0 else -30*np.expm1(-q*t)/q
        cqchecks[str(q)] = float(np.min(A @ np.r_[d, g]))
        if q >= 0:
            assert cqchecks[str(q)] > -1e-10
    checks['constant_q_minimum_slack'] = cqchecks
    coast = np.r_[t, np.ones(n)]
    wc = solve_triangular(cone.L, coast, lower=True)
    wy = solve_triangular(cone.L, y, lower=True)
    amplitude = float(wc@wy/(wc@wc))
    means = {'vertex_global_bound': np.zeros(2*n), 'coasting_pointwise': coast*amplitude,
             'fitted_null_pointwise': fitted}
    noise = rng.normal(size=(args.draws, 2*n))
    allstats, calibration = {}, {}
    for name, mean in means.items():
        values = cone.statistics(mean, noise)
        allstats[name] = values
        calibration[name] = tail(values, observed)
        print(name, calibration[name], flush=True)
    assert np.max(allstats['coasting_pointwise']-allstats['vertex_global_bound']) < 1e-7
    assert np.max(allstats['fitted_null_pointwise']-allstats['vertex_global_bound']) < 1e-7
    threshold = float(np.quantile(allstats['vertex_global_bound'], .95, method='higher'))
    power = {}
    for om in [.2, .3, .4]:
        # Free overall scale is profiled separately for each injected shape, not a sign choice.
        shape = lcdm_vector(z, om, 1)
        wshape = solve_triangular(cone.L, shape, lower=True)
        scale = float(wshape@wy/(wshape@wshape))
        mean = shape*scale
        values = cone.statistics(mean, rng.normal(size=noise.shape))
        allstats[f'power_Om_{om}'] = values
        power[str(om)] = {'injected_c_over_H0rd': scale,
                         'nominal_global_5percent_rejection': tail(values, threshold),
                         'noiseless_cone_distance': cone.fit(mean)[0]}
    # Radial-only monotonic cone in g; no transverse-distance/flatness integral.
    rrays = np.triu(np.ones((n, n)))
    ra = np.vstack([np.eye(n), np.eye(n)[:-1]-np.eye(n)[1:]])
    radial = Cone(C[n:, n:], rrays, ra)
    robs, rfit, rkkt = radial.fit(y[n:])
    rnoise = rng.normal(size=(args.draws, n))
    rv = radial.statistics(np.zeros(n), rnoise)
    allstats['radial_vertex_global_bound'] = rv
    # AP contrasts are linear; retain full covariance for interpretation.
    B = np.column_stack([np.eye(n), -np.diag(t)])
    aps = B@y; apcov = B@C@B.T
    table = pd.DataFrame({'z': z, 't': t, 'observed_DM_rd': y[:n],
                          'observed_g': y[n:], 'null_DM_rd': fitted[:n], 'null_g': fitted[n:],
                          'AP_contrast': aps, 'AP_contrast_sd': np.sqrt(np.diag(apcov)),
                          'AP_contrast_standardized': aps/np.sqrt(np.diag(apcov)),
                          'radial_null_g': rfit})
    table.to_csv(out/'endpoints.csv', index=False)
    np.savez_compressed(out/'calibration.npz', **allstats)
    np.savez_compressed(out/'geometry.npz', data=y, covariance=C, inequalities=A,
                        step_rays=rays, fitted_null=fitted, AP_covariance=apcov)
    results = {'scope': 'Flat constant-ruler BAO-only no-acceleration cone; not q0, raw-catalogue validation, or a physical correction measurement.',
               'statistic': observed, 'dimension': 2*n, 'inequalities': len(A),
               'coasting_amplitude_g': amplitude,
               'coasting_chisq': float(np.sum((wy-amplitude*wc)**2)),
               'calibration': calibration, 'global_95percent_threshold': threshold,
               'power': power, 'radial_only': {'statistic': robs, 'kkt_error': rkkt,
                                             'global_tail_bound': tail(rv, robs)},
               'validation': checks, 'constraint_labels': labels,
               'fitted_slack': (A@fitted).tolist(), 'settings': vars(args)}
    (out/'results.json').write_text(json.dumps(results, indent=2)+'\n')
    outputs = [out/'endpoints.csv', out/'calibration.npz', out/'geometry.npz', out/'results.json']
    manifest = {'created_utc': datetime.now(timezone.utc).isoformat(),
                'inputs_sha256': hashes, 'executed_source_sha256': digest(SOURCE_BYTES),
                'registered_plan_sha256': digest(PLAN_BYTES), 'argv': sys.argv,
                'outputs_sha256': {str(p.relative_to(ROOT)): digest(p.read_bytes()) for p in outputs},
                'numpy_version': np.__version__, 'source_capture': 'before data calculations; exact copies in run directory'}
    (out/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({'observed': observed, 'radial': robs, 'power': power, 'validation': checks}, indent=2))


if __name__ == '__main__':
    main()
