"""Exploratory checks selected after the registered BAO cone pilot.

Does not replace the primary statistic or its preserved calibration.
"""
from pathlib import Path
import json
import hashlib
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular
from scipy.optimize import minimize
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from bao_shape import ROOT, OUT, load_data

SOURCE = Path(__file__).read_bytes()


def predictions(par, t, evalt=None):
    """Positive g, continuous log-g, piecewise constant nonnegative q."""
    g0 = np.exp(par[0]); q = np.asarray(par[1:])
    dt = np.diff(np.r_[0, t])
    gstart = g0*np.exp(-np.r_[0, np.cumsum(q*dt)[:-1]])
    phi = np.ones_like(q)
    np.divide(-np.expm1(-q*dt), q*dt, out=phi, where=np.abs(q*dt) > 1e-10)
    d = np.cumsum(gstart*dt*phi)
    g = g0*np.exp(-np.cumsum(q*dt))
    if evalt is None:
        return np.r_[d, g]
    idx = np.minimum(np.searchsorted(t, evalt), len(t)-1)
    widths = np.asarray(evalt)-np.r_[0, t[:-1]][idx]
    qi = q[idx]
    integ = np.empty_like(widths)
    np.divide(-np.expm1(-qi*widths), qi, out=integ, where=np.abs(qi) > 1e-10)
    integ = np.where(np.abs(qi) > 1e-10, integ, widths)
    return np.r_[0, d[:-1]][idx]+gstart[idx]*integ, gstart[idx]*np.exp(-qi*widths)


def main():
    out = OUT/'post_pilot_checks'
    if out.exists():
        raise FileExistsError(out)
    out.mkdir()
    (out/'executed_source.py').write_bytes(SOURCE)
    z, y, C, hashes, order, raw, rawcov = load_data()
    n = len(z); t = np.log1p(z); L = np.linalg.cholesky(C)
    B = np.column_stack([np.eye(n), -np.diag(t)])
    ap = B@y; acov = B@C@B.T; sd = np.sqrt(np.diag(acov)); zz = ap/sd
    assert np.max(np.abs(acov-np.diag(np.diag(acov)))) < 1e-12
    # Under no acceleration each contrast mean >=0. Coasting makes all zero,
    # simultaneously, so this is a least-favourable null for this particular score.
    maxdeficit = float(max(0, -zz.min()))
    family_p = float(-np.expm1(n*np.log(norm.cdf(maxdeficit))))
    rng = np.random.default_rng(260927)
    fits = []
    objective = lambda p: float(np.sum(solve_triangular(L, y-predictions(p, t), lower=True)**2))
    bounds = [(np.log(1), np.log(200))]+[(0,5)]*n
    for i in range(8):
        start = np.r_[np.log(33), np.zeros(n) if i == 0 else rng.uniform(0,.7,n)]
        f = minimize(objective, start, bounds=bounds, method='L-BFGS-B',
                     options={'ftol': 1e-12, 'gtol': 1e-7, 'maxiter': 3000})
        fits.append({'chi2': float(f.fun), 'parameters': f.x.tolist(), 'success': bool(f.success)})
    best = min(fits, key=lambda x:x['chi2'])
    f = minimize(objective, best['parameters'], method='Powell', bounds=bounds,
                 options={'ftol': 1e-12,'xtol':1e-9, 'maxiter':3000})
    fits.append({'chi2': float(f.fun), 'parameters': f.x.tolist(), 'success': bool(f.success), 'method':'Powell'})
    best = min(fits, key=lambda x:x['chi2'])
    bgsidx = int(np.flatnonzero(raw.kind == 'DV_over_rs')[0]); bz = raw.z.iloc[bgsidx]
    bd, bg = predictions(best['parameters'],t,np.array([np.log1p(bz)]))
    bv = float((bz/(1+bz)*bd[0]**2*bg[0])**(1/3))
    bvalue = float(raw.value.iloc[bgsidx]); bsd = float(np.sqrt(rawcov[bgsidx,bgsidx]))
    results = {'status': 'exploratory secondary checks selected after pilot; primary unchanged',
               'AP_max_negative_family': {'max_deficit_sigma': maxdeficit, 'z':float(z[np.argmin(zz)]),
                    'raw_one_sided_p':float(norm.sf(maxdeficit)), 'six_contrast_global_p':family_p,
                    'least_favourable_mean':'zero contrasts, attained by coasting',
                    'interpretation':'family adjustment over six contrasts only, not all analysis choices or upstream systematic models'},
               'finite_q_history': {'best':best, 'all_starts':fits,
                    'q_intervals_z':np.r_[0,z].tolist(), 'q_bounds':[0,5],
                    'interpretation':'conditional smooth-in-log-g no-acceleration fit; no likelihood-ratio significance assigned'},
               'BGS_holdout':{'z':float(bz),'predicted_DV_rd':bv,'observed_DV_rd':bvalue,
                    'measurement_sd':bsd,'difference_in_measurement_sd':(bvalue-bv)/bsd,
                    'interpretation':'plug-in residual only; excludes expansion-fit uncertainty; not a calibrated predictive p-value'}}
    (out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    fig, axs = plt.subplots(1,2,figsize=(10.5,4.2),layout='constrained')
    axs[0].errorbar(z,y[n:],yerr=np.sqrt(np.diag(C)[n:]),fmt='o',label='DESI DR2 radial BAO')
    grid=np.linspace(.01,z.max(),600); _,gg=predictions(best['parameters'],t,np.log1p(grid))
    axs[0].plot(grid,gg,label='Nonaccelerating finite-q fit')
    axs[0].set(xlabel='Redshift z',ylabel=r'$(1+z)D_H/r_d$',title='No acceleration requires a nonincreasing curve')
    axs[0].legend(fontsize=8)
    axs[1].errorbar(z,ap,yerr=sd,fmt='o')
    axs[1].axhline(0,color='black',lw=1)
    axs[1].set(xlabel='Redshift z',ylabel=r'$D_M/r_d-\ln(1+z)(1+z)D_H/r_d$',title='Flat no-acceleration AP contrasts must be ≥ 0')
    fig.suptitle('BAO-only diagnostic: fixed released Gaussian likelihood, constant ruler')
    fig.savefig(out/'bao-shape.png',dpi=170)
    fig.savefig(out/'bao-shape.pdf')
    plt.close(fig)
    outputs=[out/'results.json',out/'bao-shape.png',out/'bao-shape.pdf']
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    (out/'manifest.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),
        'input_hashes':hashes,'executed_source_sha256':hashlib.sha256(SOURCE).hexdigest(),
        'dependency_sha256':sha(Path(__file__).with_name('bao_shape.py')),
        'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in outputs}},indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
