"""Profile the late-time CPL ridge and disclose conditional prior-cut effects.

This diagnostic follows the initial nested-sampler discrepancy. The optional
cuts are posterior conditioning exercises, not the baseline physical model or
substitutes for a measured CMB likelihood.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import argparse
import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from late_geometry import Geometry, ROOT, BAO, SN
from late_nested import sha, summary, RESULTS, WORK

DESIGN = RESULTS/'late-tail-design.json'
SEEDS = [927088, 927089]


def design():
    record = {
        'scope': 'Exploratory diagnostic following the first nested-run disagreement. A profile is not posterior probability or evidence.',
        'omega_m_grid': [.010001, .02, .05, .1, .15, .2, .25, .28, .30, .31416, .33, .35, .4, .5, .7, .9],
        'starts': [[-.85, -.5], [-.9, 1.], [-1.5, 2.]],
        'dark_energy_bounds': [[-3., 1.], [-5., 3.]],
        'scale_bounds': [5000., 15000.],
        'method': 'Profile H0*rdrag analytically using the quadratic in inverse scale, then optimize w0/wa from three declared starts at every fixed Omega_m.',
        'conditional_cuts': ['none', 'Omega_m>0.1', 'w0+wa<=0', 'Omega_m>0.1 and w0+wa<=0'],
        'cut_scope': 'Post-hoc sensitivity motivated by the initial tail. Apply separately to both refined seeds; neither cut changes the baseline run or adds observed information. The w0+wa cut alone is not a CMB likelihood.',
        'code_sha256': sha(__file__),
        'geometry_sha256': sha(Path(__file__).with_name('late_geometry.py')),
    }
    if DESIGN.exists():
        assert json.loads(DESIGN.read_text()) == record
    else:
        DESIGN.write_text(json.dumps(record, indent=2)+'\n')
    return record


def profile():
    d = design(); g = Geometry(); started = time.monotonic()
    constant = float(g.bmean@g.bprecision@g.bmean)

    def objective(wpars, om, detail=False):
        w, wa = wpars
        ova, wva, ava = np.array([om]), np.array([w]), np.array([wa])
        curve = 5*np.log10(g.integrals(g.nodes, ova, wva, ava)[0])
        sn = float(curve@g.Q@curve-2*curve@g.b+g.dd)
        integral = g.bz*g.integrals(g.bz, ova, wva, ava)[0]
        expansion = g.expansion(g.bz[None, :, None], ova, wva, ava)[0, :, 0]
        dh = 299792.458/expansion
        dm = 299792.458*integral
        dv = np.cbrt(g.bz*dh*dm*dm)
        vector = np.where(g.btype=='DH_over_rs', dh, np.where(g.btype=='DM_over_rs', dm, dv))
        aa = float(vector@g.bprecision@vector)
        bb = float(vector@g.bprecision@g.bmean)
        # y=1/(H0*rdrag); constrained positive quadratic minimizer.
        inverse = np.clip(bb/aa, 1/d['scale_bounds'][1], 1/d['scale_bounds'][0])
        scale = 1/inverse
        bao = aa*inverse**2-2*bb*inverse+constant
        if detail:
            return np.array([om, scale, w, wa]), sn+bao
        return sn+bao

    rows = []
    for om in d['omega_m_grid']:
        fits = [minimize(objective, start, args=(om,), method='Nelder-Mead',
            bounds=d['dark_energy_bounds'], options={'maxiter':10000, 'xatol':1e-9, 'fatol':1e-9})
            for start in d['starts']]
        best = min(fits, key=lambda fit:fit.fun)
        theta, chi2 = objective(best.x, om, True)
        direct = sum(float(component[0]) for component in g.components(theta))
        assert abs(chi2-direct)<1e-7
        rows.append({'omega_m':om, 'theta':theta.tolist(), 'chi2':float(chi2),
            'all_optimizer_success':[bool(f.success) for f in fits],
            'all_start_chi2':[float(f.fun) for f in fits],
            'scale_quadratic_closure':float(chi2-direct)})
    # Independent full-data distance quadrature at four ridge positions.
    quadrature = []
    for index in [0, 3, 5, 9]:
        theta = np.array(rows[index]['theta'])
        error = float(g.components(theta)[0][0]-g.direct_chi2(theta))
        assert abs(error)<1e-4
        quadrature.append({'row':index, 'SN_chi2_error':error})
    lowest = min(row['chi2'] for row in rows)
    for row in rows:
        row['delta_chi2_from_best_profile_grid'] = row['chi2']-lowest
    record = {
        'status':'completed' if all(all(row['all_optimizer_success']) for row in rows) else 'optimizer_failures',
        'rows':rows, 'direct_distance_quadrature':quadrature,
        'seconds':time.monotonic()-started,
        'code_sha256':sha(__file__),
        'dependencies_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [DESIGN, Path(__file__).with_name('late_geometry.py'), Path(__file__).with_name('late_nested.py'), SN, BAO/'desi_gaussian_bao_ALL_GCcomb_mean.txt', BAO/'desi_gaussian_bao_ALL_GCcomb_cov.txt']},
        'qualification':'Best of three starts at each fixed matter density, not a proof of global optimization or posterior mass. The lower-matter ridge reaches the declared prior boundary. No early-time data enter this calculation.',
    }
    (RESULTS/'late-tail-profile.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    print(json.dumps(record))


def conditional():
    design(); result = {}
    inputs = [DESIGN, Path(__file__).with_name('late_nested.py')]
    for seed in SEEDS:
        arrays = WORK/f'refined-seed-{seed}'/'samples.npz'
        compact = RESULTS/f'late-nested-refined-{seed}.json'
        original = json.loads(compact.read_text())
        assert original['single_run_diagnostic_passed']
        assert original['outputs_sha256'][str(arrays.relative_to(ROOT))] == sha(arrays)
        data = np.load(arrays); x, weights = data['samples'], data['weights']
        low_cut = x[:, 0]>.1
        early_cut = x[:, 2]+x[:, 3]<=0
        masks = {'none':np.ones(len(x), dtype=bool), 'Omega_m>0.1':low_cut,
                 'w0+wa<=0':early_cut, 'Omega_m>0.1 and w0+wa<=0':low_cut&early_cut}
        rows = {}
        for name, keep in masks.items():
            mass = float(weights[keep].sum()); conditioned = weights[keep]/mass
            rows[name] = {'retained_baseline_probability':mass,
                'weighted_ESS':float(1/(conditioned@conditioned)),
                'posterior':summary(x[keep], conditioned)}
        result[str(seed)] = rows
        inputs.extend([arrays, compact])
    record = {'status':'computed_conditional_sensitivities', 'seeds':result,
        'code_sha256':sha(__file__),
        'dependencies_sha256':{str(p.relative_to(ROOT)):sha(p) for p in inputs},
        'qualification':'Post-hoc conditional prior sensitivity. Each seed remains separate; use the refined comparison to assess numerical agreement. This adds no physical observation and provides no CMB measurement.'}
    (RESULTS/'late-tail-conditional.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    print(json.dumps(record))


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--freeze',action='store_true')
    parser.add_argument('--conditional',action='store_true'); args=parser.parse_args()
    if args.freeze: print(json.dumps(design()))
    elif args.conditional: conditional()
    else: profile()
