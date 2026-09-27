#!/usr/bin/env python3
"""Independent numerical checks of the parallel light-curve/cosmology study."""
from pathlib import Path
import datetime
import hashlib
import json
import sys
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor, cho_solve
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT))
from lib.ages import mu as adaptive_mu
from lib.flux_engine import Engine, BUNDLE
from lib.flux_fit import Case


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def native_flux_check():
    engine = Engine()
    records = []
    for file in sorted(BUNDLE.glob('objective_*.npz')):
        cid = file.stem.split('_')[1]
        case = Case(engine, cid)
        for name, delta in [('baseline', np.zeros(4)),
                            ('perturbed', np.array([-.09*np.log(10)/2.5, .15, .03, .15]))]:
            x = case.x0 + delta
            p = x.copy()
            p[3] += case.tref
            custom = case.flux(x, interpolation='scipy')
            native = engine.sncosmo_flux(p, case.b, case.t, case.z, case.ebv)
            diff = np.abs(custom-native)
            records.append(dict(cid=cid, branch=name,
                max_relative_flux_difference=float(np.max(diff/np.maximum(np.abs(custom), 1))),
                max_absolute_flux_difference=float(diff.max()),
                max_difference_in_sigma_flux=float(np.max(diff/np.sqrt(np.diag(case.cov))))))
    assert max(r['max_difference_in_sigma_flux'] for r in records) < .01
    return records


def cosmology_check():
    work = ROOT/'.work/age-recovery'
    d = pd.read_csv(work/'observed-design.csv')
    draws = pd.read_csv(work/'cosmology-draws.csv')
    raw = np.loadtxt(ROOT/'data/distances/Pantheon+SH0ES_STAT+SYS.cov')
    n = int(raw[0])
    c = raw[1:].reshape(n,n)
    ids = d.source_row.to_numpy()
    c = c[np.ix_(ids, ids)]
    c = (c+c.T)/2
    cf = cho_factor(c)
    l = np.linalg.cholesky(c)
    a = d.age.to_numpy()
    a = a-a.mean()
    z, zh = d.zHD.to_numpy(), d.zHEL.to_numpy()
    noise = np.random.default_rng(202609271).normal(size=(500,len(a)))
    baseline = adaptive_mu(z, zh, om=.3)
    checks = []
    for model in ['intercept_only', 'refit_width_colour_mass']:
        x = np.ones((len(a),1))
        if model != 'intercept_only':
            x = np.column_stack([x, d.x1, d.c, (d.HOST_LOGMASS>10).astype(float)])
        for slope in [0., -.03]:
            y = baseline+slope*a+l@noise[0]
            saved = draws[(draws.model==model)&(draws.injected_slope==slope)&(draws.draw==0)].iloc[0]
            for joint in [False, True]:
                design = np.column_stack([x,a]) if joint else x
                ix = cho_solve(cf, design)
                normal = np.linalg.inv(design.T@ix)
                def fit(omega):
                    r = y-adaptive_mu(z, zh, om=omega)
                    beta = normal@ix.T@r
                    residual = r-design@beta
                    return float(residual@cho_solve(cf,residual)), beta
                sol = minimize_scalar(lambda o:fit(o)[0],bounds=(.001,.999),
                    method='bounded', options={'xatol':1e-11})
                assert sol.success
                _, beta = fit(sol.x)
                target = float(saved['joint_Om' if joint else 'standard_Om'])
                checks.append(dict(model=model, slope=slope, joint_age=joint,
                    adaptive_quadrature_Omega_m=float(sol.x), reference_run_Omega_m=target,
                    absolute_Omega_m_difference=float(abs(sol.x-target)),
                    joint_age_difference=None if not joint else float(beta[-1]-saved.joint_age)))
    assert max(x['absolute_Omega_m_difference'] for x in checks)<2e-7
    assert max(abs(x['joint_age_difference'] or 0) for x in checks)<1e-9
    return checks


def main():
    output=ROOT/'studies/host_ages/results/age_recovery/independent-checks.json'
    inputs=[ROOT/'.work/age-recovery/observed-design.csv',ROOT/'.work/age-recovery/cosmology-draws.csv',
        ROOT/'data/distances/Pantheon+SH0ES_STAT+SYS.cov',BUNDLE/'exact_prior_setup.json',
        *sorted(BUNDLE.glob('objective_*.npz')),*sorted((ROOT/'data/des/model').glob('*'))]
    code=[Path(__file__),ROOT/'lib/ages.py',ROOT/'lib/flux_engine.py',ROOT/'lib/flux_fit.py']
    report=dict(checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        scope='Independent arithmetic crosscheck only. Native sncosmo Model.bandflux is compared with the custom photon integral; cosmology uses adaptive quadrature and unwhitened covariance profiling, separate from the tested cubic distance spline. Known model assumptions are unchanged.',
        native_flux=native_flux_check(), cosmology=cosmology_check(),passed=True,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},
        code_sha256={str(p.relative_to(ROOT)):sha(p) for p in code},
        command='OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/age_recovery/independent_checks.py')
    output.write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('Independent native-flux and cosmology checks passed')

if __name__=='__main__':
    main()
