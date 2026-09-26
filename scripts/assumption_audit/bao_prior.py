"""Bounded sensitivity to the CPL wa lower prior limit in BAO-only inference.

This is an alternative prior experiment, not a corrected physical posterior.
All other settings and all released data/covariances are held fixed.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone

import emcee
import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/cosmology'))
from core import BAO, Pantheon, CODE_AT_START

SCRIPT_AT_START = Path(__file__).read_bytes()
CORE_AT_START = CODE_AT_START['scripts/cosmology/core.py']


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--wa-min', type=float, required=True)
    p.add_argument('--seed', type=int, required=True)
    p.add_argument('--data', choices=['bao', 'joint', 'joint-c14fixed'], default='bao')
    p.add_argument('--steps', type=int, default=14000)
    p.add_argument('--burn', type=int, default=3000)
    a = p.parse_args()
    b = BAO()
    sn = Pantheon() if a.data != 'bao' else None
    extra_inputs = []
    if sn and a.data == 'joint-c14fixed':
        import pandas as pd
        from scipy.interpolate import PchipInterpolator
        template = ROOT / 'runs/cosmology/templates/c14-cpl63-median.csv'
        d = pd.read_csv(template)
        sn.reset_template(PchipInterpolator(d.z, d.delta_mu, extrapolate=False)(sn.z))
        extra_inputs.append(template)
    bounds = np.array([[.01, .99], [-3, 1], [a.wa_min, 2], [5000, 15000]])
    def logp(t):
        t = np.atleast_2d(t)
        ok = np.all((t > bounds[:, 0]) & (t < bounds[:, 1]), axis=1) & (t[:, 1] + t[:, 2] < 0)
        lp = np.full(len(t), -np.inf)
        if np.any(ok):
            lp[ok] = -.5 * b.chisq(t[ok, :3], t[ok, 3])
            if sn:
                lp[ok] -= .5 * sn.chisq(t[ok, :3], 'cpl', amplitude=float(a.data == 'joint-c14fixed'))
        return lp
    rng = np.random.default_rng(a.seed)
    np.random.seed(a.seed)
    opt = minimize(lambda x: -logp(x)[0], [.38, -.2, -2.7, 9130], method='Nelder-Mead',
                   options={'maxiter': 15000, 'fatol': 1e-9, 'xatol': 1e-7})
    assert opt.success
    start = []
    while len(start) < 40:
        s = opt.x + rng.normal(size=4) * .01 * np.diff(bounds, axis=1).ravel()
        if np.isfinite(logp(s)[0]):
            start.append(s)
    sampler = emcee.EnsembleSampler(40, 4, logp, vectorize=True)
    sampler.run_mcmc(np.array(start), a.steps, progress=False)
    chain = sampler.get_chain(discard=a.burn)
    flat = chain.reshape(-1, 4)
    q = .5 + 1.5 * flat[:, 1] * (1 - flat[:, 0])
    tau = sampler.get_autocorr_time(discard=a.burn, tol=0)
    qs = q.reshape(chain.shape[:2])
    blocks = np.array_split(qs, 25)
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(), settings=vars(a),
                  parameter_names=['Om', 'w0', 'wa', 'H0_rd'], prior_bounds=bounds.tolist(),
                  extra_prior='w0 + wa < 0',
                  posterior_mean=flat.mean(axis=0).tolist(), posterior_sd=flat.std(axis=0).tolist(),
                  q0_mean=float(q.mean()), q0_sd=float(q.std()),
                  p_q0_negative=float(np.mean(q < 0)),
                  q0_mean_block_mcse=float(np.std([v.mean() for v in blocks], ddof=1) / np.sqrt(25)),
                  p_q0_negative_block_mcse=float(np.std([(v < 0).mean() for v in blocks], ddof=1) / np.sqrt(25)),
                  q0_percentiles=np.quantile(q, [.025, .16, .5, .84, .975]).tolist(),
                  tau=tau.tolist(), retained_steps_per_tau=((a.steps-a.burn)/tau).tolist(),
                  acceptance_fraction=float(sampler.acceptance_fraction.mean()),
                  boundary_fraction={n:float(np.mean((flat[:, i] - bounds[i, 0] < .01 * (bounds[i,1]-bounds[i,0])) |
                                                   (bounds[i,1]-flat[:, i] < .01 * (bounds[i,1]-bounds[i,0]))))
                                     for i,n in enumerate(['Om','w0','wa','H0_rd'])},
                  posterior_wa_below_original_lower_bound=float(np.mean(flat[:, 2] < -3)),
                  mode=opt.x.tolist(), mode_chisq=float(2 * opt.fun))
    out = ROOT / 'runs/assumption_audit/bao' / f'{a.data}-wa{a.wa_min:g}-seed{a.seed}'
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / 'chains.npz', chain=chain)
    (out/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
    def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
    archive = out.parent / 'provenance'; archive.mkdir(exist_ok=True)
    for data in [SCRIPT_AT_START, CORE_AT_START]:
        (archive / (hashlib.sha256(data).hexdigest() + '.py')).write_bytes(data)
    inputs = {str(p.relative_to(ROOT)):sha(p) for p in [*b.inputs, *(sn.inputs if sn else []), *extra_inputs]}
    inputs['scripts/cosmology/core.py'] = hashlib.sha256(CORE_AT_START).hexdigest()
    (out/'manifest.json').write_text(json.dumps(dict(script_sha256=hashlib.sha256(SCRIPT_AT_START).hexdigest(),
        code_capture='At import; exact bytes in sibling provenance directory',
        inputs_sha256=inputs,
        outputs_sha256={name:sha(out/name) for name in ['summary.json','chains.npz']}), indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
