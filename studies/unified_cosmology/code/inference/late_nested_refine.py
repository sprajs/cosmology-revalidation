"""Declared refinement after the first independent seeds failed a tail check.

The original records remain unchanged. The physical likelihood and priors are
identical; only the live population and constrained sampling kernel change.
"""
from __future__ import annotations
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
import argparse
import json
import time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import dynesty
from dynesty import utils as dyfunc
from scipy.special import logsumexp
from late_nested import Target, check_likelihood, summary, quantities, sha, ROOT, WORK, RESULTS, BAO, SN

SEEDS = [927088, 927089]
DESIGN = RESULTS / 'late-nested-refinement-design.json'


def design():
    record = {
        'frozen_utc': datetime.now(timezone.utc).isoformat(),
        'reason': 'The two 1000-live random-walk runs failed the declared marginal CDF agreement gate; low-Omega_m/high-wa tail mass differed. Preserve that failure.',
        'physical_target': 'Identical baseline flat matter+CPL Dovekie+DESI DR2 distance target; no radiation, CMB or w0+wa cut.',
        'bounds': [[.01, .99], [5000., 15000.], [-3., 1.], [-5., 3.]],
        'seeds': SEEDS, 'live_points': 5000, 'sampler': 'rslice', 'slices': 5,
        'bound': 'multi', 'dlogz': .03, 'maximum_calls': 8000000,
        'gate': 'Each run finishes before call cap with weighted ESS>=2000; maximum marginal CDF distance<=.05; logZ difference<=3 combined reported standard errors. Prior-volume jitter does not include missed modes.',
        'tail_diagnostics': 'Report Omega_m probability below .1,.15,.2,.25 and w0+wa>0. These descriptive cuts follow the initial disagreement, not new physical priors.',
        'input_sha256': {str(p.relative_to(ROOT)): sha(p) for p in [Path(__file__), Path(__file__).with_name('late_nested.py'), Path(__file__).with_name('late_geometry.py'), RESULTS/'late-nested-comparison.json']},
        'method_source': 'https://dynesty.readthedocs.io/en/latest/quickstart.html',
    }
    if DESIGN.exists():
        old = json.loads(DESIGN.read_text()); record['frozen_utc'] = old['frozen_utc']; assert old == record
    else:
        DESIGN.write_text(json.dumps(record, indent=2) + '\n')
    return record


def tails(x, w):
    return {**{f'P(Omega_m<{cut})': float(w @ (x[:, 0] < cut)) for cut in [.1, .15, .2, .25]},
            'P(w0+wa>0)': float(w @ (x[:, 2] + x[:, 3] > 0))}


def run(seed):
    d = design(); out = WORK/f'refined-seed-{seed}'; out.mkdir(parents=True, exist_ok=True)
    recordfile = RESULTS/f'late-nested-refined-{seed}.json'
    if recordfile.exists(): raise RuntimeError('Completed result exists; refuse silent replacement')
    started = time.monotonic(); target = Target(); assert target.g.bounds == d['bounds']
    checks = check_likelihood(target, seed); rng = np.random.default_rng(seed)
    sampler = dynesty.NestedSampler(target.loglike, target.prior, 4, nlive=d['live_points'],
        bound=d['bound'], sample=d['sampler'], slices=d['slices'], rstate=rng)
    sampler.run_nested(dlogz=d['dlogz'], maxcall=d['maximum_calls'], print_progress=False,
        checkpoint_file=str(out/'checkpoint.pkl'), checkpoint_every=60)
    r = sampler.results; weights = np.exp(r.logwt-logsumexp(r.logwt)); ess = float(1/(weights@weights))
    arrays = out/'samples.npz'
    np.savez_compressed(arrays, samples=r.samples, weights=weights, logwt=r.logwt,
        loglikelihood=r.logl, logz=r.logz, logzerr=r.logzerr, logvol=r.logvol, ncall=r.ncall)
    jitter = []
    for _ in range(64):
        jr = dyfunc.jitter_run(r, rstate=rng); jw = np.exp(jr.logwt-logsumexp(jr.logwt))
        js = summary(jr.samples, jw)
        jitter.append([*tails(jr.samples, jw).values(), js['w0']['mean'], js['wa']['mean'], js['q0']['mean']])
    np.savez_compressed(out/'shrinkage-jitter.npz', values=jitter)
    complete = int(np.sum(r.ncall)) < d['maximum_calls']
    record = dict(status='completed_run' if complete else 'call_cap_incomplete',
        single_run_diagnostic_passed=complete and ess>=2000, seed=seed, live_points=d['live_points'],
        iterations=int(r.niter), likelihood_calls=int(np.sum(r.ncall)), weighted_ESS=ess,
        seconds_including_setup_and_checks=time.monotonic()-started,
        logZ_unnormalized=float(r.logz[-1]), logZ_error=float(r.logzerr[-1]),
        posterior=summary(r.samples, weights), tail_probabilities=tails(r.samples, weights),
        shrinkage_jitter_labels=[*tails(r.samples, weights), 'mean_w0', 'mean_wa', 'mean_q0'],
        shrinkage_jitter_sd=np.std(jitter, axis=0, ddof=1).tolist(),
        independent_likelihood_checks=checks,
        code_sha256=sha(__file__), environment={'dynesty':dynesty.__version__, 'numpy':np.__version__},
        dependencies_sha256={str(p.relative_to(ROOT)):sha(p) for p in [DESIGN, Path(__file__).with_name('late_nested.py'), Path(__file__).with_name('late_geometry.py'), SN, BAO/'desi_gaussian_bao_ALL_GCcomb_mean.txt', BAO/'desi_gaussian_bao_ALL_GCcomb_cov.txt']},
        outputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in [arrays,out/'shrinkage-jitter.npz']},
        scope=d['physical_target'], evidence_scope='Unnormalized distance likelihood; not a Bayes-factor comparison across data/covariance models.')
    recordfile.write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')
    print(json.dumps(record), flush=True)


def compare():
    design(); paths=[RESULTS/f'late-nested-refined-{s}.json' for s in SEEDS]
    records=[json.loads(p.read_text()) for p in paths]
    arrays=[np.load(WORK/f'refined-seed-{s}'/'samples.npz') for s in SEEDS]
    distances={}; qa,qb=[quantities(a['samples']) for a in arrays]
    for name in qa:
        xa,xb=qa[name],qb[name]; ia,ib=np.argsort(xa),np.argsort(xb); grid=np.unique(np.r_[xa,xb])
        ca=np.r_[0,np.cumsum(arrays[0]['weights'][ia])]; cb=np.r_[0,np.cumsum(arrays[1]['weights'][ib])]
        distances[name]=float(np.max(abs(ca[np.searchsorted(xa[ia],grid,side='right')]-cb[np.searchsorted(xb[ib],grid,side='right')])))
    difference=abs(records[0]['logZ_unnormalized']-records[1]['logZ_unnormalized'])
    error=float(np.hypot(*[r['logZ_error'] for r in records]))
    passed=all(r['single_run_diagnostic_passed'] for r in records) and max(distances.values())<=.05 and difference<=3*error
    record=dict(status='passed_independent_seed_diagnostics' if passed else 'diagnostic_failed',
        seeds=SEEDS,maximum_marginal_CDF_differences=distances, absolute_logZ_difference=difference,
        combined_reported_logZ_error=error, tail_probabilities_by_seed={str(s):r['tail_probabilities'] for s,r in zip(SEEDS,records)},
        code_sha256=sha(__file__), source_records_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths},
        qualification='Finite seed agreement tests this numerical calculation, not absent remote modes or observational/model adequacy. Initial failed records remain separate.')
    (RESULTS/'late-nested-refined-comparison.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record))


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('--seed',type=int,choices=SEEDS)
    p.add_argument('--freeze',action='store_true'); p.add_argument('--compare',action='store_true'); a=p.parse_args()
    if a.freeze: print(json.dumps(design()))
    elif a.compare: compare()
    elif a.seed is not None: run(a.seed)
    else: p.error('Specify --freeze, --seed or --compare')
