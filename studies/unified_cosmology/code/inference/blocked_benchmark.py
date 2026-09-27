"""Bounded proposal-performance benchmark; never a cosmological measurement.

Use Cobaya's unchanged BlockedProposer and unchanged model.logposterior. The
short fixed-attempt pilots retain every rejection and do not qualify a chain.
The exported blocking() is a candidate for a separately launched native Cobaya
sampler, not a modification of any running sampler.
"""
import argparse
import hashlib
import inspect
import json
from pathlib import Path
import time

import numpy as np
from scipy.linalg import solve_triangular

from modern_fast import configuration, identify
from late_geometry import sample_path
from modern_run import ROOT, HERE

COSMOLOGY = ['H0', 'ombh2', 'omch2', 'logA', 'ns', 'tau', 'w', 'wa']
FAST = ['A_planck', 'P_act', 'Tcal', 'Ecal', 'A_fg', 'epsilon']
DESIGN = HERE/'blocked-benchmark-design.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def blocking(names, factor):
    """Slow first: triangular covariance preserves fast-only cache reuse."""
    assert set(names) <= set(COSMOLOGY+FAST)
    if factor == 0:
        return [[1, list(names)]]
    assert int(factor) == factor and factor >= 1
    return [[1, [n for n in names if n in COSMOLOGY]],
            [factor, [n for n in names if n in FAST]]]


def metrics(draws, names, seconds, slow_seconds=None):
    import arviz as az
    # Diagnostic only: one short chain cannot establish stationarity or mixing.
    draws = np.asarray(draws)
    tail = draws[len(draws)//4:]
    result = {}
    for i, name in enumerate(names):
        value = float(az.ess(tail[:, i][None, :], method='bulk'))
        result[name] = {'rank_bulk_ESS': value, 'ESS_per_CPU_second': value/seconds}
        if slow_seconds is not None:
            result[name]['ESS_per_modeled_second'] = value/slow_seconds
    return {'post_25_percent_rows': len(tail), 'parameters': result,
            'minimum_cosmology_ESS': min(result[n]['rank_bulk_ESS'] for n in names if n in COSMOLOGY),
            'minimum_fast_ESS': min(result[n]['rank_bulk_ESS'] for n in names if n in FAST)}


def kernel(evaluate, start, covariance, names, factor, steps, seed):
    from cobaya.samplers.mcmc.proposal import BlockedProposer
    rng = np.random.default_rng(seed)
    blocks = blocking(names, factor)
    proposer = BlockedProposer([[names.index(n) for n in b] for _, b in blocks], rng,
                               oversampling_factors=[f for f, _ in blocks], proposal_scale=1.6)
    proposer.set_covariance(covariance)
    current = np.array(start, copy=True)
    old = float(evaluate(current))
    assert np.isfinite(old), 'Nonfinite benchmark start.'
    draws = []; accepted = 0; slow = 0; fast = 0; nonfinite = 0
    wall = time.monotonic(); cpu = time.process_time()
    for _ in range(steps):
        proposed = current.copy()
        proposer.get_proposal(proposed)
        changed_slow = any(proposed[i] != current[i] for i, n in enumerate(names) if n in COSMOLOGY)
        slow += int(changed_slow); fast += int(not changed_slow)
        new = float(evaluate(proposed))
        assert not np.isnan(new), 'NaN benchmark density.'
        nonfinite += int(not np.isfinite(new))
        if np.log(rng.uniform()) < new-old:
            current = proposed; old = new; accepted += 1
        draws.append(current.copy())
    return np.array(draws), {'attempts': steps, 'accepted': accepted,
        'acceptance': accepted/steps, 'slow_proposals': slow, 'fast_proposals': fast,
        'nonfinite_density_proposals': nonfinite, 'wall_seconds': time.monotonic()-wall,
        'CPU_seconds': time.process_time()-cpu, 'blocking': blocks, 'seed': seed}


def proposal(names, info):
    source = ROOT/'studies/unified_cosmology/results/external_probes/author-configuration.json'
    record = json.loads(source.read_text())['models']['cpl']['proposal_only']
    path = ROOT/record['path']; assert digest(path) == record['sha256']
    oldnames = path.open().readline().lstrip('#').split(); oldcov = np.loadtxt(path)
    covariance = np.zeros((len(names), len(names)))
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if a in oldnames and b in oldnames:
                covariance[i, j] = oldcov[oldnames.index(a), oldnames.index(b)]
        if names[i] not in oldnames:
            covariance[i, i] = info['params'][names[i]]['proposal']**2
    start = np.array([record['transformed_mean'].get(n, info['params'][n]['ref']) for n in names])
    return start, covariance, {'author_record_sha256': digest(source), 'covariance_path': str(path.relative_to(ROOT)),
                              'covariance_sha256': digest(path), 'epsilon_covariance': 'Declared .02 proposal, no learned correlations.'}


def reference_check(model, start, names, covariance):
    theory = model.theory['spectral_surrogate']
    dependency = {n: [c.get_name() for c in components] for n, components in model.sampled_dependence.items()}
    base = model.logposterior(start, return_derived=False)
    assert np.isfinite(base.logpost)
    records = []
    for i, name in enumerate(names):
        p = start.copy(); p[i] += .15*np.sqrt(covariance[i, i])
        before = theory.surrogate_calls+theory.exact_calls
        begin = time.process_time()
        result = model.logposterior(p, return_derived=False)
        elapsed = time.process_time()-begin
        calls = theory.surrogate_calls+theory.exact_calls-before
        assert np.isfinite(result.logpost)
        if name in FAST:
            assert calls == 0, f'{name} unexpectedly recalculated cosmology.'
        # Uncached recomputation compares every likelihood and prior, not just sum.
        raw = model.logposterior(p, return_derived=False, cached=False)
        difference = max(float(np.max(abs(np.array(raw.loglikes)-result.loglikes))),
                         float(np.max(abs(np.array(raw.logpriors)-result.logpriors))))
        assert difference < 1e-9
        records.append({'parameter': name, 'CPU_seconds': elapsed, 'theory_calls': calls,
                        'cached_uncached_max_logdensity_difference': difference})
        model.logposterior(start, return_derived=False)
    return {'dependency': dependency, 'evaluations': records}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--surrogate', type=Path, required=True)
    p.add_argument('--work', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    a.work = a.work.resolve(); a.output = a.output.resolve(); a.surrogate = a.surrogate.resolve()
    a.work.mkdir(parents=True, exist_ok=True)
    assert a.work.resolve().is_relative_to(ROOT/'.work')
    assert not (a.work/'design.json').exists(), 'Preserve each benchmark; choose a fresh directory.'
    design = json.loads(DESIGN.read_text())
    (a.work/'design.json').write_text(json.dumps(design, indent=2)+'\n')
    import camb
    from cobaya.model import get_model
    from cobaya.samplers.mcmc.proposal import BlockedProposer
    # Abort, never reject, on a cap breach: truncating support would change target.
    original_get_results = camb.get_results
    calls = {'native_spectra': 0}
    def capped(*args, **kwargs):
        if calls['native_spectra'] >= design['maximum_native_spectrum_calls']:
            raise RuntimeError('Native spectrum cap exceeded; abort entire pilot.')
        calls['native_spectra'] += 1
        return original_get_results(*args, **kwargs)
    camb.get_results = capped
    sources = [Path(__file__), DESIGN, Path(inspect.getfile(BlockedProposer))]
    report = {'status': 'in_progress', 'scope': design['scope'], 'design': design,
              'source_sha256': {str(x): digest(x) for x in sources}, 'reference': {}, 'pilots': [], 'gaussian': []}
    try:
        for evolution in ['linear', 'smooth01']:
            info = configuration('cpl', evolution, 'dovekie', 'official_planck', a.surrogate)
            target = identify(info, sample_path('dovekie'), a.surrogate)
            info['debug'] = 40
            with get_model(info) as model:
                names = list(model.parameterization.sampled_params())
                start, covariance, inputs = proposal(names, info)
                report['reference'][evolution] = {'target_identity': target, 'proposal': inputs,
                                                  **reference_check(model, start, names, covariance)}
                if evolution != 'linear':
                    continue
                evaluate = lambda x: model.logposterior(x, return_derived=False).logpost
                theory = model.theory['spectral_surrogate']
                for seed in design['pilot_seeds']:
                    for factor in design['oversampling_factors']:
                        before = theory.surrogate_calls+theory.exact_calls
                        rows, result = kernel(evaluate, start, covariance, names, factor,
                                              design['pilot_attempts'], seed)
                        result['theory_calls'] = theory.surrogate_calls+theory.exact_calls-before
                        result['metrics'] = metrics(rows, names, result['CPU_seconds'])
                        result['evolution'] = evolution
                        rawpath = a.work/f'linear-factor{factor}-seed{seed}.npz'
                        np.savez_compressed(rawpath, rows=rows, names=names)
                        result['draws_path'] = str(rawpath.relative_to(ROOT)); result['draws_sha256'] = digest(rawpath)
                        report['pilots'].append(result)
                        (a.work/'progress.json').write_text(json.dumps(report, indent=2)+'\n')
                        print(json.dumps({'pilot': seed, 'factor': factor, 'seconds': result['wall_seconds'],
                                          'theory_calls': result['theory_calls']}), flush=True)
                # Known correlated Gaussian with the same dimensionality/covariance.
                chol = np.linalg.cholesky(covariance)
                slowcost = np.median([r['CPU_seconds'] for r in report['reference']['linear']['evaluations'] if r['parameter'] in COSMOLOGY])
                fastcost = np.median([r['CPU_seconds'] for r in report['reference']['linear']['evaluations'] if r['parameter'] in FAST])
                def gaussian(x):
                    u = solve_triangular(chol, x, lower=True, check_finite=False)
                    return -.5*u@u
                for seed in design['gaussian_seeds']:
                    initial = chol@np.random.default_rng(seed+1).normal(size=len(names))
                    for factor in design['oversampling_factors']:
                        rows, result = kernel(gaussian, initial, covariance, names, factor,
                                              design['gaussian_attempts'], seed)
                        modeled = result['slow_proposals']*slowcost+result['fast_proposals']*fastcost
                        result['modeled_seconds'] = modeled
                        result['cost_model'] = {'slow_CPU_seconds': float(slowcost), 'fast_CPU_seconds': float(fastcost)}
                        result['metrics'] = metrics(rows, names, result['CPU_seconds'], modeled)
                        result['maximum_mean_error_in_known_marginal_SD'] = float(np.max(abs(rows[len(rows)//4:].mean(axis=0))/np.sqrt(np.diag(covariance))))
                        rawpath = a.work/f'gaussian-factor{factor}-seed{seed}.npz'
                        np.savez_compressed(rawpath, rows=rows, names=names)
                        result['draws_path'] = str(rawpath.relative_to(ROOT)); result['draws_sha256'] = digest(rawpath)
                        report['gaussian'].append(result)
        report['status'] = 'completed_performance_only'
    except Exception as error:
        report['status'] = 'aborted_benchmark'; report['error'] = repr(error)
        raise
    finally:
        camb.get_results = original_get_results
        report['native_spectrum_calls'] = calls['native_spectra']
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
