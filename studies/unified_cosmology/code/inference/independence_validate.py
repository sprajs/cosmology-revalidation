"""Independent density, actual four-rank sampler and completed-hold validation."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
from scipy.special import logsumexp
from scipy.stats import multivariate_normal, multivariate_t, truncnorm
from independence_proposal import DESIGN, ROOT, FrozenMixture, digest, relative

HERE = Path(__file__).resolve().parent
WORK = ROOT/'.work/unified-cosmology/inference/independence-validation'


def gaussian_loglike(x, y):
    return float(-.5*((x-.1)/.4)**2-.5*((y+.1)/.6)**2)


def synthetic(folder):
    from mpi4py import MPI
    from cobaya.run import run
    from independence_sampler import IndependenceMCMC
    from mpi_metadata import install
    rank = MPI.COMM_WORLD.rank
    assert MPI.COMM_WORLD.size == 4
    if rank == 0:
        assert not folder.exists()
        folder.mkdir(parents=True)
        np.savez(folder/'proposal.npz', mean=[.2, -.3], cov=[[.64, .1], [.1, .81]], names=['x', 'y'])
    MPI.COMM_WORLD.Barrier()
    install(folder/'chain')
    info = {'params': {'x': {'prior': {'min': -2., 'max': 2.}, 'ref': [-.6, .7, -.2, .4][rank]},
                       'y': {'prior': {'min': -3., 'max': 3.}, 'ref': [.8, -.7, .3, -.4][rank]}},
            'likelihood': {'known_gaussian': {'external': gaussian_loglike}},
            'sampler': {'independence_sampler.IndependenceMCMC': {              'proposal_file': str(folder/'proposal.npz'), 'proposal_sha256': digest(folder/'proposal.npz'),
              'covmat': [[.64, .1], [.1, .81]], 'covmat_params': ['x', 'y'], 'blocking': [[1, ['x', 'y']]],
              'learn_proposal': False, 'drag': False, 'measure_speeds': False,
              'oversample_power': 0, 'oversample_thin': False, 'temperature': 1,
              'seed': 273242, 'burn_in': 20, 'learn_every': 200, 'Rminus1_stop': .005,
              'Rminus1_cl_stop': .05, 'Rminus1_cl_level': .95, 'max_tries': 2000,
              'max_samples': 10000, 'output_every': '1s'}}, 'output': str(folder/'chain')}
    # Minimal synthetic manifest only exercises the unchanged diagnostics parser.
    manifest = {'rank': rank, 'target_identity': {'configuration': {'params': info['params']}}}
    (folder/f'run-{rank}.json').write_text(json.dumps(manifest, indent=2)+'\n')
    run(info)


def unit_checks():
    from independence_sampler import rank_gate, IndependenceMCMC
    from independence_run import prior_start, runtime
    rng = np.random.default_rng(273243)
    q = FrozenMixture([.2, -.3], [[1., .2], [.2, .7]], ['x', 'y'])
    x = rng.normal(size=(1000, 2))*3
    independently = np.logaddexp(np.log(.9)+multivariate_normal.logpdf(x, q.mean, q.covariance*1.05**2),
                                 np.log(.1)+multivariate_t.logpdf(x, q.mean, q.covariance*2.4, df=5))
    density_error = float(max(abs(q.logpdf(x)-independently)))
    assert density_error < 1e-12
    logp = multivariate_normal.logpdf(x, [-.1, .1], [[.5, -.1], [-.1, .6]])
    logq = q.logpdf(x)
    forward = logp[:-1]+logq[1:]+np.minimum(0, np.diff(logp)-np.diff(logq))
    reverse = logp[1:]+logq[:-1]+np.minimum(0, -np.diff(logp)+np.diff(logq))
    balance_error = float(max(abs(forward-reverse))); assert balance_error < 1e-11
    parameter = {'x': {'prior': {'min': -2, 'max': 2}}, 'y': {'prior': {'min': -3, 'max': 3}}}
    starts = [prior_start(q, parameter, 273241, rank) for rank in range(4)]
    assert len({tuple(s['point'].values()) for s in starts}) == 4
    assert starts == [prior_start(q, parameter, 273241, rank) for rank in range(4)]
    # Dense tails retain support rather than silently truncating the proposal.
    assert np.isfinite(q.logpdf([1e6, -1e6]))
    good = [np.column_stack([np.ones(3000), rng.normal(size=(3000, 2))]) for _ in range(4)]
    assert rank_gate(good, ['x', 'y'], ['x', 'y'])['passed']
    shifted = [r.copy() for r in good]; shifted[3][:, 1] += .5
    assert not rank_gate(shifted, ['x', 'y'], ['x', 'y'])['passed']
    repeated = [np.column_stack([np.full(20, 100), rng.normal(size=(20, 2))]) for _ in range(4)]
    assert not rank_gate(repeated, ['x', 'y'], ['x', 'y'])['passed']
    # Ensure the wrapper's guarded checkpoint never publishes a provisional true.
    dummy = object.__new__(IndependenceMCMC); dummy._defer_checkpoint = True
    dummy.write_checkpoint()  # no parent resources exist: any attempted write fails
    # Initialization guards fail before touching any model/likelihood resources.
    from types import SimpleNamespace
    for altered in [{'temperature': 2}, {'learn_proposal': True}, {'drag': True}]:
        forbidden = object.__new__(IndependenceMCMC)
        forbidden._output = SimpleNamespace(is_resuming=lambda: False)
        forbidden.temperature = 1; forbidden.learn_proposal = False
        forbidden.drag = False; forbidden.measure_speeds = False
        for key, value in altered.items(): setattr(forbidden, key, value)
        try: forbidden.initialize()
        except AssertionError: pass
        else: raise AssertionError('Forbidden adaptation/temperature accepted.')
    runtime()
    from unittest.mock import patch
    with patch.dict(os.environ, {'CLIPY_NOJAX': '0'}):
        try: runtime()
        except AssertionError: pass
        else: raise AssertionError('Wrong numerical backend accepted.')
    path = WORK/'tamper.npz'; WORK.mkdir(parents=True, exist_ok=True)
    np.savez(path, mean=q.mean, cov=q.covariance, names=q.names); expected = digest(path)
    original = path.read_bytes(); path.write_bytes(original+b'changed')
    try: FrozenMixture.read(path, expected)
    except AssertionError: pass
    else: raise AssertionError('Changed proposal accepted.')
    path.write_bytes(original)
    return {'density_max_error': density_error, 'detailed_balance_log_error': balance_error,
            'four_distinct_reproducible_initialization_streams': True,
            'shifted_and_low_effective_sample_controls_rejected': True,
            'deferred_checkpoint_cannot_publish_early_convergence': True,
            'changed_proposal_and_backend_rejected': True,
            'nonunit_temperature_and_adaptation_and_dragging_rejected': True}


def audit(folder):
    from diagnostics import check_mpi
    units = unit_checks()
    all_rows = []; all_ledgers = []; per_chain = []
    for rank in range(4):
        path = folder/f'chain.{rank+1}.txt'
        columns = path.open().readline().lstrip('#').split(); raw = np.loadtxt(path, ndmin=2)
        ledgerpath = folder/f'chain_independence-candidates-{rank}.jsonl'
        rows = [json.loads(line) for line in ledgerpath.read_text().splitlines()]
        terminal = json.loads((folder/f'chain_independence-terminal-{rank}.json').read_text())
        assert terminal['ledger_sha256'] == digest(ledgerpath)
        assert len(rows) == terminal['trials']
        expected = []
        q = FrozenMixture.read(folder/'proposal.npz', digest(folder/'proposal.npz'))
        last = None
        for index, row in enumerate(rows):
            assert row['index'] == index
            assert abs(q.logpdf(row['trial'])-row['trial_logq']) < 1e-12
            assert abs(q.logpdf(row['previous'])-row['previous_logq']) < 1e-12
            if last is not None:
                assert row['previous'] == (last['trial'] if last['accepted'] else last['previous'])
                assert row['previous_weight'] == (1 if last['accepted'] else last['previous_weight']+1)
            x, y = row['trial']
            prior_ok = -2 <= x <= 2 and -3 <= y <= 3
            assert row['prior_rejected'] == (not prior_ok)
            if prior_ok:
                independently = gaussian_loglike(x, y)-np.log(4*6)
                assert abs(independently-row['trial_logpost']) < 1e-12
                ratio = min(0, independently-row['previous_logpost']+row['previous_logq']-row['trial_logq'])
                assert abs(ratio-row['log_acceptance']) < 1e-12
                assert row['accepted'] == (row['log_uniform'] < ratio)
            else:
                assert not row['accepted'] and row['trial_logpost'] is None
            if row['accepted'] and row['burn_before'] <= 0:
                expected.append([row['previous_weight'], -row['previous_logpost'], *row['previous']])
            last = row
        expected = np.asarray(expected)
        assert len(expected) == len(raw)
        observed = raw[:, [columns.index(n) for n in ['weight', 'minuslogpost', 'x', 'y']]]
        assert np.array_equal(observed[:, 0], expected[:, 0])
        assert np.allclose(observed[:, 1:], expected[:, 1:], rtol=6e-8, atol=6e-8)
        assert terminal['terminal_weight'] == 1 and rows[-1]['accepted']
        assert terminal['converged']
        all_rows.append(raw); all_ledgers.append(rows)
        per_chain.append({'rank': rank, 'candidate_count': len(rows), 'recorded_states': len(raw),
                          'recorded_frequency_draws': int(raw[:, 0].sum()),
                          'prior_rejections': sum(r['prior_rejected'] for r in rows),
                          'density_correct_MH_decisions_checked': len(rows)})
    assert sum(c['prior_rejections'] for c in per_chain) > 0
    assert len({tuple(r[0]['trial']) for r in all_ledgers}) == 4
    check = check_mpi(folder); assert check['status'] == 'passed', check['failed_gates']
    (folder/'unchanged-diagnostics.json').write_text(json.dumps(check, indent=2)+'\n')
    history = [json.loads(line) for line in (folder/'chain_independence-convergence.jsonl').read_text().splitlines()]
    final = history[-1]; assert final['converged'] and final['rank_gate']['passed']
    assert max(final['original_multivariate_previous'], final['original_multivariate_current']) < .005
    assert final['original_bound_statistic'] < .05
    # Extra final output diagnostics must equal the in-memory gate's segments.
    assert check['equal_chain_length_for_diagnostics'] == final['rank_gate']['equal_chain_length']
    for name in ['x', 'y']:
        assert check['diagnostics'][name] == final['rank_gate']['diagnostics'][name]
    means = {}; standards = {}
    # This is a synthetic target recovery, not an observational posterior.
    for name, center, sd, bounds in [('x', .1, .4, (-2, 2)), ('y', -.1, .6, (-3, 3))]:
        truth = truncnorm((bounds[0]-center)/sd, (bounds[1]-center)/sd, loc=center, scale=sd)
        observed = check['posterior'][name]
        means[name] = float(abs(observed['mean']-truth.mean()))
        standards[name] = float(abs(observed['sd']-truth.std()))
        assert means[name] < .04 and standards[name] < .04
    source_paths = [HERE/name for name in ['independence_proposal.py', 'independence_sampler.py',
                    'independence_run.py', 'independence_runtime.py', 'independence_validate.py', 'independence-sampling-design.json', 'diagnostics.py']]
    import importlib
    inherited = importlib.import_module('cobaya.samplers.mcmc.mcmc')
    source_paths.append(Path(inherited.__file__))
    result = {'status': 'passed_synthetic_frozen_independence_validation', 'observational_measurement': False,
              'unit_checks': units, 'MPI_chains': 4, 'per_chain': per_chain,
              'synthetic_mean_absolute_errors': means, 'synthetic_sd_absolute_errors': standards,
              'final_original_and_additional_convergence': final,
              'unchanged_diagnostics_passed': True, 'rank_gate_equals_flushed_final_diagnostics': True,
              'source_sha256': {str(p.resolve()): digest(p) for p in source_paths},
              'synthetic_artifacts_sha256': {relative(p): digest(p) for p in sorted(folder.iterdir()) if p.is_file()},
              'scope': 'Actual synthetic sampling, rejection weights, unchanged diagnostics and density controls; no real cosmology measurement.'}
    target = ROOT/'studies/unified_cosmology/results/inference/independence-sampler-validation.json'
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'per_chain': per_chain, 'mean_errors': means, 'sd_errors': standards}, indent=2))


def main():
    parser = argparse.ArgumentParser(); parser.add_argument('mode', choices=['synthetic', 'audit'])
    parser.add_argument('folder', type=Path); args = parser.parse_args()
    if args.mode == 'synthetic': synthetic(args.folder.resolve())
    else: audit(args.folder.resolve())


if __name__ == '__main__':
    main()
