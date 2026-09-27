"""Frozen independence MH with Cobaya's original RLE and convergence equations."""
import json
from pathlib import Path
import time

import numpy as np
from cobaya import mpi
from cobaya.samplers.mcmc.mcmc import MCMC

from independence_proposal import FrozenMixture, digest


def rank_gate(collections, names, sampled):
    """Match diagnostics.py on complete collections, without producing estimates."""
    from diagnostics import metrics
    chains = []
    for rows in collections:
        weights = rows[:, 0]
        assert (weights > 0).all() and np.array_equal(weights, weights.astype(int))
        expanded = np.repeat(rows[:, 1:], weights.astype(int), axis=0)
        chains.append(expanded[int(.3*len(expanded)):])
    length = min(map(len, chains))
    if len(chains) != 4 or length < 100:
        return {'passed': False, 'reason': 'four chains and at least100 retained draws required',
                'equal_chain_length': length}
    values = np.array([c[-length:] for c in chains])
    checks = {}; failed = {}
    for index, name in enumerate(names):
        if np.ptp(values[:, :, index]) == 0:
            if name in sampled:
                failed[name] = {'reason': 'sampled coordinate is constant'}
            continue
        m = metrics(values[:, :, index]); checks[name] = m
        if not (m['rhat'] <= 1.01 and min(m['bulk_ess'], m['tail_ess']) >= 400):
            failed[name] = m
    return {'passed': not failed, 'failed': failed, 'diagnostics': checks,
            'equal_chain_length': length, 'discard_fraction_after_sampler_burnin': .3}


class IndependenceMCMC(MCMC):
    """Replace only the proposal transition; no target or covariance adaptation."""
    proposal_file: str = ''
    proposal_sha256: str = ''
    file_base_name = 'mcmc'

    def initialize(self):
        assert not self.output.is_resuming(), 'Fresh runs only; never resume an unrecorded RNG state.'
        assert self.temperature in (None, 1)
        assert not self.learn_proposal and not self.drag and not self.measure_speeds
        assert self.oversample_power == 0 and not self.oversample_thin
        assert self.Rminus1_stop == .005 and self.Rminus1_cl_stop == .05
        assert self.Rminus1_cl_level == .95
        assert mpi.size() == 4, 'Four independent MPI chains required.'
        self.mixture = FrozenMixture.read(self.proposal_file, self.proposal_sha256)
        assert list(self.model.parameterization.sampled_params()) == self.mixture.names
        assert not self.model.prior._periodic_bounds, 'Periodic coordinates require a different density.'
        self._defer_checkpoint = False
        self._candidate_index = 0
        self._ledger_path = Path(self.output.add_suffix(f'independence-candidates-{mpi.rank()}.jsonl'))
        assert not self._ledger_path.exists(), 'Never append to an existing candidate ledger.'
        super().initialize()
        assert self.current_point.output_thin == 1
        self._current_logq = self.mixture.logpdf(self.current_point.values)
        from independence_runtime import numerical_runtime
        Path(self.output.add_suffix(f'independence-runtime-{mpi.rank()}.json')).write_text(
            json.dumps(numerical_runtime(), indent=2)+'\n')

    def get_new_sample_metropolis(self):
        previous = self.current_point.values.copy()
        previous_logp = float(self.current_point.logpost)
        previous_logq = self._current_logq
        previous_weight = int(self.current_point.weight)
        burn_before = int(self.burn_in_left)
        trial, component = self.mixture.draw(self._rng)
        logq = self.mixture.logpdf(trial)
        started = time.perf_counter()
        result = self.model.logposterior(trial)
        seconds = time.perf_counter()-started
        logp = float(result.logpost)
        logu = float(np.log(self._rng.random()))
        logalpha = min(0., logp-previous_logp+previous_logq-logq)
        accept = bool(np.isfinite(logp) and logu < logalpha)
        self.process_accept_or_reject(accept, trial, result)
        if accept:
            self._current_logq = logq
        row = {'index': self._candidate_index, 'component': component, 'trial': trial.tolist(),
               'trial_logq': logq, 'trial_logpost': logp if np.isfinite(logp) else None,
               'previous': previous.tolist(), 'previous_logq': previous_logq,
               'previous_logpost': previous_logp, 'previous_weight': previous_weight,
               'log_uniform': logu, 'log_acceptance': logalpha if np.isfinite(logalpha) else None,
               'accepted': accept, 'prior_rejected': bool(result.logprior == -np.inf),
               'burn_before': burn_before, 'collection_rows_after': len(self.collection),
               'seconds': seconds}
        with self._ledger_path.open('a') as stream:
            stream.write(json.dumps(row, allow_nan=False)+'\n')
        self._candidate_index += 1
        return accept

    def write_checkpoint(self):
        # The parent method normally writes a checkpoint before our extra gate.
        # Suppress that transient status; publish only the final joint decision.
        if not getattr(self, '_defer_checkpoint', False):
            super().write_checkpoint()

    def check_convergence_and_learn_proposal(self):
        assert digest(self.proposal_file) == self.proposal_sha256, 'Frozen proposal changed during sampling.'
        before = self.Rminus1_last
        self._defer_checkpoint = True
        try:
            super().check_convergence_and_learn_proposal()
        finally:
            self._defer_checkpoint = False
        original_passed = bool(self.converged)
        gate = None
        if original_passed:
            from diagnostics import PARAMETERS
            self.collection.out_update()
            # Read the fully flushed numeric text, including its released precision,
            # so the eventual unchanged consumer sees exactly these segments.
            path = Path(self.collection.file_name)
            headers = path.open().readline().lstrip('#').split()
            names = sorted((set(self.mixture.names) | set(PARAMETERS)) & set(headers))
            saved = np.loadtxt(path, ndmin=2)
            assert len(saved) == len(self.collection)
            rows = saved[:, [headers.index(n) for n in ['weight', *names]]]
            gathered = mpi.gather((names, rows))
            if mpi.is_main_process():
                assert all(n == names for n, _ in gathered)
                gate = rank_gate([r for _, r in gathered], names, self.mixture.names)
            gate = mpi.share(gate)
            self.converged = bool(gate['passed'])
        if mpi.is_main_process():
            record = {'check': self.i_learn, 'original_multivariate_previous': float(before) if np.isfinite(before) else None,
                      'original_multivariate_current': float(self.Rminus1_last) if np.isfinite(self.Rminus1_last) else None,
                      'original_convergence_passed': original_passed, 'rank_gate': gate,
                      'converged': bool(self.converged), 'proposal_sha256': self.proposal_sha256}
            bound = self.progress.loc[self.i_learn].get('Rminus1_cl', np.nan)
            record['original_bound_statistic'] = float(bound) if np.isfinite(bound) else None
            path = Path(self.output.add_suffix('independence-convergence.jsonl'))
            with path.open('a') as stream:
                stream.write(json.dumps(record, allow_nan=False)+'\n')
        self.write_checkpoint()

    def run(self):
        super().run()
        # Inherited stopping occurs on an accepted transition. All completed
        # post-burn holds/rejections have been written; only new terminal state
        # (weight1) is outside Cobaya's completed-hold sample collection.
        assert self.current_point.weight == 1
        assert digest(self.proposal_file) == self.proposal_sha256
        final = {'converged': bool(self.converged), 'trials': self._candidate_index,
                 'terminal_state': self.current_point.values.tolist(),
                 'terminal_weight': int(self.current_point.weight),
                 'terminal_state_in_collection': False,
                 'scope': 'Standard completed-hold Cobaya collection; no post-burn rejected hold omitted.',
                 'ledger_sha256': digest(self._ledger_path)}
        Path(self.output.add_suffix(f'independence-terminal-{mpi.rank()}.json')).write_text(json.dumps(final, indent=2)+'\n')
