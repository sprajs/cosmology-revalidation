"""No-likelihood checks of blocked production configuration and initialization."""
import hashlib
import inspect
import json
import logging
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from modern_blocked_sample import (ROOT, HERE, configuration, parameter_blocks,
    sampler_options, initial_reference, COSMOLOGY, FAST)
from blocked_benchmark import blocking
from target_identity import canonical


def main():
    from cobaya.model import Model
    from cobaya.samplers.mcmc.proposal import BlockedProposer
    fixture = ROOT/'.work/unified-cosmology/inference/surrogate/cubic-0509.npz'
    author = json.loads((ROOT/'studies/unified_cosmology/results/external_probes/author-configuration.json').read_text())
    rows = []; initial = []
    for family in ['lcdm', 'cpl']:
        for evolution in ['none', 'linear', 'smooth01', 'smooth03']:
            info = configuration(family, evolution, 'dovekie', 'official_planck', fixture)
            before = canonical(info)
            names = [k for k, v in info['params'].items() if isinstance(v, dict) and 'prior' in v]
            for factor in [1, 4]:
                blocks = parameter_blocks(info, factor)
                assert blocks == blocking(names, factor)
                dummy = SimpleNamespace(sampled_dependence=dict.fromkeys(names), log=logging.getLogger('validation'))
                checked, checked_factors = Model.check_blocking(dummy, blocks)
                assert list(checked) == [x[1] for x in blocks]
                assert checked_factors.tolist() == [1, factor]
                options = sampler_options(info, 'synthetic-covariance-path', 272812, fast_factor=factor)
                assert options['Rminus1_stop'] == .005 and options['Rminus1_cl_stop'] == .05
                assert options['oversample_thin'] is False and options['measure_speeds'] is False
                assert options['drag'] is False and options['proposal_scale'] == 1.6
                rng = np.random.default_rng(272830)
                matrix = rng.normal(size=(len(names), len(names)))
                covariance = matrix@matrix.T+np.eye(len(names))
                proposer = BlockedProposer([[names.index(k) for k in b] for _, b in blocks], rng,
                    oversampling_factors=[1, factor], proposal_scale=1.6)
                proposer.set_covariance(covariance)
                fast_checks = 0
                for _ in range(100):
                    vector = np.zeros(len(names)); proposer.get_proposal(vector)
                    if proposer.current_iblock == 1:
                        assert all(vector[i] == 0 for i, n in enumerate(names) if n in COSMOLOGY)
                        fast_checks += 1
                assert fast_checks > 0
                assert canonical(info) == before, 'Sampler configuration mutated the target.'
                rows.append({'model': family, 'evolution': evolution, 'factor': factor,
                             'slow': len(blocks[0][1]), 'fast': len(blocks[1][1]),
                             'zero_slow_change_fast_proposals': fast_checks})
            record = author['models'][family]['proposal_only']
            path = ROOT/record['path']; covnames = path.open().readline().lstrip('#').split()
            if family == 'cpl' and evolution in ['linear', 'smooth01']:
                seed = {'linear': 272812, 'smooth01': 272813}[evolution]
                for rank in range(4):
                    outcome = initial_reference(info, covnames, record['transformed_mean'], np.loadtxt(path), seed+rank)
                    assert outcome['point']['w']+outcome['point']['wa'] <= 0
                    initial.append({'evolution': evolution, 'rank': rank, 'seed': seed+rank,
                                    'accepted_attempt': outcome['accepted_attempt']})
    info['params']['unreviewed_new_parameter'] = {'prior': {'min': 0, 'max': 1}}
    try:
        parameter_blocks(info)
    except AssertionError:
        pass
    else:
        raise AssertionError('Unaudited parameter dependency was silently allowed.')
    paths = [Path(__file__), HERE/'modern_blocked_sample.py', HERE/'modern_sample.py',
             HERE/'mpi_metadata.py', HERE/'blocked_benchmark.py', Path(inspect.getfile(BlockedProposer))]
    result = {'status': 'passed_no_likelihood_sampler_validation', 'CMB_spectrum_calls': 0,
              'likelihood_evaluations': 0, 'configurations': rows, 'planned_initializations': initial,
              'unreviewed_parameter_dependency_rejected': True,
              'source_sha256': {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
              'scope': 'Configuration, actual native proposer transforms and finite-support initialization only; no chain was launched or qualified.'}
    (ROOT/'studies/unified_cosmology/results/inference/blocked-sampler-validation.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'configurations': len(rows), 'planned_starts': len(initial)}))


if __name__ == '__main__':
    main()
