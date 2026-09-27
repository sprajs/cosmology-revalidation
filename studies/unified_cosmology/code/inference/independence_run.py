"""Freeze a global proposal or start fresh baseline-compatible independent chains."""
import argparse
import datetime
import importlib.metadata
import json
import os
from pathlib import Path

import numpy as np
from independence_proposal import DESIGN, ROOT, FrozenMixture, digest, freeze, relative

HERE = Path(__file__).resolve().parent


def runtime():
    keys = ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'CLIPY_NOJAX']
    settings = {k: os.environ.get(k) for k in keys}
    assert all(v == '1' for v in settings.values()), 'Require baseline one-thread NumPy environment.'
    assert np.dtype(float).itemsize == 8
    return {'environment': settings, 'coordinate_numpy_dtype': 'float64',
            'third_party_internal_dtypes': 'Not inferred; preserve the existing target runtime settings.',
            'versions': {name: importlib.metadata.version(name) for name in ['numpy', 'scipy', 'cobaya', 'mpi4py']}}


def prior_start(mixture, parameters, seed, rank):
    """Initialization only: explicit independent stream, no state clipping."""
    rng = np.random.default_rng(np.random.SeedSequence(seed).spawn(4)[rank])
    rejected = []
    for index in range(10000):
        values, component = mixture.draw(rng)
        point = dict(zip(mixture.names, map(float, values)))
        outside = []
        for name, value in point.items():
            p = parameters[name]['prior']
            if not np.isfinite(value) or ('min' in p and value <= p['min']) or ('max' in p and value >= p['max']):
                outside.append(name)
        if point.get('w', -1.)+point.get('wa', 0.) > 0:
            outside.append('native_CAMB_CPL_domain')
        if not outside:
            return {'point': point, 'component': component, 'attempt': index+1,
                    'previous_rejections': rejected, 'seed_root': seed, 'spawn_key': [rank]}
        rejected.append(outside)
    raise RuntimeError('No valid independent initialization in declared10000-attempt limit.')


def sample(proposal_folder, output, seed=273240, initialization_seed=273241, max_samples=None):
    from mpi4py import MPI
    rank, size = MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size
    assert size == 4
    environment = runtime()
    proposal_folder, output = Path(proposal_folder).resolve(), Path(output).resolve()
    training = json.loads((proposal_folder/'proposal.json').read_text())
    assert digest(DESIGN) == training['design_sha256']
    assert digest(HERE/'independence_proposal.py') == training['source_sha256']
    for filename, item in training['snapshot_files'].items():
        assert digest(proposal_folder/filename) == item['sha256']
    parent = json.loads((proposal_folder/'run-0.json').read_text())
    for i in range(4):
        assert json.loads((proposal_folder/f'run-{i}.json').read_text())['target_identity'] == parent['target_identity']
    from measurement_summary import verify_current_target
    verify_current_target(parent)
    original = parent['arguments']
    assert original.get('fast_lensing') and not original.get('gpu'), 'CPU modern_fast snapshots only.'
    from modern_fast import configuration, identify
    from late_geometry import sample_path
    scientific = {k: original[k] for k in ['model', 'evolution', 'sample', 'calibration', 'surrogate']}
    info = configuration(scientific['model'], scientific['evolution'], scientific['sample'],
                         scientific['calibration'], Path(scientific['surrogate']))
    identity = identify(info, sample_path(scientific['sample']), Path(scientific['surrogate'])) if rank == 0 else None
    identity = MPI.COMM_WORLD.bcast(identity, root=0)
    assert identity == parent['target_identity']
    mixture = FrozenMixture.read(proposal_folder/'proposal.npz', training['proposal_sha256'])
    assert mixture.names == [k for k, v in info['params'].items() if isinstance(v, dict) and 'prior' in v]
    initial = prior_start(mixture, info['params'], initialization_seed, rank)
    if rank == 0:
        assert not output.exists(), 'Fresh output only; no overwrite or implicit resumption.'
        output.mkdir(parents=True)
    MPI.COMM_WORLD.Barrier()
    settings = dict(scientific, fast_lensing=True, gpu=False, resume=False,
                    sampler='frozen_independence', seed=seed, initialization_seed=initialization_seed,
                    proposal_folder=str(proposal_folder), proposal_sha256=training['proposal_sha256'],
                    runtime=environment)
    sources = [Path(__file__), HERE/'independence_sampler.py', HERE/'independence_proposal.py', DESIGN,
               HERE/'diagnostics.py', HERE/'mpi_metadata.py', HERE/'independence_runtime.py']
    import cobaya.samplers.mcmc.mcmc as inherited
    sources.append(Path(inherited.__file__))
    code = {str(p.resolve()): digest(p) for p in sources}
    from independence_runtime import native_initialize
    initialization = native_initialize(info)
    manifest = {'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'arguments': settings, 'rank': rank, 'MPI_size': size, 'OMP_NUM_THREADS': '1',
                'target_identity': identity, 'surrogate_sha256': identity['surrogate_sha256'],
                'initialization': initial, 'native_numerical_initialization': initialization,
                'sampler_source_sha256': code,
                'proposal_record_sha256': digest(proposal_folder/'proposal.json'),
                'seed_sequence': {'root_entropy': seed, 'spawn_key': [rank]},
                'qualification': 'numerical_proposal_only_requires_unchanged_diagnostics_and_native_correction'}
    from mpi_metadata import install
    manifest['metadata_policy'] = install(output/'chain')
    (output/f'run-{rank}.json').write_text(json.dumps(manifest, indent=2)+'\n')
    # Scientific identity was recorded before initialization/sampler/output aids,
    # matching the established immutable-manifest contract.
    for name, value in initial['point'].items():
        info['params'][name]['ref'] = value
    from independence_sampler import IndependenceMCMC
    options = {'proposal_file': str(proposal_folder/'proposal.npz'),
               'proposal_sha256': training['proposal_sha256'], 'covmat': mixture.covariance,
               'covmat_params': mixture.names, 'blocking': [[1, mixture.names]],
               'learn_proposal': False, 'drag': False, 'measure_speeds': False,
               'oversample_power': 0, 'oversample_thin': False, 'temperature': 1,
               'seed': seed, 'burn_in': '10d', 'learn_every': '20d',
               'Rminus1_stop': .005, 'Rminus1_cl_stop': .05, 'Rminus1_cl_level': .95,
               'max_tries': '200d', 'output_every': '30s'}
    if max_samples is not None:
        options['max_samples'] = max_samples
    info['sampler'] = {'independence_sampler.IndependenceMCMC': options}; info['output'] = str(output/'chain')
    info['timing'] = True
    from cobaya.run import run
    run(info)
    for path, expected in code.items():
        assert digest(path) == expected, 'Sampler or gate code changed during execution.'
    verify_current_target(parent)


def main():
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest='action', required=True)
    p = commands.add_parser('freeze')
    p.add_argument('folder', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--manifest-folder', type=Path)
    p = commands.add_parser('sample')
    p.add_argument('proposal_folder', type=Path); p.add_argument('output', type=Path)
    p.add_argument('--seed', type=int, default=273240)
    p.add_argument('--initialization-seed', type=int, default=273241)
    p.add_argument('--max-samples', type=int)
    args = parser.parse_args()
    if args.action == 'freeze':
        result = freeze(args.folder, args.output, args.manifest_folder)
        print(json.dumps({'proposal_sha256': result['proposal_sha256'], 'names': result['names']}))
    else:
        sample(args.proposal_folder, args.output, args.seed, args.initialization_seed, args.max_samples)


if __name__ == '__main__':
    main()
