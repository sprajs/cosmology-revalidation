"""Fresh speed-blocked proposal chains; never modifies an existing live run.

Uses the unchanged modern target and Cobaya's native correlated block proposer.
All proposal chains still require the registered convergence and exact-CAMB
correction checks before scientific interpretation.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from modern_fast import configuration, identify
from modern_run import ROOT
from modern_sample import initial_reference
from late_geometry import sample_path

HERE = Path(__file__).resolve().parent
COSMOLOGY = ['H0', 'ombh2', 'omch2', 'logA', 'ns', 'tau', 'w', 'wa']
FAST = ['A_planck', 'P_act', 'Tcal', 'Ecal', 'A_fg', 'epsilon']


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parameter_blocks(info, fast_factor=1):
    if not isinstance(fast_factor, int) or fast_factor < 1:
        raise ValueError('Fast factor must be a positive integer.')
    names = [k for k, v in info['params'].items() if isinstance(v, dict) and 'prior' in v]
    assert set(names) <= set(COSMOLOGY+FAST), 'Unknown dependency: audit before adding a parameter.'
    slow = [n for n in names if n in COSMOLOGY]
    fast = [n for n in names if n in FAST]
    assert slow and fast and len(slow)+len(fast) == len(names)
    return [[1, slow], [fast_factor, fast]]


def sampler_options(info, covariance, seed, scale=1.6, fast_factor=1):
    return {'covmat': str(covariance), 'drag': False,
            'blocking': parameter_blocks(info, fast_factor),
            'proposal_scale': scale, 'measure_speeds': False,
            'oversample_thin': False, 'seed': seed,
            'learn_every': '20d', 'burn_in': '10d',
            'Rminus1_stop': .005, 'Rminus1_cl_stop': .05,
            'max_tries': '200d', 'output_every': '30s'}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model', choices=['lcdm', 'cpl'], default='cpl')
    p.add_argument('--evolution', choices=['none', 'linear', 'smooth01', 'smooth03'], default='none')
    p.add_argument('--sample', choices=['dovekie', 'pantheon', 'des3yr'], default='dovekie')
    p.add_argument('--calibration', choices=['official_planck', 'paper_literal'], default='official_planck')
    p.add_argument('--surrogate', type=Path, required=True)
    p.add_argument('--resume', action='store_true')
    p.add_argument('--proposal-scale', type=float, default=1.6)
    p.add_argument('--fast-factor', type=int, default=1)
    p.add_argument('--gpu', action='store_true', help='Use the separately validated exact GPU dispatch wrapper.')
    p.add_argument('--seed', type=int, default=272812)
    p.add_argument('--describe', action='store_true', help='Print sampler settings without launching or writing a run.')
    args = p.parse_args()
    if args.fast_factor < 1 or args.proposal_scale <= 0:
        p.error('Use a positive fast factor and proposal scale.')
    # Existing correction/measurement consumers identify the exact fast-lensing
    # factory from this flag. Blocking never enters the scientific target.
    args.fast_lensing = True
    args.surrogate = args.surrogate.resolve()
    target_configuration, target_identify = configuration, identify
    if args.gpu:
        from modern_gpu import configuration as target_configuration, identify as target_identify
    info = target_configuration(args.model, args.evolution, args.sample, args.calibration, args.surrogate)
    proposal_record_path = ROOT/'studies/unified_cosmology/results/external_probes/author-configuration.json'
    proposal_record = json.loads(proposal_record_path.read_text())['models'][args.model]['proposal_only']
    covariance = ROOT/proposal_record['path']
    assert digest(covariance) == proposal_record['sha256']
    options = sampler_options(info, covariance, args.seed, args.proposal_scale, args.fast_factor)
    if args.describe:
        print(json.dumps(options, indent=2)); return
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    modelhash = digest(args.surrogate)
    method = ('fast-gpu-blocked' if args.gpu else 'fast-blocked')
    method += f'-f{args.fast_factor}-scale{str(args.proposal_scale).replace(".", "p")}-seed{args.seed}'
    out = ROOT/f'.work/unified-cosmology/inference/modern-{args.model}-{args.evolution}-{args.sample}-{args.calibration}-{method}-{modelhash[:12]}'
    out.mkdir(parents=True, exist_ok=True)
    target = target_identify(info, sample_path(args.sample), args.surrogate) if rank == 0 else None
    target = MPI.COMM_WORLD.bcast(target, root=0)
    source_paths = [Path(__file__), HERE/'modern_sample.py', HERE/'mpi_metadata.py']
    sources = {str(path.relative_to(ROOT)): digest(path) for path in source_paths}
    settings = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items() if k != 'describe'}
    manifest = {'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'arguments': settings, 'rank': rank, 'MPI_size': MPI.COMM_WORLD.size,
                'OMP_NUM_THREADS': os.environ.get('OMP_NUM_THREADS'),
                'target_identity': target, 'surrogate_sha256': modelhash,
                'sampler_source_sha256': digest(Path(__file__)), 'sampler_dependency_sha256': sources,
                'proposal_only_sha256': digest(covariance),
                'author_starting_mean': proposal_record['transformed_mean'],
                'sampler_settings': options,
                'method': method, 'qualification': 'numerical_proposal_only_requires_exact_correction'}
    from mpi_metadata import install
    manifest['metadata_policy'] = install(out/'chain')
    manifest['metadata_guard_sha256'] = digest(HERE/'mpi_metadata.py')
    names = covariance.open().readline().lstrip('#').split()
    if not args.resume:
        manifest['initialization'] = initial_reference(info, names, proposal_record['transformed_mean'],
                                                       np.loadtxt(covariance), args.seed+rank)
        for key, value in manifest['initialization']['point'].items():
            info['params'][key]['ref'] = value
    path = out/f'run-{rank}.json'
    if path.exists():
        assert args.resume, 'Existing run: use --resume; never overwrite it.'
        previous = json.loads(path.read_text())
        assert previous['target_identity'] == target
        assert previous['sampler_dependency_sha256'] == sources, 'Sampler/helper source changed.'
        assert previous['sampler_settings'] == options, 'Resume would change sampler settings.'
        assert previous['proposal_only_sha256'] == manifest['proposal_only_sha256'], 'Proposal covariance bytes changed.'
        assert previous['MPI_size'] == MPI.COMM_WORLD.size, 'Resume MPI size changed.'
        with (out/f'run-history-{rank}.jsonl').open('a') as stream:
            stream.write(json.dumps(manifest)+'\n')
    else:
        assert not args.resume, 'Cannot resume a run without its original manifest.'
        path.write_text(json.dumps(manifest, indent=2)+'\n')
    info['sampler'] = {'mcmc': options}
    info['output'] = str(out/'chain'); info['timing'] = True
    from cobaya.run import run
    run(info, resume=args.resume)


if __name__ == '__main__':
    main()
