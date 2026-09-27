"""Fresh GPU-target independence chains, reusing the frozen CPU sampler.

freeze/describe only inspect identities and complete source rows. Native
initialization and sampling are confined to the explicit sample action.
"""
import argparse
import datetime
import importlib
import json
from pathlib import Path

from independence_proposal import ROOT, FrozenMixture, digest, freeze
from independence_run import prior_start, runtime

HERE = Path(__file__).resolve().parent
DESIGN = HERE/'independence-gpu-design.json'
VALIDATION = ROOT/'studies/unified_cosmology/results/inference/independence-gpu-validation.json'


def dependencies():
    design = json.loads(DESIGN.read_text())
    for name, expected in design['frozen_dependency_sha256'].items():
        assert digest(ROOT/name) == expected, f'Frozen dependency changed: {name}'
    source = dict(design['frozen_dependency_sha256'])
    for path in [Path(__file__), DESIGN]:
        source[str(path.relative_to(ROOT))] = digest(path)
    inherited = importlib.import_module('cobaya.samplers.mcmc.mcmc')
    source[str(Path(inherited.__file__).resolve())] = digest(inherited.__file__)
    return source


def configuration_from_parent(parent, identify=True):
    """GPU-first factory and identity; no model construction or density call."""
    settings = parent['arguments']
    assert settings.get('gpu') is True and settings.get('fast_lensing') is True
    assert settings['model'] == 'cpl' and settings['evolution'] in ('linear', 'smooth01')
    assert settings['sample'] == 'dovekie' and settings['calibration'] == 'official_planck'
    from modern_gpu import configuration, identify as identify_gpu
    from late_geometry import sample_path
    scientific = {k: settings[k] for k in ['model', 'evolution', 'sample', 'calibration', 'surrogate']}
    info = configuration(scientific['model'], scientific['evolution'], scientific['sample'],
                         scientific['calibration'], Path(scientific['surrogate']))
    if identify:
        identity = identify_gpu(info, sample_path(scientific['sample']), Path(scientific['surrogate']))
        assert identity == parent['target_identity'], 'Current GPU target differs from frozen source.'
    else:
        identity = None
    return info, identity, scientific


def base_proposal(folder):
    from measurement_summary import verify_current_target
    from independence_proposal import DESIGN as BASE_DESIGN
    folder = Path(folder).resolve()
    training = json.loads((folder/'proposal.json').read_text())
    assert digest(BASE_DESIGN) == training['design_sha256']
    assert digest(HERE/'independence_proposal.py') == training['source_sha256']
    for name, item in training['snapshot_files'].items():
        assert digest(folder/name) == item['sha256'], f'Snapshot changed: {name}'
    parents = [json.loads((folder/f'run-{rank}.json').read_text()) for rank in range(4)]
    parent = parents[0]
    assert all(p['target_identity'] == parent['target_identity'] for p in parents)
    assert parent['target_identity']['identity'] == training['parent_target_identity']
    assert parent['arguments'] == training['parent_settings']
    verify_current_target(parent)
    mixture = FrozenMixture.read(folder/'proposal.npz', training['proposal_sha256'])
    return training, parent, mixture


def freeze_gpu(source, output):
    source, output = Path(source).resolve(), Path(output).resolve()
    bound = dependencies()
    parent = json.loads((source/'run-0.json').read_text())
    configuration_from_parent(parent)
    freeze(source, output)
    training, copied_parent, mixture = base_proposal(output)
    info, identity, scientific = configuration_from_parent(copied_parent)
    assert mixture.names == [n for n,v in info['params'].items() if isinstance(v,dict) and 'prior' in v]
    sidecar = {'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'scope': 'Frozen GPU proposal only; no density evaluation or scientific qualification.',
               'source_folder': str(source.relative_to(ROOT)), 'target_identity': identity['identity'],
               'scientific_settings': scientific, 'gpu': True, 'fast_lensing': True,
               'proposal_record_sha256': digest(output/'proposal.json'),
               'proposal_sha256': training['proposal_sha256'], 'dependency_sha256': bound}
    assert dependencies() == bound
    (output/'gpu-proposal.json').write_text(json.dumps(sidecar,indent=2)+'\n')
    return sidecar


def load_gpu_proposal(folder, identify=True):
    folder = Path(folder).resolve()
    bound = dependencies(); sidecar = json.loads((folder/'gpu-proposal.json').read_text())
    assert sidecar['dependency_sha256'] == bound
    assert sidecar['gpu'] is True and sidecar['fast_lensing'] is True
    assert digest(folder/'proposal.json') == sidecar['proposal_record_sha256']
    training, parent, mixture = base_proposal(folder)
    assert training['proposal_sha256'] == sidecar['proposal_sha256']
    info, identity, scientific = configuration_from_parent(parent, identify)
    assert parent['target_identity']['identity'] == sidecar['target_identity']
    assert scientific == sidecar['scientific_settings']
    assert mixture.names == [n for n,v in info['params'].items() if isinstance(v,dict) and 'prior' in v]
    return info, identity, scientific, parent, mixture, sidecar


def sampler_options(mixture, proposal_file, proposal_sha256, seed, max_samples=None):
    options = {'proposal_file': str(proposal_file), 'proposal_sha256': proposal_sha256,
               'covmat': mixture.covariance, 'covmat_params': mixture.names, 'blocking': [[1, mixture.names]],
               'learn_proposal': False, 'drag': False, 'measure_speeds': False,
               'oversample_power': 0, 'oversample_thin': False, 'temperature': 1,
               'seed': seed, 'burn_in': '10d', 'learn_every': '20d',
               'Rminus1_stop': .005, 'Rminus1_cl_stop': .05, 'Rminus1_cl_level': .95,
               'max_tries': '200d', 'output_every': '30s'}
    if max_samples is not None:
        assert max_samples > 0
        options['max_samples'] = max_samples
    return options


def sample(proposal_folder, output, seed=None, initialization_seed=None, max_samples=None):
    from mpi4py import MPI
    rank, size = MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size
    assert size == 4, 'Four independent chains required.'
    environment = runtime(); bound = dependencies()
    checked = json.loads(VALIDATION.read_text())
    assert checked['status'] == 'passed_configuration_only_no_likelihood_calls'
    assert checked['driver_dependency_sha256'] == bound
    folder, output = Path(proposal_folder).resolve(), Path(output).resolve()
    info, identity, scientific, parent, mixture, sidecar = load_gpu_proposal(folder, identify=rank==0)
    identity = MPI.COMM_WORLD.bcast(identity, root=0)
    assert identity == parent['target_identity']
    defaults = json.loads(DESIGN.read_text())['default_seeds'][scientific['evolution']]
    seed = defaults['sampling'] if seed is None else seed
    initialization_seed = defaults['initialization'] if initialization_seed is None else initialization_seed
    assert seed != initialization_seed
    initial = prior_start(mixture,info['params'],initialization_seed,rank)
    if rank == 0:
        assert not output.exists(), 'Fresh output only; no resumption or overwrite.'
        output.mkdir(parents=True)
    MPI.COMM_WORLD.Barrier()
    from independence_runtime import native_initialize
    initialization = native_initialize(info)
    options = sampler_options(mixture,folder/'proposal.npz',sidecar['proposal_sha256'],seed,max_samples)
    settings = dict(scientific,fast_lensing=True,gpu=True,resume=False,sampler='frozen_independence_gpu',
                    seed=seed,initialization_seed=initialization_seed,proposal_folder=str(folder),
                    proposal_sha256=sidecar['proposal_sha256'],runtime=environment)
    from mpi_metadata import install
    manifest = {'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'arguments':settings,'rank':rank,'MPI_size':size,'OMP_NUM_THREADS':'1',
                'target_identity':identity,'surrogate_sha256':identity['surrogate_sha256'],
                'initialization':initial,'native_numerical_initialization':initialization,
                'sampler_source_sha256':bound,'proposal_record_sha256':digest(folder/'proposal.json'),
                'gpu_proposal_record_sha256':digest(folder/'gpu-proposal.json'),
                'configuration_validation_sha256':digest(VALIDATION),
                'sampler_settings':{k:(v.tolist() if hasattr(v,'tolist') else v) for k,v in options.items()},
                'seed_sequence':{'root_entropy':seed,'spawn_key':[rank]},
                'metadata_policy':install(output/'chain'),
                'qualification':'numerical_proposal_only_requires_unchanged_diagnostics_and_native_correction'}
    (output/f'run-{rank}.json').write_text(json.dumps(manifest,indent=2)+'\n')
    for name,value in initial['point'].items():info['params'][name]['ref']=value
    from independence_sampler import IndependenceMCMC
    info['sampler']={'independence_sampler.IndependenceMCMC':options}
    info['output']=str(output/'chain');info['timing']=True
    from cobaya.run import run
    run(info)
    assert dependencies() == bound
    load_gpu_proposal(folder)


def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='action',required=True)
    f=sub.add_parser('freeze');f.add_argument('source',type=Path);f.add_argument('output',type=Path)
    d=sub.add_parser('describe');d.add_argument('proposal_folder',type=Path)
    s=sub.add_parser('sample');s.add_argument('proposal_folder',type=Path);s.add_argument('output',type=Path)
    s.add_argument('--seed',type=int);s.add_argument('--initialization-seed',type=int);s.add_argument('--max-samples',type=int)
    a=p.parse_args()
    if a.action=='freeze':print(json.dumps(freeze_gpu(a.source,a.output),indent=2))
    elif a.action=='describe':
        _,identity,scientific,_,mixture,sidecar=load_gpu_proposal(a.proposal_folder)
        print(json.dumps({'target_identity':identity['identity'],'scientific_settings':scientific,
                          'parameters':mixture.names,'proposal_sha256':sidecar['proposal_sha256'],'sampling_started':False},indent=2))
    else:sample(a.proposal_folder,a.output,a.seed,a.initialization_seed,a.max_samples)


if __name__=='__main__':main()
