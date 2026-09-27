"""Dedicated anchored CPU proposal sampling; no original factory/manifest aliases."""
import argparse
import copy
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.metadata
import json
from pathlib import Path
import numpy as np

import anchored_adapter as adapter
from independence_proposal import FrozenMixture
from independence_run import prior_start, runtime
from independence_gpu_run import sampler_options
from target_identity import canonical, digest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
DESIGN=HERE/'anchored-sampling-design.json'
VALIDATION=ROOT/'studies/unified_cosmology/results/inference/anchored-sampling-validation.json'


def identity(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def relative(path):return str(Path(path).resolve().relative_to(ROOT))


def dump_new(path,value):
    path=Path(path)
    if path.exists():
        assert json.loads(path.read_text())==value,'Previously frozen payload changed.'
    else:
        path.parent.mkdir(parents=True,exist_ok=True)
        temporary=path.with_suffix('.part')
        temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');temporary.replace(path)


def verify_hashes(mapping):
    for path,expected in mapping.items():
        assert digest(ROOT/path)==expected,'Source/input identity changed: '+path


def dependencies():
    names=['anchored_sampling.py','anchored-sampling-design.json','anchored_adapter.py',
           'anchored-design.json','anchored_proposal.py','anchored-proposal-design.json','independence_sampler.py','independence_proposal.py',
           'independence-sampling-design.json','independence_run.py','independence_gpu_run.py',
           'independence_runtime.py','diagnostics.py','mpi_metadata.py']
    result={relative(HERE/name):digest(HERE/name) for name in names}
    inherited=importlib.import_module('cobaya.samplers.mcmc.mcmc')
    result[str(Path(inherited.__file__).resolve())]=digest(inherited.__file__)
    return result


def settings(model,surrogate,native_accuracy=1):
    design=json.loads(DESIGN.read_text())
    assert model in design['models'] and native_accuracy in design['native_accuracy_choices']
    return dict(model=model,evolution='none',sample=adapter.SAMPLE,calibration='official_planck',
                surrogate=str(Path(surrogate).resolve()),native_accuracy=int(native_accuracy),
                proposal_accuracy=1,fast_lensing=True,gpu=False,target_family='dedicated_anchored_v1')


def configuration(arguments,native=False):
    design=json.loads(DESIGN.read_text())
    expected=settings(arguments['model'],arguments['surrogate'],arguments['native_accuracy'])
    assert all(arguments[k]==v for k,v in expected.items()),'Unsupported anchored settings.'
    info=adapter.configuration(model=arguments['model'],surrogate=None if native else Path(arguments['surrogate']))
    if native:
        extra=info['theory']['camb']['extra_args']
        for key in design['accuracy_keys']:extra[key]=arguments['native_accuracy']
        assert extra['lens_potential_accuracy']==4
    return info


def identify(arguments):
    proposal=configuration(arguments)
    native=configuration(arguments,native=True)
    # Start with the complete original asset/version/calibration audit, then
    # bind the separately declared native numerical target and this factory.
    record=adapter.identify(proposal,adapter.DATA,Path(arguments['surrogate']))
    record.pop('identity')
    record['native_configuration']=canonical(native)
    record['native_accuracy']=arguments['native_accuracy']
    record['proposal_accuracy']=1
    record['target_family']='dedicated_anchored_v1'
    for path in [Path(__file__),DESIGN]:record['source_sha256'][relative(path)]=digest(path)
    record['identity']=identity(record)
    return record


def verify_manifest(manifest):
    assert manifest['schema']=='dedicated-anchored-chain-v1'
    frozen=manifest['target_identity'];arguments=manifest['arguments']
    assert manifest['validation']==validation_guard(),'Sampler prerequisite evidence changed.'
    verify_hashes(frozen['source_sha256'])
    assert manifest['sampler_source_sha256']==dependencies(),'Sampler dependency mapping changed or is incomplete.'
    verify_hashes(manifest['sampler_source_sha256'])
    _,proposal=load_proposal(Path(arguments['proposal_folder']),configuration(arguments))
    assert proposal==manifest['proposal'],'Frozen normalized proposal identity changed.'
    for name,version in frozen['versions'].items():assert importlib.metadata.version(name)==version
    assert identify(arguments)==frozen,'Current anchored target differs from sampled target.'
    inputs=dict(frozen['source_sha256'])
    inputs.update(frozen['anchored_calibration']['input_and_audit_sha256'])
    inputs.update(manifest['sampler_source_sha256'])
    inputs[relative(Path(arguments['surrogate']))]=frozen['surrogate_sha256']
    for filename in ['asset-file-inventory.json','modern-file-inventory.json']:
        p=ROOT/'.work/unified-cosmology/external-probes'/filename;inputs[relative(p)]=digest(p)
    return inputs


def load_proposal(folder,info):
    """Training is proposal evidence only; it is never a target likelihood."""
    folder=Path(folder).resolve();record=json.loads((folder/'proposal.json').read_text())
    if record.get('schema')=='anchored-bridge-weighted-proposal-v1':
        from anchored_proposal import load
        return load(folder,info)
    assert record.get('schema') in (None,'frozen-independence-proposal-v1'),'Unknown proposal schema.'
    assert digest(HERE/'independence-sampling-design.json')==record['design_sha256']
    assert digest(HERE/'independence_proposal.py')==record['source_sha256']
    for name,item in record['snapshot_files'].items():assert digest(folder/name)==item['sha256']
    mixture=FrozenMixture.read(folder/'proposal.npz',record['proposal_sha256'])
    names=[k for k,v in info['params'].items() if isinstance(v,dict) and 'prior'in v]
    assert mixture.names==names,'Proposal dimensions/order must match anchored target.'
    return mixture,{'proposal_file':str(folder/'proposal.npz'),'proposal_sha256':record['proposal_sha256'],
                    'proposal_record_path':relative(folder/'proposal.json'),
                    'proposal_record_sha256':digest(folder/'proposal.json'),
                    'training_target_identity':record['parent_target_identity'],
                    'training_settings':record['parent_settings'],
                    'role':'Normalized proposal only; the anchored density is evaluated at every proposal.'}


def validation_guard():
    record=json.loads(VALIDATION.read_text())
    assert record['status']=='passed_anchored_scaffold_synthetic_checks_no_inference'
    assert len(set(record['full_live_target_identities_checked']))==4
    assert record['sampling_dependency_sha256']==dependencies()
    verify_hashes(record['source_sha256'])
    return {'path':relative(VALIDATION),'sha256':digest(VALIDATION)}


def sample(proposal_folder,output,arguments,seed=None,initialization_seed=None,max_samples=None):
    from mpi4py import MPI
    rank,size=MPI.COMM_WORLD.rank,MPI.COMM_WORLD.size
    assert size==4,'Exactly four independent chains required.'
    environment=runtime();checked=validation_guard();bound=dependencies()
    info=configuration(arguments)
    target=identify(arguments) if rank==0 else None
    target=MPI.COMM_WORLD.bcast(target,root=0)
    mixture,proposal=load_proposal(proposal_folder,info)
    defaults=json.loads(DESIGN.read_text())['default_seeds'][arguments['model']]
    seed=defaults['sampling'] if seed is None else seed
    initialization_seed=defaults['initialization'] if initialization_seed is None else initialization_seed
    assert seed!=initialization_seed
    initial=prior_start(mixture,info['params'],initialization_seed,rank)
    output=Path(output).resolve();assert output.is_relative_to(ROOT/'.work')
    if rank==0:
        assert not output.exists(),'Fresh output only; no restart or overwrite.'
        output.mkdir(parents=True)
    MPI.COMM_WORLD.Barrier()
    from independence_runtime import native_initialize
    initialization=native_initialize(info)
    run_arguments=dict(arguments,seed=seed,initialization_seed=initialization_seed,
                       sampler='frozen_independence',resume=False,proposal_folder=str(Path(proposal_folder).resolve()))
    manifest={'schema':'dedicated-anchored-chain-v1','created_utc':datetime.now(timezone.utc).isoformat(),
              'rank':rank,'MPI_size':size,'arguments':run_arguments,'target_identity':target,
              'runtime':environment,'sampler_source_sha256':bound,'proposal':proposal,
              'initialization':initial,'native_numerical_initialization':initialization,
              'validation':checked,'seed_sequence':{'root_entropy':seed,'spawn_key':[rank]},
              'qualification':'Proposal sampling only; four-chain diagnostics and dedicated native correction required.'}
    from mpi_metadata import install
    manifest['metadata_policy']=install(output/'chain')
    dump_new(output/f'run-{rank}.json',manifest)
    for name,value in initial['point'].items():info['params'][name]['ref']=value
    info['sampler']={'independence_sampler.IndependenceMCMC':sampler_options(
        mixture,proposal['proposal_file'],proposal['proposal_sha256'],seed,max_samples)}
    info['output']=str(output/'chain');info['timing']=True
    from cobaya.run import run
    run(info)
    assert dependencies()==bound
    verify_manifest(manifest)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['describe','sample'])
    p.add_argument('--model',choices=['lcdm','cpl'],required=True)
    p.add_argument('--surrogate',type=Path,required=True)
    p.add_argument('--native-accuracy',type=int,choices=[1,2],default=1)
    p.add_argument('--proposal-folder',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--seed',type=int);p.add_argument('--initialization-seed',type=int)
    p.add_argument('--max-samples',type=int)
    a=p.parse_args();arguments=settings(a.model,a.surrogate,a.native_accuracy)
    if a.action=='describe':
        target=identify(arguments)
        dump_new(a.output,{'status':'described_only_no_native_or_background_calls','arguments':arguments,'target_identity':target})
        print(target['identity'])
    else:
        assert a.proposal_folder is not None
        sample(a.proposal_folder,a.output,arguments,a.seed,a.initialization_seed,a.max_samples)

if __name__=='__main__':main()
