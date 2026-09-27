"""Sample the unchanged modern target with a single correlated Metropolis block.

This avoids extrapolation-triggered startup timing measurements selecting an
excessively expensive dragging schedule for a usually fast spectrum proposal.
Published-chain covariance is an initialization aid only, never a prior.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import numpy as np
from modern_run import configuration,ROOT
from late_geometry import sample_path
from target_identity import identify


def initial_reference(info, names, mean, covariance, seed):
    """Draw a finite-support starting point without clipping or changing priors."""
    rng = np.random.default_rng(seed)
    rejected = []
    sampled = [name for name in names if isinstance(info['params'].get(name), dict)
               and 'prior' in info['params'][name]]
    for attempt in range(10000):
        delta = rng.multivariate_normal(np.zeros(len(names)), covariance)
        candidate = {key: float(mean[key] + d) for key, d in zip(names, delta)
                     if key in sampled}
        outside = []
        for key, value in candidate.items():
            prior = info['params'][key]['prior']
            if not np.isfinite(value) or ('min' in prior and value <= prior['min']) \
                    or ('max' in prior and value >= prior['max']):
                outside.append(key)
        if candidate.get('w', -1.) + candidate.get('wa', 0.) > 0:
            outside.append('CAMB_w_plus_wa')
        if not outside:
            return {'accepted_attempt': attempt + 1, 'point': candidate,
                    'rejected_support_constraints': rejected,
                    'scope': 'Initialization distribution only; physical prior and target unchanged.'}
        rejected.append(outside)
    raise RuntimeError('Could not draw an initial point inside the declared support.')


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--model',choices=['lcdm','cpl'],default='cpl')
    p.add_argument('--evolution',choices=['none','linear','smooth01','smooth03'],default='none')
    p.add_argument('--sample',choices=['dovekie','pantheon','des3yr'],default='dovekie')
    p.add_argument('--calibration',choices=['official_planck','paper_literal'],default='official_planck')
    p.add_argument('--surrogate',type=Path,required=True)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--fast-lensing',action='store_true')
    p.add_argument('--proposal-scale',type=float,default=2.4)
    p.add_argument('--seed',type=int,default=272652)
    a=p.parse_args()
    if a.fast_lensing:
        from modern_fast import configuration as target_configuration, identify as target_identify
    else:
        target_configuration, target_identify = configuration, identify
    from mpi4py import MPI
    rank=MPI.COMM_WORLD.rank
    modelhash=hashlib.sha256(a.surrogate.read_bytes()).hexdigest()
    method='fast-metropolis' if a.fast_lensing else 'metropolis'
    if a.proposal_scale!=2.4:
        method+='-scale'+str(a.proposal_scale).replace('.','p')
    if a.seed!=272652:
        method+='-seed'+str(a.seed)
    out=ROOT/f'.work/unified-cosmology/inference/modern-{a.model}-{a.evolution}-{a.sample}-{a.calibration}-{method}-{modelhash[:12]}'
    out.mkdir(parents=True,exist_ok=True)
    info=target_configuration(a.model,a.evolution,a.sample,a.calibration,a.surrogate)
    identity=target_identify(info,sample_path(a.sample),a.surrogate) if rank==0 else None
    identity=MPI.COMM_WORLD.bcast(identity,root=0)
    source=ROOT/'studies/unified_cosmology/results/external_probes/author-configuration.json'
    proposal=json.loads(source.read_text())['models'][a.model]['proposal_only']
    covariance=ROOT/proposal['path']
    assert hashlib.sha256(covariance.read_bytes()).hexdigest()==proposal['sha256']
    settings={k:str(v.resolve()) if isinstance(v,Path) else v for k,v in vars(a).items()}
    manifest={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'arguments':settings,'rank':rank,'MPI_size':MPI.COMM_WORLD.size,
        'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS'),
        'target_identity':identity,'surrogate_sha256':modelhash,
        'sampler_source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'proposal_only_sha256':hashlib.sha256(covariance.read_bytes()).hexdigest(),
        'author_starting_mean':proposal['transformed_mean'],
        'qualification':'numerical_proposal_only_requires_exact_correction'}
    from mpi_metadata import install
    manifest['metadata_policy']=install(out/'chain')
    manifest['metadata_guard_sha256']=hashlib.sha256((Path(__file__).parent/'mpi_metadata.py').read_bytes()).hexdigest()
    names=covariance.open().readline().lstrip('#').split()
    if not a.resume:
        manifest['initialization']=initial_reference(info,names,proposal['transformed_mean'],
                                                     np.loadtxt(covariance),a.seed+rank)
        for key,value in manifest['initialization']['point'].items():
            info['params'][key]['ref']=value
    mpath=out/f'run-{rank}.json'
    if mpath.exists():
        assert a.resume,'Existing run: use --resume.'
        assert json.loads(mpath.read_text())['target_identity']==identity
        with (out/f'run-history-{rank}.jsonl').open('a') as f:f.write(json.dumps(manifest)+'\n')
    else:mpath.write_text(json.dumps(manifest,indent=2)+'\n')
    sampled=[k for k,v in info['params'].items() if isinstance(v,dict) and 'prior' in v]
    info['sampler']={'mcmc':{'covmat':str(covariance),'drag':False,'blocking':[[1,sampled]],
        'proposal_scale':a.proposal_scale,
        'measure_speeds':False,'seed':a.seed,'learn_every':'20d','burn_in':'10d',
        'Rminus1_stop':.005,'Rminus1_cl_stop':.05,'max_tries':'200d','output_every':'30s'}}
    info['output']=str(out/'chain');info['timing']=True
    from cobaya.run import run
    run(info,resume=a.resume)


if __name__=='__main__':main()
