"""Sample an explicit contemporary CMB+BAO+SN target or a numerical proposal.

Surrogate output is provisional until held-out checks, independent-chain
diagnostics and exact-CAMB importance correction have all been reported.
"""
import argparse
import datetime
import hashlib
import json
import os
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0,str(HERE.parent/'external_probes'))
from modern_adapter import modern_info
from likelihood import ReleasedDistances,ExpansionDiagnostics
from late_geometry import sample_path
from spectral_surrogate import replace_theory


def configuration(model='cpl',evolution='none',sample='dovekie',calibration='official_planck',surrogate=None):
    info = modern_info(model,calibration=calibration)
    info['likelihood']['released_sn'] = {'external':ReleasedDistances,
        'data_file':str(sample_path(sample)),
        'smooth_sigma':{'none':0.,'linear':0.,'smooth01':.1,'smooth03':.3}[evolution]}
    info['likelihood']['expansion_diagnostics'] = {'external':ExpansionDiagnostics}
    info['params']['epsilon'] = ({'prior':{'min':-.5,'max':.5},'ref':0.,'proposal':.02}
                                if evolution=='linear' else 0.)
    if model=='cpl':
        info['params']['w']['ref'] = -.85
        info['params']['wa']['ref'] = -.55
    if surrogate:
        info = replace_theory(info,surrogate)
    return info


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--model',choices=['lcdm','cpl'],default='cpl')
    p.add_argument('--evolution',choices=['none','linear','smooth01','smooth03'],default='none')
    p.add_argument('--sample',choices=['dovekie','pantheon','des3yr'],default='dovekie')
    p.add_argument('--calibration',choices=['official_planck','paper_literal'],default='official_planck')
    p.add_argument('--surrogate',type=Path)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--evaluate',action='store_true')
    p.add_argument('--seed',type=int,default=272627)
    p.add_argument('--max-samples',type=int)
    a = p.parse_args()
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    name = f'modern-{a.model}-{a.evolution}-{a.sample}-{a.calibration}'
    if a.surrogate:name += '-proposal-'+hashlib.sha256(a.surrogate.read_bytes()).hexdigest()[:12]
    out = ROOT/'.work/unified-cosmology/inference'/name
    out.mkdir(parents=True,exist_ok=True)
    info = configuration(a.model,a.evolution,a.sample,a.calibration,a.surrogate)
    from target_identity import identify
    identity = identify(info,sample_path(a.sample),a.surrogate) if rank==0 else None
    identity = MPI.COMM_WORLD.bcast(identity,root=0)
    manifest = {'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'arguments':{k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},
        'rank':rank,'MPI_size':MPI.COMM_WORLD.size,'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS'),
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                        for p in list(HERE.glob('*.py'))+[HERE.parent/'external_probes/modern_adapter.py',HERE.parent/'external_probes/adapter.py']},
        'qualification':'numerical_proposal_only' if a.surrogate else 'exact_declared_likelihood'}
    manifest['target_identity'] = identity
    if a.surrogate:
        manifest['surrogate_sha256'] = hashlib.sha256(a.surrogate.read_bytes()).hexdigest()
    mpath = out/f'run-{rank}.json'
    if mpath.exists():
        if not a.resume and not a.evaluate:raise FileExistsError('Use --resume; do not overwrite a run.')
        with (out/f'run-history-{rank}.jsonl').open('a') as f:f.write(json.dumps(manifest)+'\n')
    else:mpath.write_text(json.dumps(manifest,indent=2)+'\n')
    if a.evaluate:
        from cobaya.model import get_model
        from adapter import reference_point
        with get_model(info) as model:
            point = reference_point(model);r = model.logposterior(point)
            result = {'parameters':point,'logpost':r.logpost,
                'loglikes':dict(zip(model.likelihood,map(float,r.loglikes))),
                'derived':dict(zip(model.parameterization.derived_params(),map(float,r.derived)))}
        (out/'reference.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2));return
    covariance = ROOT/f'.work/unified-cosmology/external-probes/proposal-modern-{a.model}-dovekie.covmat'
    if not a.resume:
        names = covariance.open().readline().lstrip('#').split()
        sampled = [k for k,v in info['params'].items() if isinstance(v,dict) and 'prior' in v and k in names]
        index = [names.index(k) for k in sampled]
        delta = np.random.default_rng(a.seed+rank).multivariate_normal(np.zeros(len(index)),np.loadtxt(covariance)[np.ix_(index,index)])
        for k,d in zip(sampled,delta):info['params'][k]['ref'] += float(d)
    info['sampler'] = {'mcmc':{'covmat':str(covariance),'drag':True,'seed':a.seed,
        'learn_every':'30d','burn_in':'10d','Rminus1_stop':.005,'Rminus1_cl_stop':.05,
        'max_tries':'200d','output_every':'30s'}}
    if a.max_samples:info['sampler']['mcmc']['max_samples'] = a.max_samples
    info['output'] = str(out/'chain');info['timing'] = True
    from cobaya.run import run
    run(info,resume=a.resume)


if __name__=='__main__':main()
