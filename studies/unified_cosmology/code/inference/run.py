"""Joint Planck spectra + DESI DR2 + one SN compilation, using a shared CAMB model."""
import argparse
import json
from pathlib import Path
import sys
import os
import hashlib
import datetime

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parent/'external_probes'))
sys.path.insert(0, str(HERE))
from adapter import external_info
from likelihood import ReleasedDistances, ExpansionDiagnostics


def configuration(model='lcdm', cmb='full', evolution='none', sn=True, seed=272609, sample='dovekie'):
    info = external_info(cmb, 'cpl' if model in ['wcdm','cpl'] else 'lcdm')
    if model == 'wcdm':
        info['params']['wa'] = 0.
    info['likelihood']['expansion_diagnostics'] = {'external': ExpansionDiagnostics}
    if sn:
        from late_geometry import sample_path
        info['likelihood']['released_sn'] = {
            'external': ReleasedDistances,
            'data_file':str(sample_path(sample)),
            'smooth_sigma': {'none':0., 'linear':0., 'smooth01':.1, 'smooth03':.3}[evolution],
        }
        info['params']['epsilon'] = ({'prior': {'min':-.5, 'max':.5}, 'ref':0., 'proposal':.02}
                                     if evolution == 'linear' else 0.)
    proposal_name = 'proposal-lcdm' if model=='lcdm' else 'proposal-cpl'
    if sn and evolution=='none':
        proposal_name += '-dovekie'
    covariance = ROOT/'.work/unified-cosmology/external-probes'/(proposal_name+'.covmat')
    info['sampler'] = {'mcmc': {'covmat':str(covariance), 'drag':True,
                               'Rminus1_stop':.01, 'Rminus1_cl_stop':.1,
                               'learn_every':'40d', 'burn_in':'20d',
                               'max_tries':'100d', 'seed':seed,
                               'output_every':'30s'}}
    info['timing'] = True
    return info


def audit_camb_failures(path):
    """Record every failed CAMB calculation without changing its returned value."""
    import functools
    from cobaya.theories.camb.camb import CAMB, CambTransfers
    for cls in [CAMB, CambTransfers]:
        original = cls.calculate
        @functools.wraps(original)
        def wrapped(self, state, *args, _original=original, _label=cls.__name__, **kwargs):
            try:
                result = _original(self,state,*args,**kwargs)
            except Exception as error:
                with path.open('a') as file:
                    file.write(json.dumps({'component':_label,'params':state.get('params',kwargs),
                                           'exception':repr(error)},default=float)+'\n')
                raise
            if result is False:
                with path.open('a') as file:
                    file.write(json.dumps({'component':_label,'params':state.get('params',kwargs),
                                           'returned':False},default=float)+'\n')
            return result
        cls.calculate = wrapped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model', choices=['lcdm','wcdm','cpl'], default='lcdm')
    parser.add_argument('--cmb', choices=['full','lite'], default='full')
    parser.add_argument('--evolution', choices=['none','linear','smooth01','smooth03'], default='none')
    parser.add_argument('--no-sn', action='store_true')
    parser.add_argument('--sample', choices=['dovekie','pantheon','des3yr'], default='dovekie')
    parser.add_argument('--seed', type=int, default=272609)
    parser.add_argument('--evaluate', action='store_true')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--max-samples', type=int)
    args = parser.parse_args()
    info = configuration(args.model,args.cmb,args.evolution,not args.no_sn,args.seed,args.sample)
    name = f'{args.cmb}-{args.model}-{args.evolution}-'+('nonsn' if args.no_sn else args.sample)
    work = ROOT/'.work/unified-cosmology/inference'/name
    work.mkdir(parents=True,exist_ok=True)
    from mpi4py import MPI
    rank = MPI.COMM_WORLD.rank
    audit_camb_failures(work/f'camb-failures-{rank}.jsonl')
    # Starting points are dispersed in a correlated cosmological proposal, not
    # independent axis draws that destroy the CMB acoustic-scale relation.
    if not args.evaluate and not args.resume:
        import numpy as np
        proposal = Path(info['sampler']['mcmc']['covmat'])
        names = proposal.open().readline().lstrip('#').split()
        eligible = [key for key,value in info['params'].items()
                    if isinstance(value,dict) and 'prior' in value and key in names]
        indices = [names.index(key) for key in eligible]
        cov = np.loadtxt(proposal)[np.ix_(indices,indices)]
        delta = np.random.default_rng(args.seed+rank).multivariate_normal(np.zeros(len(indices)),cov)
        for key,change in zip(eligible,delta):
            info['params'][key]['ref'] += float(change)
    manifest = {'start_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
                'arguments':vars(args),'MPI_size':MPI.COMM_WORLD.size,'rank':rank,
                'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS'),
                'code_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in list(HERE.glob('*.py'))+[HERE/'design.json',HERE.parent/'external_probes/adapter.py']}}
    manifest_path = work/f'run-{rank}.json'
    if manifest_path.exists():
        if not args.resume and not args.evaluate:
            raise FileExistsError('Existing run manifest: use --resume explicitly.')
        with (work/f'run-history-{rank}.jsonl').open('a') as file:
            file.write(json.dumps(manifest)+'\n')
    else:
        manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    if args.evaluate:
        from cobaya.model import get_model
        from adapter import reference_point
        with get_model(info) as model:
            point = reference_point(model)
            posterior = model.logposterior(point)
            result = {'parameters':point, 'logpost':posterior.logpost,
                      'loglikes':dict(zip(model.likelihood,posterior.loglikes)),
                      'derived':dict(zip(model.parameterization.derived_params(),posterior.derived))}
        (work/'reference-evaluation.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result,indent=2))
    else:
        from cobaya.run import run
        info['output'] = str(work/'chain')
        if args.max_samples:
            info['sampler']['mcmc']['max_samples'] = args.max_samples
        run(info,resume=args.resume)


if __name__ == '__main__':
    main()
