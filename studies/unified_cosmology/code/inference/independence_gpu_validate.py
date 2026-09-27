"""Configuration-only GPU independence validation: no model/density/physics calls."""
import argparse
import ast
from contextlib import ExitStack
import copy
import json
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace

import numpy as np
from independence_gpu_run import (HERE,ROOT,DESIGN,VALIDATION,configuration_from_parent,
    dependencies,load_gpu_proposal,sampler_options)
from independence_proposal import FrozenMixture,digest
from independence_run import prior_start,runtime
from target_identity import canonical


def reject(function,*args,**kwargs):
    try:function(*args,**kwargs)
    except AssertionError:return
    raise AssertionError('Invalid configuration unexpectedly accepted.')


def cpu_sampler_options(mixture,folder,training,seed):
    """Evaluate only the original literal options dictionary, without its driver."""
    tree=ast.parse((HERE/'independence_run.py').read_text())
    function=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='sample')
    nodes=[n for n in ast.walk(function) if isinstance(n,ast.Assign)
           and any(isinstance(t,ast.Name) and t.id=='options' for t in n.targets)]
    assert len(nodes)==1 and isinstance(nodes[0].value,ast.Dict)
    expression=ast.Expression(body=nodes[0].value);ast.fix_missing_locations(expression)
    return eval(compile(expression,'<frozen-CPU-options-only>','eval'),
                {'mixture':mixture,'proposal_folder':folder,'training':training,'seed':seed})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('parents',type=Path,nargs=2)
    p.add_argument('--proposals',type=Path,nargs=2)
    p.add_argument('--output',type=Path,default=VALIDATION)
    a=p.parse_args();assert not a.output.exists(),'Preserve earlier validation attempts.'
    bound=dependencies();environment=runtime();tests=[];identities=[]
    def forbidden(*args,**kwargs):raise AssertionError('Forbidden native/model/density call in configuration validator.')
    with ExitStack() as guard:
        for name in ['camb.get_results','camb.get_background','camb.get_transfer_functions',
                     'cobaya.model.get_model','cobaya.model.Model.logposterior']:
            guard.enter_context(patch(name,side_effect=forbidden))
        for folder in a.parents:
            folder=folder.resolve();parent=json.loads((folder/'run-0.json').read_text())
            source_before=copy.deepcopy(parent)
            info,identity,settings=configuration_from_parent(parent)
            assert parent==source_before
            assert canonical(info)==parent['target_identity']['configuration']
            from modern_fast import configuration as cpu_config
            from modern_gpu import GPUSpectralSurrogate
            from spectral_surrogate import SpectralSurrogate
            cpu=cpu_config(settings['model'],settings['evolution'],settings['sample'],settings['calibration'],Path(settings['surrogate']))
            assert info['theory']['spectral_surrogate']['external'] is GPUSpectralSurrogate
            assert cpu['theory']['spectral_surrogate']['external'] is SpectralSurrogate
            adjusted=copy.deepcopy(info);adjusted['theory']['spectral_surrogate']['external']=SpectralSurrogate
            assert canonical(adjusted)==canonical(cpu)
            names=[n for n,v in info['params'].items() if isinstance(v,dict) and 'prior'in v]
            assert len(names)==(14 if settings['evolution']=='linear' else 13)
            mean=[]
            for name in names:
                parameter=info['params'][name];ref=parameter['ref']
                mean.append(float(ref if isinstance(ref,(int,float)) else ref.get('loc',np.mean([parameter['prior'].get('min',0),parameter['prior'].get('max',1)]))))
            # This synthetic tiny covariance is an initialization-API control,
            # not the frozen proposal used by production.
            q=FrozenMixture(mean,np.eye(len(names))*1e-10,names)
            starts=[prior_start(q,info['params'],273261,r) for r in range(4)]
            assert len({tuple(x['point'].values()) for x in starts})==4
            assert starts==[prior_start(q,info['params'],273261,r) for r in range(4)]
            assert canonical(info)==parent['target_identity']['configuration']
            for seed in [273260,273270]:
                new=sampler_options(q,folder/'proposal.npz','test',seed)
                old=cpu_sampler_options(q,folder,{'proposal_sha256':'test'},seed)
                assert new.keys()==old.keys()
                for key in old:
                    if isinstance(old[key],np.ndarray):assert np.array_equal(old[key],new[key])
                    else:assert old[key]==new[key],key
            for change in [dict(gpu=False),dict(fast_lensing=False),dict(evolution='none'),dict(model='lcdm'),dict(calibration='paper_literal'),dict(sample='pantheon')]:
                invalid=copy.deepcopy(parent);invalid['arguments'].update(change);reject(configuration_from_parent,invalid)
            invalid=copy.deepcopy(parent);invalid['target_identity']['identity']='wrong';reject(configuration_from_parent,invalid)
            tests.append({'evolution':settings['evolution'],'GPU_source_identity_exact':True,
                          'CPU_configuration_differences_only_GPU_external':True,'all_scalar_priors_unchanged':True,
                          'four_distinct_deterministic_starts':True,'two_CPU_option_dictionary_comparisons':True,
                          'seven_invalid_parent_controls_rejected':True})
            identities.append({'parent':str(folder.relative_to(ROOT)),'run0_sha256':digest(folder/'run-0.json'),'target_identity':identity['identity']})
        assert {r['evolution'] for r in tests}=={'linear','smooth01'}
        # Reject wrong rank count before any model/identity work.
        from independence_gpu_run import sample
        with patch('mpi4py.MPI.COMM_WORLD',SimpleNamespace(rank=0,size=3)):
            reject(sample,Path('unused'),Path('unused'))
        # Actual frozen-dependency checksum gate, without mutating active bytes.
        victim=HERE/'independence_sampler.py'
        with patch('independence_gpu_run.digest',side_effect=lambda path:'wrong' if Path(path)==victim else digest(path)):
            reject(dependencies)
        proposal_records=[]
        if a.proposals:
            for folder in a.proposals:
                folder=folder.resolve();info,identity,settings,parent,q,sidecar=load_gpu_proposal(folder)
                reject(FrozenMixture.read,folder/'proposal.npz','wrong')
                defaults=json.loads(DESIGN.read_text())['default_seeds'][settings['evolution']]
                starts=[prior_start(q,info['params'],defaults['initialization'],rank) for rank in range(4)]
                assert len({tuple(s['point'].values()) for s in starts})==4
                proposal_records.append({'folder':str(folder.relative_to(ROOT)),
                    'proposal_sha256':sidecar['proposal_sha256'],'gpu_proposal_sha256':digest(folder/'gpu-proposal.json'),
                    'target_identity':identity['identity'],'dimension':q.d,'prior_initialization_attempts':[s['attempt'] for s in starts]})
    assert dependencies()==bound
    result={'status':'passed_configuration_only_no_likelihood_calls','driver_dependency_sha256':bound,
            'validator_sha256':digest(__file__),'environment':environment,'target_checks':tests,
            'nonfour_MPI_and_changed_helper_checksum_rejected':True,
            'parent_manifests':identities,'proposals':proposal_records,
            'native_spectra_calls':0,'background_calls':0,'likelihood_calls':0,'GPU_spectrum_matrix_calls':0,
            'not_tested':'GPU target density/throughput and native initialization at these targets; no production or pilot launched. Existing GPU dispatch and frozen sampler mathematical validations are inherited, not rerun.'}
    a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':result['status'],'targets':len(tests),'proposals':len(proposal_records)},indent=2))


if __name__=='__main__':main()
