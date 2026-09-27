"""Independent configuration and cumulative-weight stratification audit."""
import copy
import hashlib
import json
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch
import numpy as np
import yaml

ROOT=Path(__file__).resolve().parents[4]
HERE=ROOT/'studies/unified_cosmology/code/inference'
sys.path.insert(0,str(HERE))
import anchored_sampling as sampling
import anchored_correction as correction
import camb
import cobaya.model


def forbidden(*args,**kwargs):raise AssertionError('No model/background/spectrum allowed.')


with patch.object(camb,'get_results',forbidden),patch.object(camb,'get_background',forbidden),patch.object(camb,'get_transfer_functions',forbidden),patch.object(cobaya.model,'get_model',forbidden):
    configs=[]
    model_file=ROOT/'.work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz'
    for model in ['lcdm','cpl']:
        reference=sampling.adapter.configuration(model=model)
        original_proposal=sampling.adapter.configuration(model=model,surrogate=model_file)
        for accuracy in [1,2]:
            arguments=sampling.settings(model,model_file,accuracy)
            native=sampling.configuration(arguments,native=True)
            proposal=sampling.configuration(arguments)
            assert sampling.canonical(proposal)==sampling.canonical(original_proposal)
            control=copy.deepcopy(native)
            for key in ['AccuracyBoost','lAccuracyBoost','lSampleBoost']:
                assert control['theory']['camb']['extra_args'][key]==accuracy
                control['theory']['camb']['extra_args'][key]=1
            assert sampling.canonical(control)==sampling.canonical(reference)
            assert native['likelihood']['released_sn']['external'] is sampling.adapter.AnchoredReleasedSN
            assert native['params']==proposal['params']
            configs.append({'model':model,'accuracy':accuracy,'nonaccuracy_target_fields_identical':True})

    rng=np.random.default_rng(8274351)
    with tempfile.TemporaryDirectory(dir=ROOT/'.work') as temporary:
        folder=Path(temporary);paths=[];raws=[];expanded_lengths=[]
        for rank,n in enumerate([1200,1300,1400,1500]):
            # Shuffled parameter columns ensure selection does not use position
            # assumptions. Monotone row labels distinguish repeated holds.
            weights=rng.integers(1,7,n)
            raw=np.column_stack([weights,np.arange(n)/37.+rank,np.arange(n)/100.,60+np.arange(n)/200.,np.arange(n)/1000.-1])
            path=folder/f'chain.{rank+1}.txt'
            np.savetxt(path,raw,fmt='%.12g',header='weight minuslogpost label H0 w')
            paths.append(path);raws.append(np.loadtxt(path));expanded_lengths.append(int(weights.sum()))
        common=min(n-int(.3*n) for n in expanded_lengths)
        check={'equal_chain_length_for_diagnostics':common}
        selected=correction.select_points(paths,['w','H0'],check,2000)
        assert selected==correction.select_points(paths,['w','H0'],check,2000)
        assert np.array_equal(np.bincount(selected['groups']),[500]*4)
        maximum_error=0.
        for i,(point,loc,stored) in enumerate(zip(selected['points'],selected['locations'],selected['stored_logposts'])):
            rank=i//500;j=i%500;raw=raws[rank]
            assert selected['groups'][i]==rank and loc['chain']==paths[rank].name
            discard=expanded_lengths[rank]-common
            offset=loc['expanded_index']-discard
            # An independent cumulative-frequency lookup avoids np.repeat and
            # checks the original full-chain index, rather than mirroring it.
            row=int(np.searchsorted(np.cumsum(raw[:,0]),loc['expanded_index'],side='right'))
            assert row==loc['row']
            assert j*common/500-1 < offset < (j+1)*common/500
            expected={'w':float(raw[row,4]),'H0':float(raw[row,3])}
            assert point==expected and stored==-float(raw[row,1])
            maximum_error=max(maximum_error,abs(stored+raw[row,1]))

    # Verify compatibility with an existing completed sampler's file format;
    # this reads completion evidence only and does not qualify its posterior.
    oldfolder=ROOT/'.work/unified-cosmology/inference/modern-cpl-none-dovekie-official_planck-independence-seed273262-7316f91342ee'
    oldmanifest=json.loads((oldfolder/'run-0.json').read_text())
    minimal={'proposal':{'proposal_sha256':oldmanifest['arguments']['proposal_sha256']}}
    completed=correction.completion_evidence(oldfolder,minimal)
    assert completed['status']=='all_original_sampler_and_terminal_gates_passed'
    denied=[]
    with tempfile.TemporaryDirectory(dir=ROOT/'.work') as temporary:
        folder=Path(temporary)
        checkpoint=yaml.safe_load((oldfolder/'chain.checkpoint').read_text())
        last=copy.deepcopy(completed['final_convergence'])
        # Synthetic candidate ledgers/terminals exercise status/hash gates. No
        # such files are substituted into a scientific or sampled manifest.
        blobs={'chain.checkpoint':yaml.safe_dump(checkpoint),
               'chain_independence-convergence.jsonl':json.dumps(last)+'\n'}
        for rank in range(4):
            data=json.dumps({'synthetic':rank})+'\n'
            term={'converged':True,'terminal_weight':1,'terminal_state_in_collection':False,
                  'ledger_sha256':hashlib.sha256(data.encode()).hexdigest()}
            blobs[f'chain_independence-candidates-{rank}.jsonl']=data
            blobs[f'chain_independence-terminal-{rank}.json']=json.dumps(term)
        def reset():
            for name,data in blobs.items():(folder/name).write_text(data)
        reset();assert correction.completion_evidence(folder,minimal)['status']==completed['status']
        changes=[]
        for key,value in [('original_multivariate_previous',.005),('original_multivariate_current',.005),
                          ('original_bound_statistic',.05),('converged',False)]:
            bad=copy.deepcopy(last);bad[key]=value
            changes.append((key,'chain_independence-convergence.jsonl',json.dumps(bad)+'\n'))
        bad=copy.deepcopy(last);bad['rank_gate']['passed']=False
        changes.append(('rank_gate','chain_independence-convergence.jsonl',json.dumps(bad)+'\n'))
        bad=copy.deepcopy(checkpoint);next(iter(bad['sampler'].values()))['converged']=False
        changes.append(('checkpoint','chain.checkpoint',yaml.safe_dump(bad)))
        bad=json.loads(blobs['chain_independence-terminal-3.json']);bad['converged']=False
        changes.append(('terminal','chain_independence-terminal-3.json',json.dumps(bad)))
        changes.append(('candidate_ledger','chain_independence-candidates-2.jsonl','altered\n'))
        for label,name,data in changes:
            reset();(folder/name).write_text(data)
            try:correction.completion_evidence(folder,minimal)
            except (AssertionError,ValueError,KeyError):denied.append(label)
            else:raise AssertionError('Accepted invalid completion: '+label)
        assert len(denied)==8

paths=[Path(__file__)]+[HERE/name for name in ['anchored_sampling.py','anchored_correction.py','anchored_measurement.py','anchored-sampling-design.json','independence_sampler.py','independence_proposal.py','diagnostics.py']]
result={'status':'passed_independent_configuration_RLE_and_completion_checks','configurations':configs,
        'printed_RLE_selected_points':2000,'points_per_chain':[500]*4,'cumulative_weight_lookup_exact':True,
        'stored_logposterior_error':maximum_error,'native_or_background_or_model_calls':0,
        'completed_sampler_file_format_compatible':True,'completion_mutations_rejected':denied,
        'completed_sampler_input_sha256':completed['input_sha256'],
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'scope':'Independent configuration and printed-chain selection arithmetic only. Sampler completion/identity/consumer source review is separate; this record alone is not production approval.'}
out=ROOT/'studies/unified_cosmology/results/inference/anchored-scaffold-review-controls.json';out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'selected_points':2000}))
