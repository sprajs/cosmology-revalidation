"""Configuration, proposal and sealed-consumer tests; no cosmological evaluation."""
import argparse
import copy
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
from unittest.mock import patch
import numpy as np

import anchored_sampling as s
import anchored_correction as c
import anchored_measurement as m
from exact_correction import record_digest
from diagnostics import check_mpi


def rejected(call):
    try:call()
    except (AssertionError,ValueError,KeyError):return
    raise AssertionError('Invalid fixture accepted.')


def validate(full_identity=False):
    import camb
    import cobaya.model
    def forbidden(*args,**kwargs):raise AssertionError('No cosmological calls in scaffold validation.')
    model_file=s.ROOT/'.work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz'
    bound=s.dependencies();cases=[];native_ids=[];starts=[]
    with patch.object(camb,'get_results',forbidden),patch.object(camb,'get_background',forbidden),patch.object(camb,'get_transfer_functions',forbidden),patch.object(cobaya.model,'get_model',forbidden):
        for model in ['lcdm','cpl']:
            nominal=None
            for accuracy in [1,2]:
                settings=s.settings(model,model_file,accuracy)
                proposal=s.configuration(settings);native=s.configuration(settings,native=True)
                standard=s.canonical(s.adapter.configuration(model=model))
                actual=s.canonical(native)
                aextra=actual['theory']['camb']['extra_args'];bextra=standard['theory']['camb']['extra_args']
                for key in json.loads(s.DESIGN.read_text())['accuracy_keys']:
                    assert aextra.pop(key)==accuracy
                    assert bextra.pop(key,1)==1
                assert actual==standard
                assert native['theory']['camb']['extra_args']['lens_potential_accuracy']==4
                if nominal is None:nominal=s.canonical(proposal)
                else:assert nominal==s.canonical(proposal),'Native accuracy must not silently alter the sampled proposal.'
                if full_identity:
                    identity=s.identify(settings)
                    assert identity['native_configuration']==s.canonical(native)
                    assert identity['configuration']==s.canonical(proposal)
                    assert identity['native_accuracy']==accuracy
                    native_ids.append(identity['identity'])
                cases.append({'model':model,'native_accuracy':accuracy,'only_three_native_accuracy_changes':True,'proposal_unchanged':True})
            folder=s.ROOT/f'.work/unified-cosmology/inference/independence-proposal-{model}-v1'
            mixture,record=s.load_proposal(folder,proposal)
            seeds=json.loads(s.DESIGN.read_text())['default_seeds'][model]
            for rank in range(4):
                initial=s.prior_start(mixture,proposal['params'],seeds['initialization'],rank)
                starts.append({'model':model,'rank':rank,'attempt':initial['attempt'],'point':initial['point']})
                assert initial==s.prior_start(mixture,proposal['params'],seeds['initialization'],rank)
            options=s.sampler_options(mixture,record['proposal_file'],record['proposal_sha256'],seeds['sampling'])
            assert [options[k] for k in ['Rminus1_stop','Rminus1_cl_stop','Rminus1_cl_level']]==[.005,.05,.95]
            assert options['learn_proposal'] is False and options['temperature']==1
            assert options['blocking']==[[1,mixture.names]]
        if full_identity:assert len(set(native_ids))==4
        for invalid in [dict(s.settings('lcdm',model_file),gpu=True),dict(s.settings('lcdm',model_file),evolution='linear'),dict(s.settings('lcdm',model_file),sample='dovekie')]:
            rejected(lambda invalid=invalid:s.configuration(invalid))
        rejected(lambda:s.settings('lcdm',model_file,3))

        # Schema dispatch is separate from the proposal helper's full training
        # and ancestry validation, independently tested by its own validator.
        import anchored_proposal
        with tempfile.TemporaryDirectory(dir=s.ROOT/'.work') as td:
            dispatch=Path(td);(dispatch/'proposal.json').write_text(json.dumps({'schema':'anchored-bridge-weighted-proposal-v1'}))
            marker=(object(),{'proposal_kind':'anchored_bridge_weighted_v1'})
            with patch.object(anchored_proposal,'load',return_value=marker) as loader:
                assert s.load_proposal(dispatch,proposal) is marker
                loader.assert_called_once_with(dispatch.resolve(),proposal)
            (dispatch/'proposal.json').write_text(json.dumps({'schema':'unknown-proposal'}))
            rejected(lambda:s.load_proposal(dispatch,proposal))

        # Explicit synthetic fixture: target reidentification alone is mocked.
        # Chain diagnostics, RLE selection, native records, all raw weights,
        # recomputation, file seals and consumer gates are real.
        rng=np.random.default_rng(8274301)
        with tempfile.TemporaryDirectory(dir=s.ROOT/'.work') as td:
            folder=Path(td);target={'identity':'synthetic-anchored-only','configuration':{'params':{'H0':{'prior':{'min':50.,'max':90.}}}},
                                   'native_configuration':{'likelihood':{'released_sn':{},'synthetic_CMB_BAO':{}}},
                                   'native_accuracy':2}
            arguments={'synthetic_fixture':True,'native_accuracy':2,'model':'lcdm','evolution':'none','sample':s.adapter.SAMPLE,'calibration':'official_planck'}
            for rank in range(4):
                x=rng.normal(67.,1.,4000)
                weight=rng.integers(1,4,len(x));logpost=-.5*(x-67.)**2-2
                columns=np.column_stack([weight,-logpost,x,-.5+.03*(x-67.),-.1+.03*(x-67.),.15+.03*(x-67.),1.+.03*(x-67.),1.+.01*(x-67.)**2])
                np.savetxt(folder/f'chain.{rank+1}.txt',columns,header='weight minuslogpost H0 q0 q05 q1 j0 sn_chi2')
                manifest={'schema':'dedicated-anchored-chain-v1','rank':rank,'MPI_size':4,'arguments':arguments,'target_identity':target,'sampler_source_sha256':{},'proposal':{'proposal_sha256':'synthetic-proposal'}}
                s.dump_new(folder/f'run-{rank}.json',manifest)
            checkpoint={'sampler':{'independence_sampler.IndependenceMCMC':{'converged':True,'mpi_size':4,'burn_in':0,'Rminus1_last':.001}}}
            s.dump_new(folder/'chain.checkpoint',checkpoint)
            convergence={'original_multivariate_previous':.002,'original_multivariate_current':.001,
                         'original_bound_statistic':.03,'converged':True,'original_convergence_passed':True,
                         'proposal_sha256':'synthetic-proposal','rank_gate':{'passed':True,'failed':{},'discard_fraction_after_sampler_burnin':.3}}
            (folder/'chain_independence-convergence.jsonl').write_text(json.dumps(convergence)+'\n')
            for rank in range(4):
                ledger=folder/f'chain_independence-candidates-{rank}.jsonl';ledger.write_text('{"synthetic_only": true}\n')
                s.dump_new(folder/f'chain_independence-terminal-{rank}.json',{'converged':True,'terminal_weight':1,
                           'terminal_state_in_collection':False,'ledger_sha256':s.digest(ledger)})
            stop=c.completion_evidence(folder,manifest)
            raw=(folder/'chain.checkpoint').read_bytes();bad=copy.deepcopy(checkpoint);bad['sampler']['independence_sampler.IndependenceMCMC']['converged']=False
            (folder/'chain.checkpoint').write_text(json.dumps(bad));rejected(lambda:c.completion_evidence(folder,manifest));(folder/'chain.checkpoint').write_bytes(raw)
            invalid_stops=0
            convergence_path=folder/'chain_independence-convergence.jsonl';original_convergence=convergence_path.read_bytes()
            for key,value in [('original_multivariate_previous',.005),('original_multivariate_current',.006),
                              ('original_bound_statistic',.05),('original_convergence_passed',False),('converged',False),
                              ('proposal_sha256','wrong')]:
                changed=copy.deepcopy(convergence);changed[key]=value
                convergence_path.write_text(json.dumps(changed)+'\n')
                rejected(lambda:c.completion_evidence(folder,manifest));invalid_stops+=1
            convergence_path.write_bytes(original_convergence)
            terminal_path=folder/'chain_independence-terminal-3.json';original_terminal=terminal_path.read_bytes()
            changed=json.loads(original_terminal);changed['converged']=False;terminal_path.write_text(json.dumps(changed))
            rejected(lambda:c.completion_evidence(folder,manifest));invalid_stops+=1;terminal_path.write_bytes(original_terminal)
            ledger_path=folder/'chain_independence-candidates-3.jsonl';original_ledger=ledger_path.read_bytes();ledger_path.write_bytes(original_ledger+b' ')
            rejected(lambda:c.completion_evidence(folder,manifest));invalid_stops+=1;ledger_path.write_bytes(original_ledger)
            check=check_mpi(folder,.3);assert check['status']=='passed',check['failed_gates']
            diagnostics=folder/'diagnostics.json';s.dump_new(diagnostics,check)
            with patch.object(s,'verify_manifest',return_value={}):
                observed_manifest,actual_check,actual_stop=c.qualified_chain_inputs(folder,diagnostics)
                assert actual_stop==stop
                assert actual_check==check
            selected=c.select_points(sorted(folder.glob('chain.*.txt')),['H0'],check,2000)
            correction_bound=c.dependencies();cid=s.identity({'target_identity':target['identity'],'correction_dependencies':correction_bound})
            selected.update(schema='dedicated-anchored-native-selection-v1',settings=arguments,correction_identity=cid,
                correction_dependency_sha256=correction_bound,proposal_target_identity=target['identity'],native_accuracy=2,
                sampler_completion=stop,chain_sha256=check['input_sha256'],diagnostics_path=s.relative(diagnostics),diagnostics_sha256=s.digest(diagnostics),
                rank_manifest_sha256={s.relative(folder/f'run-{i}.json'):s.digest(folder/f'run-{i}.json') for i in range(4)})
            output=folder/'native';selection_path=output/'selection.json';s.dump_new(selection_path,selected)
            class FakeModel:
                likelihood={'released_sn':None,'synthetic_CMB_BAO':None}
                parameterization=SimpleNamespace(derived_params=lambda:['q0','q05','q1','j0','sn_chi2'])
                def logposterior(self,point):
                    t=point['H0']-67.;ll=-.5*t*t-2
                    return SimpleNamespace(logpost=ll,loglikes=[ll+1,-1.],derived=[-.5+.03*t,-.1+.03*t,.15+.03*t,1.+.03*t,1.+.01*t*t])
            records=[]
            with patch.object(c,'_bound',correction_bound),patch.object(c,'_initialization',{}),patch.object(c,'_exact',FakeModel()),patch.object(c,'_proposal',FakeModel()):
                for i,point in enumerate(selected['points']):records.append(c.evaluate((i,point,str(output),cid,selected['stored_logposts'][i])))
                assert c.evaluate((0,selected['points'][0],str(output),cid,selected['stored_logposts'][0]))==records[0]
                corrupt=output/'00000.json';before=corrupt.read_bytes();row=json.loads(before);row['log_weight']+=.1;corrupt.write_text(json.dumps(row))
                rejected(lambda:c.evaluate((0,selected['points'][0],str(output),cid,selected['stored_logposts'][0])))
                corrupt.write_bytes(before)
            summary=c.summarize(records,np.array(selected['groups']));assert summary['status']=='passed_importance_weight_gates',summary.get('failed_gates')
            summary.update(schema='dedicated-anchored-native-summary-v1',selection_path=s.relative(selection_path),selection_sha256=s.digest(selection_path),
                code_sha256=s.digest(c.__file__),correction_dependency_sha256=correction_bound,correction_identity=cid,
                proposal_target_identity=target['identity'],native_accuracy=2,diagnostics_path=selected['diagnostics_path'],diagnostics_sha256=selected['diagnostics_sha256'],
                initialization_record_sha256={},native_record_sha256={s.relative(output/f'{i:05d}.json'):s.digest(output/f'{i:05d}.json') for i in range(2000)})
            summary_path=folder/'summary.json';s.dump_new(summary_path,summary)
            with patch.object(s,'verify_manifest',return_value={}):
                measurement=m.summarize_run(folder,summary_path)
                assert measurement['qualified_under_declared_numerical_gates'] and measurement['settings']['native_accuracy']==2
                assert abs(measurement['native_correction']['raw_weight_ess']-2000)<1e-9
                expected=np.var([p['H0'] for p in selected['points']])
                got=measurement['weighted_covariance']['matrix'][0][0];assert abs(got-expected)<1e-12
                raw=summary_path.read_bytes();changed=json.loads(raw);changed['posterior']['H0']['mean']+=.1;summary_path.write_text(json.dumps(changed))
                rejected(lambda:m.summarize_run(folder,summary_path));summary_path.write_bytes(raw)
                raw=summary_path.read_bytes();changed=json.loads(raw);changed['status']='insufficient_importance_overlap';summary_path.write_text(json.dumps(changed))
                rejected(lambda:m.summarize_run(folder,summary_path));summary_path.write_bytes(raw)
                raw=(output/'00000.json').read_bytes();changed=json.loads(raw);changed['log_weight']+=.01;(output/'00000.json').write_text(json.dumps(changed))
                rejected(lambda:m.summarize_run(folder,summary_path));(output/'00000.json').write_bytes(raw)
            from measurement_summary import summarize_run as old_consumer
            rejected(lambda:old_consumer(folder,summary_path))
            # Qualified native summaries remain conditional on weight support.
            bad=copy.deepcopy(records);bad[0]['log_weight']=1000
            failed=c.summarize(bad,np.array(selected['groups']))
            assert failed['status']=='insufficient_importance_overlap' and 'raw_weight_ESS'in failed['failed_gates']
    assert s.dependencies()==bound
    paths=[Path(__file__),Path(s.__file__),Path(c.__file__),Path(m.__file__),s.DESIGN,
           s.HERE/'anchored_adapter.py',s.HERE/'exact_correction.py',s.HERE/'measurement_summary.py',
           s.HERE/'luminosity_sensitivity.py',s.HERE/'luminosity-sensitivity-design.json']
    return {'status':'passed_anchored_scaffold_synthetic_checks_no_inference',
            'configuration_cases':cases,'full_live_target_identities_checked':native_ids,
            'initializations':starts,'synthetic_native_records':2000,
            'synthetic_covariance_error':abs(got-expected),'synthetic_cache_tamper_rejected':True,'invalid_sampler_completion_controls_rejected':invalid_stops+1,
            'synthetic_consumer_recompute_and_failure_checks':True,'generic_old_consumer_rejected':True,'genuine_anchored_proposal_schema_dispatch_checked':True,
            'native_CMB_or_background_calls':0,'sampling_dependency_sha256':bound,
            'source_sha256':{s.relative(p):s.digest(p) for p in paths},
            'limitations':'Synthetic pipeline/identity validation only. No cosmological density was evaluated, and no sampler efficiency, posterior or unseen-mode coverage is established.'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--full-identity',action='store_true');p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=validate(a.full_identity);s.dump_new(a.output,result)
    print(json.dumps({'status':result['status'],'full_identities':len(result['full_live_target_identities_checked'])}))

if __name__=='__main__':main()
