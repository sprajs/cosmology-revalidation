"""Synthetic-only validation of preparatory native-accuracy contracts."""
import argparse
from contextlib import ExitStack
from dataclasses import FrozenInstanceError
import copy
import json
from pathlib import Path
from types import SimpleNamespace,MethodType
import tempfile
from unittest.mock import patch

import numpy as np

import native_accuracy_correction as contracts
from native_accuracy_measurement import target_descriptor

ROOT=contracts.ROOT
OUT=ROOT/'studies/unified_cosmology/results/inference/native-accuracy-contract-validation.json'


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=OUT)
    args=parser.parse_args()
    import camb
    import cobaya.model
    from cobaya.theories.camb.camb import CAMB
    def forbidden(*args,**kwargs):
        raise AssertionError('No model/background/spectrum calls authorized in preparatory validation.')
    source=contracts.sources()
    rng=np.random.default_rng(273781)
    refused=[]
    def refuse(name,call):
        try:call()
        except (AssertionError,ValueError,FrozenInstanceError):refused.append(name)
        else:raise AssertionError('Did not refuse '+name)
    with ExitStack() as guard:
        for owner,name in [(camb,'get_results'),(camb,'get_background'),(camb,'get_transfer_functions'),
                           (cobaya.model,'get_model'),(cobaya.model.Model,'__init__'),(CAMB,'__init__')]:
            guard.enter_context(patch.object(owner,name,forbidden))
        configuration_count=0
        for model in ('lcdm','cpl'):
            for evolution in ('none','linear','smooth01','smooth03'):
                for gpu in (False,True):
                    settings={'model':model,'evolution':evolution,'sample':'dovekie',
                              'calibration':'official_planck','fast_lensing':True,'gpu':gpu}
                    one,two=contracts.configuration_pair(settings)
                    target=target_descriptor(one,two,'synthetic-parent',{'data':'0'*64},
                                             {'synthetic':'1'},source)
                    assert target.numerical_accuracy==2 and target.scope.endswith('not_posterior_qualification')
                    initial=target.configuration
                    initial['theory']['camb']['extra_args']['AccuracyBoost']=99
                    assert target.configuration['theory']['camb']['extra_args']['AccuracyBoost']==2
                    refuse('frozen_descriptor_'+str(configuration_count),lambda:setattr(target,'numerical_accuracy',1))
                    other=target_descriptor(one,two,'synthetic-parent',{'data':'1'*64},{'synthetic':'1'},source)
                    assert other.target_identity!=target.target_identity
                    bad=copy.deepcopy(two);bad['params']['H0']['prior']['max']+=1
                    refuse('changed_prior_'+str(configuration_count),lambda:contracts.assert_only_boost_change(one,bad))
                    bad=copy.deepcopy(two);bad['theory']['camb']['extra_args']['mnu']=.12
                    refuse('changed_physics_'+str(configuration_count),lambda:contracts.assert_only_boost_change(one,bad))
                    configuration_count+=1
        max_conversion=0.;max_roundtrip=0.
        for n in (8,31,127,9002):
            for tcmb in (2.7,2.7255):
                # Fake state only; the actual installed provider conversion
                # method is invoked without constructing CAMB or a model.
                total=rng.normal(size=(n,4))*1e-12;total[0]=0
                potential=rng.normal(size=(n,1))*1e-9
                fake=SimpleNamespace(current_state={'Cl':{'total':total,'lens_potential':potential},
                                                     'derived_extra':{'TCMB':tcmb}},
                                     _needs_lensing_cross=False)
                fake._cmb_unit_factor=MethodType(CAMB._cmb_unit_factor,fake)
                dl=CAMB._get_Cl(fake,ell_factor=True,units='FIRASmuK2')
                expected=CAMB._get_Cl(fake,ell_factor=False,units='FIRASmuK2')
                archived={key:dl[key] for key in ('ell',*contracts.SPECTRA)}
                raw=contracts.raw_from_dl(archived);roundtrip=contracts.dl_from_raw(raw)
                for key in contracts.SPECTRA:
                    err=float(np.max(abs(raw[key]-expected[key]))/max(float(np.max(abs(expected[key]))),1e-300))
                    back=float(np.max(abs(roundtrip[key]-archived[key]))/max(float(np.max(abs(archived[key]))),1e-300))
                    max_conversion=max(max_conversion,err);max_roundtrip=max(max_roundtrip,back)
                assert raw['pp'][0]==archived['pp'][0]
        assert max_conversion<1e-13 and max_roundtrip<1e-13
        for mode in ('nonzero_ell0','nonfinite','wrong_shape','missing_component','noncontiguous_ell'):
            bad={key:np.copy(value) for key,value in archived.items()}
            if mode=='nonzero_ell0':bad['tt'][0]=1
            elif mode=='nonfinite':bad['ee'][2]=np.nan
            elif mode=='wrong_shape':bad['pp']=bad['pp'][:-1]
            elif mode=='missing_component':del bad['bb']
            else:bad['ell'][3]=99
            refuse(mode,lambda:contracts.raw_from_dl(bad))
        # Independent exact Gaussian target shift, retaining the same proposal.
        n=40000;x=rng.normal(size=(n,2));groups=np.repeat(np.arange(4),n//4)
        mu1=np.array([.1,0]);mu2=np.array([.25,-.15]);prior=-.2
        q=-.5*np.sum(x*x,axis=1)-np.log(2*np.pi)
        p1=q+x@mu1-.5*(mu1@mu1);p2=q+x@mu2-.5*(mu2@mu2)
        records=[];closure=0.
        for i in range(n):
            parent={'status':'finite','exact_logpost':float(p1[i]),'proposal_logpost':float(q[i]),
                    'log_weight':float(p1[i]-q[i]),'exact_loglikes':{'toy':float(p1[i]-prior)},
                    'derived':{'z':float(2*x[i,0]+x[i,1])}}
            derived={'z':float(3*x[i,0]+x[i,1])}
            row=contracts.accuracy2_density(parent,{'toy':float(p2[i]-prior)},[prior],float(p2[i]),derived)
            assert row['observational_qualification'] is False
            closure=max(closure,max(abs(v) for v in row['closure'].values()))
            records.append({'status':'finite','point':{'x':float(x[i,0]),'y':float(x[i,1])},
                            'derived':derived,'log_weight':row['accuracy2_logweight']})
        from exact_correction import summarize
        result=summarize(records,groups)
        assert result['status']=='passed_importance_weight_gates'
        expected={'x':(.25,1),'y':(-.15,1),'z':(.6,np.sqrt(10))}
        gaussian_errors={}
        for name,(mean,sd) in expected.items():
            value=result['posterior'][name]
            gaussian_errors[name]={'mean':abs(value['mean']-mean),'sd':abs(value['sd']-sd)}
            assert gaussian_errors[name]['mean']<.04 and gaussian_errors[name]['sd']<.04
        assert abs(result['posterior']['z']['mean']-.35)>.15, 'Accidentally used the old derived value.'
        # Reuse the identical gate implementation; a new target receives no
        # relaxed thresholds when its shifted weights lose support.
        chosen=np.concatenate([np.flatnonzero(groups==g)[:1000] for g in range(4)])
        subset=[copy.deepcopy(records[int(i)]) for i in chosen];sg=groups[chosen]
        for row,g in zip(subset,sg):row['log_weight']+=8*(g==3)
        failed=summarize(subset,sg)
        assert failed['status']=='insufficient_importance_overlap' and failed['failed_gates']
        prior_changed=[prior+.01]
        refuse('prior_change',lambda:contracts.accuracy2_density(parent,{'toy':float(p2[-1]-prior)},prior_changed,float(p2[-1]+.01),derived))
        corrupted=copy.deepcopy(parent);corrupted['log_weight']+=.02
        refuse('parent_weight_accounting',lambda:contracts.accuracy2_density(corrupted,{'toy':float(p2[-1]-prior)},[prior],float(p2[-1]),derived))
        refuse('missing_native2_derived',lambda:contracts.accuracy2_density(parent,{'toy':float(p2[-1]-prior)},[prior],float(p2[-1]),{}))
        refuse('nonfinite_native2',lambda:contracts.accuracy2_density(parent,{'toy':float('nan')},[prior],float(p2[-1]),derived))
        # Explicitly mocked evidence boundary: tests gate semantics, not the
        # real consumer's independent provenance validation already upstream.
        import native_precision_review_consumer as consumer
        import native_precision_refinement as refinement
        receipt={'reviewed_screen_status':'numerical_sensitivity_requires_followup',
                 'reviewed_diagnostic_flags':['total_centered_RMS'],
                 'parent_target_identity':'synthetic-parent','original_screen_status':'incomplete_or_failed_numerical_screen'}
        plan={'verified_original_receipt_identity':refinement.original.identity(receipt),
              'parent_target_identity':'synthetic-parent','identity':'synthetic-refinement',
              'points':[{'point':{'x':float(i)}} for i in range(32)]}
        report={'status':'no_large_variation_detected_on_fixed32','diagnostic_flags':[],
                'known_native_logposterior_invocations':32,'unknown_count_attempts':0,
                'cumulative_known_native_logposterior_invocations':33,'replay_gate_sha256':'synthetic-replay',
                'points':[{'point':{'x':float(i)}} for i in range(32)]}
        with tempfile.TemporaryDirectory(prefix='native-accuracy-contract-',dir=ROOT/'.work/unified-cosmology') as tmp:
            work=Path(tmp);(work/'refine-summary.json').write_text(json.dumps(report))
            def call_gate():
                with ExitStack() as stack:
                    stack.enter_context(patch.object(consumer,'verify',lambda *a:copy.deepcopy(receipt)))
                    stack.enter_context(patch.object(refinement,'verify_plan',lambda *a:copy.deepcopy(plan)))
                    stack.enter_context(patch.object(refinement,'require_replay',lambda *a:'synthetic-replay'))
                    stack.enter_context(patch.object(refinement,'run_stage',lambda *a:copy.deepcopy(report)))
                    stack.enter_context(patch.object(refinement.original,'sealed_read',lambda *a:dict(copy.deepcopy(report),payload_sha256='synthetic')))
                    return contracts.prerequisite_snapshot('synthetic-screen','synthetic-review',work)
            good=call_gate();assert good['posterior_qualified'] is False and good['physical_execution_authorized'] is False
            original_receipt=copy.deepcopy(receipt);original_report=copy.deepcopy(report)
            for key,value in [('status','numerical_sensitivity_requires_followup'),
                              ('unknown_count_attempts',1),('known_native_logposterior_invocations',31),
                              ('cumulative_known_native_logposterior_invocations',32),
                              ('replay_gate_sha256','wrong'),('diagnostic_flags',['total_centered_RMS'])]:
                report.clear();report.update(copy.deepcopy(original_report));report[key]=value
                refuse('refinement_'+key,call_gate)
            report.clear();report.update(copy.deepcopy(original_report))
            report['points'][31]['point']['x']=-1
            refuse('changed_refinement_point',call_gate)
            report.clear();report.update(copy.deepcopy(original_report))
            receipt['reviewed_diagnostic_flags']=['unresolved_native_failure']
            refuse('nongenuine_density_variation',call_gate)
            receipt.clear();receipt.update(original_receipt)
        assert source==contracts.sources()
    output={'status':'passed_preparatory_native_accuracy_contract_validation',
            'implementation_stage':'preparation_only','physical_executor_implemented':False,
            'observational_qualifier_implemented':False,'physical_execution_authorized':False,
            'model_background_spectrum_calls':0,'configuration_pairs':configuration_count,
            'spectrum_conversion':dict(contracts.converter_contract(),fake_installed_provider_states=8,
                                       maximum_relative_error=max_conversion,maximum_roundtrip_relative_error=max_roundtrip),
            'gaussian_density_cases':n,'maximum_density_accounting_error':closure,
            'known_gaussian_errors':gaussian_errors,'shared_gate_failure_flags':failed['failed_gates'],
            'refused_cases':refused,'mocked_evidence_boundary_explicit':True,
            'source_sha256':source,'seed':273781,
            'limits':'No observational posterior, persistent-worker equivalence, throughput or physical execution is qualified. The actual2000record identity consumer and processpool remain to be implemented and reviewed.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:output[k] for k in ('status','configuration_pairs','gaussian_density_cases','model_background_spectrum_calls')}))


if __name__=='__main__':main()
