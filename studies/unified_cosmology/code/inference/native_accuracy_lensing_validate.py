"""Synthetic density, typed-boundary and stored-spectrum tests; no physics calls."""
import argparse
import copy
from dataclasses import replace
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

import numpy as np
from scipy.special import logsumexp

import native_accuracy_lensing as c
import native_accuracy_qualify as qualifier


def reject(call):
    try: call()
    except (AssertionError, ValueError, KeyError, FileNotFoundError, RuntimeError): return True
    raise AssertionError('Invalid synthetic input accepted.')


def pure():
    rng=np.random.default_rng(274032); errors=[]
    inherited=json.loads(c.bridge.DESIGN.read_text()); old=inherited['old_factors']
    for _ in range(100):
        likes=dict(zip(old+['unchanged'],rng.normal(size=3)))
        prior=float(rng.normal()); post=sum(likes.values())+prior; q=float(rng.normal())
        new={'kernel_loglike':float(rng.normal()), 'normalized_gaussian_constant':float(rng.normal())}
        constants=dict(zip(old,rng.normal(size=2)))
        row={'native2_loglikes':likes,'native2_logpost':post,'proposal_logpost':q,
             'native2_logpriors':[prior],'log_weight':post-q}
        value=c.density(row,new,constants,old,likes,1e-9)
        expected=likes['unchanged']+prior+new['kernel_loglike']-q
        errors.append(abs(value['target_logweight']-expected))
        assert errors[-1]<1e-12
        expected_normalized=new['kernel_loglike']+new['normalized_gaussian_constant']-sum(likes[k]+constants[k] for k in old)
        assert abs(value['normalized_density_logratio']-expected_normalized)<1e-12
        changed=copy.deepcopy(row); changed['native2_logpost']+=1
        reject(lambda:c.density(changed,new,constants,old,likes,1e-9))
    return {'independent_full_density_cases':100,'maximum_absolute_error':max(errors)}


def archived_spectra():
    """Compare array-only compression with the unmodified released evaluator."""
    source=c.ROOT/'.work/unified-cosmology/external-probes/spectral-training/train'
    rows=[]; files={}
    for name in ('0000.npz','0001.npz'):
        path=source/name; files[c.rt.relative(path)]=c.rt.digest(path)
        with np.load(path,allow_pickle=False) as data:
            ell=data['ell']; values=data['spectra']; assert values.shape[0]==5
            factor=np.where(ell>=2,ell*(ell+1)/(2*np.pi),1.)
            raw={k:values[j]/(factor if k!='pp' else np.where(ell>=2,ell**2*(ell+1)**2/(2*np.pi),1.))
                 for j,k in enumerate(('tt','ee','bb','te','pp'))}
            raw['ell']=ell
            for k in ('tt','ee','bb','te','pp'):raw[k][:2]=0.
        c.capture.validate_spectra(raw)
        for variant in json.loads(c.DESIGN.read_text())['variants']:
            response=c.bridge.JointResponse(variant)
            value=response.evaluate(raw,native_check=True)
            assert value['native_loglike_absolute_difference']<1e-10
            assert value['native_bandpower_max_absolute_difference']<1e-16
            # Independent full covariance solve uses effective covariance C/h.
            n=response.native; kk=n.pp_to_kk(raw['pp'],ell)
            _,bp=n.generic_lnlike(response.data,ell,kk,ell,raw['tt'],raw['ee'],raw['te'],raw['bb'],return_theory=True)
            residual=response.data['data_binned_clkk']-bp
            effective_cov=np.linalg.inv(response.data['cinv'])
            independent=-.5*residual@np.linalg.solve(effective_cov,residual)
            constant=-.5*np.linalg.slogdet(2*np.pi*effective_cov)[1]
            assert abs(independent-value['kernel_loglike'])<1e-9
            assert abs(constant-value['normalized_gaussian_constant'])<1e-10
            rows.append({'archive':name,'variant':variant,
                         'native_loglike_error':value['native_loglike_absolute_difference'],
                         'independent_covariance_error':float(abs(independent-value['kernel_loglike'])),
                         'normalization_error':float(abs(constant-value['normalized_gaussian_constant']))})
    return rows,files


def fixture(folder):
    rng=np.random.default_rng(274033); n=2000
    x=rng.normal(size=n); a=rng.uniform(0,2,n)
    for i in (167,1229):x[i]=x[i-1];a[i]=a[i-1]
    inherited=json.loads(c.bridge.DESIGN.read_text()); old=inherited['old_factors']
    config={'theory':{'camb':{'extra_args':{k:2 for k in ('AccuracyBoost','lAccuracyBoost','lSampleBoost')}}},
        'params':{'A_fg':{'prior':{'min':0.,'max':2.}}},
        'likelihood':{old[0]:{'variant':'actplanck_baseline','lens_only':False,'apply_hartlap':True,'nsims_planck':400},
                      old[1]:{'clear_internal_priors':True},'unchanged':{}}}
    spectrum=folder/'native2.npz'; ell=np.arange(3102)
    np.savez(spectrum,ell=ell,**{k:np.ones(len(ell))*.01 for k in ('tt','ee','te','bb','pp')})
    data_hash=c.rt.digest(spectrum); inputs={c.rt.relative(spectrum):data_hash}
    rows=[]; locations=[]; constants={old[0]:1.25,old[1]:-.7}
    originals=folder/'native1'; originals.mkdir(); (originals/'spectra').mkdir()
    for i in range(n):
        point={'x':float(x[i]),'A_fg':float(a[i])}
        derived={k:float(s*x[i]+b) for k,s,b in [('q0',.1,-.5),('q05',.1,-.1),('q1',.1,.2),('j0',.2,1.)]}
        derived['sn_chi2']=float(x[i]**2)
        likes={old[0]:float(-.1*x[i]**2),old[1]:float(-.2*a[i]),'unchanged':float(-.4*x[i]**2+.2*a[i])}
        q=sum(likes.values())-np.log(2.)
        old_likes={k:v+1 for k,v in likes.items()}; old_post=q+3.
        parent={'index':i,'point':point,'status':'finite','target_identity':'c'*64,
                'exact_logpost':old_post,'exact_loglikes':old_likes,'log_weight':3.,'derived':dict(derived,q0=99.)}
        parentpath=originals/f'{i:05d}.json'; c.rt.write_new(parentpath,parent)
        sealed=c.rt.read(parentpath)
        sidecar={'index':i,'point':point,'status':'captured_from_parent_native_evaluation',
                 'parent_native_record_sha256':c.rt.digest(parentpath),'parent_payload_sha256':sealed['payload_sha256'],
                 'correction_identity':'c'*64,'source_sha256':c.capture.dependencies(), 'CLIPY_NOJAX':'1',
                 'additional_native_evaluations':0,'source_normalized_gaussian_constants':constants,
                 'spectral_units':'raw C_l: TT/EE/TE/BB FIRASmuK2, pp dimensionless potential',
                 'spectra_path':c.rt.relative(folder/'old-spectrum-does-not-exist.npz'),'spectra_sha256':'f'*64}
        c.rt.write_new(c.capture.sidecar_path(parentpath),sidecar)
        row={'index':i,'point':point,'group':i//500,'location':None,'status':'finite_native_accuracy2',
            'numerical_target_identity':'a'*64,'parent_native_path':c.rt.relative(parentpath),
            'parent_native_sha256':c.rt.digest(parentpath),'native1_logpost':old_post,
            'native1_loglikes':old_likes,'native1_logweight':3.,'native1_derived':parent['derived'],
            'native2_logpost':q,'native2_loglikes':likes,'native2_logpriors':[-float(np.log(2.))],
            'native2_derived':derived,'proposal_logpost':q,'log_weight':0.,
            'spectrum_path':c.rt.relative(spectrum),'spectrum_sha256':data_hash}
        location={'chain':str(i//500),'expanded_index':i%500,'row':i%500}
        if i in (167,1229): location=copy.deepcopy(locations[-1])
        locations.append(location); row['location']=location; rows.append(row)
    summary=folder/'summary.json'; summary.write_text('{}\n')
    inputs[c.rt.relative(summary)]=c.rt.digest(summary)
    target=qualifier.QualifiedNativeAccuracy2(2,'a'*64,config,
        {'model':'lcdm','evolution':'none','sample':'synthetic','calibration':'synthetic'},'b'*64,
        tuple(rows),tuple(r['point'] for r in rows),np.repeat(np.arange(4),500),tuple(locations),
        np.array([r['proposal_logpost'] for r in rows]),np.zeros(n),np.ones(n)/n,{},inputs,
        {'evidence':{'versions':{'numpy':np.__version__}}},folder/'parent-plan.json',summary)
    return target,summary,spectrum


def boundary():
    outcomes=[]
    with tempfile.TemporaryDirectory(dir=c.ROOT/'.work',prefix='native2-lensing-validation-') as tmp:
        folder=Path(tmp); target,summary,spectrum=fixture(folder)
        expected_inputs=dict(target.input_sha256)
        audit=folder/'audit.json';audit.write_text(json.dumps({'status':'released_likelihood_evaluation_reproduced_no_inference','input_sha256':{},'source_sha256':{}}))
        validation=folder/'validation.json'
        class Response:
            def __init__(self,variant):self.variant=variant
            def evaluate(self,spectra):
                assert float(spectra['tt'][2])==.01 # specifically the native2 fixture
                return {'kernel_loglike':0.,'normalized_gaussian_constant':2.,'normalized_gaussian_loglike':2.,'bandpowers':[.1]}
        with patch.object(c.bridge,'AUDIT',audit),patch.object(c,'VALIDATION',validation):
            validation.write_text(json.dumps({'status':'passed_native_accuracy2_lensing_validation','source_sha256':c.sources(),'physical_calls':0}))
            with patch.object(qualifier,'qualify_run',return_value=target),patch.object(c.bridge,'JointResponse',Response):
                result=c.actual(folder,summary,folder/'good-cache')
            for value in result['variants'].values():
                assert value['qualified_under_declared_numerical_gates'],value.get('failed_gates')
                assert 'A_fg' not in value['posterior']
                assert 'A_fg' in value['weight_diagnostics_including_auxiliary']['weighted_summaries_for_diagnostics']
                assert value['posterior']['q0']['mean']<0 # native1 derived q0=99 cannot leak
                assert value['weighted_covariance']['qualified']
            assert sum(x['repeated_slots'] for x in result['chronology'].values())==2
            c.rt.verify_hashes(result['input_sha256']);c.rt.verify_hashes(result['generated_sha256'])
            outcomes.append('fresh_typed_2000_slots_two_variants_native2_only_old_spectra_not_opened')
            # Direct independent weights for the nontrivial factor removal.
            expected=np.array([.1*r['point']['x']**2+.2*r['point']['A_fg'] for r in target.records])
            expected=np.exp(expected-logsumexp(expected))
            for value in result['variants'].values():
                actual_mean=value['posterior']['x']['mean']
                assert abs(actual_mean-expected@np.array([r['point']['x'] for r in target.records]))<1e-12
            outcomes.append('independent_normalized_target_weights_and_native2_derived')
            for name,altered in [('wrong_type',object()),('accuracy1',replace(target,numerical_accuracy=1)),
                ('native1_config',replace(target,native_configuration=copy.deepcopy(target.native_configuration))),
                ('wrong_slot',replace(target,records=(dict(target.records[0],index=1),)+target.records[1:])),
                ('wrong_native_identity',replace(target,records=(dict(target.records[0],numerical_target_identity='z'*64),)+target.records[1:]))]:
                if name=='native1_config':altered.native_configuration['theory']['camb']['extra_args']['AccuracyBoost']=1
                with patch.object(qualifier,'qualify_run',return_value=altered),patch.object(c.bridge,'JointResponse',side_effect=AssertionError('Must reject before response')):
                    assert reject(lambda:c.actual(folder,summary,folder/name))
                assert not (folder/name).exists();outcomes.append(name+'_rejected_before_response')
            for name,mutation in [('sidecar_wrong_point',lambda r:r.update(point={'x':-999})),
                ('sidecar_wrong_source',lambda r:r.update(source_sha256={})),
                ('sidecar_wrong_correction',lambda r:r.update(correction_identity='z'*64)),
                ('sidecar_wrong_parent',lambda r:r.update(parent_native_record_sha256='0'*64)),
                ('sidecar_changed_constants',lambda r:r['source_normalized_gaussian_constants'].update(SPT2023_lensing=3.))]:
                path=c.capture.sidecar_path(c.ROOT/target.records[0]['parent_native_path']);original=path.read_bytes()
                changed=c.rt.plain(c.rt.read(path));mutation(changed);path.unlink();c.rt.write_new(path,changed)
                with patch.object(qualifier,'qualify_run',return_value=target),patch.object(c.bridge,'JointResponse',side_effect=AssertionError('Must reject')):
                    assert reject(lambda:c.actual(folder,summary,folder/name))
                path.write_bytes(original);outcomes.append(name+'_rejected')
            bad=list(target.records);bad[0]=dict(bad[0],native2_logpost=100.)
            row={'status':'finite_bridge','density':{'target_logweight':0.,'fixed_gaussian_normalization_offset':0.}}
            failed_rows=[copy.deepcopy(row) for _ in range(2000)]
            failed_rows[0]={'status':'failed_bridge','error':'synthetic evaluation failure'}
            gates=json.loads(c.bridge.GATES.read_text())['overlap_gates']
            failed=c.summarize(target.records,failed_rows,target.groups,gates)
            assert all(failed[k] is None for k in ('posterior','weighted_covariance','conditional_sign_fractions'))
            outcomes.append('failed_evaluation_withholds_all_physical_summaries')
            for i in range(2000):
                failed_rows[i]=copy.deepcopy(row);failed_rows[i]['density']['target_logweight']=100. if i==0 else 0.
            failed=c.summarize(target.records,failed_rows,target.groups,gates)
            assert not failed['qualified_under_declared_numerical_gates']
            assert all(failed[k] is None for k in ('posterior','weighted_covariance','conditional_sign_fractions'))
            outcomes.append('failed_overlap_withholds_covariance_signs_and_posterior')
            assert target.input_sha256==expected_inputs
            # All output annotations are native2; old sidecar spectrum was absent.
            assert not (folder/'old-spectrum-does-not-exist.npz').exists()
            # End-of-run checks independently reject changed provenance after
            # successful calculation, without repeating a physical calculation.
            source=c.sources();vsha=c.rt.digest(validation)
            generated={c.rt.relative(summary):c.rt.digest(summary)}
            versions={'numpy':np.__version__}
            c.completion_guard(target.input_sha256,generated,source,versions,vsha)
            assert reject(lambda:c.completion_guard(target.input_sha256,generated,source,{'numpy':'wrong'},vsha))
            assert reject(lambda:c.completion_guard(target.input_sha256,generated,{},versions,vsha))
            assert reject(lambda:c.completion_guard(target.input_sha256,generated,source,versions,'0'*64))
            corrupted=dict(generated);corrupted[c.rt.relative(summary)]='0'*64
            assert reject(lambda:c.completion_guard(target.input_sha256,corrupted,source,versions,vsha))
            outcomes.extend(['end_version_change_rejected','end_source_change_rejected',
                             'end_validation_change_rejected','end_output_change_rejected'])
    return outcomes


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=c.VALIDATION);args=parser.parse_args()
    with c.rt.guards_without_physics():
        algebra=pure();cases=boundary();spectra,files=archived_spectra()
    audit=json.loads(c.bridge.AUDIT.read_text())
    inputs=c.rt.merge(files,audit['input_sha256'],audit['source_sha256'])
    c.rt.verify_hashes(inputs)
    result={'status':'passed_native_accuracy2_lensing_validation','physical_calls':0,'background_calls':0,
        'model_calls':0,'new_observational_posterior_evaluations':0,'synthetic_density':algebra,
        'typed_boundary_cases':cases,'stored_spectrum_closure':spectra,
        'source_sha256':c.sources(),'input_sha256':inputs,
        'scope':'Actual released array kernels on existing archived spectra plus explicit synthetic typed-parent fixtures; no cosmological posterior or physical evaluation.'}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'status':result['status'],'boundary_cases':len(cases),'sha256':c.rt.digest(args.output)}))


if __name__=='__main__':main()
