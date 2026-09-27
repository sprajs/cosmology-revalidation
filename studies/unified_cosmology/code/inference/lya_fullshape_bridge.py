"""Qualified-parent replacement of the existing Lyalpha BAO pair.

This consumes a rounded-paper Gaussian full-shape approximation. It never adds
that approximation on top of the same forest's BAO measurements. All printed
rounding corners reuse the same background predictions and are separate targets,
not posterior draws or a probability model for the rounding uncertainty.
"""
import argparse
import copy
import importlib.metadata
import json
import os
from pathlib import Path
import time
from unittest.mock import patch

import numpy as np

import anchored_bridge as summary_kernel
import lya_fullshape as kernel
from exact_correction import verify_record
from expansion_history import native_thermal_parameters
from measurement_summary import summarize_run
from probe_omission import configuration_factory,identity,payload,verify_hashes
from target_identity import canonical,digest

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
DESIGN=HERE/'lya-fullshape-design.json'
GATES=HERE/'luminosity-sensitivity-design.json'
COMPONENT='bao.desi_dr2'
VALIDATION=ROOT/'studies/unified_cosmology/results/inference/lya-fullshape-bridge-validation.json'
KERNEL_VALIDATION=ROOT/'studies/unified_cosmology/results/inference/lya-fullshape-validation.json'
SOURCES=ROOT/'studies/unified_cosmology/results/external_probes/lya-fullshape-sources.json'


def relative(path):return str(Path(path).resolve().relative_to(ROOT))


def configurations(settings,frozen,design):
    assert settings['model']in design['allowed_models']
    assert settings['evolution']in design['allowed_evolution']
    assert settings['sample']in design['allowed_samples']
    assert settings['calibration']in design['allowed_calibration']
    factory=configuration_factory(settings)
    options={k:settings[k]for k in ['model','evolution','sample','calibration']}
    assert canonical(factory(**options,surrogate=Path(settings['surrogate'])))==frozen['configuration'],'Frozen proposal configuration changed.'
    native=canonical(factory(**options))
    assert COMPONENT in native['likelihood']
    # Reject non-default or second BAO factors, instead of silently interpreting
    # an arbitrary redshift compilation as this declared thirteen-row release.
    assert native['likelihood'][COMPONENT]=={}
    assert [key for key in native['likelihood']if key.startswith('bao.')]==[COMPONENT]
    target=copy.deepcopy(native)
    target['likelihood'][COMPONENT]={'semantic_factor':'lya_fullshape.GaussianBAOReplacement',
        'replaces_existing_z2p33_pair':True,'adds_independent_lya_data':False,
        'normalization':'The inherited Cobaya quadratic -chi2/2; no evidence claim.',
        'source':'Rounded Eq26 of arXiv:2607.27410v3; not the author unrounded likelihood.'}
    left,right=copy.deepcopy(native),copy.deepcopy(target)
    left['likelihood'].pop(COMPONENT);right['likelihood'].pop(COMPONENT)
    assert left==right,'A non-BAO likelihood, prior or physical setting changed.'
    return native,target


def background_prediction(point,extra,z):
    """One thermal-compatible background, no power spectra or transfer calls."""
    import camb
    cosmology={key:point[key]for key in ['H0','ombh2','omch2','ns','tau']}
    cosmology.update(As=1e-10*np.exp(point['logA']),w=point.get('w',-1.),wa=point.get('wa',0.))
    unique,inverse=np.unique(np.asarray(z),return_inverse=True)
    def forbidden(*a,**k):raise RuntimeError('Lyalpha bridge forbids CMB spectra and transfers.')
    with patch.object(camb,'get_results',forbidden),patch.object(camb,'get_transfer_functions',forbidden):
        parameters=camb.set_params(**cosmology,**extra)
        adjusted=native_thermal_parameters(parameters)
        background=camb.get_background(adjusted)
        DM=np.asarray(background.angular_diameter_distance(unique))*(1+unique)
        H=np.asarray(background.hubble_parameter(unique))
        rd=float(background.get_derived_params()['rdrag'])
    assert DM.shape==H.shape==unique.shape and np.isfinite(DM).all()and np.isfinite(H).all()
    assert np.all(DM>0)and np.all(H>0)and np.isfinite(rd)and rd>0
    return {'DM_Mpc':DM[inverse].tolist(),'H_km_s_Mpc':H[inverse].tolist(),'rdrag_Mpc':rd,
        'unique_redshifts':unique.tolist(),'ordered_redshifts':np.asarray(z).tolist(),
        'background_calls':1,'CMB_spectrum_calls':0,
        'thermal_adapter':{'input_WantTransfer':bool(parameters.WantTransfer),'background_WantTransfer':bool(adjusted.WantTransfer)}}


def replace_BAO(record,old_loglike,new_loglike,components,design):
    """Replace only stored native BAO after an independently computed closure."""
    exact,proposal=record['exact_loglikes'],record['proposal_loglikes']
    assert record['status']=='finite'and set(exact)==set(proposal)==set(components)
    tolerance=design['density_accounting_tolerance']
    assert np.isfinite(list(exact.values())+list(proposal.values())+[record[k]for k in ['log_weight','exact_logpost','proposal_logpost']]+[old_loglike,new_loglike]).all()
    assert abs(record['log_weight']-record['exact_logpost']+record['proposal_logpost'])<=tolerance
    prior=float(record['exact_logpost']-sum(exact.values()))
    assert abs(prior-record['proposal_logpost']+sum(proposal.values()))<=tolerance,'Native/proposal priors differ.'
    stored=float(exact[COMPONENT]);closure=float(old_loglike-stored)
    ratio=float(new_loglike-stored);lw=float(record['log_weight']+ratio)
    replaced=dict(exact,**{COMPONENT:float(new_loglike)})
    logpost=float(prior+sum(replaced.values()))
    accounting=float(logpost-record['proposal_logpost']-lw)
    assert abs(accounting)<=tolerance and np.isfinite([lw,logpost,accounting]).all()
    return {'status':'finite_BAO_replacement'if abs(closure)<=design['old_BAO_loglike_closure_tolerance']else'source_BAO_density_mismatch',
        'old_BAO_stored_native_loglike':stored,'old_BAO_reconstructed_loglike':float(old_loglike),'source_BAO_loglike_closure':closure,
        'target_BAO_loglike':float(new_loglike),'BAO_log_likelihood_ratio':ratio,'target_logweight':lw,
        'parent_native_proposal_logweight':float(record['log_weight']),'unchanged_logprior':prior,
        'target_logpost_in_inherited_normalization':logpost,'target_loglikes':replaced,'density_accounting_closure':accounting}


def summarize_variant(records,rows,groups,gates):
    # Reuse all existing overlap and posterior-summary equations. The field
    # adaptation is local and explicit; native SN chi-square stays unchanged.
    adapted=[]
    for record,row in zip(records,rows):
        if row['status']!='finite_BAO_replacement':adapted.append(row);continue
        adapted.append(dict(row,status='finite_anchored_replacement',
            source_SN_loglike_closure=row['source_BAO_loglike_closure'],
            SN_log_likelihood_ratio=row['BAO_log_likelihood_ratio'],target_sn_chi2=record['derived']['sn_chi2']))
    result=summary_kernel.summarize(records,adapted,groups,gates)
    result['status']={'qualified_conditional_anchored_bridge':'qualified_conditional_rounded_lya_replacement',
        'insufficient_anchored_overlap_or_stability':'insufficient_lya_overlap_or_stability',
        'failed_anchored_bridge_evaluation':'failed_lya_bridge_evaluation'}[result['status']]
    for old,new in [('source_SN_max_absolute_loglike_closure','source_BAO_max_absolute_loglike_closure'),
                    ('SN_log_likelihood_ratio_quantiles','BAO_log_likelihood_ratio_quantiles')]:
        if old in result:result[new]=result.pop(old)
    return result


def identify_target(settings,frozen,native,semantic,source_hashes,input_hashes):
    if settings.get('gpu'):
        from modern_gpu import identify as identify_native
    elif settings.get('fast_lensing'):
        from modern_fast import identify as identify_native
    else:
        from target_identity import identify as identify_native
    sample=Path(native['likelihood']['released_sn']['data_file'])
    reference=identify_native(native,sample)
    assert reference['assets']==frozen['assets']and reference['sample_sha256']==frozen['sample_sha256']
    assert all(reference['versions'][k]==v for k,v in frozen['versions'].items())
    target=copy.deepcopy(reference);target.pop('identity');target['configuration']=semantic
    target['configuration_role']='Non-executable semantic BAO replacement; inherited executable non-BAO target and priors are unchanged.'
    target['original_native_reference_identity']=reference['identity']
    target['source_sha256'].update(source_hashes);target['replacement_input_sha256']=input_hashes
    target['replacement_variants']={name:value.tolist()for name,value in kernel.variants().items()}
    target['scientific_scope']='Rounded published Lyalpha Gaussian; no raw correlation reanalysis, no author-tail reproduction, no extra BAO pair.'
    target['identity']=identity(target)
    return target


def read_record(path,cache_identity,index,point,native_hash):
    row=json.loads(Path(path).read_text());seal=row.pop('payload_sha256')
    assert identity(row)==seal,'Cached Lyalpha payload changed.'
    assert row['identity']==cache_identity and row['index']==index and row['point']==point and row['native_record_sha256']==native_hash
    row['payload_sha256']=seal;return row


def validation_guard():
    report=json.loads(VALIDATION.read_text())
    assert report['status']=='passed_synthetic_lya_fullshape_bridge_checks_no_observational_evaluation'
    required=[Path(__file__),HERE/'lya_fullshape_bridge_validate.py',Path(kernel.__file__),DESIGN]
    assert all(report['source_sha256'].get(relative(p))==digest(p)for p in required),'Stale Lyalpha consumer validation.'
    verify_hashes(report['source_sha256']);verify_hashes(report['input_sha256'])
    assert all(report['input_sha256'].get(relative(p))==digest(p)for p in [KERNEL_VALIDATION,SOURCES])
    assert report['CMB_spectrum_calls']==report['CAMB_background_calls']==report['observational_bridge_points']==0
    return {relative(p):digest(p)for p in [VALIDATION,KERNEL_VALIDATION,SOURCES]}


def load_released():
    model=kernel.GaussianBAOReplacement.from_release()
    assert len(model.z)==13 and len(model.order)==2 and np.sum(~model.lya)==11
    inputs={relative(p):digest(p)for p in [kernel.MEAN,kernel.COVARIANCE,SOURCES,KERNEL_VALIDATION]}
    acquisition=json.loads(SOURCES.read_text())
    assert acquisition['status']=='versioned_sources_acquired_printed_Gaussian_is_explicit_approximation'
    for key in ['source_sha256','input_sha256']:
        verify_hashes(acquisition[key]);inputs.update(acquisition[key])
    for resource in acquisition['resources'].values():
        path=ROOT/resource['path'];assert path.stat().st_size==resource['bytes']
        assert digest(path)==resource['sha256'],'Published Lyalpha source bytes changed.'
        inputs[resource['path']]=resource['sha256']
    validation=json.loads(KERNEL_VALIDATION.read_text())
    assert validation['status']=='passed'
    assert all(validation['source_sha256'].get(relative(p))==digest(p)for p in [Path(kernel.__file__),DESIGN])
    for key in ['source_sha256','input_sha256']:
        verify_hashes(validation[key]);inputs.update(validation[key])
    return model,inputs


def evaluation_checks(evaluation,variant_names,design):
    assert set(evaluation['variants'])==set(variant_names)
    assert np.isfinite([evaluation[k]for k in ['old_loglike','old_chi2','retained_chi2','old_lya_chi2','old_block_closure']]).all()
    assert abs(evaluation['old_block_closure'])<=design['block_quadratic_closure_tolerance']
    for item in evaluation['variants'].values():
        assert np.isfinite(list(item.values())).all()
        assert abs(item['new_block_closure'])<=design['block_quadratic_closure_tolerance']


def actual(folder,summary_path,cache):
    parent=summarize_run(folder,summary_path)
    assert parent['qualified_under_declared_numerical_gates']is True
    prerequisites=validation_guard()
    assert all(os.environ.get(key)=='1'for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','CLIPY_NOJAX'])
    folder,summary_path,cache=[Path(p).resolve()for p in [folder,summary_path,cache]]
    assert cache.is_relative_to(ROOT/'.work')
    verify_hashes(parent['input_sha256'])
    design=json.loads(DESIGN.read_text());gates=json.loads(GATES.read_text())['overlap_gates']
    summary=json.loads(summary_path.read_text());selection_path=ROOT/summary['selection_path']
    assert selection_path.parent.parent==folder
    selection=json.loads(selection_path.read_text());settings=selection['settings']
    frozen=json.loads((folder/'run-0.json').read_text())['target_identity'];assert frozen['identity']==parent['target_identity']
    native,semantic=configurations(settings,frozen,design)
    model,input_hashes=load_released()
    source_paths=[Path(__file__),Path(kernel.__file__),DESIGN,GATES,Path(summary_kernel.__file__),
        HERE/'probe_omission.py',HERE/'measurement_summary.py',HERE/'exact_correction.py',HERE/'target_identity.py',
        HERE/'expansion_history.py',HERE/'modern_run.py',HERE/'modern_fast.py',HERE/'luminosity_sensitivity.py']
    if settings.get('gpu'):source_paths.append(HERE/'modern_gpu.py')
    sources={relative(p):digest(p)for p in source_paths}
    target=identify_target(settings,frozen,native,semantic,sources,input_hashes)
    lineage={'qualified_parent_inputs':parent['input_sha256'],'source_sha256':sources,'replacement_input_sha256':input_hashes,
        'validation_prerequisites_sha256':prerequisites,'parent_proposal_identity':frozen['identity'],
        'native_target':target,'parent_settings':settings,'all_non_BAO_likelihoods_priors_physics_identical':True}
    cache_identity=identity(lineage);cache.mkdir(parents=True,exist_ok=True)
    line=cache/'lineage.json';payload(line,lineage);ledger=cache/'record-hashes.json'
    if ledger.exists():verify_hashes(json.loads(ledger.read_text()))
    output_hashes={relative(line):digest(line)};records=[];rows=[];start=time.monotonic()
    for i,point in enumerate(selection['points']):
        native_path=selection_path.parent/f'{i:05d}.json';record=json.loads(native_path.read_text());verify_record(record)
        assert record['status']=='finite'and record['point']==point
        records.append(record);path=cache/f'{i:05d}.json'
        if path.exists():row=read_record(path,cache_identity,i,point,digest(native_path))
        else:
            try:
                prediction=background_prediction(point,native['theory']['camb']['extra_args'],model.z)
                evaluation=model.evaluate(prediction['DM_Mpc'],prediction['H_km_s_Mpc'],prediction['rdrag_Mpc'])
                evaluation_checks(evaluation,model.new,design)
                variants={name:replace_BAO(record,evaluation['old_loglike'],item['new_loglike'],native['likelihood'],design)
                          for name,item in evaluation['variants'].items()}
                row={'status':'completed_background_and_variants','background':prediction,'evaluation':evaluation,'variants':variants}
            except Exception as error:row={'status':'failed_background_or_density','error':repr(error)}
            row.update(identity=cache_identity,index=i,point=point,native_record_sha256=digest(native_path));row['payload_sha256']=identity(row);payload(path,row)
        rows.append(row);output_hashes[relative(path)]=digest(path)
        if(i+1)%100==0:print(json.dumps({'lya_bridge_points':i+1,'seconds':time.monotonic()-start}),flush=True)
    payload(ledger,output_hashes)
    names=list(model.new)
    variants={name:summarize_variant(records,[r['variants'][name]if r['status']=='completed_background_and_variants'else r for r in rows],selection['groups'],gates)for name in names}
    for name,value in variants.items():
        parameters=kernel.variants()[name].tolist()
        value['printed_parameter_scenario']=parameters
        value['scenario_target_identity']=identity({'shared_native_target_identity':target['identity'],'variant':name,'printed_parameters':parameters})
    result={'status':'completed_lya_replacement_sensitivities','nominal':variants['nominal'],'variants':variants,
        'all_rounding_variants_qualified':all(r['qualified_under_declared_numerical_gates']for r in variants.values()),
        'rounding_scope':'Thirty-two deterministic printed rounding corners; no scenario probabilities and no averaged posterior.',
        'parent_settings':settings,'native_target':target,'source_sha256':sources,'replacement_input_sha256':input_hashes,
        'validation_prerequisites_sha256':prerequisites,'bridge_identity':cache_identity,'lineage_path':relative(line),'lineage_sha256':digest(line),
        'cache_manifest_path':relative(ledger),'cache_manifest_sha256':digest(ledger),'correction_summary_path':relative(summary_path),
        'correction_summary_sha256':digest(summary_path),'parent_qualified_under_declared_numerical_gates':True,
        'background_calls_in_recorded_cohort':sum(r.get('background',{}).get('background_calls',0)for r in rows),'CMB_spectrum_calls':0,
        'native_precision_screen':'Separate pending qualification, not established by replacement support.',
        'units':{'DM':'Mpc','H':'km/s/Mpc','rdrag':'Mpc','compressed_distances':'dimensionless'},'limits':design['limits']}
    for mapping in [sources,input_hashes,parent['input_sha256'],frozen['source_sha256'],target['source_sha256']]:verify_hashes(mapping)
    assert all(importlib.metadata.version(k)==v for k,v in target['versions'].items())
    assert validation_guard()==prerequisites
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['chain-folder','correction-summary','cache','output']:parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();result=actual(args.chain_folder,args.correction_summary,args.cache)
    args.output.parent.mkdir(parents=True,exist_ok=True);payload(args.output,result)
    print(json.dumps({'status':result['status'],'nominal_status':result['nominal']['status'],'CMB_spectrum_calls':0}))


if __name__=='__main__':main()
