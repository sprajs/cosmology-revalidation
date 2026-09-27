"""Qualified-parent HF replacement with calibrators withheld from all weights.

The original anchored-bridge arithmetic is reused only as a Gaussian-factor
replacement kernel. Its input here is explicitly the noncalibrator factor.
No anchored likelihood or extra H0 factor enters the replacement weights.
"""
import argparse
import copy
import importlib.metadata
import json
import os
from pathlib import Path
import time

import numpy as np

import anchored_bridge as base
import calibration_holdout as holdout
from exact_correction import verify_record
from measurement_summary import summarize_run

HERE=Path(__file__).resolve().parent
ROOT=base.ROOT
DESIGN=HERE/'calibration-holdout-design.json'
GATES=base.GATES
VALIDATION=ROOT/'studies/unified_cosmology/results/inference/calibration-holdout-bridge-validation.json'
KERNEL_REVIEW=ROOT/'studies/unified_cosmology/results/inference/calibration-holdout-review.json'


def validation_guard():
    """Require current, source-bound numerical and consumer controls; no bypass."""
    expected=[(KERNEL_REVIEW,'passed_independent_calibration_holdout_review',
               [Path(holdout.__file__),DESIGN,HERE/'calibration_holdout_review.py']),
              (VALIDATION,'passed_synthetic_calibration_holdout_bridge_checks_no_observational_evaluation',
               [Path(__file__),Path(holdout.__file__),DESIGN,
                HERE/'calibration_holdout_bridge_validate.py',HERE/'calibration_holdout_review.py'])]
    hashes={}
    for path,status,required in expected:
        record=json.loads(path.read_text())
        assert record['status']==status,'Unqualified holdout validation prerequisite.'
        assert all(record['source_sha256'].get(base.relative(p))==base.digest(p)for p in required), 'Stale holdout validation prerequisite.'
        for field in ['source_sha256','input_sha256']:
            base.verify_hashes(record.get(field,{}))
        hashes[base.relative(path)]=base.digest(path)
    validation=json.loads(VALIDATION.read_text())
    assert validation['input_sha256'].get(base.relative(KERNEL_REVIEW))==base.digest(KERNEL_REVIEW),'Kernel review changed after bridge validation.'
    assert validation['CMB_spectrum_calls']==validation['CAMB_background_calls']==validation['observational_bridge_points']==0
    return hashes


class NoncalibratorInterface:
    """Expose only the HF factor as the numerical replacement density."""
    def __init__(self, model, design):
        self.model=model
        self.z_hd_noncalibrator=model.release.z_hd_noncalibrator
        self.design=design

    def evaluate(self, distances):
        result=self.model.evaluate(distances)
        assert result['rows_noncalibrator']==1580 and result['rows_calibrator']==77
        names=['noncalibrator_loglike','noncalibrator_chi2','conditional_calibrator_loglike',
               'conditional_calibrator_chi2','calibration_contrast_mag','contrast_sigma_mag',
               'contrast_z','full_loglike','full_loglike_closure','released_full_loglike_closure','conditional_logcdf','conditional_logsf']
        assert np.isfinite([result[k]for k in names]).all(),'Nonfinite HF/conditional calibration result.'
        assert result['contrast_sigma_mag']>0
        assert abs(result['full_loglike_closure'])<=self.design['full_loglike_closure_absolute_tolerance']
        assert abs(result['released_full_loglike_closure'])<=self.design['full_loglike_closure_absolute_tolerance']
        assert max(result['conditional_logcdf'],result['conditional_logsf'])<=1e-12
        assert abs(np.logaddexp(result['conditional_logcdf'],result['conditional_logsf']))<1e-10
        return {'loglike':result['noncalibrator_loglike'],'chi2':result['noncalibrator_chi2'],
                'withheld_calibration':result}


def background_densities(point, extra, old, new, design):
    # Exactly one background at the original Dovekie/new HF union, using the
    # already validated thermal adapter. Full row order is restored by base.
    value=base.background_densities(point,extra,old,NoncalibratorInterface(new,design))
    replacement=value.pop('anchored_SN')
    value['HF_loglike']=replacement['loglike'];value['HF_chi2']=replacement['chi2']
    value['holdout']=replacement['withheld_calibration']
    value['replacement_semantics']='1580 noncalibrator rows only;77 calibrators enter conditional diagnostics, never weights.'
    return value


def replace_SN(record, background, components, design):
    # Only these two HF fields are visible to the existing weight kernel.
    arithmetic={'old_SN_reconstructed_loglike':background['old_SN_reconstructed_loglike'],
                'anchored_SN':{'loglike':background['HF_loglike'],'chi2':background['HF_chi2']}}
    row=base.replace_SN(record,arithmetic,components,design)
    if row['status']=='finite_anchored_replacement':row['status']='finite_HF_replacement'
    row['holdout']=copy.deepcopy(background['holdout'])
    row['calibrator_observations_used_in_weights']=False
    return row


def summarize(records, rows, groups, gates):
    compatible=[dict(row,status='finite_anchored_replacement') if row['status']=='finite_HF_replacement' else dict(row) for row in rows]
    result=base.summarize(records,compatible,groups,gates)
    result['status']={'qualified_conditional_anchored_bridge':'qualified_conditional_HF_replacement',
                      'insufficient_anchored_overlap_or_stability':'insufficient_HF_overlap_or_stability',
                      'failed_anchored_bridge_evaluation':'failed_calibration_holdout_bridge_evaluation'}[result['status']]
    result['withheld_calibration_predictive']=None
    if result['qualified_under_declared_numerical_gates']:
        result['withheld_calibration_predictive']=holdout.predictive_summary(
            [row['holdout']for row in rows],np.array([row['target_logweight']for row in rows]),np.asarray(groups))
    result['calibrator_observations_used_in_weights']=False
    result['prediction_scope']='Released-likelihood Bayesian holdout conditional on unchanged HF/CMB/BAO assumptions; not a raw-blind or frequentist sigma test.'
    return result


def identify_target(anchored_config, model, design, sources):
    reference=base.anchored.identify(anchored_config)
    target=copy.deepcopy(reference);target.pop('identity')
    target['configuration']['likelihood']['released_sn']={
        'semantic_factor':'calibration_holdout.ReleasedHoldout.evaluate.noncalibrator_loglike',
        'data_file':str(base.anchored.DATA.resolve()),
        'selected_original_row_indices':model.release.original_indices[~model.release.calibrator].tolist(),
        'rows':1580,'flat_M_measure':'Same improper flat dM; normalized1580-row Gaussian with intercept determinant.',
        'calibrator_likelihood_in_target':False}
    target['configuration_role']='Canonical inherited physical/prior/probe configuration with an explicit semantic SN factor; not an executable Cobaya factory.'
    target['source_calibration_release']=target.pop('anchored_calibration')
    target['full_anchored_release_reference_identity']=reference['identity']
    target['sample_semantics']={'name':design['target_sample'],'replaces_other_SN_compilations':True,
        'absolute_Cepheid_calibrator_means_in_target':False,'extra_H0_or_Cepheid_factor':False,
        'noncalibrator_rows':1580,'withheld_calibrator_rows':77,
        'covariance':'HF marginal block of the same symmetric released total covariance; full blocks retained only for proper conditional prediction.',
        'no_evidence_claim_under_flat_M':True}
    target['source_sha256'].update(sources)
    target['identity']=base.identity(target)
    return target


def actual(folder, summary_path, cache):
    # Sole evaluation entry: fresh real parent qualification precedes every
    # background, release decomposition, predictive calculation and cache write.
    parent=summarize_run(folder,summary_path)
    assert parent['qualified_under_declared_numerical_gates'] is True
    prerequisite_hashes=validation_guard()
    assert all(os.environ.get(k)=='1'for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','CLIPY_NOJAX'])
    folder,summary_path,cache=[Path(p).resolve()for p in [folder,summary_path,cache]]
    assert cache.is_relative_to(ROOT/'.work')
    base.verify_hashes(parent['input_sha256'])
    design=json.loads(DESIGN.read_text());gates=json.loads(GATES.read_text())['overlap_gates']
    summary=json.loads(summary_path.read_text());selection_path=ROOT/summary['selection_path']
    assert selection_path.parent.parent==folder
    selection=json.loads(selection_path.read_text());settings=selection['settings']
    frozen=json.loads((folder/'run-0.json').read_text())['target_identity']
    assert parent['target_identity']==frozen['identity']
    old_config,anchored_config=base.configurations(settings,frozen,design)
    old_path=Path(old_config['likelihood']['released_sn']['data_file']).resolve()
    assert base.digest(old_path)==frozen['sample_sha256']
    source_paths=[Path(__file__),DESIGN,GATES,Path(holdout.__file__),Path(base.__file__),base.DESIGN,
        HERE/'anchored_adapter.py',HERE/'anchored-design.json',HERE/'expansion_history.py',
        HERE/'luminosity_sensitivity.py',HERE/'measurement_summary.py',HERE/'probe_omission.py',
        HERE/'exact_correction.py',HERE/'target_identity.py',HERE/'modern_fast.py',HERE/'modern_run.py',
        HERE.parent/'distance_ladder/calibration_interface.py']
    sources={base.relative(p):base.digest(p)for p in source_paths}
    with np.load(old_path,allow_pickle=False)as data:
        old=base.IntegratedLuminosity(*(data[k]for k in ['zHD','zHEL','MU','covariance']))
    new=holdout.ReleasedHoldout()
    assert len(new.release.data)==1657 and np.count_nonzero(new.release.calibrator)==77
    target=identify_target(anchored_config,new,design,sources)
    assert target['assets']==frozen['assets']
    assert all(target['versions'][k]==v for k,v in frozen['versions'].items())
    lineage={'qualified_parent_inputs':parent['input_sha256'],'source_sha256':sources,
        'validation_prerequisites_sha256':prerequisite_hashes,
        'parent_proposal_target_identity':frozen['identity'],'native_target':target,
        'parent_settings':settings,'parent_SN_sha256':base.digest(old_path),
        'all_non_SN_factors_and_priors_identical':True,'calibrator_observations_used_in_weights':False}
    cache_identity=base.identity(lineage)
    cache.mkdir(parents=True,exist_ok=True)
    lineage_path=cache/'lineage.json';base.payload(lineage_path,lineage)
    ledger=cache/'record-hashes.json'
    if ledger.exists():base.verify_hashes(json.loads(ledger.read_text()))
    records=[];rows=[];hashes_out={base.relative(lineage_path):base.digest(lineage_path)}
    start=time.monotonic()
    for index,point in enumerate(selection['points']):
        native_path=selection_path.parent/f'{index:05d}.json'
        record=json.loads(native_path.read_text());verify_record(record)
        assert record['point']==point and record['status']=='finite'
        records.append(record);path=cache/f'{index:05d}.json'
        if path.exists():row=base.read_record(path,cache_identity,index,point,base.digest(native_path))
        else:
            try:
                background=background_densities(point,old_config['theory']['camb']['extra_args'],old,new,design)
                row=replace_SN(record,background,old_config['likelihood'],design)
                row['background']={k:v for k,v in background.items()if k!='holdout'}
            except Exception as error:
                row={'status':'failed_background_or_density','error':repr(error)}
            row.update(identity=cache_identity,index=index,point=point,native_record_sha256=base.digest(native_path))
            row['payload_sha256']=base.identity(row);base.payload(path,row)
        rows.append(row);hashes_out[base.relative(path)]=base.digest(path)
        if (index+1)%100==0:print(json.dumps({'calibration_holdout_points':index+1,'seconds':time.monotonic()-start}),flush=True)
    base.payload(ledger,hashes_out)
    result=summarize(records,rows,selection['groups'],gates)
    result.update(parent_settings=parent['settings'],native_target=target,source_sha256=sources,bridge_identity=cache_identity,
        validation_prerequisites_sha256=prerequisite_hashes,
        lineage_path=base.relative(lineage_path),lineage_sha256=base.digest(lineage_path),
        cache_manifest_path=base.relative(ledger),cache_manifest_sha256=base.digest(ledger),
        correction_summary_path=base.relative(summary_path),correction_summary_sha256=base.digest(summary_path),
        parent_qualified_under_declared_numerical_gates=True,CMB_spectrum_calls=0,
        background_calls_in_recorded_cohort=sum(r.get('background',{}).get('background_calls',0)for r in rows),
        cached_or_computed_points=len(rows),limits=design['limits'],
        runtime_environment={k:os.environ.get(k)for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','CLIPY_NOJAX']},
        units={'H0':'km/s/Mpc','rdrag':'Mpc','calibration_contrast':'mag','w_wa_q_j':'dimensionless'},
        weight_rule='parent stored native/proposal + HF-only SN loglike - stored native Dovekie SN loglike',
        native_precision_screen='Separate; neither parent qualification nor this bridge establishes numerical-accuracy convergence.',
        status_of_failed_weight_summaries='Diagnostics only; no admitted posterior or withheld-calibration predictive claim.')
    for mapping in [sources,parent['input_sha256'],frozen['source_sha256'],target['source_sha256'],
                    target['source_calibration_release']['input_and_audit_sha256'],{base.relative(old_path):frozen['sample_sha256']}]:
        base.verify_hashes(mapping)
    assert all(importlib.metadata.version(k)==v for k,v in target['versions'].items())
    assert validation_guard()==prerequisite_hashes,'Holdout validation changed during execution.'
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--chain-folder',type=Path,required=True);p.add_argument('--correction-summary',type=Path,required=True)
    p.add_argument('--cache',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=actual(a.chain_folder,a.correction_summary,a.cache)
    a.output.parent.mkdir(parents=True,exist_ok=True);base.payload(a.output,result)
    print(json.dumps({'status':result['status'],'HF_qualified':result['qualified_under_declared_numerical_gates'],'CMB_spectrum_calls':0}))


if __name__=='__main__':main()
