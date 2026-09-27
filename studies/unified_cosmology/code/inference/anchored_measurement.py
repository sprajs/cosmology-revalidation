"""Admit only independently converged and native-qualified anchored measurements."""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp

import anchored_sampling as sampling
import anchored_correction as correction
from exact_correction import verify_record,summarize
from measurement_summary import PARAMETERS,weighted_fraction


def summarize_run(folder,summary_path):
    folder=Path(folder).resolve();summary_path=Path(summary_path).resolve()
    summary=json.loads(summary_path.read_text())
    assert summary['schema']=='dedicated-anchored-native-summary-v1'
    assert summary['status']=='passed_importance_weight_gates'
    selection_path=(sampling.ROOT/summary['selection_path']).resolve()
    assert selection_path.name=='selection.json' and selection_path.parent.parent==folder
    assert sampling.digest(selection_path)==summary['selection_sha256']
    selection=json.loads(selection_path.read_text())
    assert selection['schema']=='dedicated-anchored-native-selection-v1'
    diagnostics_path=sampling.ROOT/selection['diagnostics_path']
    manifest,check,completed=correction.qualified_chain_inputs(folder,diagnostics_path)
    assert completed==selection['sampler_completion']
    target=manifest['target_identity']
    assert selection['settings']==manifest['arguments']
    assert summary['native_accuracy']==selection['native_accuracy']==target['native_accuracy']
    assert summary['code_sha256']==sampling.digest(correction.__file__)
    bound=correction.dependencies()
    assert bound==selection['correction_dependency_sha256']==summary['correction_dependency_sha256']
    identity=sampling.identity({'target_identity':target['identity'],'correction_dependencies':bound})
    assert identity==selection['correction_identity']==summary['correction_identity']
    assert target['identity']==selection['proposal_target_identity']==summary['proposal_target_identity']
    assert sampling.relative(diagnostics_path)==summary['diagnostics_path']
    assert sampling.digest(diagnostics_path)==selection['diagnostics_sha256']==summary['diagnostics_sha256']
    assert selection['chain_sha256']==check['input_sha256']
    sampling.verify_hashes(selection['rank_manifest_sha256'])
    assert set(selection['rank_manifest_sha256'])=={sampling.relative(folder/f'run-{i}.json') for i in range(4)}
    names=[k for k,v in target['configuration']['params'].items() if isinstance(v,dict) and 'prior'in v]
    paths=sorted(folder.glob('chain.[0-9]*.txt'))
    selected=correction.select_points(paths,names,check,len(selection['points']))
    assert all(selection[k]==v for k,v in selected.items()),'Selected cohort differs from frozen stratification.'
    expected={sampling.relative(selection_path.parent/f'{i:05d}.json') for i in range(len(selection['points']))}
    assert set(summary['native_record_sha256'])==expected
    sampling.verify_hashes(summary['native_record_sha256'])
    sampling.verify_hashes(summary['initialization_record_sha256'])
    records=[]
    for index,point in enumerate(selection['points']):
        record=json.loads((selection_path.parent/f'{index:05d}.json').read_text());verify_record(record)
        assert record['schema']=='dedicated-anchored-native-point-v1' and record['status']=='finite'
        assert record['index']==index and record['point']==point and record['target_identity']==identity
        assert record['stored_chain_logpost']==selection['stored_logposts'][index]
        assert abs(record['log_weight']-record['exact_logpost']+record['proposal_logpost'])<1e-9
        assert abs(record['proposal_chain_logpost_difference'])<=.01
        assert abs(record['proposal_logpost']-record['stored_chain_logpost']-record['proposal_chain_logpost_difference'])<1e-9
        exact_prior=record['exact_logpost']-sum(record['exact_loglikes'].values())
        proposal_prior=record['proposal_logpost']-sum(record['proposal_loglikes'].values())
        assert abs(exact_prior-proposal_prior)<=1e-9
        assert abs(exact_prior-proposal_prior-record['native_proposal_logprior_difference'])<1e-9
        assert set(record['exact_loglikes'])==set(record['proposal_loglikes'])==set(target['native_configuration']['likelihood'])
        assert set(record['initialization_record_sha256'])<=set(summary['initialization_record_sha256'])
        for path,expected in record['initialization_record_sha256'].items():
            assert summary['initialization_record_sha256'][path]==expected
        records.append(record)
    groups=np.asarray(selection['groups'])
    assert set(groups)=={0,1,2,3} and len(records)>=2000
    recomputed=summarize(records,groups)
    assert recomputed=={key:summary[key] for key in recomputed}
    assert recomputed['status']=='passed_importance_weight_gates' and not recomputed['failed_gates']
    logweights=np.array([r['log_weight'] for r in records]);weights=np.exp(logweights-logsumexp(logweights))
    values={name:np.array([dict(r['point'],**r['derived'])[name] for r in records]) for name in recomputed['posterior']}
    published=[name for name in PARAMETERS if name in values]
    matrix=np.column_stack([values[name] for name in published]);centered=matrix-weights@matrix
    signs={name+'_below_zero':weighted_fraction(values[name]<0,weights,groups) for name in ['q0','q05','q1','j0']}
    signs['accelerating_with_decreasing_scale_factor_acceleration_today']=weighted_fraction((values['q0']<0)&(values['j0']<0),weights,groups)
    boundaries={}
    for name,definition in target['configuration']['params'].items():
        if name not in values or not isinstance(definition,dict):continue
        prior=definition.get('prior',{})
        if not isinstance(prior,dict) or not {'min','max'}<=set(prior):continue
        low,high=prior['min'],prior['max'];margin=.01*(high-low)
        boundaries[name]={'prior':[low,high],
                         'fraction_within_one_percent_of_lower_bound':float(weights@(values[name]<low+margin)),
                         'fraction_within_one_percent_of_upper_bound':float(weights@(values[name]>high-margin))}
    support=None
    if 'w'in values and 'wa'in values:
        early=values['w']+values['wa'];assert np.max(early)<=1e-10
        support={'constraint':'w0+wa<=0','diagnostic_boundary_width':.05,
                 'fraction_within_005_of_boundary':float(weights@(early>-.05)),
                 'interpretation':'Diagnostic only for the imposed CAMB support, not excluded-region coverage.'}
    inputs=sampling.verify_manifest(manifest)
    for mapping in [selection['rank_manifest_sha256'],selection['chain_sha256'],summary['native_record_sha256'],summary['initialization_record_sha256'],completed['input_sha256'],bound]:inputs.update(mapping)
    for path in [selection_path,summary_path,diagnostics_path,Path(__file__)]:inputs[sampling.relative(path)]=sampling.digest(path)
    return {'schema':'qualified-dedicated-anchored-measurement-v1',
            'settings':{k:manifest['arguments'][k] for k in ['model','evolution','sample','calibration','native_accuracy']},
            'target_identity':target['identity'],'qualified_under_declared_numerical_gates':True,
            'posterior':{name:recomputed['posterior'][name] for name in published},
            'conditional_sign_fractions':signs,
            'weighted_covariance':{'parameter_order':published,'matrix':((centered*weights[:,None]).T@centered).tolist(),
                                   'normalization':'Raw normalized importance weights; posterior covariance without Bessel correction.'},
            'scalar_box_prior_boundary_fractions':boundaries,'coupled_CAMB_support_boundary':support,
            'native_correction':{key:recomputed[key] for key in ['exact_points','raw_weight_ess','pareto_k','largest_normalized_weight','log_weight_quantiles']},
            'input_sha256':inputs,'limits':json.loads(sampling.DESIGN.read_text())['limitations'],
            'normalization':'Complete anchored Gaussian, including flat-M determinant. No normalized evidence under the improper flat dM measure.'}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('folder',type=Path);p.add_argument('--correction-summary',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args();result=summarize_run(a.folder,a.correction_summary)
    sampling.dump_new(a.output,result);print(json.dumps({'status':result['schema'],'settings':result['settings']}))

if __name__=='__main__':main()
