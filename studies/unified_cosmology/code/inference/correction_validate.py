"""Independent analytically soluble check of raw importance-weight summaries."""
import hashlib
import json
import tempfile
from pathlib import Path
import numpy as np
from exact_correction import summarize,ROOT,correction_dependencies,record_digest,evaluate


def main():
    cache_checks=[]
    with tempfile.TemporaryDirectory() as directory:
        path=Path(directory)/'00000.json'
        row={'index':0,'point':{'x':1.},'target_identity':'synthetic',
             'status':'finite','log_weight':.2}
        row['payload_sha256']=record_digest(row)
        path.write_text(json.dumps(row))
        assert evaluate((0,{'x':1.},directory,'synthetic',0.))==row
        cache_checks.append('unchanged_payload_reused_without_native_call')
        for mutation in ['numeric_change','missing_payload_hash','wrong_point','wrong_identity','wrong_index']:
            changed=dict(row)
            if mutation=='numeric_change':changed['log_weight']=.3
            elif mutation=='missing_payload_hash':del changed['payload_sha256']
            elif mutation=='wrong_point':changed['point']={'x':2.}
            elif mutation=='wrong_identity':changed['target_identity']='different'
            else:changed['index']=1
            if mutation in ['wrong_point','wrong_identity','wrong_index']:
                changed['payload_sha256']=record_digest(changed)
            path.write_text(json.dumps(changed))
            try:evaluate((0,{'x':1.},directory,'synthetic',0.))
            except AssertionError:cache_checks.append(mutation+'_rejected')
            else:raise AssertionError('Altered cached native record accepted')
    rng=np.random.default_rng(272642);n=40000;shift=np.array([.2,-.3])
    draws=rng.normal(size=(n,2))
    # Unit-covariance target N(shift,I), proposal N(0,I): exact log density ratio.
    logweights=draws@shift-.5*shift@shift
    records=[{'status':'finite','point':{'x':float(x),'y':float(y)},'derived':{},'log_weight':float(w)}
             for (x,y),w in zip(draws,logweights)]
    report=summarize(records,np.repeat(np.arange(4),n//4))
    mean_errors=np.array([report['posterior'][k]['mean'] for k in ['x','y']])-shift
    mcse=np.array([report['posterior'][k]['weighted_batch_mean_mcse'] for k in ['x','y']])
    sd_errors=np.array([report['posterior'][k]['sd'] for k in ['x','y']])-1
    expected_ess_fraction=float(np.exp(-shift@shift))
    assert np.all(abs(mean_errors)<4*mcse)
    assert np.max(abs(sd_errors))<.025
    assert abs(report['raw_weight_ess']/n-expected_ess_fraction)<.02
    assert report['status']=='passed_importance_weight_gates'
    assert not report['failed_gates']
    # Constant weights have no Pareto tail to estimate; do not manufacture k.
    groups=np.repeat(np.arange(4),n//4)
    uniform=[dict(row,log_weight=0.) for row in records]
    constant=summarize(uniform,groups)
    assert constant['status']=='passed_importance_weight_gates'
    assert constant['pareto_k'] is None
    assert constant['pareto_k_status']=='constant_weights_no_tail_to_fit'
    assert abs(constant['raw_weight_ess']-n)<1e-8
    # An entire chain can receive numerical zero mass: retain raw weights and
    # record null undefined means, never NaN or a falsely passing status.
    collapsed=[dict(row,log_weight=1000. if index==0 else 0.) for index,row in enumerate(records)]
    zero=summarize(collapsed,groups)
    assert zero['status']=='insufficient_importance_overlap'
    assert zero['per_chain_weight_diagnostics']['1']['raw_weight_ESS']==0
    assert zero['independent_chain_weighted_means']['x'][1] is None
    json.dumps(zero,allow_nan=False)
    # Global ESS/k alone cannot pass when independently sampled chains disagree.
    shifted=[]
    for row,group in zip(uniform,groups):
        point=dict(row['point']);point['x'] += {0:1.5,1:-1.5,2:0.,3:0.}[group]
        shifted.append(dict(row,point=point))
    drift=summarize(shifted,groups)
    assert drift['raw_weight_ess']>=400 and drift['pareto_k'] is None
    assert drift['status']=='insufficient_importance_overlap' and 'x_chain_means' in drift['failed_gates']
    # A sufficiently large global ESS can hide low ESS in one chain.
    small=8000;local_groups=np.repeat(np.arange(4),small//4)
    local=[dict(row,log_weight=0.) for row in records[:small]]
    for index in range(6000,small):local[index]['log_weight']=float(np.log(50. if index<6040 else 1e-100))
    local_check=summarize(local,local_groups)
    assert local_check['raw_weight_ess']>=400
    assert 'chain_3_weight_ESS' in local_check['failed_gates']
    assert local_check['status']=='insufficient_importance_overlap'
    failed=summarize([dict(records[0],status='nonfinite',log_weight=float('nan'))],[0])
    assert failed['status']=='failed_exact_evaluation'
    assert failed['failures'][0]['log_weight'] is None
    json.dumps(failed,allow_nan=False)
    # Existing downstream users retain the same posterior and MCSE keys.
    assert set(report['posterior']['x']['batch_sensitivity'])=={'10','20','40'}
    assert set(report['posterior']['x']['batch_sensitivity']['20'])=={
        'mean_mcse','minimum_points_per_batch','minimum_batch_weight','maximum_batch_weight'}
    for checked in [report,constant,drift,local_check]:json.dumps(checked,allow_nan=False)
    output={'status':'passed','draws':n,'target_mean':shift.tolist(),
        'cached_native_integrity_checks':cache_checks,
        'mean_errors':mean_errors.tolist(),'mean_mcse':mcse.tolist(),'sd_errors':sd_errors.tolist(),
        'observed_weight_ess_fraction':report['raw_weight_ess']/n,
        'analytic_asymptotic_weight_ess_fraction':expected_ess_fraction,
        'pareto_k':report['pareto_k'],
        'registered_weighted_stability_gates':report['weighted_stability_gates'],
        'all_known_gaussian_gates_pass':not report['failed_gates'],
        'constant_weights_no_fictitious_tail':True,
        'zero_chain_weights_fail_with_null_means_and_no_clipping':True,
        'chain_drift_fails_despite_full_global_ESS':True,
        'low_per_chain_ESS_fails_despite_sufficient_global_ESS':True,
        'nonfinite_exact_failure_preserved_as_strict_JSON':True,
        'legacy_posterior_and_batch_keys_preserved':True,
        'correction_dependency_sha256':correction_dependencies(),
        'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('exact_correction.py')]},
        'scope':'Algebra and numerical summary check under known Gaussian density ratios, not validation of cosmological overlap or survey coverage.'}
    (ROOT/'studies/unified_cosmology/results/inference/correction-validation.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))


if __name__=='__main__':main()
