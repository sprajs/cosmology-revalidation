"""Summarize a preserved dense-metric retry; no model or observed-data calls."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import bayesn_synthetic_geometry as geometry
import bayesn_synthetic_retune as retry
import numpy as np
from scipy.special import logsumexp


def score_draws(mean_flux, flux, errors):
    """Rebuild one complete-vector score from saved shared posterior draws."""
    assert mean_flux.shape == (4,1000,len(flux))
    assert errors.shape == flux.shape and np.all(errors>0)
    assert all(np.isfinite(x).all() for x in (mean_flux,flux,errors))
    joint=-.5*np.sum(((flux-mean_flux)/errors)**2+np.log(2*np.pi*errors**2),axis=-1)
    # A single scale for all chains preserves differences in their support.
    normalized=np.exp(joint-joint.max())
    average=normalized.mean()
    chain_se=np.std(normalized.mean(axis=1),ddof=1)/2
    block_se=np.std(np.stack([normalized[c,k:k+50].mean() for c in range(4)
                             for k in range(0,1000,50)]),ddof=1)/np.sqrt(80)
    return joint,dict(log_predictive_density=float(logsumexp(joint)-np.log(4000)),
                      logscore_MCSE_delta=float(max(chain_se,block_se)/average),
                      between_chain_MCSE_delta=float(chain_se/average),
                      block_MCSE_delta=float(block_se/average),
                      chain_log_predictive_density=(logsumexp(joint,axis=1)-np.log(1000)).tolist(),
                      score_support_passed=bool(np.all(normalized.mean(axis=1)>0)),
                      likelihood_weight_ESS=float(normalized.sum()**2/np.sum(normalized**2)))


def predictive_record(work,result):
    """Only read a synthetic held-out vector after both frozen samplers pass."""
    path=work/'case-00/synthetic-predictive-score.json'
    if not path.exists():
        return None
    assert all(row['official_optical_diagnostics']['sampler_passed'] for row in result['arms'].values())
    saved=json.loads(path.read_text())
    assert saved['source_sha256']==retry.model.source_hashes()
    assert saved['case_index']==0
    folder=path.parent
    freeze_path=folder/'pre-score-optical-freeze.json'
    assert retry.model.sha(freeze_path)==saved['optical_freeze_sha256']
    freeze=json.loads(freeze_path.read_text())
    assert freeze['source_sha256']==retry.model.source_hashes()
    for name,checksum in freeze['optical_records_sha256'].items():
        assert retry.model.sha(work/name)==checksum
    nir_path=folder/'synthetic-nir.json'
    prepared=json.loads((work/'prepared.json').read_text())
    assert retry.model.sha(nir_path)==saved['synthetic_nir_sha256']==prepared['cases'][0]['nir_sha256']
    data=json.loads(nir_path.read_text())
    assert data['scope']=='synthetic' and all(r['role']=='heldout_NIR' for r in data['rows'])
    source_paths=[path,freeze_path,nir_path,work/'before-synthetic-nir-copy.json']
    copy_freeze=json.loads(source_paths[-1].read_text())
    assert copy_freeze['source_sha256']==retry.dependencies()
    for name,checksum in copy_freeze['optical_sha256'].items():
        assert retry.model.sha(name)==checksum
    rebuilt={}
    for arm in ('LCDM','broad'):
        archive=folder/f'{arm}-synthetic-prediction.npz'
        assert retry.model.sha(archive)==saved['arms'][arm]['draw_archive_sha256']
        with np.load(archive) as arrays:
            joint,summary=score_draws(arrays['mean_flux'],np.asarray(data['flux']),np.asarray(data['errors']))
            assert np.array_equal(joint,arrays['joint_log_likelihood'])
            latent=arrays['mean_flux'].reshape(4000,-1)
            for key,value in dict(posterior_mean_flux=latent.mean(axis=0),
                                  posterior_mean_flux_covariance=np.cov(latent,rowvar=False,ddof=1),
                                  predictive_flux_covariance=np.cov(latent,rowvar=False,ddof=1)+np.diag(np.asarray(data['errors'])**2)).items():
                assert np.array_equal(np.asarray(saved['arms'][arm][key]),value)
        for key,value in summary.items():
            assert np.allclose(saved['arms'][arm][key],value,rtol=1e-13,atol=1e-13)
        assert saved['arms'][arm]['row_order']==data['rows']
        rebuilt[arm]=summary
        source_paths.append(archive)
    difference=rebuilt['LCDM']['log_predictive_density']-rebuilt['broad']['log_predictive_density']
    mcse=float(np.hypot(*(rebuilt[a]['logscore_MCSE_delta'] for a in ('LCDM','broad'))))
    assert np.isclose(difference,saved['LCDM_minus_broad_logscore'],rtol=1e-13,atol=1e-13)
    assert np.isclose(mcse,saved['independent_arm_difference_MCSE'],rtol=1e-13,atol=1e-13)
    passed=mcse<.05 and all(x['score_support_passed'] for x in rebuilt.values())
    assert passed==saved['score_precision_passed']
    for source in source_paths:
        result['input_sha256'][str(source)]=geometry.sha(source)
    return dict(arms=rebuilt,LCDM_minus_broad_logscore=difference,
                independent_arm_difference_MCSE=mcse,score_precision_passed=passed,
                joint_heldout_vector_length=len(data['flux']),
                source_path=str(path.relative_to(geometry.ROOT)),source_sha256=geometry.sha(path),
                scope='One synthetic realization only. Monte Carlo precision is not population uncertainty or predictive coverage; no observed infrared measurements.')


def summarize(work, outcome, full_path):
    plan_path=work/'retuning-plan.json'
    plan=json.loads(plan_path.read_text())
    assert plan['source_sha256']==retry.dependencies()
    assert retry.model.sha(plan['validation_path'])==plan['validation_sha256']
    assert plan['original_target_source_sha256']==retry.model.source_hashes()
    assert retry.model.sha(work/'prepared.json')==plan['original_prepared_sha256']
    result=geometry.analyze(work,0,outcome)
    old='The final adapted metric, step size and full divergent trajectories were not archived; their exact Hamiltonian error cannot be reconstructed.'
    assert result['caveats'].count(old)==1
    result['caveats'][result['caveats'].index(old)]=(
        'This retry archives the final adapted metric and step size. Full divergent leapfrog trajectories are still absent; '
        'retained states and the final metric do not reconstruct an earlier failing trajectory.')
    result['sampler_attempt']='dense_metric_target_accept_0_99'
    result['sampler_settings']=json.loads(retry.DESIGN.read_text())['sampler']
    result['retuning_source_sha256']=plan['source_sha256']
    result['source_sha256'].update({str(Path(__file__).resolve().relative_to(geometry.ROOT)):geometry.sha(__file__)})
    for path in (plan_path,work/'retuning-execution.json',Path(plan['validation_path'])):
        result['input_sha256'][str(path)]=geometry.sha(path)
    for arm,entry in result['arms'].items():
        folder=work/'case-00'/arm
        started=folder/'worker-start.json'
        if started.exists():
            record=json.loads(started.read_text())
            assert record['retuning_source_sha256']==retry.dependencies()
            assert record['optical_sha256']==plan['original_optical_sha256']
            result['input_sha256'][str(started)]=geometry.sha(started)
        metrics=[]
        for chain in entry['chains']:
            path=folder/chain['archive']
            with np.load(path) as arrays:
                assert len(chain['dense_metric'])==1
                definition=chain['dense_metric'][0]
                metric=arrays[definition['array']]
                assert metric.shape==(47,47) and list(metric.shape)==definition['shape']
                assert set(definition['latent_site_order'])=={'D','AV','RV','theta','epsilon_white','tau'}
                assert np.isfinite(metric).all()
                step=float(arrays['adapt_step_size'])
                assert step==chain['final_adapted_step_size'] and step>0
                eigenvalues=np.linalg.eigvalsh(metric)
                metrics.append(dict(chain=chain['chain'],final_adapted_step_size=step,
                                    latent_site_order=definition['latent_site_order'],
                                    maximum_asymmetry=float(np.max(abs(metric-metric.T))),
                                    minimum_eigenvalue=float(eigenvalues[0]),maximum_eigenvalue=float(eigenvalues[-1]),
                                    positive_definite=bool(eigenvalues[0]>0),
                                    eigenvalue_condition_number=float(eigenvalues[-1]/eigenvalues[0]) if eigenvalues[0]>0 else None))
        entry['adaptation_diagnostics']=metrics
        failure=folder/'worker-failure.json'
        if failure.exists():
            entry['worker_failure']=json.loads(failure.read_text())
            result['input_sha256'][str(failure)]=geometry.sha(failure)
    prediction=predictive_record(work,result)
    if prediction is not None:
        result['synthetic_prediction']=prediction
        result['caveats'][0]=('The numerical sampling gates and saved synthetic infrared-vector score concern one engineering realization; '
                              'they do not establish population coverage or an observed-data result.')
    retry.model.dump(full_path,result)
    compact=geometry.compact(result,full_path)
    compact.update(sampler_attempt=result['sampler_attempt'],sampler_settings=result['sampler_settings'],
                   retuning_source_sha256=result['retuning_source_sha256'])
    for arm in compact['arms']:
        compact['arms'][arm]['adaptation_diagnostics']=result['arms'][arm]['adaptation_diagnostics']
        if 'worker_failure' in result['arms'][arm]:
            compact['arms'][arm]['worker_failure']=result['arms'][arm]['worker_failure']
    if prediction is not None:
        compact['synthetic_prediction']=prediction
        compact['status']='synthetic_sampling_and_score_precision_passed' if prediction['score_precision_passed'] else 'synthetic_prediction_precision_gate_failed'
        compact['scope']='One synthetic optical-posterior and held-out infrared-vector test; no observed fit, population comparison or coverage claim.'
    return compact


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,required=True)
    parser.add_argument('--external-outcome',choices=('running','completed','timeout','failed'),required=True)
    parser.add_argument('--full-output',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=summarize(args.work.resolve(),args.external_outcome,args.full_output.resolve())
    retry.model.dump(args.output.resolve(),result)
    print(json.dumps({arm:{'completed_chains':row['completed_chains'],
                           'divergences_by_chain':row['divergences_by_chain']} for arm,row in result['arms'].items()}))
