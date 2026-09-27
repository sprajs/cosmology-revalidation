"""Joint synthetic NIR scoring only after optical posterior hashes/diagnostics."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.special import logsumexp
from model import ROOT, Kernel, check_parent, configure_cpu, dump, sha, source_hashes


def monte_carlo_score(loglikelihood):
    """Equal posterior draws; keep shared-vector density and chain dependence."""
    assert loglikelihood.shape == (4,1000) and np.isfinite(loglikelihood).all()
    shift = float(np.max(loglikelihood))
    scaled = np.exp(loglikelihood-shift)
    mean = float(scaled.mean())
    assert mean>0
    chain_means = scaled.mean(axis=1)
    se_chain = float(chain_means.std(ddof=1)/np.sqrt(4))
    blocks = scaled.reshape(4,20,50).mean(axis=2)
    se_block = float(blocks.std(ddof=1)/np.sqrt(80))
    se = max(se_chain,se_block)
    return dict(log_predictive_density=float(shift+np.log(mean)),
                logscore_MCSE_delta=se/mean,between_chain_MCSE_delta=se_chain/mean,
                block_MCSE_delta=se_block/mean,blocks_per_chain=20,draws_per_block=50,
                chain_log_predictive_density=(logsumexp(loglikelihood,axis=1)-np.log(1000)).tolist(),
                score_support_passed=bool(np.all(chain_means>0)),
                chain_scaled_mean_underflow=(chain_means==0).tolist(),
                likelihood_weight_ESS=float(scaled.sum()**2/np.sum(scaled**2)),
                caution='MC error only; max of between-chain and contiguous-block estimates using a delta approximation. Not calibrated population uncertainty.')


def score(work,case_index):
    configure_cpu()
    case_dir = work/f'case-{case_index:02d}'
    prepared = json.loads((work/'prepared.json').read_text())
    assert prepared['source_sha256']==source_hashes()
    check_parent(Path(prepared['archive']))
    before = {}
    manifests = {}
    for arm in ('LCDM','broad'):
        folder = case_dir/arm
        manifest = json.loads((folder/'posterior-manifest.json').read_text())
        diagnostic = json.loads((folder/'diagnostics.json').read_text())
        assert diagnostic['posterior_manifest_sha256']==sha(folder/'posterior-manifest.json')
        assert diagnostic['sampler_passed'], f'No predictive score before {arm} sampler diagnostics pass'
        assert manifest['worker_spec']['source_sha256']==source_hashes()
        for r in manifest['chains']:
            assert sha(folder/r['archive'])==r['sha256']
            before[str((folder/r['archive']).relative_to(work))] = r['sha256']
        for name in ('posterior-manifest.json','diagnostics.json'):
            before[str((folder/name).relative_to(work))] = sha(folder/name)
        manifests[arm] = manifest
    dump(case_dir/'pre-score-optical-freeze.json',{'source_sha256':source_hashes(),'optical_records_sha256':before,
                                                'scope':'All optical chains and passing diagnostics frozen before synthetic NIR read.'})
    path = case_dir/'synthetic-nir.json'
    assert sha(path)==prepared['cases'][case_index]['nir_sha256']
    data = json.loads(path.read_text())
    assert data['scope']=='synthetic' and all(r['role']=='heldout_NIR' for r in data['rows'])
    config = work/'release-filter-config.yaml'
    assert sha(config)==prepared['filter_config_sha256']
    kernel = Kernel(Path(prepared['archive']),config,data['metadata'],data['rows'])
    arms = {}
    for arm,manifest in manifests.items():
        means = []
        for r in manifest['chains']:
            x = np.load(case_dir/arm/r['archive'])
            p = np.column_stack([x['D'],x['AV'],x['RV'],x['theta'],x['epsilon_white'],x['tau']])
            assert p.shape == (1000,47)
            means.append(np.concatenate([np.asarray(kernel.batch(p[k:k+32])) for k in range(0,len(p),32)]))
        means = np.stack(means)
        y,sigma = np.asarray(data['flux']),np.asarray(data['errors'])
        # Sum within each complete shared latent draw BEFORE posterior averaging.
        joint = -.5*np.sum(((y-means)/sigma)**2+np.log(2*np.pi*sigma**2),axis=2)
        summary = monte_carlo_score(joint)
        flat = means.reshape(-1,len(y))
        summary['posterior_mean_flux'] = flat.mean(axis=0).tolist()
        summary['posterior_mean_flux_covariance'] = np.cov(flat,rowvar=False,ddof=1).tolist()
        summary['predictive_flux_covariance'] = (np.cov(flat,rowvar=False,ddof=1)+np.diag(sigma**2)).tolist()
        summary['row_order'] = data['rows']
        output = case_dir/f'{arm}-synthetic-prediction.npz'
        assert not output.exists()
        np.savez_compressed(output,mean_flux=means,joint_log_likelihood=joint)
        summary['draw_archive_sha256'] = sha(output)
        arms[arm] = summary
    difference = arms['LCDM']['log_predictive_density']-arms['broad']['log_predictive_density']
    mcse = float(np.hypot(arms['LCDM']['logscore_MCSE_delta'],arms['broad']['logscore_MCSE_delta']))
    for name,digest in before.items():
        assert sha(work/name)==digest
    check_parent(Path(prepared['archive']))
    result = dict(scope='One synthetic dataset, joint NIR-vector prediction; no observed outcome or population inference.',
                  case_index=case_index,arms=arms,LCDM_minus_broad_logscore=difference,independent_arm_difference_MCSE=mcse,
                  score_precision_passed=mcse<.05 and all(x['score_support_passed'] for x in arms.values()),synthetic_nir_sha256=sha(path),
                  optical_freeze_sha256=sha(case_dir/'pre-score-optical-freeze.json'),source_sha256=source_hashes())
    dump(case_dir/'synthetic-predictive-score.json',result)
    return result


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work',type=Path,default=ROOT/'.work/infrared/heldout-synthetic')
    parser.add_argument('--case',type=int,default=0)
    args = parser.parse_args()
    print(json.dumps(score(args.work.resolve(),args.case),indent=2),flush=True)
