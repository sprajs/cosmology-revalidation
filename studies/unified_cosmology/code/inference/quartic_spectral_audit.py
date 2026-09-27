"""Data-lineage and independent-QR audit of the declared quartic repair."""
import ast
from datetime import datetime,timezone
import hashlib
import json
from math import comb
from pathlib import Path
import numpy as np
from scipy.linalg import lstsq

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
WORK=ROOT/'.work/unified-cosmology/inference/quartic-spectral'
RESULT=ROOT/'studies/unified_cosmology/results/inference'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    design_path=HERE/'quartic-spectral-design.json';design=json.loads(design_path.read_text())
    acquisition=json.loads((RESULT/'quartic-spectral-acquisition.json').read_text())
    fit=json.loads((RESULT/'quartic-spectral-fit.json').read_text())
    holdout=json.loads((RESULT/'quartic-spectral-holdout.json').read_text())
    checked={}
    for path,digest in design['frozen_inputs_sha256'].items():
        assert sha(ROOT/path)==digest,path;checked[path]=digest
    for record in [acquisition,fit,holdout]:
        for key in ['source_sha256','dependencies_sha256']:
            for path,digest in record.get(key,{}).items():
                assert sha(ROOT/path)==digest,path;checked[path]=digest
    assert acquisition['status']=='complete_exact_acquisition' and acquisition['design_sha256']==sha(design_path)
    assert fit['training_rows']==1198>=990 and fit['training_counts']=={'original':509,'broad':689}
    assert fit['holdout_rows_used']==0 and fit['coefficients_per_output']==495
    manifests={(r['split'],r['index']):r for r in acquisition['rows']};assert len(manifests)==128
    centre=np.array(design['centre']);chol=np.array(design['coordinate_cholesky'])
    counts={};errors=[]
    for split,n in design['counts'].items():
        rng=np.random.default_rng(design['seeds'][split]);normal=rng.normal(size=(n,8))
        widths=np.where(rng.random(n)<.1,2.5,1.5)
        requested=centre+(normal*widths[:,None])@np.array(design['holdout_cholesky'][split]).T
        counts[split]={}
        for index in range(n):
            path=WORK/split/f'{index:04d}.json';meta=json.loads(path.read_text())
            assert sha(path)==manifests[(split,index)]['manifest_sha256']
            assert meta['design_sha256']==sha(design_path) and meta['source_sha256']==acquisition['source_sha256']
            assert np.array_equal(meta['requested_coordinates'],requested[index])
            status=meta['status'];counts[split][status]=counts[split].get(status,0)+1
            if status!='finite_exact':
                errors.append(dict(split=split,index=index,status=status,exception=meta.get('exception'),message=meta.get('message')));continue
            raw=ROOT/meta['file'];assert sha(raw)==meta['sha256']
            checked[str(path.relative_to(ROOT))]=sha(path);checked[str(raw.relative_to(ROOT))]=sha(raw)
            data=np.load(raw,allow_pickle=False)
            assert np.array_equal(data['coordinates'],meta['actual_coordinates'])
            assert np.array_equal(data['loglikes'],meta['loglikes'])
            assert meta['CAMB_Params_max_l']==10251 and meta['provider_Dl_length']==10152
    model_file=ROOT/fit['model_file'];assert sha(model_file)==fit['model_sha256']==holdout['model_sha256']
    model=np.load(model_file,allow_pickle=False)
    assert np.array_equal(model['centre'],centre) and np.array_equal(model['coordinate_cholesky'],chol)
    exponent=model['exponents'];assert exponent.shape==(comb(8+4,4),8)
    assert len({tuple(x) for x in exponent})==495 and set(exponent.sum(axis=1))=={0,1,2,3,4}
    length=int(model['length']);ell=np.unique(np.r_[0,1,np.linspace(2,length-1,24).astype(int)])
    output_indices=np.concatenate([ell+k*length for k in range(5)])
    points,targets=[],[]
    folders=[ROOT/'.work/unified-cosmology/external-probes/spectral-training/train',ROOT/'.work/unified-cosmology/inference/broad-spectral/train']
    for folder in folders:
        for manifest in sorted(folder.glob('*.json')):
            meta=json.loads(manifest.read_text());assert str(manifest.relative_to(ROOT)) in fit['source_sha256']
            if meta['status']!='finite_exact':continue
            path=ROOT/meta['file'];assert sha(path)==meta['sha256']
            data=np.load(path,allow_pickle=False);points.append(data['coordinates'])
            targets.append(data['spectra'].ravel()[output_indices]/model['output_scale'][output_indices])
    assert len(points)==1198
    white=np.linalg.solve(chol,(np.array(points)-centre).T).T
    features=np.ones((len(white),len(exponent)))
    for axis in range(8):features*=white[:,axis,None]**exponent[None,:,axis]
    coefficients,_,rank,_=lstsq(features,np.array(targets),lapack_driver='gelsy')
    assert rank==495
    difference=float(np.max(abs(features@(coefficients-model['coefficients'][:,output_indices]))))
    assert difference<1e-8,difference
    assert len(holdout['rows'])==128 and len({(r['split'],r['index']) for r in holdout['rows']})==128
    assert not any(r['status']=='nonfinite_proposal' for r in holdout['rows'])
    for split in list(counts)+['combined']:
        group=[r for r in holdout['rows'] if r['status']=='compared' and (r['split']==split or split=='combined')]
        delta=np.array([r['delta_loglike_total'] for r in group]);summary=holdout['summaries'][split]
        assert summary['finite_comparisons']==len(group)
        assert np.allclose(np.quantile(delta,[0,.025,.16,.5,.84,.975,1]),summary['delta_loglike_quantiles'],atol=1e-12,rtol=0)
        assert abs(np.sqrt(np.mean(delta**2))-summary['delta_loglike_rms'])<1e-12
        if split!='combined':assert len(group)==counts[split].get('finite_exact',0)
    log=WORK/'acquisition.log';warnings=log.read_text().count('WARNING: mismatch in integrated times')
    checked[str(log.relative_to(ROOT))]=sha(log)
    for path in HERE.glob('quartic_spectral*.py'):ast.parse(path.read_text())
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='passed_numerical_provenance_not_posterior_qualification',
        training_rows=1198,features=495,independent_QR_rank=int(rank),independent_outputs=len(output_indices),
        independent_QR_max_scaled_prediction_difference=difference,holdout_attempts_checked=128,holdout_counts=counts,
        native_failed_attempts=errors,native_integrated_time_warning_lines=warnings,
        warning_scope='Combined interleaved native log retained; warning assignment to individual points is unavailable. Finite likelihoods do not certify native integration accuracy.',
        scientific_prior_changed=False,holdout_values_used_for_fitting=False,
        limitations='Validates declared numerical construction, data lineage and untouched holdout accounting. Does not qualify a cosmological posterior or guarantee no surrogate-induced mode.',
        dependencies_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),design_path,
            RESULT/'quartic-spectral-acquisition.json',RESULT/'quartic-spectral-fit.json',RESULT/'quartic-spectral-holdout.json']},
        checked_sha256=checked)
    (RESULT/'quartic-spectral-validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checked_sha256'},indent=2))


if __name__=='__main__':main()
