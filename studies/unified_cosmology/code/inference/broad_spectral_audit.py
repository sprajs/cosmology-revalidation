"""Rebuild numerical-support identities and independently cross-check the fit."""
from datetime import datetime,timezone
import ast
import hashlib
import json
from pathlib import Path
from math import comb
import numpy as np
from scipy.linalg import lstsq

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
WORK=ROOT/'.work/unified-cosmology/inference/broad-spectral'
RESULT=ROOT/'studies/unified_cosmology/results/inference'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    design_path=HERE/'broad-spectral-design.json';design=json.loads(design_path.read_text())
    original_path=HERE.parent/'external_probes/spectral-training-design.json'
    original=json.loads(original_path.read_text())
    assert sha(original_path)==design['original_design_sha256']
    modern_acquisition=ROOT/'studies/unified_cosmology/results/external_probes/modern-acquisition.json'
    assert sha(modern_acquisition)==design['modern_acquisition_sha256']
    for path,digest in design['source_sha256'].items():assert sha(ROOT/path)==digest,path
    centre=np.array(design['centre']);old=np.array(original['coordinate_cholesky'])
    chol=np.array(design['coordinate_cholesky']);assert np.array_equal(chol,old@np.diag([1]*6+[3,3]))
    assert np.array_equal(centre,np.array(original['centre']))
    cov=chol@chol.T;oldcov=old@old.T
    assert np.array_equal(cov[:6,:6],oldcov[:6,:6]) and np.array_equal(cov[:6,6:],oldcov[:6,6:])
    acquired=json.loads((RESULT/'broad-spectral-training.json').read_text())
    fit=json.loads((RESULT/'broad-spectral-fit.json').read_text())
    holdout=json.loads((RESULT/'broad-spectral-holdout.json').read_text())
    assert acquired['status']=='complete_exact_acquisition'
    assert acquired['design_sha256']==sha(design_path)
    assert len(acquired['rows'])==sum(design['counts'].values())==896
    manifest_index={(r['split'],r['index']):r for r in acquired['rows']}
    assert len(manifest_index)==896
    checks={}
    for record in [acquired,fit,holdout]:
        for path,digest in record.get('source_sha256',{}).items():
            assert sha(ROOT/path)==digest,path;checks[path]=digest
    archive=np.load(WORK/'design-arrays.npz')
    source_points=[];source_spectra=[];statuses={};metadata=[]
    for split,count in design['counts'].items():
        rng=np.random.default_rng(design['seeds'][split]);normal=rng.normal(size=(count,8))
        widths=np.where(rng.random(count)<.1,2.5,1.5)
        requested=centre+(normal*widths[:,None])@chol.T
        assert np.array_equal(archive[split],requested)
        statuses[split]={}
        for index in range(count):
            path=WORK/split/f'{index:04d}.json';record=json.loads(path.read_text())
            assert sha(path)==manifest_index[(split,index)]['manifest_sha256']
            assert record['split']==split and record['index']==index
            assert np.array_equal(record['requested_coordinates'],requested[index])
            assert record['source_sha256']==acquired['source_sha256']
            assert record['design_sha256']==sha(design_path)
            checks[str(path.relative_to(ROOT))]=sha(path)
            status=record['status'];statuses[split][status]=statuses[split].get(status,0)+1
            if status!='finite_exact':continue
            file=ROOT/record['file'];assert sha(file)==record['sha256'];checks[str(file.relative_to(ROOT))]=sha(file)
            with np.load(file,allow_pickle=False) as arrays:
                assert np.array_equal(arrays['requested_coordinates'],requested[index])
                assert np.array_equal(arrays['coordinates'],record['actual_coordinates'])
                assert np.array_equal(arrays['loglikes'],record['loglikes'])
                assert np.isfinite(arrays['spectra']).all()
                assert arrays['spectra'].shape==(5,record['provider_Dl_length'])
                if split=='train':
                    source_points.append(arrays['coordinates']);source_spectra.append(arrays['spectra'])
            assert record['CAMB_Params_max_l']==10251
            extra=record['finalized_theory_extra_args']
            assert extra['lmax']==9001 and extra['lens_margin']==1250 and extra['lens_potential_accuracy']==4
            metadata.append((record['provider_Dl_length'],record['CAMB_Params_max_l']))
    for split,counts in statuses.items():
        for status,count in counts.items():assert acquired['counts'][split][status]==count
    assert fit['training_rows']==len(source_points)>=330 and fit['holdout_rows_used']==0
    assert fit['training_attempts']==768
    # Compare the producer's SVD to an independent rank-revealing QR fit on a
    # deterministic set of spectral outputs; no heldout values enter either.
    model_path=ROOT/fit['model_file'];assert sha(model_path)==fit['model_sha256']==holdout['model_sha256']
    model=np.load(model_path,allow_pickle=False)
    assert np.array_equal(model['centre'],centre) and np.array_equal(model['coordinate_cholesky'],chol)
    exponent=model['exponents'];assert exponent.shape==(comb(8+3,3),8)
    assert len({tuple(e) for e in exponent})==165
    assert set(exponent.sum(axis=1))=={0,1,2,3}
    coordinates=np.linalg.solve(chol,(np.array(source_points)-centre).T).T
    # Explicit products, independent of the production polynomial helper.
    feature=np.ones((len(coordinates),len(exponent)))
    for axis in range(8):feature*=coordinates[:,axis,None]**exponent[None,:,axis]
    spectra=np.array(source_spectra);length=spectra.shape[-1]
    selected=np.unique(np.r_[np.linspace(2,length-1,32).astype(int),[0,1]])
    outputs=np.concatenate([selected+k*length for k in range(5)])
    target=spectra.reshape(len(spectra),-1)[:,outputs]/model['output_scale'][outputs]
    coeff,_,rank,_=lstsq(feature,target,lapack_driver='gelsy')
    assert rank==165
    model_coeff=model['coefficients'][:,outputs]
    coefficient_difference=float(np.max(abs(coeff-model_coeff)))
    fitted_difference=float(np.max(abs(feature@(coeff-model_coeff))))
    assert fitted_difference<1e-8,(coefficient_difference,fitted_difference)
    assert len(holdout['rows'])==128 and sorted(r['index'] for r in holdout['rows'])==list(range(128))
    assert holdout['finite_comparisons']==statuses['holdout'].get('finite_exact',0)
    assert all(r['status']!='nonfinite_proposal' for r in holdout['rows'])
    delta=np.array([r['delta_loglike_total'] for r in holdout['rows'] if r['status']=='compared'])
    assert np.allclose(np.quantile(delta,[0,.025,.16,.5,.84,.975,1]),holdout['delta_loglike_quantiles'],atol=1e-12,rtol=0)
    assert abs(np.sqrt(np.mean(delta**2))-holdout['delta_loglike_rms'])<1e-12
    for path in HERE.glob('broad_spectral_*.py'):ast.parse(path.read_text())
    native_log=WORK/'acquisition.log'
    warning_count=native_log.read_text().count('WARNING: mismatch in integrated times')
    checks[str(native_log.relative_to(ROOT))]=sha(native_log)
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='passed_numerical_provenance_not_posterior_qualification',
        attempts_checked=896,statuses=statuses,identities_checked=len(checks),
        original_centre_unchanged=True,first_six_covariance_and_crosscovariance_unchanged=True,
        conditional_w_wa_scale_factor=3,scientific_prior_changed=False,
        source_rows_in_fit=len(source_points),holdout_rows_in_fit=0,independent_QR_rank=int(rank),
        independent_QR_outputs_checked=len(outputs),maximum_scaled_coefficient_difference=coefficient_difference,
        maximum_scaled_fitted_spectrum_difference=fitted_difference,
        theoretical_spectrum_shapes=sorted(set(metadata)),
        native_integrated_time_warning_lines=warning_count,
        native_warning_scope='Combined native worker log is retained and hashed; warnings cannot be assigned reliably to individual points from interleaved buffered output. A finite native likelihood is not an independent integration-accuracy certificate.',
        limitations='Numerical support and data lineage checks do not establish posterior convergence, exact-importance overlap, accuracy away from tested points or empirical cosmological model adequacy.',
        dependencies_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),design_path,original_path,
            RESULT/'broad-spectral-training.json',RESULT/'broad-spectral-fit.json',RESULT/'broad-spectral-holdout.json',modern_acquisition]},
        checked_sha256=checks)
    (RESULT/'broad-spectral-validation.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='checked_sha256'},indent=2))


if __name__=='__main__':main()
