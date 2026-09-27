"""Fit the predeclared cubic numerical proposal using only new train rows."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import numpy as np
from spectral_surrogate import polynomial,polynomial_exponents

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
WORK=ROOT/'.work/unified-cosmology/inference/broad-spectral'
RESULT=ROOT/'studies/unified_cosmology/results/inference'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    design_path=HERE/'broad-spectral-design.json';design=json.loads(design_path.read_text())
    acquisition_path=RESULT/'broad-spectral-training.json'
    acquisition=json.loads(acquisition_path.read_text())
    assert acquisition['status']=='complete_exact_acquisition'
    assert acquisition['design_sha256']==sha(design_path)
    output=args.output or ROOT/design['fit']['model_file']
    assert not output.exists(),'Use a new versioned path; never overwrite any active model.'
    exponents=polynomial_exponents(degree=design['fit']['degree'])
    assert len(exponents)==design['fit']['features']==165
    centre=np.array(design['centre']);chol=np.array(design['coordinate_cholesky'])
    points,spectra,paths,manifests=[],[],[],[]
    for j in range(design['counts']['train']):
        manifest=WORK/'train'/f'{j:04d}.json';row=json.loads(manifest.read_text())
        assert row['split']=='train' and row['index']==j and row['design_sha256']==sha(design_path)
        identity=next(r for r in acquisition['rows'] if r['split']=='train' and r['index']==j)
        assert sha(manifest)==identity['manifest_sha256']
        manifests.append(manifest)
        if row['status']!='finite_exact':continue
        path=ROOT/row['file'];assert sha(path)==row['sha256']
        with np.load(path,allow_pickle=False) as data:
            points.append(data['coordinates']);spectra.append(data['spectra'])
        paths.append(path)
    assert len(paths)>=design['fit']['minimum_finite_train']>=2*len(exponents)
    points=np.array(points);spectra=np.array(spectra)
    assert np.isfinite(points).all() and np.isfinite(spectra).all()
    length=spectra.shape[-1]
    white=np.linalg.solve(chol,(points-centre).T).T
    features=polynomial(white,exponents);target=spectra.reshape(len(paths),-1)
    scale=np.maximum(np.sqrt(np.mean(target**2,axis=0)),1e-30)
    u,s,vh=np.linalg.svd(features,full_matrices=False)
    assert s[-1]/s[0]>design['fit']['feature_singular_value_ratio_min']
    coefficients=((vh.T/s)@u.T)@(target/scale)
    residual=features@coefficients-target/scale
    output.parent.mkdir(parents=True,exist_ok=True)
    np.savez_compressed(output,centre=centre,coordinate_cholesky=chol,
        exponents=exponents,coefficients=coefficients,output_scale=scale,length=length)
    record=dict(status='fitted_not_posterior_qualified',created_utc=datetime.now(timezone.utc).isoformat(),
        degree=3,training_rows=len(paths),training_attempts=len(manifests),holdout_rows_used=0,
        coefficients_per_spectral_element=len(exponents),feature_condition_number=float(s[0]/s[-1]),
        training_scaled_rms=float(np.sqrt(np.mean(residual**2))),training_scaled_max=float(np.max(abs(residual))),
        model_file=str(output.relative_to(ROOT)),model_sha256=sha(output),
        source_sha256={str(p.relative_to(ROOT)):sha(p) for p in paths+manifests+
            [design_path,Path(__file__),HERE/'spectral_surrogate.py',HERE/'broad_spectral_training.py',acquisition_path]},
        scope='Numerical proposal only. The scientific priors, target and native fallback remain unchanged. No holdout values were fitted or used to choose this cubic model. Independent heldout checks and exact posterior correction are still required.')
    output.with_suffix('.json').write_text(json.dumps(record,indent=2)+'\n')
    (RESULT/'broad-spectral-fit.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':main()
