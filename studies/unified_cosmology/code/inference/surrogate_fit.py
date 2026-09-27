"""Fit a numerical spectrum proposal; independent holdout data are never fitted."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import numpy as np
from spectral_surrogate import polynomial, polynomial_exponents

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def fit(output,degree=3,limit=None):
    design_path = HERE.parent/'external_probes/spectral-training-design.json'
    design = json.loads(design_path.read_text())
    paths = sorted((ROOT/'.work/unified-cosmology/external-probes/spectral-training/train').glob('*.npz'))
    if limit: paths = paths[:limit]
    exponents = polynomial_exponents(degree=degree)
    assert len(paths)>=2*len(exponents),'Require at least twice as many training rows as coefficients.'
    centre = np.array(design['centre']);chol = np.array(design['coordinate_cholesky'])
    points,spectra = [],[]
    for path in paths:
        with np.load(path,allow_pickle=False) as row:
            points.append(row['coordinates'])
            spectra.append(row['spectra'])
    points = np.array(points);spectra = np.array(spectra)
    assert np.isfinite(spectra).all()
    length = spectra.shape[-1]
    x = np.linalg.solve(chol,(points-centre).T).T
    features = polynomial(x,exponents)
    target = spectra.reshape(len(paths),-1)
    # Scaling prevents tiny lensing numbers from being hidden by temperature
    # spectra in numerical diagnostics; each output is still ordinary least squares.
    scale = np.maximum(np.sqrt(np.mean(target**2,axis=0)),1e-30)
    u,s,vh = np.linalg.svd(features,full_matrices=False)
    assert s[-1]>1e-10*s[0], 'Ill-conditioned training features.'
    coefficients = ((vh.T/s)@u.T)@(target/scale)
    residual = features@coefficients-target/scale
    output.parent.mkdir(parents=True,exist_ok=True)
    if output.exists():raise FileExistsError('Use a new versioned model path; never overwrite an active model.')
    np.savez_compressed(output,centre=centre,coordinate_cholesky=chol,
                        exponents=exponents,coefficients=coefficients,
                        output_scale=scale,length=length)
    manifest = {'status':'fitted_not_posterior_qualified','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'degree':degree,'training_rows':len(paths),'coefficients_per_spectral_element':len(exponents),
        'feature_condition_number':float(s[0]/s[-1]),
        'training_scaled_rms':float(np.sqrt(np.mean(residual**2))),
        'training_scaled_max':float(np.max(abs(residual))),
        'model_sha256':hashlib.sha256(output.read_bytes()).hexdigest(),
        'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths+[design_path,Path(__file__),HERE/'spectral_surrogate.py']},
        'scope':'Numerical proposal only. Requires independent held-out likelihood checks and exact posterior correction. No holdout rows used in fit.'}
    output.with_suffix('.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k!='source_sha256'},indent=2))


if __name__=='__main__':
    p = argparse.ArgumentParser()
    p.add_argument('output',type=Path)
    p.add_argument('--degree',type=int,choices=[2,3],default=3)
    p.add_argument('--limit',type=int)
    a = p.parse_args();fit(a.output,a.degree,a.limit)
