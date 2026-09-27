"""Compare untouched exact-CAMB validation points with a frozen spectrum proposal."""
import argparse
import hashlib
import json
import sys
import time
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0,str(HERE.parent/'external_probes'))
from adapter import reference_point
from modern_adapter import modern_info
from spectral_surrogate import replace_theory
from cobaya.model import get_model


def main():
    p = argparse.ArgumentParser()
    p.add_argument('model_file',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--fast-lensing',action='store_true')
    a = p.parse_args()
    folder = ROOT/'.work/unified-cosmology/external-probes/spectral-training/holdout'
    paths = sorted(folder.glob('*.npz'))
    assert len(paths)>=64,'Do not assess final proposal accuracy from a tiny holdout.'
    rows = []
    info=replace_theory(modern_info(),a.model_file)
    if a.fast_lensing:
        from fast_lensing import use_fast_lensing
        info=use_fast_lensing(info)
    with get_model(info) as model:
        reference = reference_point(model)
        theory = model.theory['spectral_surrogate']
        for path in paths:
            meta = json.loads(path.with_suffix('.json').read_text())
            with np.load(path,allow_pickle=False) as data:
                point = dict(reference,**meta['physical_point'])
                previous = theory.exact_calls
                started = time.monotonic()
                result = model.logposterior(point)
                elapsed = time.monotonic()-started
                exact = dict(zip(meta['likelihood_names'],data['loglikes']))
                candidate = dict(zip(model.likelihood,result.loglikes))
                delta = {k:float(candidate[k]-v) for k,v in exact.items()}
                rows.append({'index':meta['index'],'delta_loglike_by_component':delta,
                             'delta_loglike_total':sum(delta.values()),
                             'exact_fallback':theory.exact_calls>previous,'seconds':elapsed,
                             'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest()})
    differences = np.array([r['delta_loglike_total'] for r in rows])
    assert np.isfinite(differences).all()
    out = {'status':'checked_requires_exact_posterior_correction','holdout_rows':len(rows),
           'exact_lensing_reassociation':a.fast_lensing,
           'delta_loglike_quantiles':np.quantile(differences,[0,.025,.16,.5,.84,.975,1]).tolist(),
           'delta_loglike_rms':float(np.sqrt(np.mean(differences**2))),
           'exact_fallback_count':sum(r['exact_fallback'] for r in rows),
           'rows':rows,'model_sha256':hashlib.sha256(a.model_file.read_bytes()).hexdigest(),
           'code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'spectral_surrogate.py']},
           'scope':'No model fitting or tuning uses these holdout rows. Accuracy at fixed nuisance reference values is not final posterior qualification.'}
    a.output.parent.mkdir(parents=True,exist_ok=True)
    if a.fast_lensing:
        path=HERE.parent/'external_probes/fast_lensing.py'
        out['code_sha256'][path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='rows'},indent=2))


if __name__=='__main__':main()
