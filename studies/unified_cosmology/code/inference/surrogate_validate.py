"""Close the independent theory-provider adapter against actual Cobaya/CAMB."""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0,str(HERE.parent/'external_probes'))
from modern_adapter import modern_info
from adapter import reference_point
from spectral_surrogate import replace_theory
from cobaya.model import get_model


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--model-file',type=Path)
    p.add_argument('--exact',action='store_true')
    a = p.parse_args()
    work = ROOT/'.work/unified-cosmology/inference/surrogate'
    work.mkdir(parents=True,exist_ok=True)
    pilot = ROOT/'.work/unified-cosmology/external-probes/spectral-training/train/0000.npz'
    design = json.loads((HERE.parent/'external_probes/spectral-training-design.json').read_text())
    if not a.model_file:
        # A placeholder spectrum is never used for scientific inference: this
        # validates only the exact-mode provider interface, with no fitted model.
        with np.load(pilot) as f:
            length = len(f['ell']);spectra = f['spectra']
        a.model_file = work/'interface-placeholder.npz'
        np.savez_compressed(a.model_file,centre=design['centre'],
            coordinate_cholesky=design['coordinate_cholesky'],exponents=np.zeros((1,8),int),
            coefficients=spectra.reshape(1,-1),output_scale=np.ones(spectra.size),length=length)
        assert a.exact,'Placeholder permitted only with --exact.'
    info = modern_info()
    with get_model(info) as native:
        point = reference_point(native)
        result = native.logposterior(point)
        exact_likes = dict(zip(native.likelihood,map(float,result.loglikes)))
        exact_cls = {k:v.copy() for k,v in native.provider.get_Cl(ell_factor=True).items()}
        exact_derived = dict(zip(native.parameterization.derived_params(),map(float,result.derived)))
    with get_model(replace_theory(modern_info(),a.model_file,exact=a.exact)) as candidate:
        tested = candidate.logposterior(point)
        candidate_likes = dict(zip(candidate.likelihood,map(float,tested.loglikes)))
        candidate_cls = candidate.provider.get_Cl(ell_factor=True)
        candidate_derived = dict(zip(candidate.parameterization.derived_params(),map(float,tested.derived)))
    delta = {k:candidate_likes[k]-value for k,value in exact_likes.items()}
    cls_errors = {k:float(np.max(abs(candidate_cls[k]-v))/max(np.max(abs(v)),1e-100))
                  for k,v in exact_cls.items() if k in ['tt','ee','bb','te','pp']}
    result = {'status':'passed' if max(abs(x) for x in delta.values())<1e-5 else 'failed',
              'mode':'exact_provider_interface' if a.exact else 'surrogate_reference_only',
              'loglike_differences':delta,'spectrum_peak_normalized_max_errors':cls_errors,
              'derived_differences':{k:candidate_derived[k]-v for k,v in exact_derived.items() if k in ['rdrag','omegam']},
              'code_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),HERE/'spectral_surrogate.py']},
              'scope':'One reference-point adapter closure only; not heldout surrogate accuracy or posterior qualification.'}
    (work/('interface-validation.json' if a.exact else 'reference-validation.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    assert result['status']=='passed'


if __name__=='__main__':
    main()
