"""One untouched broad-design holdout evaluation; never training feedback."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[3]
WORK=ROOT/'.work/unified-cosmology/inference/broad-spectral'
RESULT=ROOT/'studies/unified_cosmology/results/inference'
sys.path.insert(0,str(HERE.parent/'external_probes'))
from adapter import reference_point
from modern_adapter import modern_info
from fast_lensing import use_fast_lensing
from spectral_surrogate import replace_theory
from cobaya.model import get_model


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--model-file',type=Path)
    args=parser.parse_args()
    design_path=HERE/'broad-spectral-design.json';design=json.loads(design_path.read_text())
    fit_path=RESULT/'broad-spectral-fit.json';fit=json.loads(fit_path.read_text())
    model_path=args.model_file or ROOT/design['fit']['model_file']
    assert sha(model_path)==fit['model_sha256'] and fit['holdout_rows_used']==0
    acquisition=json.loads((RESULT/'broad-spectral-training.json').read_text())
    assert acquisition['status']=='complete_exact_acquisition'
    sources={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'spectral_surrogate.py',
        HERE.parent/'external_probes/modern_adapter.py',HERE.parent/'external_probes/fast_lensing.py',
        design_path,fit_path,model_path]}
    rows=[]
    info=use_fast_lensing(replace_theory(modern_info(),model_path))
    with get_model(info) as model:
        reference=reference_point(model);theory=model.theory['spectral_surrogate']
        for j in range(design['counts']['holdout']):
            manifest=WORK/'holdout'/f'{j:04d}.json';meta=json.loads(manifest.read_text())
            identity=next(r for r in acquisition['rows'] if r['split']=='holdout' and r['index']==j)
            assert sha(manifest)==identity['manifest_sha256'];sources[str(manifest.relative_to(ROOT))]=sha(manifest)
            assert meta['index']==j and meta['split']=='holdout' and meta['design_sha256']==sha(design_path)
            if meta['status']!='finite_exact':
                rows.append(dict(index=j,status='not_evaluated_native_'+meta['status'],native_status=meta['status']));continue
            path=ROOT/meta['file'];assert sha(path)==meta['sha256'];sources[str(path.relative_to(ROOT))]=sha(path)
            with np.load(path,allow_pickle=False) as data:
                previous=theory.exact_calls;point=dict(reference,**meta['physical_point'])
                started=time.monotonic();proposal=model.logposterior(point);elapsed=time.monotonic()-started
                native=dict(zip(meta['likelihood_names'],data['loglikes']))
                candidate=dict(zip(model.likelihood,proposal.loglikes))
                assert set(candidate)==set(native)
                if not np.isfinite(proposal.logpost):
                    rows.append(dict(index=j,status='nonfinite_proposal',seconds=elapsed));continue
                delta={k:float(candidate[k]-v) for k,v in native.items()}
                rows.append(dict(index=j,status='compared',delta_loglike_by_component=delta,
                    delta_loglike_total=sum(delta.values()),exact_fallback=theory.exact_calls>previous,seconds=elapsed))
    good=[r for r in rows if r['status']=='compared'];differences=np.array([r['delta_loglike_total'] for r in good])
    assert np.isfinite(differences).all()
    record=dict(created_utc=datetime.now(timezone.utc).isoformat(),
        status='checked_requires_exact_posterior_correction' if not any(r['status']=='nonfinite_proposal' for r in rows) else 'nonfinite_proposal_failures',
        requested_holdout=len(rows),finite_comparisons=len(good),model_sha256=sha(model_path),
        delta_loglike_quantiles=np.quantile(differences,[0,.025,.16,.5,.84,.975,1]).tolist() if len(good) else None,
        delta_loglike_rms=float(np.sqrt(np.mean(differences**2))) if len(good) else None,
        exact_fallback_count=sum(r['exact_fallback'] for r in good),rows=rows,source_sha256=sources,
        scope='Untouched broad numerical holdout, nuisance values fixed at the native acquisition reference. No hyperparameter/model selection uses these outcomes; no cosmological posterior qualification or prior inference follows.')
    destination=RESULT/'broad-spectral-holdout.json'
    assert not destination.exists(),'Preserve the first untouched holdout; use a new explicitly declared version for any repeat.'
    destination.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:v for k,v in record.items() if k not in ['source_sha256','rows']},indent=2))


if __name__=='__main__':main()
