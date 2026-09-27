"""Post-holdout diagnosis, explicitly exploratory and not fresh validation."""
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import sys
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


def summary(rows,key='delta_loglike'):
    values=np.array([r[key] for r in rows])
    return dict(n=len(rows),rms=float(np.sqrt(np.mean(values**2))) if len(rows) else None,
        quantiles=np.quantile(values,[0,.16,.5,.84,1]).tolist() if len(rows) else None)


def main():
    design_path=HERE/'broad-spectral-design.json';design=json.loads(design_path.read_text())
    old_design_path=HERE.parent/'external_probes/spectral-training-design.json'
    old=json.loads(old_design_path.read_text());old_chol=np.array(old['coordinate_cholesky'])
    model_path=ROOT/design['fit']['model_file']
    broad_path=RESULT/'broad-spectral-holdout.json';broad=json.loads(broad_path.read_text())
    old_check_path=RESULT/'surrogate-holdout.json';old_check=json.loads(old_check_path.read_text())
    dependencies={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),design_path,old_design_path,
        broad_path,old_check_path,model_path,HERE/'spectral_surrogate.py',
        HERE.parent/'external_probes/modern_adapter.py',HERE.parent/'external_probes/fast_lensing.py']}
    rows=[]
    for check in broad['rows']:
        if check['status']!='compared':continue
        manifest=WORK/'holdout'/f"{check['index']:04d}.json"
        meta=json.loads(manifest.read_text());point=np.array(meta['actual_coordinates'])
        dependencies[str(manifest.relative_to(ROOT))]=sha(manifest)
        rows.append(dict(index=check['index'],delta_loglike=check['delta_loglike_total'],
            native_fixed_nuisance_CMB_BAO_loglike=sum(meta['loglikes']),
            old_envelope_coordinate_max=float(np.max(abs(np.linalg.solve(old_chol,point-np.array(old['centre']))))),
            w=float(point[-2]),wa=float(point[-1]),w_plus_wa=float(point[-2]+point[-1]),
            exact_fallback=check['exact_fallback'],component_delta=check['delta_loglike_by_component']))
    best=max(r['native_fixed_nuisance_CMB_BAO_loglike'] for r in rows)
    for r in rows:r['native_loglike_minus_best_observed']=r['native_fixed_nuisance_CMB_BAO_loglike']-best
    groups={'all':rows,'old_envelope_max4':[r for r in rows if r['old_envelope_coordinate_max']<=4],
        'old_envelope_max2':[r for r in rows if r['old_envelope_coordinate_max']<=2],
        'native_loglike_within10':[r for r in rows if r['native_loglike_minus_best_observed']>=-10],
        'native_loglike_within30':[r for r in rows if r['native_loglike_minus_best_observed']>=-30],
        'w_plus_wa_greater_minus02':[r for r in rows if r['w_plus_wa']>-.2],
        'w_plus_wa_atmost_minus02':[r for r in rows if r['w_plus_wa']<=-.2]}
    # The original holdout was already exposed in earlier model development.
    # It is a common-point diagnostic only, never a new untouched holdout.
    old_rows={r['index']:r for r in old_check['rows']};common=[]
    folder=ROOT/'.work/unified-cosmology/external-probes/spectral-training/holdout'
    with get_model(use_fast_lensing(replace_theory(modern_info(),model_path))) as model:
        reference=reference_point(model);theory=model.theory['spectral_surrogate']
        for path in sorted(folder.glob('*.npz')):
            meta=json.loads(path.with_suffix('.json').read_text());index=meta['index']
            assert sha(path)==meta['sha256']
            dependencies[str(path.relative_to(ROOT))]=sha(path)
            dependencies[str(path.with_suffix('.json').relative_to(ROOT))]=sha(path.with_suffix('.json'))
            before=theory.exact_calls
            result=model.logposterior(dict(reference,**meta['physical_point']))
            assert np.isfinite(result.logpost)
            delta=float(np.sum(result.loglikes)-sum(meta['loglikes']))
            common.append(dict(index=index,old_delta_loglike=old_rows[index]['delta_loglike_total'],
                broad_delta_loglike=delta,old_fallback=old_rows[index]['exact_fallback'],
                broad_fallback=theory.exact_calls>before,
                native_fixed_nuisance_CMB_BAO_loglike=sum(meta['loglikes'])))
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='exploratory_diagnosis_complete',
        selected_after='Broad cubic heldout RMS17.26/max139.90 was inspected; cuts and common-point comparison were then requested. These diagnostics are not preregistered accuracy gates.',
        broad_holdout_strata={name:summary(group) for name,group in groups.items()},
        broad_holdout_rows=rows,
        common_original_holdout=dict(old=summary(common,'old_delta_loglike'),broad=summary(common,'broad_delta_loglike'),
            old_fallbacks=sum(r['old_fallback'] for r in common),broad_fallbacks=sum(r['broad_fallback'] for r in common),rows=common),
        interpretation='Fixed-reference-nuisance CMB+BAO likelihood errors, not complete SN posterior errors or an importance-ESS forecast. Large errors in low-native-likelihood regions remain real numerical defects; stratification is not permission to silently remove those points or narrow scientific priors.',
        dependencies_sha256=dependencies)
    (RESULT/'broad-spectral-diagnosis.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['dependencies_sha256','broad_holdout_rows','common_original_holdout']},indent=2))
    print(json.dumps({k:v for k,v in result['common_original_holdout'].items() if k!='rows'},indent=2))


if __name__=='__main__':main()
