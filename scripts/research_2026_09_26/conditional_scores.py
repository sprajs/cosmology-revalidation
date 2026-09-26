"""Frozen paired selected-cohort score comparisons; see conditional-repair.md."""
import hashlib,json
from pathlib import Path
import numpy as np,pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/research_2026_09_26/conditional'
DEST=OUT/'scoring';DEST.mkdir(parents=True,exist_ok=True)
MODELS=['none','stretch','colour','tripp','host','host_colour','broken_colour','evolution','flexible']
REPRESENTATIVE=['tripp','host','flexible']
VARIANTS=['recovered_mask','omit_1307748','covariance_1p20','eds_reference','desitter_reference']
RUNS=[f'main_{noise}_{model}' for noise in ('gaussian','student4') for model in MODELS]
RUNS += [f'sensitivity_{variant}_{model}' for variant in VARIANTS for model in REPRESENTATIVE]
cohort=json.loads((ROOT/'phase2/hierarchy/conditional-cohort.json').read_text())['cid_to_fold']
expected={cid for cid,fold in cohort.items() if fold==0}
metadata={};manifests={}
for run in RUNS:
    path=OUT/run
    gate=json.loads((path/'convergence.json').read_text())
    if not gate['valid_for_scoring'] or gate['n_chains']<4 or gate['max_r_hat']>1.01 or gate['min_n_eff']<400:
        raise RuntimeError(f'Invalid convergence: {run}')
    manifest=json.loads((path/'manifest.json').read_text())
    for relative,digest in {**manifest['inputs_sha256'],**manifest['outputs_sha256']}.items():
        if hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()!=digest:
            raise RuntimeError(f'Hash mismatch: {run} {relative}')
    table=pd.read_csv(path/'heldout.csv',dtype={'CID':str},usecols=['CID','field','fold','observed_mB'])
    if len(table)!=213 or table.CID.duplicated().any() or set(table.CID)!=expected or not (table.fold==0).all():
        raise RuntimeError(f'Invalid held-out cohort: {run}')
    metadata[run]=table.set_index('CID').sort_index()
    manifests[run]=hashlib.sha256((path/'manifest.json').read_bytes()).hexdigest()
base=metadata['main_gaussian_tripp']
ids=base.index
for run,table in metadata.items():
    if not table.index.equals(ids) or not table.field.equals(base.field):
        raise RuntimeError(f'Field or ID mismatch: {run}')
    target_ref=metadata['sensitivity_recovered_mask_tripp'] if run.startswith('sensitivity_recovered_mask_') else base
    if not np.allclose(table.observed_mB,target_ref.observed_mB,atol=1e-10,rtol=0):
        raise RuntimeError(f'Target mismatch within refit arm: {run}')
if base.field.nunique()!=10:raise RuntimeError('Expected ten DES fields')

# Only after every run passes the checks above do we read score values.
scores={};summaries={};batch={}
for run in RUNS:
    path=OUT/run
    table=pd.read_csv(path/'heldout.csv',dtype={'CID':str},usecols=['CID','log_predictive_density','log_score_batch_mcse']).set_index('CID').loc[ids]
    summary=json.loads((path/'summary.json').read_text())
    values=table.log_predictive_density.to_numpy()
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(table.log_score_batch_mcse)):
        raise RuntimeError(f'Nonfinite score or MCSE: {run}')
    if not np.isclose(values.sum(),summary['heldout_log_score_sum'],atol=1e-8,rtol=0) or not np.isclose(values.mean(),summary['heldout_log_score_mean'],atol=1e-10,rtol=0):
        raise RuntimeError(f'Summary score mismatch: {run}')
    arr=np.load(path/'predictive_batch_scores.npz')['batch_scores']
    if arr.shape!=(40,213) or not np.all(np.isfinite(arr)):
        raise RuntimeError(f'Invalid score batches: {run}')
    batch_mcse=float(arr.sum(axis=1).std(ddof=1)/np.sqrt(40))
    if not np.isclose(batch_mcse,summary['heldout_log_score_batch_mcse'],atol=1e-8,rtol=0):
        raise RuntimeError(f'Batch MCSE mismatch: {run}')
    scores[run]=values;summaries[run]=summary;batch[run]=batch_mcse

seed=20260926;rng=np.random.default_rng(seed);n=len(ids)
object_draws=rng.integers(0,n,size=(10000,n))
fields=base.field.to_numpy();unique_fields=np.sort(np.unique(fields));nfields=len(unique_fields)
field_draws=rng.integers(0,nfields,size=(10000,nfields))
contrasts=[];paired_rows=[]
def compare(kind,variant,label,terms):
    # Terms are signed score vectors; shared CIDs make every difference paired.
    weights={}
    for sign,run in terms:weights[run]=weights.get(run,0)+sign
    weights={run:weight for run,weight in weights.items() if weight}
    delta=sum((weight*scores[run] for run,weight in weights.items()),np.zeros(n))
    obj_means=delta[object_draws].mean(axis=1)
    field_sums=np.array([delta[fields==f].sum() for f in unique_fields]);field_counts=np.array([(fields==f).sum() for f in unique_fields])
    cluster_means=field_sums[field_draws].sum(axis=1)/field_counts[field_draws].sum(axis=1)
    leave_one=np.array([(delta.sum()-field_sums[j])/(n-field_counts[j]) for j in range(nfields)])
    # Independent-chain quadrature is a conservative approximation to the
    # MCSE of a sum of separately fit posterior-integrated score estimates.
    mcse=float(np.sqrt(sum((weight*batch[run])**2 for run,weight in weights.items())))
    paired_rows.extend({'kind':kind,'variant':variant,'contrast':label,'CID':cid,'field':field,'delta_log_score':float(value)}
                       for cid,field,value in zip(ids,fields,delta))
    contrasts.append({'kind':kind,'variant':variant,'contrast':label,'n_objects':n,'n_fields':nfields,
                      'sum_delta_log_score':float(delta.sum()),'mean_delta_log_score':float(delta.mean()),
                      'object_bootstrap_mean_ci_low':float(np.quantile(obj_means,.025)),'object_bootstrap_mean_ci_high':float(np.quantile(obj_means,.975)),
                      'field_bootstrap_mean_ci_low':float(np.quantile(cluster_means,.025)),'field_bootstrap_mean_ci_high':float(np.quantile(cluster_means,.975)),
                      'leave_one_field_mean_min':float(leave_one.min()),'leave_one_field_mean_max':float(leave_one.max()),
                      'contrast_sum_mcse_approx':mcse})
for noise in ('gaussian','student4'):
    tripp=f'main_{noise}_tripp'
    for model in MODELS:
        run=f'main_{noise}_{model}'
        compare('within_noise_vs_tripp',noise,model+' minus tripp',[(1,run),(-1,tripp)])
for model in MODELS:
    compare('student4_vs_gaussian','primary',model+' student4 minus gaussian',[(1,f'main_student4_{model}'),(-1,f'main_gaussian_{model}')])
for variant in VARIANTS:
    prefix=f'sensitivity_{variant}_'
    for model in ('host','flexible'):
        compare('within_sensitivity_vs_tripp',variant,model+' minus tripp',[(1,prefix+model),(-1,prefix+'tripp')])
    if variant=='recovered_mask':
        for model in ('host','flexible'):
            compare('sensitivity_contrast_change',variant,model+' vs tripp contrast minus primary contrast',[(1,prefix+model),(-1,prefix+'tripp'),(-1,f'main_gaussian_{model}'),(1,'main_gaussian_tripp')])
    else:
        for model in REPRESENTATIVE:
            compare('sensitivity_vs_primary',variant,model+' sensitivity minus primary',[(1,prefix+model),(-1,f'main_gaussian_{model}')])

frame=pd.DataFrame(contrasts)
frame.to_csv(DEST/'paired_comparisons.csv',index=False)
pd.DataFrame(paired_rows).to_csv(DEST/'paired_object_differences.csv',index=False)
absolute=pd.DataFrame([{'run':run,'sum_log_score':summaries[run]['heldout_log_score_sum'],'mean_log_score':summaries[run]['heldout_log_score_mean'],'sum_score_batch_mcse':batch[run]} for run in RUNS])
absolute.to_csv(DEST/'run_scores.csv',index=False)
verification={'cohort_n':n,'fields':unique_fields.tolist(),'primary_runs':18,'sensitivity_runs':15,'contrasts':len(contrasts),'bootstrap_seed':seed,'bootstrap_replicates':10000,'max_primary_sum_score_mcse':max(batch[r] for r in RUNS[:18]),'max_primary_object_score_mcse':max(float(pd.read_csv(OUT/r/'heldout.csv',usecols=['log_score_batch_mcse']).log_score_batch_mcse.max()) for r in RUNS[:18]),'input_manifest_sha256':manifests,'protocol':'docs/research-2026-09-26/conditional-repair.md'}
(DEST/'verification.json').write_text(json.dumps(verification,indent=2)+'\n')
outputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in DEST.iterdir() if p.name!='manifest.json'}
(DEST/'manifest.json').write_text(json.dumps({'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'outputs_sha256':outputs},indent=2)+'\n')
print(json.dumps({k:verification[k] for k in ('cohort_n','fields','primary_runs','sensitivity_runs','contrasts','max_primary_sum_score_mcse')},indent=2))
