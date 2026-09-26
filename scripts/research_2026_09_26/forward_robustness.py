"""Run frozen forward sensitivities and compare only common-CID outcomes."""
import argparse,hashlib,json,os,sys
from pathlib import Path
import numpy as np,pandas as pd

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'scripts/phase2/forward_discrimination'
sys.path.insert(0,str(SOURCE))
import compare

OUT=ROOT/'runs/research_2026_09_26/forward_robustness'
PREREG=ROOT/'phase2/forward_discrimination/preregistration.json'
PRIMARY=ROOT/'phase2/forward_discrimination/primary'
SETTINGS={
 'all-des':{'all_des':True},
 'recovered-mask':{'arm':'recovered-mask'},
 'strict-cuts':{'strict':True},
 'bandwidth075':{'bw':.075},
 'bandwidth150':{'bw':.15},
 'reference-qplus05':{'q':.5},
 'reference-qminus1':{'q':-1},
}
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def run_variant(name):
    if name not in SETTINGS:raise ValueError(name)
    dest=OUT/name
    if dest.exists():raise RuntimeError(f'Refusing to overwrite {dest}')
    OUT.mkdir(parents=True,exist_ok=True)
    compare.OUT=OUT
    summary=compare.run(name,nboot=300,**SETTINGS[name])
    files=set(compare.INPUTS)|{Path(__file__),SOURCE/'compare.py',PREREG,ROOT/'phase2/hierarchy/forward/manifest.json',ROOT/'phase2/forward_discrimination/manifest.json'}
    manifest={'variant':name,'configuration':summary['configuration'],
      'inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(files)},
      'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(dest.iterdir()) if p.is_file() and p.name!='manifest.json'}}
    (dest/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'variant':name,'cohort':summary['cohort'],'manifest':str((dest/'manifest.json').relative_to(ROOT))}),flush=True)

def verify_run(path):
    manifest=json.loads((path/'manifest.json').read_text())
    for rel,digest in {**manifest['inputs_sha256'],**manifest['outputs_sha256']}.items():
        if sha(ROOT/rel)!=digest:raise RuntimeError(f'Hash mismatch: {rel}')
    summary=json.loads((path/'summary.json').read_text())
    test=pd.read_csv(path/'real-cohort.csv',dtype={'CID':str})
    test=test.loc[test.fold==0]
    if len(test)!=summary['cohort']['n_test_common'] or test.CID.duplicated().any():
        raise RuntimeError(f'Invalid test membership: {path}')
    return summary,test

def comparisons():
    # Primary was saved before this continuation; do not alter it.
    original=json.loads((ROOT/'phase2/forward_discrimination/manifest.json').read_text())
    for rel,digest in original['outputs_sha256'].items():
        if rel.startswith('phase2/forward_discrimination/primary/') and sha(ROOT/rel)!=digest:
            raise RuntimeError(f'Primary changed: {rel}')
    base_summary,base_test=verify_run(PRIMARY) if (PRIMARY/'manifest.json').exists() else (json.loads((PRIMARY/'summary.json').read_text()),pd.read_csv(PRIMARY/'real-cohort.csv',dtype={'CID':str}).query('fold==0'))
    base_ids=set(base_test.CID)
    variant_data={}
    for name in SETTINGS:
        summary,test=verify_run(OUT/name)
        variant_data[name]=(summary,test)
    # Only after all manifests/cohorts pass do we read per-object scores.
    score_keys=('energy_cx','energy_mcx')
    models=('BS21','G10')
    rng=np.random.default_rng(20260926)
    outputs=[]
    for name,(summary,test) in variant_data.items():
        common=sorted(base_ids&set(test.CID));lost=sorted(base_ids-set(test.CID));gained=sorted(set(test.CID)-base_ids)
        if len(common)<10:raise RuntimeError(f'Too few common held-out CIDs in {name}')
        base_scores={m:pd.read_csv(PRIMARY/(m+'-test-scores.csv'),dtype={'CID':str}).set_index('CID').loc[common] for m in ('P21',*models)}
        variant_scores={m:pd.read_csv(OUT/name/(m+'-test-scores.csv'),dtype={'CID':str}).set_index('CID').loc[common] for m in ('P21',*models)}
        draws=rng.integers(0,len(common),size=(10000,len(common)))
        for model in models:
            for metric in score_keys:
                baseline=(base_scores[model][metric]-base_scores['P21'][metric]).to_numpy()
                alternative=(variant_scores[model][metric]-variant_scores['P21'][metric]).to_numpy()
                change=alternative-baseline
                vals=change[draws].mean(axis=1)
                within=summary['models'][model]['metrics'][metric]
                outputs.append({'variant':name,'model':model,'metric':metric,'baseline_n_test':len(base_ids),'variant_n_test':len(test),'common_n_test':len(common),'lost_baseline_n':len(lost),'gained_variant_n':len(gained),'variant_n_train':summary['cohort']['n_train_common'],'variant_n_fields':summary['cohort']['n_fields_test'],
                    'variant_full_cohort_delta_vs_P21':within['delta_vs_P21'],'variant_full_cohort_delta_ci_low':within['delta_vs_P21_ci95'][0],'variant_full_cohort_delta_ci_high':within['delta_vs_P21_ci95'][1],
                    'common_primary_contrast_mean':float(baseline.mean()),'common_variant_contrast_mean':float(alternative.mean()),'common_change_mean':float(change.mean()),'common_change_object_bootstrap_ci_low':float(np.quantile(vals,.025)),'common_change_object_bootstrap_ci_high':float(np.quantile(vals,.975))})
    score_path=OUT/'comparison';score_path.mkdir(exist_ok=False)
    pd.DataFrame(outputs).to_csv(score_path/'common_cid_comparisons.csv',index=False)
    support=[]
    for name,(_,test) in variant_data.items():
        support.append({'variant':name,'n_test':len(test),'n_common':len(base_ids&set(test.CID)),'lost_baseline_cids':sorted(base_ids-set(test.CID)),'gained_variant_cids':sorted(set(test.CID)-base_ids),'fields':sorted(test.field.unique().tolist())})
    (score_path/'support.json').write_text(json.dumps(support,indent=2)+'\n')
    manifest={'script_sha256':sha(Path(__file__)),'primary_manifest_sha256':sha(ROOT/'phase2/forward_discrimination/manifest.json'),
              'variant_manifest_sha256':{name:sha(OUT/name/'manifest.json') for name in SETTINGS},
              'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in score_path.iterdir() if p.is_file() and p.name!='manifest.json'}}
    (score_path/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'variants':len(variant_data),'comparisons':len(outputs),'support':[{k:v for k,v in item.items() if k in ('variant','n_test','n_common')} for item in support]}),flush=True)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--variant',choices=list(SETTINGS));ap.add_argument('--compare',action='store_true');a=ap.parse_args()
    if (a.variant is None)==(not a.compare):ap.error('Choose exactly one --variant or --compare')
    if a.variant:run_variant(a.variant)
    else:comparisons()
if __name__=='__main__':main()
