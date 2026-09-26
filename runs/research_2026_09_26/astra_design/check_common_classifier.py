"""Independent common-clump/inference provenance and membership checks."""
from pathlib import Path
import hashlib,json,re
import numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[3];P=Path(__file__).resolve().parent;D=ROOT/'runs/research_2026_09_26/common_classifier_residual';C=ROOT/'phase2/classification/reconstruction_20260926';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
count=0
for man in [P/'common-classifier-manifest.json',D/'final_source_rerun/inference/manifest.json']:
 m=json.loads(man.read_text())
 for key in ['sha256','inputs_sha256','outputs_sha256']:
  for rel,h in m.get(key,{}).items():assert sha(ROOT/rel)==h,rel;count+=1
for name in ['simulation-probabilities.csv','real-probabilities.csv','peak-ledger.csv','single-batch-check.csv','inference-gate.json']:
 assert (D/'inference'/name).read_bytes()==(D/'final_source_rerun/inference'/name).read_bytes(),name
base=(C/'clump_run/clump.nml').read_text();binary=sha(C/'clump_run/snana.exe');clumps=[]
for arm in ['P21','G10']:
 d=D/'clump'/arm;p=json.loads((d/'prepared.json').read_text());g=json.loads((d/'run-gate.json').read_text());assert p['binary_sha256']==binary and p['base_nml_sha256']==sha(C/'clump_run/clump.nml');assert g['returncode']==0 and g['graceful']
 x=(d/'clump.nml').read_text();x=re.sub(r"PRIVATE_DATA_PATH\s*=.*",next(l for l in base.splitlines() if 'PRIVATE_DATA_PATH' in l).strip(),x);x=x.replace(f"'PH2_pilot02_{arm}'","'DES-SN5YR_DES'");assert x==base
 assert sha(d/'clump.nml')==p['nml_sha256'];assert sha(d/f'PH2_pilot02_{arm}.SNANA.TEXT')==g['output_sha256'];clumps.append({'arm':arm,'only_declared_three_fields_changed':True,'same_binary':True})
sim=pd.read_csv(D/'inference/simulation-probabilities.csv',dtype={'CID':str});real=pd.read_csv(D/'inference/real-probabilities.csv',dtype={'CID':str});checks=pd.read_csv(D/'inference/single-batch-check.csv');assert len(sim)==512 and len(real)==1020 and len(checks)==10;assert checks.abs_delta.max()==0
ref=pd.read_csv(C/'des_clump_diagnostic.csv',dtype={'CID':str}).set_index('CID');delta=real.pIa.to_numpy()-ref.loc[real.CID,'pIa'].to_numpy();assert abs(delta).max()<1e-6
old=pd.read_csv(ROOT/'runs/research_2026_09_26/classifier_residual_sensitivity/membership_and_scores.csv',dtype={'CID':str});assert set(real.loc[real.pIa>.999,'CID'])==set(old.loc[old.retained_both_gt_0p999,'CID'])
for arm in ['P21','G10']:
 ids=set((P/f'simulation_design/{arm}-cids.txt').read_text().split());q=sim[sim.arm==arm];assert set(q.CID)==ids and len(q)==256;assert np.array_equal(q.gt999,q.pIa>.999);assert (q.n_raw-q.n_window_excluded-q.n_flag_excluded==q.n_retained).all();assert q.n_retained.gt(0).all()
assert (real.n_raw-real.n_window_excluded-real.n_flag_excluded==real.n_retained).all()
paths=[]
for arm in ['P21','G10']:
 v=f'PH2_pilot02_{arm}';b=ROOT/'phase2/literature/simulations/outputs'/v;paths += [b/f'{v}_HEAD.FITS',b/f'{v}_PHOT.FITS']
report={'status':'PASS: frozen inputs/outputs, same binary/three-field clump change, measured window count closure and exact declared memberships','hash_checks':count,'clump_checks':clumps,'sim_counts':sim.groupby('arm').gt999.sum().to_dict(),'real_count':int((real.pIa>.999).sum()),'max_real_probability_repeat_difference':float(abs(delta).max()),'single_batch_difference':float(checks.abs_delta.max()),'near_threshold_count':int((abs(sim.pIa-.999)<=1e-6).sum()),'historical_runner_note':'Original inference runner gained score() after its manifest; final-source full rerun has byte-identical five principal inference outputs and closes current source hash; original manifest preserved.', 'scope':'Confirms executable common measured-input pipeline; original published-probability reproduction remains failed, no calibrated purity or truth peak substitution','additional_raw_inputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in paths},'checker_sha256':sha(__file__)};(P/'common-classifier-independent.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:v for k,v in report.items() if 'sha256' not in k},indent=2))
