"""Read-only diagnosis of already completed initialization branches; no native calls."""
from pathlib import Path
import json,hashlib,csv,re
Q=Path(__file__).resolve().parent;P=Q.parent.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
paths={n:P/'fits-restricted'/n for n in ['joint12','joint9','joint_minus']}
def parse(p):
 out={};csp=[];adjust=[];support={};science=[]
 for line in (p/'native.log').read_text().splitlines():
  v=line.split()
  if not v:continue
  if v[0]=='PROSP_PEAK_BOUNDS':out.setdefault(v[1],{}).setdefault('bounds',[]).append(v)
  if v[0]=='CSP_ENTRY:':out.setdefault(v[1],{}).setdefault('entries',[]).append(v)
  if v[0]=='CSP_OBJECTIVE:':out.setdefault(v[1],{}).setdefault('objectives',[]).append(v)
  if v[0].startswith('CSP_'):csp.append(v)
  if v[0]=='PROSP_SUPPORT':support[v[1]]=v
  if 'Adjust PKMJD=' in line:adjust.append(line.strip())
 for line in (p/'fit.FITRES.TEXT').read_text().splitlines():
  if line.startswith(('SN:','VARNAMES:')):science.append(line)
 return out,csp,adjust,support,science
p={n:parse(d) for n,d in paths.items()};rows=[]
for cid in sorted(p['joint12'][0],key=int):
 nom=p['joint12'][0][cid];alt=p['joint_minus'][0][cid];nine=p['joint9'][0][cid]
 assert all(int(x['bounds'][0][2])==int(x['entries'][0][2])==1 for x in [nom,alt,nine])
 first=lambda b:float(b['bounds'][0][16]);entry=lambda b:float(b['entries'][0][9]);Dentry=lambda b:float(b['entries'][0][6]);final=lambda b:b['objectives'][-1]
 a,b,c=final(nom),final(alt),final(nine)
 # Every saved callback explicitly prints the native prior center.
 prior_equal=all(x[12]==y[12] for x,y in zip(nom['objectives'],alt['objectives']))
 rows.append(dict(CID=cid,nominal_pregrid_peak=first(nom),minus_pregrid_peak=first(alt),pregrid_difference=first(alt)-first(nom),nominal_optimizer_entry=entry(nom),minus_optimizer_entry=entry(alt),optimizer_entry_difference=entry(alt)-entry(nom),minus_native_grid_adjustment=entry(alt)-first(alt),nominal_Dentry=Dentry(nom),minus_Dentry=Dentry(alt),Dentry_difference=Dentry(alt)-Dentry(nom),all_saved_prior_centers_exact=prior_equal,final_D_difference=float(b[8])-float(a[8]),final_peak_difference=float(b[11])-float(a[11]),nine_vs_twelve_final_D=float(c[8])-float(a[8]),nine_vs_twelve_final_peak=float(c[11])-float(a[11]),nominal_mean_calls=int(p['joint12'][3][cid][2]),minus_mean_calls=int(p['joint_minus'][3][cid][2]),nominal_unsupported_calls=int(p['joint12'][3][cid][3]),minus_unsupported_calls=int(p['joint_minus'][3][cid][3])))
assert len(rows)==8 and all(r['pregrid_difference']==-2 and r['optimizer_entry_difference']==0 for r in rows)
assert all(r['all_saved_prior_centers_exact'] for r in rows)
with (Q/'per-object.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
result={'status':'Initializer intervention applied; intended optimizer-entry displacement NOT achieved. Frozen gate correctly fails; no tolerance change.','native_calls':0,'frozen_inputs_modified':False,'objects':8,'all_pregrid_differences_days':sorted(set(x['pregrid_difference'] for x in rows)),'all_optimizer_entry_differences_days':sorted(set(x['optimizer_entry_difference'] for x in rows)),'all_minus_native_grid_adjustments_days':sorted(set(x['minus_native_grid_adjustment'] for x in rows)),'all_CSP_tokens_nominal_vs_minus_exact':p['joint12'][1]==p['joint_minus'][1],'all_FITRES_science_nominal_vs_minus_exact':p['joint12'][4]==p['joint_minus'][4],'nominal_CSP_record_count':len(p['joint12'][1]),'minus_CSP_record_count':len(p['joint_minus'][1]),'all_saved_prior_centers_exact':all(r['all_saved_prior_centers_exact'] for r in rows),'maximum_absolute_final_D_difference':max(abs(r['final_D_difference']) for r in rows),'maximum_absolute_final_peak_difference':max(abs(r['final_peak_difference']) for r in rows),'nominal_grid_MJD':[57703.80078125,57705.80078125,57707.80078125,57709.80078125,57711.80078125],'minus_grid_MJD':[57701.80078125,57703.80078125,57705.80078125,57707.80078125,57709.80078125],'grid_winner_all_objects_both_cases':57707.80078125,'nominal_adjust_prints':p['joint12'][2],'minus_adjust_prints':p['joint_minus'][2],'existing_new_estimator_spent_seconds':sum(x['wall_seconds'] for x in json.loads((Q.parent/'native-activity.json').read_text()))}
(Q/'result.json').write_text(json.dumps(result,indent=2)+'\n')
files=[Q/'diagnose.py',Q/'per-object.csv',Q/'result.json',Q.parent/'noisy-failure.json',Q.parent/'native-activity.json',Q.parent/'execution-protocol.json',Q.parent/'execution-freeze.json',Q.parent/'run_restricted.py',Q.parent/'build/src/snlc_fit.car',Q.parent/'build/src/snana.car']
for d in paths.values():files += [d/x for x in ['native.log','fit.FITRES.TEXT','fit.nml','execution.json']]
(Q/'evidence-manifest.json').write_text(json.dumps({'files':{str(p.relative_to(R)):sha(p) for p in files}},indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if not k.endswith('_prints')},indent=2))
