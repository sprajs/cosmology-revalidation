"""Frozen arithmetic and convergence audit; no fitting or outcome-dependent choices."""
from pathlib import Path
import json, hashlib, importlib.util, csv
from collections import defaultdict
import numpy as np
from astropy.io import fits
O=Path(__file__).resolve().parent;ROOT=Path.cwd();V=O.parent/'convergence-diagnostic'
def module(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
A=module('nominal_audit',V/'check.py');R=module('response_runner',O/'run_v2.py')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def dump(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def lines(path,tag,cid):return [l for l in path.read_text().splitlines() if l.startswith(tag+' ') and l.split()[1]==cid]
def main():
 p=json.loads((O/'protocol.json').read_text());freeze=json.loads((O/'postprocessor-freeze.json').read_text())
 for n,h in {**p['inputs_sha256'],**freeze['files_sha256']}.items():assert sha(ROOT/n)==h,n
 runs={j['name']:A.groups(A.parse(O/'fits'/j['name']/'fit.log')) for j in p['jobs']}
 nominal={n:A.groups(A.parse(V/'fits'/n/'fit.log')) for n in ['iter12_default','iter12_minus','iter12_plus']}
 rawmap=defaultdict(list)
 for r in json.loads((O/'input-map.json').read_text()):rawmap[(r['CID'],tuple(r['key_MJD_dataF_dataE_band']))].append(r['new_band'])
 with fits.open(ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.fits') as h:
  t=h['FilterTrans'].data;w=t['wavelength (A)'];bands={r['band'] for data in runs.values() for bs in data.values() for b in bs for r in b['rows']};bounds={c:(float(w[t['CSP-'+c]!=0].min()),float(w[t['CSP-'+c]!=0].max())) for c in bands}
 for j in p['jobs']:assert set(runs[j['name']])==set(p['membership']) and all(len(bs)==j['iterations'] for bs in runs[j['name']].values())
 histories={};cases=[];all_ar=[];all_maps=[];all_steps=[];controls=[];structural=[]
 for cid in p['membership']:
  case={'CID':cid,'metadata_affected':cid in p['affected_32'],'runs':{}}
  for j in p['jobs']:
   name=j['name'];bs=runs[name][cid];refname='iter12_minus' if name.endswith('minus') else 'iter12_plus' if name.endswith('plus') else 'iter12_default';ref=nominal[refname][cid]
   ar=[A.arithmetic(b) for b in bs];all_ar.extend(ar);steps=[A.step(a,b) for a,b in zip(bs[1:-1],bs[2:])];all_steps.extend(steps[-2:]);cm=[]
   for i in range(2,len(bs)):
    a,b,c=bs[i-2:i+1]
    if A.ident(b,c):
     s=np.exp(-A.K*(b['objective']['D']-a['objective']['D']));pred=A.E(c)+(b['C']-A.E(b))*s*s;cm.append(float(np.linalg.norm(c['C']-pred)/np.linalg.norm(c['C'])))
    else:cm.append(None)
   all_maps.extend(v for v in cm if v is not None)
   physical=[R.physical(ref[i],b,rawmap) if j['arm']=='known_Jdw_to_j' else A.ident(ref[i],b) for i,b in enumerate(bs)]
   control_exact=all(R.exact(ref[i],b) for i,b in enumerate(bs)) if cid in p['controls_10'] or name=='nominal12_copy' else True
   row_support=[]
   for b in bs:
    row_support.append(all(bounds[r['band']][0]/(1+r['z'])>=r['rest_lambda_fit_min'] and bounds[r['band']][1]/(1+r['z'])<=r['rest_lambda_fit_max'] for r in b['rows']))
   gates={'physical_rows':all(physical),'control_states_exact':control_exact,'last_two_D':all(abs(s['delta_D'])<=.001 for s in steps[-2:]),'last_two_C':all(s['C_whitened_operator_norm']<=.001 for s in steps[-2:]),'same_final_mask':all(s['identity'] for s in steps),'final_frozen_C_D':ar[-1]['frozen_C_delta_D'] is not None and abs(ar[-1]['frozen_C_delta_D'])<=.001,'positive_means_covariance':all(a['positive_means'] and a['C_min_eigenvalue']>0 for a in ar),'all_objectives_close':all(abs(a['Q_closure'])<=1e-8 for a in ar),'Cmap':all(v is not None and v<=5e-6 for v in cm),'final_support':all(-20<=a['phase_min'] and a['phase_max']<=70 and .7<=a['shape']<=1.3 for a in ar[1:]),'throughput_support':all(row_support)}
   structural.extend(physical);structural.append(control_exact)
   science=True;lc=True
   if cid in p['controls_10'] or name=='nominal12_copy':
    science=lines(O/'fits'/name/'fit.FITRES.TEXT','SN:',cid)==lines(V/'fits'/refname/'fit.FITRES.TEXT','SN:',cid)
    lc=lines(O/'fits'/name/'fit.LCPLOT.TEXT','OBS:',cid)==lines(V/'fits'/refname/'fit.LCPLOT.TEXT','OBS:',cid);controls.append({'CID':cid,'job':name,'states':control_exact,'FITRES':science,'LCPLOT':lc});gates.update(control_FITRES=science,control_LCPLOT=lc)
   case['runs'][name]={'final_D':ar[-1]['D'],'final_Q':ar[-1]['Q_native'],'last_two_steps':steps[-2:],'initial_empirical_grid_support':-20<=ar[0]['phase_min'] and ar[0]['phase_max']<=70,'gates':gates}
   histories[name+'__'+cid]={'arithmetic':ar,'steps':steps,'covariance_map_relative_errors':cm,'throughput_support_by_callback':row_support}
  nine=runs['changed09_default'][cid];twelve=runs['changed12_default'][cid];nom=nominal['iter12_default'][cid];vals=[runs[n][cid][-1]['objective']['D'] for n in ['changed12_default','changed12_minus','changed12_plus']];actual=[runs[n][cid][0]['entry']['D_entry']-twelve[0]['entry']['D_entry'] for n in ['changed12_minus','changed12_plus']]
  cross={'nine_to_twelve':abs(vals[0]-nine[-1]['objective']['D'])<=.001,'exact_nine_prefix':all(R.exact(a,b) for a,b in zip(nine,twelve[:9],strict=True)),'multistart':max(vals)-min(vals)<=.001,'actual_starts':bool(np.max(np.abs(np.array(actual)-[-.2,.2]))<=1e-12),'final_mask_across_starts':all(A.ident(runs[n][cid][-1],twelve[-1]) for n in ['changed12_minus','changed12_plus'])}
  delta12=vals[0]-nom[-1]['objective']['D'];delta3=twelve[2]['objective']['D']-nom[2]['objective']['D']
  case.update(nominal_D12=nom[-1]['objective']['D'],changed_D12=vals[0],delta_D12=delta12,delta_D3=delta3,secondary_minus_primary=delta3-delta12,changed_D12_minus_D9=vals[0]-nine[-1]['objective']['D'],final_D_multistart_range=max(vals)-min(vals),actual_first_entry_offsets=actual,cross_run_gates=cross)
  case['numeric_pass']=all(cross.values()) and all(all(v['gates'].values()) for n,v in case['runs'].items() if n!='changed09_default');cases.append(case)
 ledger=json.loads((O/'execution-ledger.json').read_text());recorded_structural=all(json.loads((O/'fits'/j['name']/'structural-gates.json').read_text())['pass'] for j in p['jobs']);prefix=json.loads((O/'fits/changed12_default/prefix-gate.json').read_text())['changed09_to_changed12_exact'];copy=json.loads((O/'fits/nominal12_copy/copy-gate.json').read_text())['science_and_LCPLOT_exact']
 passed=all(c['numeric_pass'] for c in cases) and recorded_structural and all(structural) and prefix and copy and len(ledger['jobs'])==5 and all(j['returncode']==0 for j in ledger['jobs']) and ledger['additional_seconds']<=120 and ledger['overall_seconds']<=600
 def aggregates(field):
  return {'mean_all42_mag':float(np.mean([c[field] for c in cases])),'mean_affected32_mag':float(np.mean([c[field] for c in cases if c['metadata_affected']])),'high_minus_low_change_mag':float(-np.mean([c[field] for c in cases])),'min_mag':min(c[field] for c in cases),'max_mag':max(c[field] for c in cases),'maximum_abs_control_mag':max(abs(c[field]) for c in cases if c['CID'] in p['controls_10'])}
 summary={'certified_all42_numeric_pass':bool(passed),'interpretation':'Conditional native current-input RC1-to-RC2 processing response; no BBC regeneration, physical WIRC correction or cosmology. Aggregates are descriptive and uncertified if any gate fails.','failed_CIDs':[c['CID'] for c in cases if not c['numeric_pass']],'prefix_exact':prefix,'nominal_copy_exact':copy,'all_controls_science_LCPLOT_states_exact':all(all(x[k] for k in ['states','FITRES','LCPLOT']) for x in controls),'max_abs_Q_closure':max(abs(a['Q_closure']) for a in all_ar),'max_Cmap_relative_error':max(all_maps),'max_last_two_D_increment':max(abs(s['delta_D']) for s in all_steps),'max_last_two_C_whitened_operator_norm':max(s['C_whitened_operator_norm'] for s in all_steps),'max_final_frozen_C_D_gap':max(abs(histories[n+'__'+c]['arithmetic'][-1]['frozen_C_delta_D']) for n in runs for c in p['membership']),'max_final_D_multistart_range':max(c['final_D_multistart_range'] for c in cases),'max_changed_D12_minus_D9':max(abs(c['changed_D12_minus_D9']) for c in cases),'initial_grid_unsupported':[c['CID'] for c in cases if not c['runs']['changed12_default']['initial_empirical_grid_support']],'primary':aggregates('delta_D12'),'secondary':aggregates('delta_D3'),'secondary_minus_primary':aggregates('secondary_minus_primary'),'execution':ledger}
 dump(O/'result.json',{**summary,'cases':cases});dump(O/'summary.json',summary);dump(O/'state-histories.json',histories);dump(O/'control-gates.json',controls)
 fields=['CID','metadata_affected','numeric_pass','nominal_D12','changed_D12','delta_D12','delta_D3','secondary_minus_primary','changed_D12_minus_D9','final_D_multistart_range']
 with (O/'responses.csv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(cases)
 print(json.dumps(summary,indent=2,allow_nan=False))
if __name__=='__main__':main()
