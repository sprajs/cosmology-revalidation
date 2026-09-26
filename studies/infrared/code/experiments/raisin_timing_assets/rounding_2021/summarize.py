from pathlib import Path
import json,hashlib,csv
P=Path(__file__).resolve().parent
load=lambda p:json.loads(p.read_text())
e=load(P/'endpoint-result.json');c=load(P/'corner-result.json');allc=e['cases']+c['cases'];bad=[x for x in allc if not all(x[k] for k in ['distance_cap_pass','data_Q_cap_pass','native_gate','clone_gate'])]
rem=[]
for x in c['cases']:
 if x['kind']!='directed_corner':continue
 actual=x['delta_D'] if x['statistic']=='D' else x['delta_data_Q'];rem.append({'CID':x['original_CID'],'statistic':x['statistic'],'direction':x['direction'],'linear_prediction':x['linear_prediction'],'actual_response':actual,'nonlinear_remainder':actual-x['linear_prediction']})
source=[x for x in c['cases'] if x['kind']=='candidate_SIMLIB_times'];states=e['native_state']+c['native_state']
result={'status':'STOP_precision_Q_gate_failed','endpoint_pass':e['all_gates_pass'],'corner_pass':c['all_gates_pass'],'native_cases':len(allc),'iteration_blocks':len(states),'all_native_mask_count_fixed_coordinate_gates_pass':all(x['native_gate'] for x in allc),'all_distance_caps_pass':all(x['distance_cap_pass'] for x in allc),'max_abs_distance_response_mag':max(abs(x['delta_D']) for x in allc),'max_abs_dataQ_response':max(abs(x['delta_data_Q']) for x in allc),'maximum_objective_reconstruction_error':max(abs(x['objective_closure']) for x in states),'maximum_nonlinear_remainder':{s:max(abs(x['nonlinear_remainder']) for x in rem if x['statistic']==s) for s in ['D','Q']},'failing_cases':[{k:x.get(k) for k in ['new_CID','original_CID','kind','statistic','direction','delta_D','delta_data_Q','data_Q_cap','native_gate']} for x in bad],'source_time_branch':{'all_original_caps_pass':all(x['distance_cap_pass'] and x['data_Q_cap_pass'] for x in source),'max_abs_distance_response_mag':max(abs(x['delta_D']) for x in source),'max_abs_dataQ_response':max(abs(x['delta_data_Q']) for x in source),'scope':'Unique source-compatible timestamps, not selected as a replacement; no original flags/full-precision flux recovery.'},'native_seconds':sum(load(x)['wall_seconds'] for x in (P/'fits').glob('*/execution.json')),'paired_timing_executed':False,'scope':'Finite deterministic precision probes, not a rigorous global bound or posterior distribution. Fixed original gate preserved; no altered cohort or tolerance.'}
(P/'result.json').write_text(json.dumps(result,indent=2)+'\n');(P/'nonlinear-remainders.json').write_text(json.dumps(rem,indent=2)+'\n')
with (P/'response-ledger.csv').open('w',newline='') as f:
 keys=['new_CID','original_CID','kind','statistic','direction','coordinate_index','side','delta_D','delta_data_Q','data_Q_cap','distance_cap_pass','data_Q_cap_pass','native_gate','clone_gate','native_NFITDATA','FITRES_NDOF','fixed_gate'];w=csv.DictWriter(f,fieldnames=keys,extrasaction='ignore');w.writeheader();w.writerows(allc)
print(json.dumps(result,indent=2))
