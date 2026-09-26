"""Read-only occurrence-aware diagnosis of saved native final retries.
Does not alter any frozen checker, release, input, source, fit, or threshold.
"""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='1'
import json, hashlib, re, csv
from collections import Counter, defaultdict
import numpy as np
O=Path(__file__).resolve().parent;P=O.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova')
S=P.parent/'restricted-peak-engineering/start-only-hook/build/src';L=P/'fits/joint12/native.log'
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while x:=f.read(1<<20):h.update(x)
 return h.hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
rows=[];by=defaultdict(list);occ=Counter();bounds=defaultdict(list);summary={};support={};bad=[];b=None
retained={};maxq=0.;mineig=float('inf');allfinite=True

def finish():
 global b,maxq,mineig,allfinite
 if b is None:return
 assert len(b['rows'])==len(b['weights'])==117
 x=np.array([list(map(float,t[6:])) for t in b['rows']]);ww=np.array(b['weights']);W=ww if b['usecov'] else np.diag(ww[:,0]);C=np.linalg.inv(W)
 assert W.shape==(117,117)
 ob=b['objective'];r=x[:,4]-x[:,2];q=float(r@W@r);err=q+ob[1]+ob[2]-ob[0];eig=float(np.linalg.eigvalsh(C).min());amp=float(x[:,2]@W@x[:,4]/(x[:,2]@W@x[:,2]))
 assert np.isfinite(x).all() and np.isfinite(C).all() and eig>0 and abs(err)<1e-7
 maxq=max(maxq,abs(err));mineig=min(mineig,eig)
 d={k:v for k,v in b.items() if k not in ['rows','weights']};d.update(Q_closure=err,C_min_eigenvalue=eig,fixed_C_delta_D=float(-2.5*np.log10(amp)),phase_min=float(x[:,1].min()),phase_max=float(x[:,1].max()))
 rows.append(d);by[b['cid']].append(d)
 if b['iteration']>=11:retained[b['cid'],b['iteration'],b['occurrence']]={'x':x,'W':W,'C':C,'rawrows':b['rows'],'entry':b['entry'],'objective':ob}
 b=None

with L.open() as f:
 for lineno,l in enumerate(f,1):
  t=l.split()
  if not t:continue
  if t[0]=='CSP_ENTRY:':
   finish();key=t[1],int(t[2]);occ[key]+=1;b={'cid':key[0],'iteration':key[1],'occurrence':occ[key],'repeat':t[5]=='T','entry_line':lineno,'entry':list(map(float,t[6:])),'rows':[],'weights':[]}
  elif t[0]=='CSP_ROW:':assert b is not None;b['rows'].append(t)
  elif t[0]=='CSP_OBJECTIVE:':
   assert (t[1],int(t[2]))==(b['cid'],b['iteration']);b['objective']=list(map(float,t[5:]));b['objective_line']=lineno;b['usecov']=t[4]=='T'
  elif t[0] in ['CSP_WROW:','CSP_WDIAG:']:b['weights'].append(list(map(float,t[4:])))
  elif 'MNFIT_DRIVER: MIGRAD returns IERR' in l:b['MINIMIZE_return']=int(t[-1]);b['minimize_line']=lineno
  elif 'MNFIT_DRIVER: EXIT returns IERR' in l:b['EXIT_return']=int(t[-1])
  elif 'Bad COV matrix (MNSTAT_COV=' in l:
   m=re.search(r'MNSTAT_COV=\s*(\d+)\).*CID=(\S+)',l);bad.append({'cid':m[2],'MNSTAT_COV':int(m[1]),'line':lineno})
  elif t[0]=='PROSP_PEAK_BOUNDS':bounds[t[1]].append({'line':lineno,'tokens':t})
  elif t[0]=='PROSP_DOMAIN_SUMMARY':assert t[1] not in summary;summary[t[1]]=t
  elif t[0]=='PROSP_SUPPORT':assert t[1] not in support;support[t[1]]=t
finish()
assert len(rows)==770 and len(by)==64
comparisons=[]
for c,rr in by.items():
 iterations=[x['iteration'] for x in rr]
 expected=list(range(1,13))+([12] if c in ['22','64'] else [])
 assert iterations==expected
 assert [x['repeat'] for x in rr]==[False]*12+([True] if c in ['22','64'] else [])
 assert len(bounds[c])==len(rr)
 for bd,st in zip(bounds[c],rr):assert int(bd['tokens'][2])==st['iteration'] and bd['line']<st['entry_line']
 if c not in ['22','64']:continue
 a=retained[c,12,1];z=retained[c,12,2];prev=retained[c,11,1]
 d={'cid':c,'final12_repeated12_row_tokens_equal':a['rawrows']==z['rawrows'],'full_numeric_rows_exact':np.array_equal(a['x'],z['x']),'W_exact':np.array_equal(a['W'],z['W']),'C_exact':np.array_equal(a['C'],z['C']),'objective_values_exact':a['objective']==z['objective'],'entry_numeric_values_exact':a['entry']==z['entry'],'entry_flags':[False,True],'final_objective':z['objective'],'repeat_minus_original_D':z['objective'][3]-a['objective'][3],'repeat_minus_original_peak':z['objective'][6]-a['objective'][6],'repeat_minus_iter11_D':z['objective'][3]-prev['objective'][3],'repeat_minus_iter11_peak':z['objective'][6]-prev['objective'][6],'bounds_tokens_exact':bounds[c][-1]['tokens']==bounds[c][-2]['tokens'],'MINIMIZE_returns':[rr[-2]['MINIMIZE_return'],rr[-1]['MINIMIZE_return']]}
 ll=np.linalg.cholesky(prev['C']);dc=np.linalg.solve(ll,z['C']-prev['C']);dc=np.linalg.solve(ll,dc.T).T;d['C_change_from_iter11_whitened_norm']=float(np.linalg.norm(dc,2))
 comparisons.append(d)
 np.savez_compressed(O/f'CID{c}-repeated-final-states.npz',first_rows=a['x'],repeat_rows=z['x'],first_W=a['W'],repeat_W=z['W'],first_C=a['C'],repeat_C=z['C'])

schema=json.loads((P/'build-protocol.json').read_text())['schema'];summ={c:dict(zip(schema['PROSP_DOMAIN_SUMMARY'],t[1:])) for c,t in summary.items()}
assert set(summ)==set(by)==set(support)
for c,x in summ.items():
 assert int(x['bounds_assertion_count'])==len(by[c]);assert int(x['physical_mean_calls'])==int(support[c][2]);assert int(support[c][3])==int(support[c][4])==0
 assert float(x['min_phase'])>=-20 and float(x['max_phase'])<=70 and float(x['minimum_phase_edge_distance_rest_days'])>=0
 assert float(x['min_FCN_absolute_peak'])>=57671.342999365806 and float(x['max_FCN_absolute_peak'])<=57715.067000181196 and float(x['minimum_peak_edge_distance_observer_days'])>=0
save(O/'callback-ledger.json',rows);save(O/'domain-ledger.json',summ)
result={'diagnosis':'Native one-time final-iteration retry after MINUIT parameter-covariance status1; not a duplicate printer call. Frozen domain checker assumes exactly one native visit per nominal iteration.','saved_native_process_return':json.loads((P/'fits/joint12/execution.json').read_text())['returncode'],'fit_native_seconds':json.loads((P/'fits/joint12/execution.json').read_text())['wall_seconds'],'callbacks':len(rows),'objects':len(by),'callback_counts':dict(Counter(len(x) for x in by.values())),'MINIMIZE_return_counts':dict(Counter(x['MINIMIZE_return'] for x in rows)),'final12_MINIMIZE_return_counts':dict(Counter(x['MINIMIZE_return'] for x in rows if x['iteration']==12)),'bad_parameter_covariance_retries':bad,'repeated_states':comparisons,'max_abs_Q_closure_all770':maxq,'minimum_observation_C_eigenvalue_all770':mineig,'all_support_and_domain_call_counters_closed_to_occurrences':True,'physical_mean_calls':sum(int(x['physical_mean_calls']) for x in summ.values()),'final_parameter_covariance_status':'Not directly exported in this FITRES. Initial final12 status1 explicitly printed; repeated final12 MINIMIZE also returns4. No accurate final parameter-Hessian certification.','convergence_status':'Point-state arithmetic/support checks pass. Native MINIMIZE abnormal termination remains material; full predefined9-vs12 and true-start comparisons have not run. Exact repetition is not proof of peak stationarity or global optimality.','frozen_checker_and_failure_unchanged':True,'native_calls_in_this_review':0}
save(O/'result.json',result)
source_ranges={'snlc_fit.car':[(1848,1870),(2038,2072),(2318,2336),(7168,7177),(8085,8133)],'snana.car':[(9180,9272),(9363,9372),(36675,36778)],'minuit.F':[(3690,3707),(3778,3802),(3970,3984),(12522,12544)]}
with (O/'source-evidence.txt').open('w') as out:
 for name,rr in source_ranges.items():
  p=S/name;lines=p.read_text().splitlines();out.write(f'FILE {p}\nSHA256 {sha(p)}\n')
  for lo,hi in rr:
   for n in range(lo,hi+1):out.write(f'{n}: {lines[n-1]}\n')
inputs=[L,P/'fits/joint12/execution.json',P/'fits/joint12/fit.FITRES.TEXT',P/'fits/joint12/fit.nml',P/'joint-failure.json',P/'validators.py',P/'run.py',P/'protocol.json',P/'freeze.json']+[S/k for k in source_ranges]
save(O/'input-hashes.json',{str(p.relative_to(R)):sha(p) for p in inputs})
print(json.dumps(result,indent=2))
