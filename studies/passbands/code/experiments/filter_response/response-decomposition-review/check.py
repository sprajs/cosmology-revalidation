"""Independent Cholesky/alternative-multiset review, no native model calls."""
from pathlib import Path
from collections import defaultdict,Counter
import sys,json,hashlib
import numpy as np
R=Path.cwd();O=Path(__file__).resolve().parent;N=O.parent;sys.path.insert(0,str(N/'instrumentation'));from parse_native import parse
K=np.log(10)/2.5
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def groups(p):
 out=defaultdict(list)
 for b in parse(p):out[b['CID']].append(b)
 return out
def key(r):return (r['MJD'],r['dataF'],r['data_error'],r['band'])
def main():
 F=N/'stable-filter-response';V=N/'convergence-diagnostic';P=R/'runs/research_2026_09_26/csp_native_filter_response/root-state-review';nom=groups(V/'fits/iter12_default/fit.log');changed=groups(F/'fits/changed12_default/fit.log');p=json.loads((F/'protocol.json').read_text());target=json.loads((P/'response-decomposition/result.json').read_text());expected={r['CID']:r for r in target['objects']};rawmap=defaultdict(list)
 for r in json.loads((F/'input-map.json').read_text()):rawmap[(r['CID'],tuple(r['key_MJD_dataF_dataE_band']))].append(r['new_band'])
 rows=[];maxerr=0.;ambiguous=[];dups=[]
 for cid in p['membership']:
  a=nom[cid][-1];b=changed[cid][-1];idx=defaultdict(list);newidx=defaultdict(list)
  for i,r in enumerate(a['rows']):idx[key(r)].append(i)
  for i,r in enumerate(b['rows']):newidx[key(r)].append(i)
  order=np.full(a['NFITDATA'],-1,dtype=int)
  for k,ids in idx.items():
   targets=rawmap[(cid,k)]
   if len(set(targets))==1:targets=[targets[0]]*len(ids)
   else:assert len(targets)==len(ids),'Partial ambiguous copies cannot be aligned'
   # Canonical target order with reversed old occurrence order, distinct from root deque alignment.
   if len(ids)>1:
    dup={'CID':cid,'key':k,'multiplicity':len(ids),'targets':targets};dups.append(dup)
    for i in ids[1:]:
     perm=np.arange(a['NFITDATA']);perm[ids[0]],perm[i]=perm[i],perm[ids[0]]
     err=float(np.max(np.abs(a['W'][np.ix_(perm,perm)]-a['W']))/np.max(np.abs(a['W'])));assert err<1e-12
     dup.setdefault('old_covariance_swap_relative_errors',[]).append(err)
    if len(set(targets))>1:ambiguous.append(dup)
   for i,band in zip(sorted(ids,reverse=True),sorted(targets),strict=True):
    nk=(*k[:3],band);assert newidx[nk];order[i]=newidx[nk].pop() # reverse queue, not root popleft
  assert not any(newidx.values()) and np.array_equal(np.sort(order),np.arange(b['NFITDATA']))
  x=a['array'];z=b['array'][order];assert np.array_equal(x[:,[0,4,5]],z[:,[0,4,5]])
  d=b['objective']['D']-a['objective']['D'];h=z[:,2]*np.exp(K*d);y=x[:,4]
  L=np.linalg.cholesky(a['C']);u=np.linalg.solve(L,h);v=np.linalg.solve(L,y);amp=float(np.dot(u,v)/np.dot(u,u));assert amp>0
  mean=-np.log1p(amp-1)/K;remaining=d-mean
  c1=b['C'][np.ix_(order,order)];L1=np.linalg.cholesky(c1);u1=np.linalg.solve(L1,z[:,2]);v1=np.linalg.solve(L1,y);amp1=float(np.dot(u1,v1)/np.dot(u1,u1));assert amp1>0;gap=-np.log1p(amp1-1)/K;assert abs(gap)<=.001
  # The old numerical fit need not be an exactly analytic optimum, even for identical operators.
  u0=np.linalg.solve(L,x[:,2]);a0=float(np.dot(u0,v)/np.dot(u0,u0));gap0=-np.log1p(a0-1)/K
  r={'CID':cid,'affected':cid in p['affected_32'],'total_delta_D':d,'fixed_old_C_mean_response':float(mean),'remaining_weight_update_response':float(remaining),'final_C_optimum_gap':float(gap),'nominal_C_optimum_gap':float(gap0),'oldC_analytic_to_analytic_response':float(mean-gap0),'permutation_changed_positions':int(np.sum(order!=np.arange(len(order))))}
  errs={k:abs(r[k]-expected[cid][k]) for k in ['total_delta_D','fixed_old_C_mean_response','remaining_weight_update_response','final_C_optimum_gap']};r['max_root_difference']=max(errs.values());maxerr=max(maxerr,r['max_root_difference']);assert r['max_root_difference']<1e-12;rows.append(r)
 means={str(n):{k:float(np.mean([r[k] for r in rows if n==42 or r['affected']])) for k in ['total_delta_D','fixed_old_C_mean_response','remaining_weight_update_response']} for n in [32,42]}
 for n in means:
  for k,v in means[n].items():assert abs(v-target['means'][n][k])<1e-12
 controls=[r for r in rows if not r['affected']];assert all(r['total_delta_D']==0 for r in controls)
 out={'pass':True,'scope':'Independent arithmetic only; root prescribed order-dependent bookkeeping, no physical attribution/likelihood ratio.','method':'Own native parser, frozen pre-fit raw-map multiset, reverse occurrence mapping with duplicate covariance invariance, covariance Cholesky whitened dot products.','max_abs_root_difference_mag':maxerr,'means':means,'duplicate_groups':dups,'mixed_target_duplicate_groups':ambiguous,'control_total_exact_zero':True,'max_abs_control_component_mag':max(abs(r['fixed_old_C_mean_response']) for r in controls),'max_abs_control_mean_component_minus_nominal_optimum_gap':max(abs(r['fixed_old_C_mean_response']-r['nominal_C_optimum_gap']) for r in controls),'max_abs_control_analytic_to_analytic_response':max(abs(r['oldC_analytic_to_analytic_response']) for r in controls),'objects':rows}
 (O/'result.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['objects','duplicate_groups','mixed_target_duplicate_groups']},indent=2));print('duplicate groups',len(dups),'mixed',len(ambiguous))
 paths=[O/'check.py',O/'result.json',N/'instrumentation/parse_native.py',F/'protocol.json',F/'input-map.json',V/'fits/iter12_default/fit.log',F/'fits/changed12_default/fit.log',P/'response-decomposition-protocol.json',P/'response-decomposition/result.json',R/'scripts/research_2026_09_26/csp_filter_response_decomposition.py']
 (O/'manifest.json').write_text(json.dumps({'sha256':{str(q.relative_to(R)):sha(q) for q in paths}},indent=2)+'\n')
if __name__=='__main__':main()
