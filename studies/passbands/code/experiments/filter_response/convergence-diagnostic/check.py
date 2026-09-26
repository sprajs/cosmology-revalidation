"""Audit frozen nominal iteration experiment without changing fits or gates."""
from pathlib import Path
import json,hashlib,sys
import numpy as np
from astropy.io import fits
ROOT=Path.cwd();O=Path(__file__).resolve().parent;I=O.parent/'instrumentation';sys.path.insert(0,str(I));from parse_native import parse
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest();K=np.log(10)/2.5

def groups(bs):
 d={}
 for b in bs:d.setdefault(b['CID'],[]).append(b)
 return d

def ident(a,b):return [(r['band'],r['MJD'],r['dataF'],r['data_error'],r['source_epoch']) for r in a['rows']]==[(r['band'],r['MJD'],r['dataF'],r['data_error'],r['source_epoch']) for r in b['rows']]
def E(b):
 x=b['array'];return np.diag(np.asarray(np.asarray(x[:,5],dtype=np.float32)**2+np.asarray(x[:,9],dtype=np.float32)**2,dtype=float))
def step(a,b):
 x,y=a['C'],b['C'];L=np.linalg.cholesky(x);t=np.linalg.solve(L,y-x);wh=np.linalg.solve(L,t.T).T
 return dict(iterations=[a['ITER'],b['ITER']],delta_D=b['objective']['D']-a['objective']['D'],identity=ident(a,b),C_relative_Frobenius=float(np.linalg.norm(y-x)/np.linalg.norm(x)),W_relative_Frobenius=float(np.linalg.norm(b['W']-a['W'])/np.linalg.norm(a['W'])),C_whitened_operator_norm=float(np.max(np.abs(np.linalg.eigvalsh((wh+wh.T)/2)))))
def arithmetic(b):
 x=b['array'];f=x[:,2];y=x[:,4];W=b['W'];r=y-f;o=b['objective'];Q=float(r@W@r);amp=float(f@W@y/(f@W@f));delta=float(-np.log(amp)/K) if amp>0 else None
 return dict(ITER=b['ITER'],D=o['D'],n=b['NFITDATA'],Q_native=o['totalQ'],Q_data=Q,Q_closure=Q+o['priorQ']+o['sigmaQ']-o['totalQ'],priorQ=o['priorQ'],sigmaQ=o['sigmaQ'],frozen_C_amplitude=amp,frozen_C_delta_D=delta,positive_means=bool(np.all(f>0)),C_min_eigenvalue=float(np.linalg.eigvalsh(b['C']).min()),phase_min=float(x[:,1].min()),phase_max=float(x[:,1].max()),model_mag_error_max=float(x[:,3].max()),shape=o['shape'])

def main():
 p=json.loads((O/'protocol.json').read_text());runs={j['name']:groups(parse(O/'fits'/j['name']/'fit.log')) for j in p['jobs']};old=groups(parse(ROOT/'runs/research_2026_09_26/csp_native_filter_response/cohort-execution/full/fits/nominal/fit.log'));assert len(old)==42
 source=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b';
 with fits.open(source/'kcor/kcor_CSPDR3_BD17.fits') as h:
  t=h['FilterTrans'].data;w=t['wavelength (A)'];bands={r['band'] for bs in old.values() for b in bs for r in b['rows']};bound={c:(float(w[t['CSP-'+c]!=0].min()),float(w[t['CSP-'+c]!=0].max())) for c in bands}
 results=[];histories={};maxq=0;allCmap=[]
 for cid in old:
  case={'CID':cid,'original_D3':old[cid][-1]['objective']['D'],'runs':{}}
  for name,data in runs.items():
   bs=data[cid];assert len(bs)==int(next(j for j in p['jobs'] if j['name']==name)['iterations']);ars=[arithmetic(b) for b in bs];maxq=max(maxq,max(abs(x['Q_closure']) for x in ars));steps=[step(a,b) for a,b in zip(bs[1:-1],bs[2:])];cm=[]
   for i in range(2,len(bs)):
    a,b,c=bs[i-2:i+1];assert ident(b,c);s=np.exp(-K*(b['objective']['D']-a['objective']['D']));pred=E(c)+(b['C']-E(b))*s*s;err=float(np.linalg.norm(c['C']-pred)/np.linalg.norm(c['C']));cm.append(err);allCmap.append(err)
   ranges=[(min(bound[r['band']][0]/(1+r['z']) for r in b['rows']),max(bound[r['band']][1]/(1+r['z']) for r in b['rows'])) for b in bs]
   final_support=all(-20<=z['phase_min'] and z['phase_max']<=70 and .7<=z['shape']<=1.3 for z in ars[1:]);initial_support=-20<=ars[0]['phase_min'] and ars[0]['phase_max']<=70
   throughput=all(lo>=b['array'][:,11].min() and hi<=b['array'][:,12].max() for b,(lo,hi) in zip(bs,ranges))
   gates={'last_two_D':all(abs(s['delta_D'])<=.001 for s in steps[-2:]),'last_two_C':all(s['C_whitened_operator_norm']<=.001 for s in steps[-2:]),'same_final_mask':all(s['identity'] for s in steps),'final_frozen_C_D':abs(ars[-1]['frozen_C_delta_D'])<=.001,'positive_means_covariance':all(a['positive_means'] and a['C_min_eigenvalue']>0 for a in ars),'all_objectives_close':all(abs(a['Q_closure'])<=1e-8 for a in ars),'Cmap':max(cm)<=5e-6,'final_support':final_support,'throughput_support':bool(throughput)}
   case['runs'][name]=dict(final_D=ars[-1]['D'],final_Q=ars[-1]['Q_native'],last_two_steps=steps[-2:],initial_empirical_grid_support=initial_support,gates=gates)
   histories[name+'__'+cid]=dict(arithmetic=ars,steps=steps,covariance_map_relative_errors=cm,rest_throughput_bounds=ranges)
  d9=runs['iter09_default'][cid][-1]['objective']['D'];d12=runs['iter12_default'][cid][-1]['objective']['D'];vals=[runs[n][cid][-1]['objective']['D'] for n in ['iter12_default','iter12_minus','iter12_plus']]
  actual=[runs[n][cid][0]['entry']['D_entry']-runs['iter12_default'][cid][0]['entry']['D_entry'] for n in ['iter12_minus','iter12_plus']];case.update(D9=d9,D12=d12,delta_D12_D3=d12-case['original_D3'],delta_D12_D9=d12-d9,final_D_multistart_range=max(vals)-min(vals),actual_first_entry_offsets=actual)
  case['cross_run_gates']={'nine_to_twelve':abs(d12-d9)<=.001,'multistart':max(vals)-min(vals)<=.001,'actual_starts':np.max(np.abs(np.array(actual)-[-.2,.2]))<=1e-12,'final_mask_across_starts':all(ident(runs[n][cid][-1],runs['iter12_default'][cid][-1]) for n in ['iter12_minus','iter12_plus'])};case['cross_run_gates']={k:bool(v) for k,v in case['cross_run_gates'].items()}
  # A formal Gaussian score decomposition, not an alternative estimator or bias claim.
  bs=runs['iter12_default'][cid];b=bs[-1];prev=bs[-2];f=b['array'][:,2];r=b['array'][:,4]-f;s=np.exp(-K*(d12-prev['objective']['D']));ee=E(b);m=s*s*(b['C']-ee);cc=ee+m;ww=np.linalg.inv(cc);wres=ww@r;sg,ld=np.linalg.slogdet(cc);assert sg>0
  mean=float(2*K*f@wres);quad=float(2*K*wres@m@wres);det=float(-2*K*np.trace(ww@m))
  case['formal_Gaussian_score_per_mag']=dict(mean_term=mean,covariance_quadratic_term=quad,logdet_term=det,sum=mean+quad+det,logdetC=ld,native_frozen_C_mean_term=float(2*K*f@b['W']@r),M_min_eigenvalue=float(np.linalg.eigvalsh(m).min()),Cmap_gate_pass=bool(all(v['gates']['Cmap'] for v in case['runs'].values())),qualification='Conditional sourceC(D) algebra, no alternativeMLE, preference, unbiasedness claim or selection-normalized physical likelihood. C12 usesD11; here C(D12) reconstructed onlyafterCmapgates.')
  case['numeric_pass']=bool(all(all(v['gates'].values()) for name,v in case['runs'].items() if name.startswith('iter12')) and all(case['cross_run_gates'].values()));results.append(case)
 prefix_gates=[json.loads((O/'fits'/n/'prefix-gate.json').read_text())['pass'] for n in ['iter09_default','iter12_default']]
 out={'scope':'All42 nominal numerical-recipe convergence; archived2007A header/Qfailure and9initialtail extrapolations remain. No changed-filter response or cosmology.','all42_numeric_pass':all(r['numeric_pass'] for r in results) and all(prefix_gates),'prefix_pass':all(prefix_gates),'max_abs_Q_closure':maxq,'max_Cmap_relative_error':max(allCmap),'initial_grid_unsupported':[r['CID'] for r in results if not r['runs']['iter12_default']['initial_empirical_grid_support']],'failures':[r['CID'] for r in results if not r['numeric_pass']],'mean_all42_D12_minus_D3':float(np.mean([r['delta_D12_D3'] for r in results])),'max_abs_D12_minus_D3':max(abs(r['delta_D12_D3']) for r in results),'cases':results}
 (O/'result.json').write_text(json.dumps(out,indent=2,allow_nan=False)+'\n');(O/'state-histories.json').write_text(json.dumps(histories,indent=2,allow_nan=False)+'\n');print(json.dumps({k:v for k,v in out.items() if k!='cases'},indent=2))
if __name__=='__main__':main()
