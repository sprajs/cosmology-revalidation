"""Targeted native-tangent / marginal-row sensitivity; no primary artifact edits."""
from pathlib import Path
import sys,os,re,json,hashlib,subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
ROOT=Path(__file__).resolve().parents[3];P=Path(__file__).resolve().parent;O=P/'negmean';V=P/'validation1020';S=ROOT/'runs/research_2026_09_26/sed_identification/nonpositive-model-epochs.csv';EXE=ROOT/'phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe';PY=ROOT/'phase2/env-official/bin/python';EXPORT=ROOT/'scripts/research_2026_09_26/export_flux_objectives.py';K=np.log(10)/2.5
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
ledger=pd.read_csv(S,dtype={'CID':str});ids=sorted(ledger.CID.unique());assert len(ids)==8 and len(ledger)==10;steps={1:.001,2:.0001,3:.01}
assert not O.exists();O.mkdir();protocol={'scope':'Post-result implementation-gate sensitivity selected solely by nonpositive integrated model means; no residual-based outlier rejection','affected_CIDs':ids,'rows':10,'selection':'Native Python mean<=0 OR exported official mean<=0 in unchanged SED ledger','primary_sensitivity':'Full rows/exact native mean and C, replace approximate nuisance J with direct native central differences at saved coordinates','secondary_sensitivity':'Remove only10 flagged rows; use marginal C submatrix and refit nuisance projection, with cached J and direct-native J separately','native_steps':steps,'half_steps':'All increments halved for numerical convergence; preserve both outputs','amplitude':'Exact43 discovery vector unchanged; no refits to real residual coefficients or scientific data','unchanged_rest':'Other1012 validation objects untouched; report full1020 aggregate M,I,gain and chi-square with differences from original','hashes':{str(p.relative_to(ROOT)):sha(p) for p in [S,Path(__file__),V/'fit.nml',V/'frozen-discovery-coefficients.npz',EXE,EXPORT]}};(O/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
variants=[('base',0,0.)]+[(f'{stage}_p{j}_{"plus" if sign>0 else "minus"}',j,sign*h*factor) for stage,factor in [('full',1),('half',.5)] for j,h in steps.items() for sign in [-1,1]];jobs=[]
for name,j,delta in variants:
 d=O/name;d.mkdir();seed=d/'seed.FITRES';lines=['VARNAMES: CID PKMJD x0 x1 c']
 for cid in ids:
  pars=np.load(V/f'objectives/objective_{cid}.npz')['parameters_x0_x1_c_t0'].copy();pars[j]+=delta;x0,x1,c,t0=pars;lines.append(f'SN: {cid} {t0:.17g} {x0:.17g} {x1:.17g} {c:.17g}')
 seed.write_text('\n'.join(lines)+'\n');nml=(V/'fit.nml').read_text();nml=re.sub(r"(?m)^\s*SNCID_LIST_FILE\s*=.*$",f" SNCID_LIST_FILE = '{seed}'",nml);nml=re.sub(r"(?m)^\s*TEXTFILE_PREFIX\s*=.*$",f" TEXTFILE_PREFIX = '{d/'fit'}'",nml);nml=nml.replace('&FITINP','&FITINP\n LFIXPAR_ALL = T',1);(d/'fit.nml').write_text(nml);(d/'cids.txt').write_text('\n'.join(ids)+'\n');jobs.append(d)
env=os.environ.copy();env.update(SNANA_DIR=str(ROOT/'phase2/official/build/SNANA-audit-v3'),SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
def run(d):
 with (d/'fit.log').open('x') as f:r=subprocess.run([str(EXE),str(d/'fit.nml')],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 assert r.returncode==0,d
 r=subprocess.run([str(PY),str(EXPORT),'--log',str(d/'fit.log'),'--output',str(d/'objectives'),'--expected-cids',str(d/'cids.txt')],cwd=ROOT,env=env,capture_output=True,text=True);(d/'export.log').write_text(r.stdout+r.stderr);assert r.returncode==0,(d,r.stderr)
 return d.name
with ThreadPoolExecutor(max_workers=2) as pool:
 for d in pool.map(run,jobs):print('finished',d,flush=True)
basis=np.load(V/'frozen-discovery-coefficients.npz');c=basis['basis_mean'];griz=basis['gauge_griz']@c;original=pd.read_csv(V/'analysis/object-scores.csv',dtype={'CID':str});original=original[original.arm=='published_mask'].set_index('CID');rows=[];jacobians={};checks=[]
def score(C,J,f,y,b):
 L=np.linalg.cholesky(C);jw=solve_triangular(L,J,lower=True);U,s,_=np.linalg.svd(jw,full_matrices=False);assert np.sum(s>s[0]*1e-10)==4
 raw=solve_triangular(L,-K*f*np.array([griz['griz'.index(k)] for k in b]),lower=True);v=raw-U@(U.T@raw);r=solve_triangular(L,y-f,lower=True);rp=r-U@(U.T@r);a=float(v@r);i=float(v@v)
 return dict(M=a,I=i,gain=a-i/2,chi2=float(rp@rp),dimension=len(f)-4)
for cid in ids:
 a=np.load(V/f'objectives/objective_{cid}.npz');cache=np.load(V/f'analysis/objects/{cid}.npz');order=pd.DataFrame({'MJD':a['MJD'],'band':a['band']}).sort_values(['MJD','band'],kind='stable').index.to_numpy()
 for field,cached in [('MJD','MJD'),('band','band'),('model_flux','official_flux_model'),('data_flux','observed_flux'),('data_fluxerr','quoted_error')]:assert np.array_equal(a[field][order],cache[cached]),(cid,field)
 zero=np.load(O/f'base/objectives/objective_{cid}.npz')
 for field in ['MJD','band','model_flux','data_flux','data_fluxerr','parameters_x0_x1_c_t0']:assert np.array_equal(a[field],zero[field]),(cid,'base',field)
 Js={'cached':cache['jacobian_flux']}
 for stage,factor in [('full',1.),('half',.5)]:
  J=np.zeros((len(order),4));J[:,0]=-K*a['model_flux']
  for j,h in steps.items():
   pp=[]
   for sign in [-1,1]:
    q=np.load(O/f'{stage}_p{j}_{"plus" if sign>0 else "minus"}'/f'objectives/objective_{cid}.npz')
    for key in ['MJD','band','data_flux','data_fluxerr']:assert np.array_equal(a[key],q[key]),(cid,stage,j,sign,key)
    par=a['parameters_x0_x1_c_t0'].copy();par[j]+=sign*h*factor;assert np.max(abs(par-q['parameters_x0_x1_c_t0']))<1e-12
    pp.append(q['model_flux'])
   J[:,j]=(pp[1]-pp[0])/(2*h*factor)
  Js['native_'+stage]=J[order];jacobians[cid+'__'+stage]=J[order]
 bad=(cache['native_flux_model']<=0)|(cache['official_flux_model']<=0);assert bad.sum()==sum(ledger.CID==cid)
 check=score(cache['exact_covariance'],Js['cached'],cache['official_flux_model'],cache['observed_flux'],cache['band']);r0=original.loc[cid];assert max(abs(check['M']-r0.matched_filter),abs(check['I']-r0.information),abs(check['gain']-r0.fixed_prediction_gain))<1e-10
 for label,J in Js.items():
  for mask,keep in [('full_mask',np.arange(len(bad))),('drop_nonpositive',np.flatnonzero(~bad))]:
   sc=score(cache['exact_covariance'][np.ix_(keep,keep)],J[keep],cache['official_flux_model'][keep],cache['observed_flux'][keep],cache['band'][keep]);rows.append(dict(CID=cid,tangent=label,mask=mask,removed=int(len(bad)-len(keep)),**sc))
 checks.append(dict(CID=cid,nonpositive=int(bad.sum()),base_mean_identity=True,all_derivative_coordinate_epoch_identity=True))
f=pd.DataFrame(rows);f.to_csv(O/'object-scores.csv',index=False);np.savez_compressed(O/'native-jacobians.npz',**jacobians)
base={'M':float(original.matched_filter.sum()),'I':float(original.information.sum()),'gain':float(original.fixed_prediction_gain.sum()),'chi2':float(original.projected_chi2.sum()),'dimension':int(original.projected_dimension.sum())};old=f[(f.tangent=='cached')&(f['mask']=='full_mask')];summary=[]
for (label,mask),g in f.groupby(['tangent','mask']):
 delta={k:float(g[k].sum()-old[k].sum()) for k in base};summary.append(dict(tangent=label,mask=mask,affected_objects=8,removed=int(g.removed.sum()),aggregate={k:base[k]+delta[k] for k in base},change=delta))
result={'status':'PASS: native base identity, exact derivative coordinate/mask identity, independent original score reconstruction','primary_original_1020':base,'sensitivities':summary,'checks':checks,'source_sha256':sha(__file__),'protocol_sha256':sha(O/'protocol.json'),'scope':'Targeted8-object/10-epoch sensitivity; neither all1020 native-tangent closure nor calibrated new significance'};(O/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
