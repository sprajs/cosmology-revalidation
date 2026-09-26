"""Independent read-only parsing/arithmetic. Does not import either executor/checker."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import hashlib,json,csv,collections,re
import numpy as np
from astropy.io import fits
O=Path(__file__).resolve().parent;H=O.parent;Q=H.parent;P=Q.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
freeze=json.loads((H/'execution-freeze.json').read_text())
for p,h in freeze['files'].items():assert sha(R/p)==h,p
paths={'nominal':P/'fits-restricted/joint12','nominal9':P/'fits-restricted/joint9','minus':P/'fits-start-only/joint_start_minus','plus':P/'fits-start-only/joint_start_plus'}
heads=list((P/'readme-adapted/ledger/PTE').glob('*HEAD.FITS*'));phots=list((P/'readme-adapted/ledger/PTE').glob('*PHOT.FITS*'));assert len(heads)==len(phots)==1
with fits.open(heads[0]) as f:head=f[1].data.copy()
with fits.open(phots[0]) as f:phot=f[1].data.copy()
CID=[str(x['SNID']).strip() for x in head];assert CID==[str(x) for x in range(1,9)]
data={}
for h in head:
 c=str(h['SNID']).strip();pp=phot[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])];assert len(pp)==117
 data[c]=collections.Counter((str(x['BAND']).strip(),float(x['MJD']),float(x['FLUXCAL']),float(x['FLUXCALERR'])) for x in pp)

def parse(d):
 blocks={};entries={};bounds={};mn={};readback={};support={};domain={}
 for line in (d/'native.log').read_text().splitlines():
  v=line.split()
  if not v:continue
  tag=v[0]
  if tag.startswith(('CSP_','PROSP_PEAK_BOUNDS','PROSP_MNPARM','PROSP_MNPOUT')):
   assert tag not in ['PROSP_MNPARM_ABORT']
  if tag in ['CSP_ROW:','CSP_OBJECTIVE:','CSP_WROW:','CSP_WDIAG:']:
   key=(v[1],int(v[2]));b=blocks.setdefault(key,{'rows':[],'Wrows':[]})
   if tag=='CSP_ROW:':assert len(v)==19;b['rows'].append(v)
   elif tag=='CSP_OBJECTIVE:':assert len(v)==15 and 'o' not in b;b['o']=v
   else:b['Wrows'].append(v)
  elif tag=='CSP_ENTRY:':assert len(v)==18;key=(v[1],int(v[2]));assert key not in entries;entries[key]=v
  elif tag=='PROSP_PEAK_BOUNDS':assert len(v)==19;key=(v[1],int(v[2]));assert key not in bounds;bounds[key]=v
  elif tag=='PROSP_MNPARM_PEAK':assert len(v)==12;key=(v[1],int(v[2]));assert key not in mn;mn[key]=v
  elif tag=='PROSP_MNPOUT_PEAK':assert len(v)==13;key=(v[1],int(v[2]));assert key not in readback;readback[key]=v
  elif tag=='PROSP_SUPPORT':assert len(v)==11 and v[1] not in support;support[v[1]]=v
  elif tag=='PROSP_DOMAIN_SUMMARY':assert len(v)==12 and v[1] not in domain;domain[v[1]]=v
  elif tag.startswith(('PROSP_SUPPORT_','PROSP_DOMAIN_')):raise AssertionError(line)
 nit=int(re.search(r'NFIT_ITERATION\s*=\s*(\d+)',(d/'fit.nml').read_text()).group(1));expected={(c,i) for c in CID for i in range(1,nit+1)}
 assert set(blocks)==set(entries)==set(bounds)==expected;assert set(support)==set(domain)==set(CID)
 maxq=0;maxasym=0;mineig=np.inf
 for (c,it),b in blocks.items():
  o=b['o'];n=int(o[3]);assert n==117 and len(b['rows'])==len(b['Wrows'])==n
  assert [int(x[3]) for x in b['rows']]==list(range(1,n+1)) and [int(x[3]) for x in b['Wrows']]==list(range(1,n+1))
  X=np.array([[float(q) for q in x[6:]] for x in b['rows']]);assert X.shape==(117,13)
  W=np.array([[float(q) for q in x[4:]] for x in b['Wrows']]);W=W if o[4]=='T' else np.diag(W[:,0]);assert W.shape==(117,117)
  assert np.isfinite(X).all() and np.isfinite(W).all()
  C=np.linalg.inv(W);ev=np.linalg.eigvalsh(C);assert ev[0]>0;mineig=min(mineig,float(ev[0]));maxasym=max(maxasym,float(np.max(np.abs(W-W.T))))
  y=X[:,4];f=X[:,2];res=y-f;err=abs(float(res@W@res)+float(o[6])+float(o[7])-float(o[5]));maxq=max(maxq,err);assert err<1e-7
  assert collections.Counter((x[5],float(x[6]),float(x[10]),float(x[11])) for x in b['rows'])==data[c]
  assert np.all(X[:,1]>=-20)&np.all(X[:,1]<=70);assert float(o[9])==1 and float(o[10])==0
  prior=float(o[12]);entry=float(entries[c,it][9]);assert prior==entry
  bb=bounds[c,it];assert int(bb[3])==117 and float(bb[4])==float(np.float32(.453)) and float(bb[5])==0
  assert list(map(float,bb[6:10]))==[-20,70,-20,85]
  h=next(x for x in head if str(x['SNID']).strip()==c);t=phot[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])]['MJD'];z=float(h['REDSHIFT_HELIO'])
  lo=float(np.max(t-(1+z)*70));hi=float(np.min(t+(1+z)*20));assert list(map(float,bb[10:16]))==[lo,hi,lo,hi,lo,hi]
  assert lo<=float(o[11])<=hi
  b.update(X=X,W=W,C=C,identity=[(x[4],x[5],x[6],x[10],x[11]) for x in b['rows']],lo=lo,hi=hi)
 for c in CID:
  s=support[c];g=domain[c];assert int(s[2])>0 and int(s[3])==int(s[4])==0
  assert int(g[2])==1 and int(g[3])==nit and int(g[4])==int(s[2]) and int(g[5])>0
  assert float(g[6])>=-20 and float(g[7])<=70 and float(g[8])>=0 and float(g[11])>=0
  final=blocks[c,nit];prev=blocks[c,nit-1]
  assert all(blocks[c,i]['identity']==final['identity'] for i in range(1,nit+1))
  D,T=float(final['o'][8]),float(final['o'][11]);assert abs(D-float(prev['o'][8]))<=.001 and abs(T-float(prev['o'][11]))<=.01 and abs(T-57707.80078125)<=4
  L=np.linalg.cholesky(prev['C']);B=np.linalg.solve(L,final['C']-prev['C']);B=np.linalg.solve(L,B.T).T;norm=float(np.linalg.norm(B,2));assert norm<=.001
  ff=final['X'][:,2];yy=final['X'][:,4];aa=float(ff@final['W']@yy/(ff@final['W']@ff));assert aa>0
  delta=-2.5*np.log10(aa);assert abs(delta)<=.001
  final.update(last_C_change=norm,fixedC_D_stationarity=float(delta),last_D_change=D-float(prev['o'][8]),last_T_change=T-float(prev['o'][11]))
 names=None;fitrows=[]
 for line in (d/'fit.FITRES.TEXT').read_text().splitlines():
  v=line.split()
  if not v:continue
  if v[0]=='VARNAMES:':names=v[1:]
  elif v[0]=='SN:':assert len(v)-1==len(names);fitrows.append(dict(zip(names,v[1:])))
 assert len(fitrows)==8 and {x['CID'] for x in fitrows}==set(CID)
 assert all(int(x['ERRFLAG_FIT'])==0 and int(x['CUTFLAG_SNANA'])==3 for x in fitrows)
 return dict(blocks=blocks,entries=entries,bounds=bounds,mn=mn,readback=readback,support=support,domain=domain,nit=nit,maxq=maxq,maxW_asymmetry=maxasym,minC_eigenvalue=mineig)

allp={k:parse(d) for k,d in paths.items()};rows=[]
for label,shift in [('minus',-2),('plus',2)]:
 p=allp[label];base=allp['nominal'];assert (paths[label]/'fit.nml').read_bytes()==(paths['nominal']/'fit.nml').read_bytes()
 assert set(p['mn'])==set(p['readback'])==set(p['blocks'])
 for c in CID:
  assert p['entries'][c,1]==base['entries'][c,1] and p['bounds'][c,1]==base['bounds'][c,1]
  for it in range(1,13):
   m=p['mn'][c,it];r=p['readback'][c,it];src=float(m[5]);proposed=float(m[6]);stored=float(r[7]);delta=shift if it==1 else 0
   assert int(m[3])==shift and int(m[4])==delta
   assert src==float(p['entries'][c,it][9])==float(p['blocks'][c,it]['o'][12])==float(m[7])==float(r[5])==float(r[8])
   assert float(r[6])==proposed==src+delta and abs(stored-src-delta)<.004 and abs(stored-proposed)<.004
   if it>1:assert stored==src
   assert r[4]=='PKMJD' and int(r[3])>0
   b=p['blocks'][c,it];assert [float(m[8]),float(m[9])]==[b['lo'],b['hi']]
   assert list(map(float,r[9:13]))==[b['lo'],b['hi'],b['lo'],b['hi']]
  f=p['blocks'][c,12];n=base['blocks'][c,12];Ddiff=float(f['o'][8])-float(n['o'][8]);Tdiff=float(f['o'][11])-float(n['o'][11]);assert abs(Ddiff)<=.001 and abs(Tdiff)<=.01
  rows.append({'branch':label,'CID':c,'actual_first_shift_days':float(p['readback'][c,1][7])-float(p['mn'][c,1][5]),'first_prior_center_exact':True,'final_D_minus_nominal':Ddiff,'final_peak_minus_nominal_days':Tdiff,'last_D_change':f['last_D_change'],'last_peak_change':f['last_T_change'],'last_C_whitened_change':f['last_C_change'],'fixedC_D_stationarity':f['fixedC_D_stationarity'],'final_lower_margin_days':float(f['o'][11])-f['lo'],'final_upper_margin_days':f['hi']-float(f['o'][11]),'mean_calls':int(p['support'][c][2]),'unsupported_calls':int(p['support'][c][3]),'first_iteration_W_max_difference_at_different_final_parameters':float(np.max(np.abs(p['blocks'][c,1]['W']-base['blocks'][c,1]['W'])))})
for c in CID:
 a=allp['nominal']['blocks'][c,12]['o'];b=allp['nominal9']['blocks'][c,9]['o'];assert abs(float(a[8])-float(b[8]))<=.001 and abs(float(a[11])-float(b[11]))<=.01
with (O/'per-object.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
result={'independent_gate_pass':True,'native_calls':0,'frozen_input_hashes_verified':len(freeze['files']),'callbacks_checked':sum(len(p['blocks']) for p in allp.values()),'actual_start_pairs_checked':sum(len(allp[k]['mn']) for k in ['minus','plus']),'actual_first_shifts':sorted(set(x['actual_first_shift_days'] for x in rows)),'first_common_prior_entry_and_pregrid_records_exact_all8_both_signs':True,'all117FITS_signed_measurements_exact_each_callback':True,'max_objective_closure_error':max(p['maxq'] for p in allp.values()),'max_absolute_D_vs_nominal':max(abs(x['final_D_minus_nominal']) for x in rows),'max_absolute_peak_vs_nominal_days':max(abs(x['final_peak_minus_nominal_days']) for x in rows),'max_last_C_whitened_change':max(x['last_C_whitened_change'] for x in rows),'max_abs_fixedC_D_stationarity':max(abs(x['fixedC_D_stationarity']) for x in rows),'minimum_final_boundary_margin_days':min(min(x['final_lower_margin_days'],x['final_upper_margin_days']) for x in rows),'total_start_branch_mean_calls':sum(x['mean_calls'] for x in rows),'unsupported_start_branch_mean_calls':sum(x['unsupported_calls'] for x in rows),'weight_qualification':'First common initialization and prior state are exact. ITER1 W depends on parameter values, so final W equality across displaced starts is neither required nor claimed; source function is unchanged. Later covariance/prior histories may differ.','scope':'Supports numerical path stability of this conditional restricted estimator on the same8 photons; no NIR paired response or population timing/cosmology correction is established here.'}
(O/'result.json').write_text(json.dumps(result,indent=2)+'\n')
files=[O/'review.py',O/'per-object.csv',O/'result.json',H/'execution-protocol.json',H/'execution-freeze.json',H/'starts-gate.json']
for d in paths.values():files += [d/x for x in ['native.log','fit.FITRES.TEXT','fit.nml','execution.json']]
(O/'manifest.json').write_text(json.dumps({'files':{str(p.relative_to(R)):sha(p) for p in files}},indent=2)+'\n');print(json.dumps(result,indent=2))
