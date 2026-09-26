"""Independent saved-output NIR verification. No fitter/executor/parser import."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,hashlib,re,collections,csv
import numpy as np
from astropy.io import fits
O=Path(__file__).resolve().parent;H=O.parent;P=H.parent.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
CID=[str(i) for i in range(1,9)]
paths={name:P/'fits-start-only'/name for name in ['NIR_true12','NIR_true9','NIR_trueminus','NIR_trueplus','NIR_estimated12','NIR_estimated9','NIR_estimatedminus','NIR_estimatedplus','NIR_null12']}
for ff in [H/'execution-freeze.json',H/'head-byte-adapter/freeze.json']:
 for f,h in json.loads(ff.read_text())['files'].items():assert sha(R/f)==h,f
with fits.open(P/'readme-adapted/ledger/PTE/PTE_HEAD.FITS') as hd:sourcehead=hd[1].data.copy()
with fits.open(P/'readme-adapted/ledger/PTE/PTE_PHOT.FITS') as hd:phot=hd[1].data.copy()
wanted={};domain={}
for h in sourcehead:
 c=str(h['SNID']).strip();pp=phot[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])];assert len(pp)==117
 wanted[c]=collections.Counter((str(x['BAND']).strip(),float(x['MJD']),float(x['FLUXCAL']),float(x['FLUXCALERR'])) for x in pp if str(x['BAND']).strip() in ['J','H']);assert sum(wanted[c].values())==6
 z=float(h['REDSHIFT_HELIO']);domain[c]=(float(np.max(pp['MJD']-(1+z)*70)),float(np.min(pp['MJD']+(1+z)*20)))

def parse(w):
 nml=(w/'fit.nml').read_text();nit=int(re.search(r'NFIT_ITERATION\s*=\s*(\d+)',nml).group(1));ds=w/re.search(r"PRIVATE_DATA_PATH\s*=\s*'([^']+)'",nml).group(1)/'PTE'
 assert 'INISTP_PEAKMJD = 0' in nml
 with fits.open(ds/'PTE_HEAD.FITS') as hd:head=hd[1].data.copy()
 assert (ds/'PTE_PHOT.FITS').resolve()==(P/'readme-adapted/ledger/PTE/PTE_PHOT.FITS').resolve()
 peaks={str(x['SNID']).strip():float(x['PEAKMJD']) for x in head};assert set(peaks)==set(CID)
 entries={};blocks={};sup={};dsum={};bounds={};allcsp=[];allguard=[];shift=[]
 for line in (w/'native.log').read_text().splitlines():
  v=line.split()
  if not v:continue
  tag=v[0]
  if tag.startswith('CSP_'):allcsp.append(v)
  if tag.startswith(('PROSP_DOMAIN','PROSP_PEAK_BOUNDS','PROSP_SUPPORT')):allguard.append(v)
  assert not tag.startswith(('PROSP_MNPARM','PROSP_MNPOUT')),'NIR peak-offset environment active'
  if tag=='CSP_ENTRY:':
   assert len(v)==18;key=(v[1],int(v[2]));assert key not in entries;entries[key]=v
  elif tag in ['CSP_ROW:','CSP_OBJECTIVE:','CSP_WROW:','CSP_WDIAG:']:
   key=(v[1],int(v[2]));b=blocks.setdefault(key,{'rows':[],'w':[]})
   if tag=='CSP_ROW:':assert len(v)==19;b['rows'].append(v)
   elif tag=='CSP_OBJECTIVE:':assert len(v)==15 and 'obj' not in b;b['obj']=v
   else:b['w'].append(v)
  elif tag=='CSP_START_SHIFT:':assert len(v)==6;shift.append(v)
  elif tag=='PROSP_PEAK_BOUNDS':assert len(v)==19;key=(v[1],int(v[2]));assert key not in bounds;bounds[key]=v
  elif tag=='PROSP_SUPPORT':assert len(v)==11 and v[1] not in sup;sup[v[1]]=v
  elif tag=='PROSP_DOMAIN_SUMMARY':assert len(v)==12 and v[1] not in dsum;dsum[v[1]]=v
  elif tag.startswith(('PROSP_SUPPORT','PROSP_DOMAIN')):raise AssertionError(line)
 expected={(c,it) for c in CID for it in range(1,nit+1)};assert set(blocks)==set(entries)==set(bounds)==expected and set(sup)==set(dsum)==set(CID)
 maxq=0
 for (c,it),b in blocks.items():
  rows=b['rows'];o=b['obj'];assert len(rows)==len(b['w'])==int(o[3])==6
  assert [int(x[3]) for x in rows]==[int(x[3]) for x in b['w']]==list(range(1,7))
  X=np.array([[float(q) for q in x[6:]] for x in rows]);W=np.array([[float(q) for q in x[4:]] for x in b['w']]);W=W if o[4]=='T' else np.diag(W[:,0]);assert X.shape==(6,13) and W.shape==(6,6)
  assert np.isfinite(X).all() and np.isfinite(W).all();assert np.linalg.eigvalsh(W).min()>0;C=np.linalg.inv(W)
  assert collections.Counter((x[5],float(x[6]),float(x[10]),float(x[11])) for x in rows)==wanted[c]
  res=X[:,4]-X[:,2];err=abs(float(res@W@res)+float(o[6])+float(o[7])-float(o[5]));maxq=max(maxq,err);assert err<1e-7
  assert float(o[9])==1 and float(o[10])==0 and float(o[11])==float(entries[c,it][9])==peaks[c]
  assert np.all(X[:,1]>=-20) and np.all(X[:,1]<=70);assert abs(peaks[c]-57707.80078125)<=4
  bd=bounds[c,it];lo,hi=domain[c];assert list(map(float,bd[10:16]))==[lo,hi,lo,hi,lo,hi] and int(bd[3])==117 and float(bd[5])==0
  assert list(map(float,bd[6:10]))==[-20,70,-20,85]
  b.update(X=X,W=W,C=C,D=float(o[8]),peak=float(o[11]),identity=[(x[4],x[5],x[6],x[10],x[11]) for x in rows])
 for c in CID:
  ss=sup[c];g=dsum[c];assert int(ss[2])>0 and int(ss[3])==int(ss[4])==0
  assert int(g[2])==1 and int(g[3])==nit and int(g[4])==int(ss[2]) and int(g[5])>0
  assert float(g[6])>=-20 and float(g[7])<=70 and float(g[8])>=0 and float(g[11])>=0
  lo,hi=domain[c];assert lo<=float(g[9])<=float(g[10])<=hi
  f=blocks[c,nit];prev=blocks[c,nit-1];assert all(blocks[c,it]['identity']==f['identity'] for it in range(1,nit+1))
  assert abs(f['D']-prev['D'])<=.001 and f['peak']==prev['peak']
  L=np.linalg.cholesky(prev['C']);cw=np.linalg.solve(L,f['C']-prev['C']);cw=np.linalg.solve(L,cw.T).T;cn=float(np.linalg.norm(cw,2));assert cn<=.001
  model=f['X'][:,2];y=f['X'][:,4];amp=float(model@f['W']@y/(model@f['W']@model));assert amp>0;dd=float(-2.5*np.log10(amp));assert abs(dd)<=.001
  f.update(last_C_norm=cn,last_D_change=f['D']-prev['D'],fixedC_delta_D=dd)
 names=None;science=[];fitrows=[]
 for line in (w/'fit.FITRES.TEXT').read_text().splitlines():
  if line.startswith(('VARNAMES:','SN:')):
   science.append(line);v=line.split()
   if v[0]=='VARNAMES:':names=v[1:]
   else:assert len(v)-1==len(names);fitrows.append(dict(zip(names,v[1:])))
 assert len(fitrows)==8 and {x['CID'] for x in fitrows}==set(CID) and all(x['ERRFLAG_FIT']=='0' and x['CUTFLAG_SNANA']=='3' for x in fitrows)
 ex=json.loads((w/'execution.json').read_text());assert ex['minuit_peak_shift'] is None and ex['hard_domain_mode']==1 and ex['returncode']==0
 return dict(blocks=blocks,entries=entries,bounds=bounds,sup=sup,dsum=dsum,allcsp=allcsp,allguard=allguard,nit=nit,science=science,head=head,headpath=ds/'PTE_HEAD.FITS',maxq=maxq,execution=ex)

cases={n:parse(w) for n,w in paths.items()}
rows=[];all_numerical=[]
for arm in ['true','estimated']:
 base=cases['NIR_'+arm+'12']
 for suffix in ['9','minus','plus']:
  alt=cases['NIR_'+arm+suffix]
  for c in CID:
   a=base['blocks'][c,12];b=alt['blocks'][c,alt['nit']];delta=b['D']-a['D'];assert abs(delta)<=.001 and a['peak']==b['peak'] and a['identity']==b['identity']
   if suffix in ['minus','plus']:
    ds=-.2 if suffix=='minus' else .2;actual=float(alt['entries'][c,1][6])-float(base['entries'][c,1][6]);assert abs(actual-ds)<1e-4
   all_numerical.append(abs(delta))
null=cases['NIR_null12'];true=cases['NIR_true12']
assert null['allcsp']==true['allcsp'] and null['allguard']==true['allguard'] and null['science']==true['science']
assert sha(null['headpath'])==sha(true['headpath'])
for n in true['head'].dtype.names:assert np.array_equal(true['head'][n],null['head'][n])
# Actual assigned header intervention equals frozen nominal joint peak rounded to native R4.
joint={}
for line in (P/'fits-restricted/joint12/native.log').read_text().splitlines():
 v=line.split()
 if v and v[0]=='CSP_OBJECTIVE:' and v[2]=='12':joint[v[1]]=float(v[11])
reported=json.loads((H/'engineering-result.json').read_text());assert reported['gate_pass']
truth={x['CID']:x['DLMU_true'] for x in json.loads((P/'generation-gate.json').read_text())['branches']['ledger']}
for c in CID:
 t=cases['NIR_true12']['blocks'][c,12];e=cases['NIR_estimated12']['blocks'][c,12]
 assert t['peak']==57707.80078125 and e['peak']==float(np.float32(joint[c])) and t['identity']==e['identity']
 row={'CID':c,'peak_error_days':e['peak']-t['peak'],'D_truth':truth[c],'D_true_peak':t['D'],'D_fitted_peak':e['D'],'paired_D':e['D']-t['D'],'error_true_peak':t['D']-truth[c],'error_fitted_peak':e['D']-truth[c]};rows.append(row)
 original=next(x for x in reported['rows'] if x['CID']==c)
 for a,b in [('D_NIR_truepeak','D_true_peak'),('D_NIR_fittedpeak','D_fitted_peak'),('paired_delta_D','paired_D'),('error_truepeak','error_true_peak'),('error_fittedpeak','error_fitted_peak')]:assert original[a]==row[b]
with (O/'per-object.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
activity=json.loads((H/'native-activity.json').read_text());spent=sum(x['wall_seconds'] for x in activity);assert spent<=120
result={'independent_gate_pass':True,'native_calls_by_review':0,'NIR_jobs':9,'callbacks_checked':sum(len(x['blocks']) for x in cases.values()),'all_frozen_hashes_verified':True,'all6_original_FITS_measurements_exact_each_callback':True,'all_actual_fixed_peaks_equal_assigned_HEAD':True,'full_null_model_error_phase_entry_objective_covariance_science_FITRES_exact':True,'null_HEAD_file_hash_exact':True,'max_objective_closure_error':max(x['maxq'] for x in cases.values()),'max_D_numerical_control_difference':max(all_numerical),'max_last_iteration_C_whitened_change':max(b['last_C_norm'] for x in cases.values() for k,b in x['blocks'].items() if k[1]==x['nit']),'max_fixedC_stationarity_abs_D':max(abs(b['fixedC_delta_D']) for x in cases.values() for k,b in x['blocks'].items() if k[1]==x['nit']),'NIR_mean_calls':sum(int(s[2]) for x in cases.values() for s in x['sup'].values()),'unsupported_mean_calls':0,'all_cumulative_native_seconds':spent,'rows':rows,'engineering_only_summary':{'mean_paired_D_mag':float(np.mean([x['paired_D'] for x in rows])),'max_abs_paired_D_mag':float(np.max(np.abs([x['paired_D'] for x in rows])))},'scope':'Eight engineering draws on a single fixed cadence/truth/model/noise setting. Same photons with only assigned NIR peak changed. Numerical pass does not estimate survey population timing bias or cosmology correction.'}
(O/'result.json').write_text(json.dumps(result,indent=2)+'\n')
fs=[O/'review.py',O/'per-object.csv',O/'result.json',H/'engineering-result.json',H/'head-byte-adapter/resume-result.json',H/'head-byte-adapter/protocol.json',H/'head-byte-adapter/freeze.json',H/'native-activity.json']
for w in paths.values():fs += [w/x for x in ['native.log','fit.nml','fit.FITRES.TEXT','execution.json']]
(O/'manifest.json').write_text(json.dumps({'files':{str(p.relative_to(R)):sha(p) for p in fs}},indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))
