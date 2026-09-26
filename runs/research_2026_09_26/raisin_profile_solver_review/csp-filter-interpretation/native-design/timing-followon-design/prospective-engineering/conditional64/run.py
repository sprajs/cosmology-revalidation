"""Fixed64 conditional study; metadata/imports do not invoke native code.
Every stage requires an explicit root release bound to the frozen protocol.
"""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,hashlib,subprocess,time,argparse,csv,shutil,resource,gc
import numpy as np
from astropy.io import fits
import validators as v
import full_identity as identity
from head_adapter import patch_head
O=Path(__file__).resolve().parent;P=O.parent;H=P/'restricted-peak-engineering/start-only-hook';R=Path('/home/szymon/Documents/ChatGPT/supernova');SHORT=R/'phase2/pt64';OLD=R/'phase2/pte';CIDS=[str(x) for x in range(1,65)]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
check=v.check;save=v.save

def size():return sum(p.stat().st_size for p in O.rglob('*') if p.is_file() and not p.is_symlink())
def resources():
 rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024;child=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*1024;disk=size();check(max(rss,child)<=2_000_000_000,'2GB process RSS cap');check(disk<=2_000_000_000,'2GB output cap');return {'max_RSS_bytes':rss,'max_native_RSS_bytes':child,'output_bytes':disk}
def release(stage,path):
 x=json.loads(path.read_text());check(x['protocol_sha256']==sha(O/'protocol.json') and x['freeze_sha256']==sha(O/'freeze.json'),'release hash mismatch');check(stage in x['stages'],'stage not released')
 for f,h in json.loads((O/'freeze.json').read_text())['files'].items():check(sha(R/f)==h,'frozen hash mismatch '+f)
 check((R/'phase2/pt64').resolve()==O,'short alias changed');return sha(path)

def native(name,w,kind,branch=None,peak_shift=None,dshift=None):
 act=O/(kind+'-native-activity.json');rows=json.loads(act.read_text()) if act.exists() else []
 cap=120. if kind=='generation' else 600.;remaining=cap-sum(x['wall_seconds'] for x in rows);limit=min(90.,remaining);check(limit>0,'native cap exhausted');resources()
 env=os.environ.copy()
 for key in ['CSP_POSTINIT_DSHIFT','PROSP_LEDGER','PROSP_FIT_SUPPORT','PROSP_HARD_PEAK_DOMAIN','PROSP_MINUIT_PEAK_SHIFT']:env.pop(key,None)
 if kind=='generation':
  snana=OLD/'repair-v5'/('original-build' if branch=='original' else 'ledger-build');binary=snana/'bin/snlc_sim.exe';args=['sim.input'];sndata=SHORT/'generator-sndata'
  if branch=='ledger':env['PROSP_LEDGER']='1'
 else:
  snana=OLD/'restricted-peak-engineering/start-only-hook/build';binary=snana/'bin/snlc_fit.exe';args=['fit.nml'];sndata=SHORT/'fit-sndata';env.update(PROSP_FIT_SUPPORT='1',PROSP_HARD_PEAK_DOMAIN='1')
  if peak_shift is not None:env['PROSP_MINUIT_PEAK_SHIFT']=str(peak_shift)
  if dshift is not None:env['CSP_POSTINIT_DSHIFT']=str(dshift)
 env.update(SNANA_DIR=str(snana),SNDATA_ROOT=str(sndata),LD_LIBRARY_PATH=str(R/'phase2/official/build/sysroot/usr/lib'))
 for k in ['SNANA_DIR','SNDATA_ROOT']:check(len(env[k])<120,k+' native buffer')
 t=time.monotonic()
 with (w/'native.log').open('x') as f:
  try:rc=subprocess.run([str(binary)]+args,cwd=w,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=limit).returncode
  except subprocess.TimeoutExpired:rc='timeout'
 row={'label':name,'returncode':rc,'wall_seconds':time.monotonic()-t,'timeout_seconds':limit,'binary_sha256':sha(binary),'kind':kind,'hard_domain_mode':1 if kind=='fits' else None,'minuit_peak_shift':peak_shift,'dshift':dshift,'work':str(w.relative_to(R)),'SNANA_DIR':env['SNANA_DIR'],'SNDATA_ROOT':env['SNDATA_ROOT']};rows.append(row);save(act,rows);save(w/'execution.json',row)
 check(rc==0,'native process failed '+name);resources()

def records(p):
 schema=json.loads((P/'instrumentation/schema-v2.json').read_text());out={x:[] for x in schema}
 for line in p.read_text().splitlines():
  x=line.split()
  if x and x[0] in schema:check(len(x)-1==len(schema[x[0]]),'generation ledger schema');out[x[0]].append(dict(zip(schema[x[0]],x[1:])))
 return out

def audit_generation():
 src=O/'generation/ledger/output/PTE';orig=O/'generation/original/output/PTE'
 with fits.open(src/'PTE_HEAD.FITS') as f:h=f[1].data.copy()
 with fits.open(src/'PTE_PHOT.FITS') as f:p=f[1].data.copy()
 with fits.open(orig/'PTE_HEAD.FITS') as f:h0=f[1].data.copy()
 with fits.open(orig/'PTE_PHOT.FITS') as f:p0=f[1].data.copy()
 for x,y in [(h,h0),(p,p0)]:
  check(x.dtype==y.dtype and x.shape==y.shape,'generation replay table layout')
  for k in x.dtype.names:check(np.array_equal(x[k],y[k],equal_nan=True) if x[k].dtype.kind in 'fc' else np.array_equal(x[k],y[k]),'generation hook changed column '+k)
 reference=records(P/'generation-v5/ledger/native.log');refnoise={(float(x['mjd']),x['band']):x for x in reference['PROSP_NOISE'] if x['cid']=='1'};check(len(refnoise)==117,'engineering reference mean cadence')
 rec=records(O/'generation/ledger/native.log');check(len(rec['PROSP_ATTEMPT'])==64 and not rec['PROSP_REJECT'],'attempt count/hidden rejection');check([int(x['actual_attempt']) for x in rec['PROSP_ATTEMPT']]==list(range(1,65)),'attempt sequence');check(len(rec['PROSP_EVENT'])==len(h)==64,'written event count')
 events={x['cid']:x for x in rec['PROSP_EVENT']};check(len(events)==64 and set(events)==set(CIDS),'unique CID coverage');check([str(x['SNID']).strip() for x in h]==CIDS,'native event ordering')
 cad=list(csv.DictReader((O/'inputs-v2/cadence-ledger.csv').open()));physical=sorted((float(x['mjd']),x['band']) for x in cad if x['keep']=='True');expected=json.loads((O/'inputs-v2/manifest.json').read_text());rows={c:[] for c in CIDS};noise={c:[] for c in CIDS}
 for x in rec['PROSP_ROW']:rows[x['cid']].append(x)
 for x in rec['PROSP_NOISE']:noise[x['cid']].append(x)
 summary=[];full=[]
 for hdr in h:
  c=str(hdr['SNID']).strip();ev=events[c];check(int(ev['attempt'])==int(c) and int(ev['libid'])==11 and int(ev['forced_accept'])==0,'event provenance')
  check(float(ev['peak_true'])==float(ev['peak_header'])==float(hdr['PEAKMJD'])==expected['generation_peak'],'truth peak');check(float(ev['zhel_true'])==float(ev['zcmb_true'])==expected['generation_redshift'],'truth redshift');check(float(ev['shape'])==1 and float(ev['AV'])==0 and float(ev['RV'])==expected['RV'],'truth nuisance');check(float(ev['MWEBV_map'])==float(ev['MWEBV_true']) and float(ev['MWEBV_error'])==0 and float(np.float32(ev['MWEBV_map']))==float(hdr['MWEBV']),'MW truth/header')
  pp=p[int(hdr['PTROBS_MIN'])-1:int(hdr['PTROBS_MAX'])];check(len(pp)==117 and sorted((float(x['MJD']),str(x['BAND']).strip()) for x in pp)==physical,'117 exposure identity')
  written=[x for x in rows[c] if int(x['obsflag_write'])];nn={int(x['epoch']):x for x in noise[c]};check(len(nn)==len(noise[c])==117 and len(written)==117,'complete unique noise ledger');check(set(nn)=={int(x['epoch']) for x in written},'exact written/noise epoch match')
  for x,y in zip(written,pp):
   ep=int(x['epoch']);n=nn[ep];check(float(x['mjd'])==float(y['MJD'])==float(n['mjd']),'native R8 time');check(np.float32(x['fluxcal_native_R4'])==y['FLUXCAL'] and np.float32(x['fluxcal_error_native_R4'])==y['FLUXCALERR'] and int(x['photflag'])==int(y['PHOTFLAG']),'native R4 flux/error/flag');check(n['band']==str(y['BAND']).strip(),'noise band join');check(int(n['attempt'])==int(c),'noise attempt join')
   reference_mean=refnoise[float(n['mjd']),n['band']]
   for field in ['flux_true_pe','src_variance_pe2','sky_variance_pe2','zp_variance_pe2','host_variance_pe2','template_source_pe','template_sky_variance_pe2','true_S_variance_pe2','true_SZ_variance_pe2','true_T_variance_pe2','true_F_variance_pe2','Npe_per_fluxcal']:
    check(n[field]==reference_mean[field],'reused native mean/true-variance reference differs '+field)
   q={k:float(t) for k,t in n.items() if k not in ['band','field']};check(q['saturation_excess_pe']<=0 and q['random_template_option']==0,'noise stage branch');check(q['template_source_pe']==q['true_T_variance_pe2']==q['true_F_variance_pe2']==0,'unexpected template/fudge variance')
   calc=q['flux_true_pe']+q['shift_SZ_pe']+q['shift_T_pe']+q['shift_F_pe'];check(abs(calc-q['flux_observed_pe'])<1e-12*max(1,abs(calc)),'true+noise identity');rv=q['reported_variance_before_realization_pe2']+max(q['flux_observed_pe'],0)-q['flux_true_pe'];check(abs(rv-q['reported_variance_final_pe2'])<1e-10*max(1,abs(rv)) and rv>0,'reported realized variance');check(q['smearflag_flux']==1,'noise switched off')
  summary.append({'CID':c,'DLMU_true':float(ev['DLMU_true']),'negative_rows':int(np.sum(pp['FLUXCAL']<0)),'zero_rows':int(np.sum(pp['FLUXCAL']==0)),'row_count':117});full.append({'event':ev,'rows_including_nonwritten':rows[c],'noise_by_written_epoch':noise[c]})
 save(O/'truth-measurement-ledger.json',{'native_internal_numbers':'17digit original strings; FITS input flux/error R4 and MJD R8 preserved','records':full})
 save(O/'generation-gate.json',{'gate_pass':True,'original_vs_ledger_all_native_columns_exact':True,'attempts':64,'rejects':0,'events':64,'all64_R8_mean_and_true_variance_exact_engineering_reference':True,'noiseless_reused':str((P/'restricted-peak-engineering/noiseless-gate.json').relative_to(R)),'rows':summary,'FITS_HEAD_sha256':sha(src/'PTE_HEAD.FITS'),'FITS_PHOT_sha256':sha(src/'PTE_PHOT.FITS')})
 return src

def adapt_readme(src):
 dst=O/'datasets/ledger/PTE';dst.mkdir(parents=True,exist_ok=False)
 for f in src.iterdir():
  if f.name=='PTE.README':
   old=f.read_bytes();check(old.splitlines()[0].strip()==b'DOCUMENTATION:','README first line');i=old.index(b'\n')+1;new=old[:i]+b'\n'+old[i:];check(old.split()==new.split(),'README tokens changed');check(max(map(len,b'\n'.join(new.splitlines()[:6]).split()))<60,'README six-line buffer');(dst/f.name).write_bytes(new)
   import yaml
   check(yaml.safe_load(old.split(b'DOCUMENTATION_END:')[0])==yaml.safe_load(new.split(b'DOCUMENTATION_END:')[0]),'README YAML changed')
   save(O/'README-adapter.json',{'source_sha256':sha(f),'target_sha256':sha(dst/f.name),'one_blank_line_only':True,'all_tokens_and_DOCANA_equal':True})
  elif f.is_file():(dst/f.name).symlink_to(f.resolve())
 for f in ['PTE_HEAD.FITS','PTE_PHOT.FITS','PTE.LIST']:check((dst/f).resolve()==(src/f).resolve(),'data adapter identity')
 check((O/'fit-sndata/SIM/PATH_SNDATA_SIM.LIST').read_text()=='','fit registry not empty')

def generation():
 for b in ['original','ledger']:
  w=O/'generation'/b;check((w/'output').is_dir() and os.access(w/'output',os.W_OK),'missing/unwritable output parent');native('generate_'+b,w,'generation',branch=b)
 adapt_readme(audit_generation());resources()

def compact(by,detail):
 return {c:{'final':bb[-1]['objective'],'first_entry':bb[0]['entry'],'measurement':[(r['source_epoch'],r['band'],r['MJD'],r['dataF'],r['data_error']) for r in bb[-1]['rows']]} for c,bb in by.items()}

def review(w,n):
 by,detail=v.base_block_review(w,n,CIDS,True);domain=v.domain_review(w,by,CIDS);start=v.shift_review(w,by,CIDS,domain)
 save(w/'gate.json',{'gate_pass':True,'details':detail,'domain':domain,'starts':start});small=compact(by,detail);del by;gc.collect();resources();return small

def run_fit(name,ds='../../datasets/ledger',bands='JH',it=12,peakstep=0,pshift=None,dshift=None):
 w=O/'fits'/name;w.mkdir(parents=True,exist_ok=False);tpl=(O/'inputs-v2/fit-template-v3.nml').read_text();s=tpl.format(data_path=ds,iterations=it,peak_step=peakstep,peak_initializer='INIVAL_PEAKMJD = 57707.80078125' if peakstep else '! fixed peak uses native HEAD PEAKMJD',bands=bands);(w/'fit.nml').write_text(s)
 if pshift is not None:check((w/'fit.nml').read_bytes()==(O/'fits/joint12/fit.nml').read_bytes(),'true start nominal NML mismatch')
 save(w/'input-freeze.json',{'NML_sha256':sha(w/'fit.nml'),'protocol_sha256':sha(O/'protocol.json'),'dataset':ds});native(name,w,'fits',peak_shift=pshift,dshift=dshift);return review(w,117 if bands=='grizJH' else 6)

def compare(a,b,shift=None):
 for c in CIDS:
  x,y=a[c],b[c];check(x['measurement']==y['measurement'],'same signed measurement mask');check(abs(x['final']['D']-y['final']['D'])<=.001,'D numerical agreement');check(abs(x['final']['peak_absolute']-y['final']['peak_absolute'])<=.01,'peak numerical agreement')
  if shift is not None:check(abs(y['first_entry']['D_entry']-x['first_entry']['D_entry']-shift)<1e-4,'D start not applied')

def joint():
 check(json.loads((O/'generation-gate.json').read_text())['gate_pass'],'generation gate');a=run_fit('joint12',bands='grizJH',peakstep=2)
 for name,it,ps in [('joint9',9,None),('joint_minus',12,-2),('joint_plus',12,2)]:compare(a,run_fit(name,bands='grizJH',peakstep=2,it=it,pshift=ps))
 save(O/'joint-gate.json',{'gate_pass':True,'nominal_final':a,'all64_in_all_four_jobs':True})

def adapt_head(kind,peaks):
 src=O/'datasets/ledger/PTE';dst=O/'datasets'/kind/'PTE';dst.mkdir(parents=True,exist_ok=False)
 for f in src.iterdir():
  if f.name!='PTE_HEAD.FITS' and f.is_file():(dst/f.name).symlink_to(f.resolve())
 z=patch_head(src/'PTE_HEAD.FITS',dst/'PTE_HEAD.FITS',peaks,require_null=kind=='null');save(dst.parent/'adapter-ledger.json',z);return '../../datasets/'+kind

def nir():
 j=json.loads((O/'joint-gate.json').read_text());check(j['gate_pass'],'joint gates');peaks={c:j['nominal_final'][c]['final']['peak_absolute'] for c in CIDS}
 ds=adapt_head('estimated',peaks);null=adapt_head('null',{c:57707.80078125 for c in CIDS});final={}
 for arm,data in [('true','../../datasets/ledger'),('estimated',ds)]:
  a=run_fit('NIR_'+arm+'12',ds=data);final[arm]=a
  for suffix,it,dshift in [('9',9,None),('minus',12,-.2),('plus',12,.2)]:compare(a,run_fit('NIR_'+arm+suffix,ds=data,it=it,dshift=dshift),dshift)
 z=run_fit('NIR_null12',ds=null);compare(final['true'],z)
 ex=identity.compare(O/'fits/NIR_true12',O/'fits/NIR_null12',CIDS,6,True,12)
 for tag in ['PROSP_PEAK_BOUNDS','PROSP_DOMAIN_SUMMARY']:check(v.records(O/'fits/NIR_true12/native.log',tag)==v.records(O/'fits/NIR_null12/native.log',tag),'full null domain identity')
 check(sha(O/'datasets/ledger/PTE/PTE_HEAD.FITS')==sha(O/'datasets/null/PTE/PTE_HEAD.FITS'),'null HEAD bytes');save(O/'full-null-identity.json',ex)
 truth={x['CID']:x['DLMU_true'] for x in json.loads((O/'generation-gate.json').read_text())['rows']};rows=[]
 for c in CIDS:
  d0=final['true'][c]['final']['D'];d1=final['estimated'][c]['final']['D'];eta=float(np.float32(peaks[c]))-57707.80078125
  rows.append({'CID':c,'eta_days':eta,'eta_internal_days':peaks[c]-57707.80078125,'peak_quantization_days':float(np.float32(peaks[c]))-peaks[c],'D_truth':truth[c],'D_truepeak':d0,'D_fittedpeak':d1,'error_truepeak':d0-truth[c],'error_fittedpeak':d1-truth[c],'delta_D':d1-d0})
 save(O/'fit-gate.json',{'gate_pass':True,'all64':True,'rows':rows});summarize(rows);resources()

def summarize(rows):
 check(len(rows)==64 and {x['CID'] for x in rows}==set(CIDS),'64 unconditional complete draws')
 keys=['eta_days','error_truepeak','error_fittedpeak','delta_D'];X=np.array([[r[k] for k in keys] for r in rows]);mean=X.mean(0);sd=X.std(0,ddof=1);cov=np.cov(X,rowvar=False,ddof=1)
 scalars={k:{'mean':float(mean[i]),'SD_ddof1':float(sd[i]),'MCSE_mean':float(sd[i]/8),'mean_square':float(np.mean(X[:,i]**2)),'MCSE_mean_square':float(np.std(X[:,i]**2,ddof=1)/8)} for i,k in enumerate(keys)}
 cross=[]
 for i in range(len(keys)):
  for k in range(i+1,len(keys)):
   product=X[:,i]*X[:,k];leave=np.array([np.cov(np.delete(X[:,[i,k]],r,axis=0),rowvar=False,ddof=1)[0,1] for r in range(64)]);se=np.sqrt(63/64*np.sum((leave-leave.mean())**2))
   cross.append({'a':keys[i],'b':keys[k],'raw_crossmoment':float(product.mean()),'MCSE_crossmoment':float(product.std(ddof=1)/8),'covariance_ddof1':float(cov[i,k]),'jackknife_MCSE_covariance':float(se)})
 with (O/'per-draw.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 save(O/'result.json',{'gate_pass':True,'seed':26092672,'n':64,'engineering8_pooled':False,'estimands':scalars,'covariance_order':keys,'covariance_ddof1':cov.tolist(),'crossmoments':cross,'scope':'Conditional native timing/noise/weight effect at one fixed truth/cadence. No population/selection-normalized survey correction or cosmology inference.','native_resource':resources()})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['generation','joint','nir']);ap.add_argument('--release',type=Path,required=True);a=ap.parse_args();rel=release(a.stage,a.release)
 try:globals()[a.stage]()
 except Exception as err:
  save(O/(a.stage+'-failure.json'),{'exception':type(err).__name__,'message':str(err),'release_sha256':rel,'protocol_sha256':sha(O/'protocol.json'),'action':'stop; retain all64 membership, no refill/reseed/effect-based selection','completed_outputs':sorted(str(p.relative_to(O)) for p in (O/'fits').glob('*/fit.FITRES.TEXT')) if (O/'fits').exists() else []});raise
