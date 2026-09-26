"""Frozen staged native engineering driver. Requires a separate root release file.
No import, preparation, or inspection action generates photons or invokes a fitter.
"""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,hashlib,subprocess,time,shutil,argparse,importlib.util
import numpy as np
from astropy.io import fits
P=Path(__file__).resolve().parent;R=Path('/home/szymon/Documents/ChatGPT/supernova');SHORT=R/'phase2/pte'
A=R/'runs/research_2026_09_26/astra_design/raisin_timing_assets/instrumentation_2021'
spec=importlib.util.spec_from_file_location('native_parser',A/'parse_native.py');parser=importlib.util.module_from_spec(spec);spec.loader.exec_module(parser)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def release(stage,path):
 d=json.loads(path.read_text());check(d['protocol_sha256']==sha(P/'protocol-v5.json'),'release protocol mismatch');check(d['freeze_sha256']==sha(P/'freeze-fit-private.json'),'release freeze mismatch');check(stage in d['stages'],'stage not released')
 for f,h in json.loads((P/'freeze-fit-private.json').read_text())['files'].items():check(sha(R/f)==h,'input hash mismatch '+f)
 return sha(path)
def native(label,work,binary,snana,args,ledger=False,support=False,dshift=None):
 activity=P/'native-activity-v5.json';rows=json.loads(activity.read_text()) if activity.exists() else []
 spent=sum(x['wall_seconds'] for x in rows);allowed=min(40.,120.-spent);check(allowed>0,'120s stage budget exhausted')
 env=os.environ.copy();env.pop('CSP_POSTINIT_DSHIFT',None);env.pop('PROSP_LEDGER',None);env.pop('PROSP_FIT_SUPPORT',None)
 env.update(SNANA_DIR=str(snana),SNDATA_ROOT=str(SHORT/'fit-private-lookup/SNDATA_ROOT') if binary.name=='snlc_fit.exe' else str(R/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(R/'phase2/official/build/sysroot/usr/lib'))
 for k in ['SNANA_DIR','SNDATA_ROOT']:check(len(env[k])<120,k+' too long')
 if ledger:env['PROSP_LEDGER']='1'
 if support:env['PROSP_FIT_SUPPORT']='1'
 if dshift is not None:env['CSP_POSTINIT_DSHIFT']=str(dshift)
 t=time.monotonic()
 with (work/'native.log').open('x') as f:
  try:rc=subprocess.run([str(binary)]+args,cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=allowed).returncode
  except subprocess.TimeoutExpired:rc='timeout'
 row=dict(label=label,returncode=rc,wall_seconds=time.monotonic()-t,timeout_seconds=allowed,binary_sha256=sha(binary),args=args,work=str(work.relative_to(R)),SNANA_DIR=env['SNANA_DIR'],ledger=ledger,support=support,dshift=dshift)
 rows.append(row);save(activity,rows);save(work/'execution.json',row);check(rc==0,'native failed '+label)
def records(path):
 schema=json.loads((P/'instrumentation/schema-v2.json').read_text());out={k:[] for k in schema}
 for line in path.read_text().splitlines():
  v=line.split()
  if not v or v[0] not in schema:continue
  check(len(v)-1==len(schema[v[0]]),'ledger schema '+v[0]);out[v[0]].append(dict(zip(schema[v[0]],v[1:])))
 return out
def files(branch):
 d=P/'generation-v5'/branch/'output/PTE';heads=list(d.glob('*HEAD.FITS*'));phots=list(d.glob('*PHOT.FITS*'));check(len(heads)==len(phots)==1,'one native FITS file pair expected');return d,heads[0],phots[0]
def read(branch):
 d,h,p=files(branch)
 with fits.open(h) as f:head=f[1].data.copy();hc=f[1].columns
 with fits.open(p) as f:phot=f[1].data.copy();pc=f[1].columns
 return head,phot,hc,pc

def audit_generation():
 h0,p0,_,_=read('original');h,p,hc,pc=read('ledger');hn,pn,_,_=read('noiseless');check(h0.dtype.names==h.dtype.names and p0.dtype.names==p.dtype.names,'FITS schema equality')
 for x,y in [(h0,h),(p0,p)]:
  for name in x.dtype.names:check(np.array_equal(x[name],y[name],equal_nan=True) if x[name].dtype.kind in 'fc' else np.array_equal(x[name],y[name]),'instrumentation changed '+name)
 report={'original_vs_instrumented_all_native_table_columns_exact':True,'branches':{},'FITS_HEAD_formats':dict(zip(hc.names,hc.formats)),'FITS_PHOT_formats':dict(zip(pc.names,pc.formats))}
 expected=json.loads((P/'inputs-v2/manifest.json').read_text());cad=list(__import__('csv').DictReader((P/'inputs-v2/cadence-ledger.csv').open()));keep=[x for x in cad if x['keep']=='True']
 physical=sorted((float(x['mjd']),x['band']) for x in keep)
 for branch,head,phot,n in [('ledger',h,p,8),('noiseless',hn,pn,1)]:
  rec=records(P/'generation-v5'/branch/'native.log');check(len(rec['PROSP_ATTEMPT'])==n,'actual attempt count');check(not rec['PROSP_REJECT'],'hidden generation rejection');check(len(rec['PROSP_EVENT'])==len(head)==n,'event count')
  events={x['cid']:x for x in rec['PROSP_EVENT']};rr={};noises={}
  for x in rec['PROSP_ROW']:rr.setdefault(x['cid'],[]).append(x)
  for x in rec['PROSP_NOISE']:noises.setdefault(x['cid'],[]).append(x)
  out=[]
  for hdr in head:
   cid=str(hdr['SNID']).strip();ev=events[cid];check(int(ev['forced_accept'])==0,'forcedaccept');check(float(ev['peak_true'])==float(ev['peak_header'])==expected['generation_peak'],'truth/header peak mismatch');check(float(hdr['PEAKMJD'])==expected['generation_peak'],'native R4 peak mismatch');check(float(ev['zhel_true'])==float(ev['zcmb_true'])==expected['generation_redshift'],'frame/redshift mismatch');check(float(ev['shape'])==1 and float(ev['AV'])==0,'wrongtruth nuisance')
   rows=phot[int(hdr['PTROBS_MIN'])-1:int(hdr['PTROBS_MAX'])];check(len(rows)==117,'cadence count');pairs=sorted((float(x['MJD']),str(x['BAND']).strip()) for x in rows);check(pairs==physical,'exposure identity')
   written=[x for x in rr[cid] if int(x['obsflag_write'])];check(len(written)==117,'written ledger count')
   epkeys=[int(x['epoch']) for x in noises[cid]];check(len(epkeys)==len(set(epkeys)),'duplicate noise record');check(set(int(x['epoch']) for x in written).issubset(set(epkeys)),'missing written noise record')
   check(float(ev['RV'])==expected['RV'],'RV truth mismatch');check(float(ev['MWEBV_map'])==float(ev['MWEBV_true']) and float(ev['MWEBV_error'])==0,'MW truth/map/error mismatch');check(float(np.float32(ev['MWEBV_map']))==float(hdr['MWEBV']),'MW header rounding mismatch')
   # Native epoch sequence is expected to agree with native PHOT output, including all signs.
   for x,y in zip(written,rows):
    check(float(x['mjd'])==float(y['MJD']),'MJD serialization mismatch');check(np.float32(x['fluxcal_native_R4'])==y['FLUXCAL'] and np.float32(x['fluxcal_error_native_R4'])==y['FLUXCALERR'],'native R4 data mismatch');check(int(x['photflag'])==int(y['PHOTFLAG']),'native flag mismatch')
   maximum=0.
   for x in noises[cid]:
    v={k:float(q) for k,q in x.items() if k not in ['band','field']};check(v['saturation_excess_pe']<=0 and v['random_template_option']==0,'unplanned noise branch')
    check(v['template_source_pe']==v['true_T_variance_pe2']==v['true_F_variance_pe2']==0,'unexpected template/fudge noise')
    calc=v['flux_true_pe']+v['shift_SZ_pe']+v['shift_T_pe']+v['shift_F_pe'];er=abs(calc-v['flux_observed_pe'])/max(1,abs(calc));maximum=max(maximum,er);check(er<1e-12,'noise shift equation')
    expected_var=v['reported_variance_before_realization_pe2']+max(v['flux_observed_pe'],0)-v['flux_true_pe'];check(abs(expected_var-v['reported_variance_final_pe2'])<1e-10*max(1,abs(expected_var)),'reported variance equation')
    check(v['reported_variance_final_pe2']>0,'nonpositive reported variance')
    if branch=='noiseless':check(v['shift_SZ_pe']==v['shift_T_pe']==v['shift_F_pe']==0,'noiseless draw nonzero')
   out.append({'CID':cid,'negative_rows':int(np.sum(rows['FLUXCAL']<0)),'zero_rows':int(np.sum(rows['FLUXCAL']==0)),'row_count':len(rows),'noise_mean_equation_maxrelative':maximum,'DLMU_true':float(ev['DLMU_true'])})
  report['branches'][branch]=out
 report['gate_pass']=True;save(P/'generation-gate.json',report);return report

def run_generation():
 for branch in ['original','ledger','noiseless']:
  w=P/'generation-v5'/branch;snana=SHORT/'repair-v5'/('original-build' if branch=='original' else 'ledger-build');binary=snana/'bin/snlc_sim.exe'
  native('generate_'+branch,w,binary,snana,['sim.input'],ledger=branch!='original')
 audit_generation()

def prepare_fit(name,dataset,bands='JH',iterations=12,peak_step=0,peak_start=0):
 w=P/'fits-private'/name;w.mkdir(parents=True,exist_ok=False);tpl=(P/'inputs-v2/fit-template-v3.nml').read_text();(w/'fit.nml').write_text(tpl.format(data_path=dataset,iterations=iterations,peak_step=peak_step,peak_initializer=('INIVAL_PEAKMJD = '+format(57707.80078125+peak_start,'.17g')) if peak_step else '! fixed peak uses native HEAD PEAKMJD',bands=bands));save(w/'input-freeze.json',{'NML_sha256':sha(w/'fit.nml'),'dataset':dataset,'protocol_sha256':sha(P/'protocol-v5.json')});return w

def run_fit(name,dataset,bands='JH',iterations=12,peak_step=0,peak_start=0,dshift=None,original=False):
 w=prepare_fit(name,dataset,bands,iterations,peak_step,peak_start)
 snana=SHORT/'fit-original' if original else SHORT/'fit-support/build';native(name,w,snana/'bin/snlc_fit.exe',snana,['fit.nml'],support=not original,dshift=dshift)
 return w

def block_review(w,expected_n,expected_cids,support=True):
 blocks=parser.parse(w/'native.log');by={};detail=[]
 table=[];names=None
 for line in (w/'fit.FITRES.TEXT').read_text().splitlines():
  v=line.split()
  if v and v[0]=='VARNAMES:':names=v[1:]
  elif v and v[0]=='SN:':check(names is not None and len(v)-1==len(names),'FITRES schema');table.append(dict(zip(names,v[1:])))
 check({x['CID'] for x in table}==set(expected_cids) and len(table)==len(expected_cids),'FITRES complete membership');check(all(int(x['ERRFLAG_FIT'])==0 and int(x['CUTFLAG_SNANA'])==0 for x in table),'native error/cut flag')
 nml=(w/'fit.nml').read_text();fixed='INISTP_PEAKMJD = 0' in nml;arm_peaks={}
 if fixed:
  datapath=__import__('re').search(r"PRIVATE_DATA_PATH = '([^']+)'",nml).group(1);headers=list((w/datapath/'PTE').glob('*HEAD.FITS*'));check(len(headers)==1,'arm HEAD')
  with fits.open(headers[0]) as hh:arm_peaks={str(x['SNID']).strip():float(x['PEAKMJD']) for x in hh[1].data}
 for b in blocks:
  by.setdefault(b['CID'],[]).append(b);check(b['NFITDATA']==expected_n,'accepted cadence size');ob=b['objective'];
  if fixed:check(ob['peak_absolute']==arm_peaks[b['CID']],'fixed peak differs from assigned native HEAD')
  x=b['array'];W=b['W'];C=b['C'];r=x[:,4]-x[:,2];q=float(r@W@r)+ob['priorQ']+ob['sigmaQ'];check(abs(q-ob['totalQ'])<1e-7,'objective arithmetic');check(np.linalg.eigvalsh(C).min()>0,'C not SPD');check(np.isfinite(x).all() and np.isfinite(C).all(),'nonfinite state');check(x[:,1].min()>=-20 and x[:,1].max()<=70,'finalcallback phase');check(ob['shape']==1 and ob['AV']==0,'fixed nuisance changed')
 check(set(by)==set(expected_cids),'complete membership')
 branch='noiseless' if w.name.startswith('noiseless') else 'ledger';head,phot,_,_=read(branch)
 for hdr in head:
  c=str(hdr['SNID']).strip();native_rows=phot[int(hdr['PTROBS_MIN'])-1:int(hdr['PTROBS_MAX'])];wanted=[x for x in native_rows if expected_n==117 or str(x['BAND']).strip() in 'JH'];observed=by[c][-1]['rows'];check(len(wanted)==len(observed),'input row count')
  remaining=[(str(x['BAND']).strip(),float(np.float32(x['MJD'])),float(x['FLUXCAL']),float(x['FLUXCALERR'])) for x in wanted]
  for x in observed:
   key=(x['band'],x['MJD'],x['dataF'],x['data_error']);check(key in remaining,'fitter native R4 measurement mismatch '+repr(key));remaining.remove(key)
  check(not remaining,'missing fitter input exposure')
 sup={}
 if support:
  for line in (w/'native.log').read_text().splitlines():
   v=line.split()
   if not v:continue
   check(v[0] not in ['PROSP_SUPPORT_FIRST_BAD','PROSP_SUPPORT_OVERFLOW'],'all-call support violation')
   if v[0]=='PROSP_SUPPORT':check(len(v)==11,'supportschema');check(v[1] not in sup,'duplicate support summary');sup[v[1]]=v[2:];check(int(v[2])>0 and int(v[3])==int(v[4])==0,'support counters')
  check(set(sup)==set(expected_cids),'missing allcall summary')
 for cid,bb in by.items():
  b=bb[-1];prev=next(x for x in reversed(bb[:-1]) if x['ITER']<b['ITER']);ob=b['objective'];last=prev['objective'];check(abs(ob['D']-last['D'])<=.001,'last iterations D');check(abs(ob['peak_absolute']-last['peak_absolute'])<=.01,'last iterations peak');check(abs(ob['peak_absolute']-57707.80078125)<=4,'final peak guard')
  # same native data/order in every callback; native fitter may R4-round exposure MJD.
  identity=lambda z:[(r['source_epoch'],r['band'],r['MJD'],r['dataF'],r['data_error']) for r in z['rows']]
  check(all(identity(x)==identity(b) for x in bb),'state accepted rows changed');L=np.linalg.cholesky(prev['C']);d=np.linalg.solve(L,b['C']-prev['C']);d=np.linalg.solve(L,d.T).T;metric=float(np.linalg.norm(d,2));check(metric<=.001,'last iterations whitened C')
  x=b['array'];f=x[:,2];y=x[:,4];a=float(f@b['W']@y/(f@b['W']@f));check(a>0,'nonpositive amplitude');da=-2.5*np.log10(a);check(abs(da)<=.001,'fixedC amplitude stationarity');detail.append({'CID':cid,'D':ob['D'],'peak':ob['peak_absolute'],'dataQ':ob['totalQ']-ob['priorQ']-ob['sigmaQ'],'priorQ':ob['priorQ'],'sigmaQ':ob['sigmaQ'],'fixedC_delta_D':da,'last_C_whitened_norm':metric,'max_standardized_model_residual':float(np.max(np.abs((x[:,4]-x[:,2])/x[:,5]))),'mean_calls':int(sup[cid][0]) if support else None})
 save(w/'gate.json',{'gate_pass':True,'details':detail});return by,detail

def compare(w1,w2,n,cids,tolD=.001,tolT=.01):
 b1,d1=block_review(w1,n,cids,support=w1.name!='noiseless_joint_original');b2,d2=block_review(w2,n,cids)
 for c in cids:
  a,b=b1[c][-1],b2[c][-1];check(abs(a['objective']['D']-b['objective']['D'])<=tolD,'paired numerical D');check(abs(a['objective']['peak_absolute']-b['objective']['peak_absolute'])<=tolT,'paired numerical peak');check(np.array_equal(a['array'][:,[0,4,5]],b['array'][:,[0,4,5]]),'paired measurement identity')
 return b1,b2

def run_noiseless():
 check(json.loads((P/'generation-gate.json').read_text())['gate_pass'],'generation gate');cid=[json.loads((P/'generation-gate.json').read_text())['branches']['noiseless'][0]['CID']];data='../../generation-v5/noiseless/output'
 original=run_fit('noiseless_joint_original',data,'grizJH',peak_step=2,original=True);instrumented=run_fit('noiseless_joint',data,'grizJH',peak_step=2)
 a,b=compare(original,instrumented,117,cid,0,0)
 for c in cid:
  check(len(a[c])==len(b[c]),'callback count original/instrumented')
  for x,y in zip(a[c],b[c]):check(np.array_equal(x['array'],y['array']) and np.array_equal(x['W'],y['W']) and x['objective']==y['objective'],'instrumentation changed native mean/C/objective')
 def science(p):return [x for x in p.read_text().splitlines() if x.startswith(('SN:','VARNAMES:'))]
 check(science(original/'fit.FITRES.TEXT')==science(instrumented/'fit.FITRES.TEXT'),'instrumentation FITRES mismatch')
 nir=run_fit('noiseless_NIR',data);_,d1=block_review(instrumented,117,cid);_,d2=block_review(nir,6,cid);truth=json.loads((P/'generation-gate.json').read_text())['branches']['noiseless'][0]['DLMU_true']
 for d in d1+d2:check(abs(d['D']-truth)<=.001,'noiseless D recovery');check(abs(d['peak']-57707.80078125)<=.01,'noiseless peak recovery');check(d['max_standardized_model_residual']<=.02,'noiseless native mean mismatch')
 save(P/'noiseless-gate.json',{'gate_pass':True,'truth_D':truth,'joint':d1,'NIR':d2,'default_instrumentation_exact':True})

def adapt_head(kind,peaks):
 out=P/'derived'/kind;out.mkdir(parents=True,exist_ok=False);src,h,p=files('ledger');dst=out/'PTE';dst.mkdir()
 for f in src.iterdir():
  if f==h:continue
  if f.is_file():(dst/f.name).symlink_to(f.resolve())
 with fits.open(h) as hd:
  before=hd[1].data.copy()
  for row in hd[1].data:row['PEAKMJD']=np.float32(peaks[str(row['SNID']).strip()])
  for key in before.dtype.names:
   if key!='PEAKMJD':check(np.array_equal(before[key],hd[1].data[key]),'adapter touched '+key)
  target=dst/h.name;hd.writeto(target,overwrite=False)
  if kind=='null':check(np.array_equal(before,hd[1].data),'null adapter changed HEADtable')
 save(out/'adapter-ledger.json',{'source_HEAD_sha256':sha(h),'source_PHOT_sha256':sha(p),'output_HEAD_sha256':sha(target),'peaks':peaks,'changed_column_only':'PEAKMJD','source_PHOT_same_file':(dst/p.name).resolve()==p.resolve()});return '../../derived/'+kind

def run_noisy():
 check(json.loads((P/'noiseless-gate.json').read_text())['gate_pass'],'noiseless gate');cids=[x['CID'] for x in json.loads((P/'generation-gate.json').read_text())['branches']['ledger']];data='../../generation-v5/ledger/output';joint=run_fit('joint12',data,'grizJH',peak_step=2);jb,jd=block_review(joint,117,cids)
 for name,it,start in [('joint9',9,0),('joint_minus',12,-2),('joint_plus',12,2)]:
  w=run_fit(name,data,'grizJH',it,2,start);aa,bb=compare(joint,w,117,cids)
  for c in cids:check(abs((bb[c][0]['entry']['peak_entry_absolute']-aa[c][0]['entry']['peak_entry_absolute'])-start)<.004,'peak start was not applied')
 peaks={c:jb[c][-1]['objective']['peak_absolute'] for c in cids};ndata=adapt_head('estimated',peaks);nulldata=adapt_head('null',{c:57707.80078125 for c in cids})
 finals={}
 for arm,ds in [('true',data),('estimated',ndata)]:
  w=run_fit('NIR_'+arm+'12',ds);_,dd=block_review(w,6,cids);finals[arm]={x['CID']:x for x in dd}
  for suffix,it,shift in [('9',9,None),('minus',12,-.2),('plus',12,.2)]:
   q=run_fit('NIR_'+arm+suffix,ds,iterations=it,dshift=shift);aa,bb=compare(w,q,6,cids)
   if shift is not None:
    for c in cids:check(abs((bb[c][0]['entry']['D_entry']-aa[c][0]['entry']['D_entry'])-shift)<1e-4,'amplitude start was not applied')
 null=run_fit('NIR_null12',nulldata);compare(P/'fits-private/NIR_true12',null,6,cids,0,0)
 rows=[]
 for c in cids:
  truth=next(x['DLMU_true'] for x in json.loads((P/'generation-gate.json').read_text())['branches']['ledger'] if x['CID']==c);d0=finals['true'][c]['D'];d1=finals['estimated'][c]['D'];rows.append({'CID':c,'peak_true':57707.80078125,'peak_fitted_native_R4':float(np.float32(peaks[c])),'D_truth':truth,'D_NIR_truepeak':d0,'D_NIR_fittedpeak':d1,'paired_delta_D':d1-d0,'error_truepeak':d0-truth,'error_fittedpeak':d1-truth})
 save(P/'engineering-result.json',{'gate_pass':True,'scope':'8 engineering draws, conditional native model/noise/weights and fixed cadence; no population or survey bias estimate','rows':rows,'next64':'not frozen or authorized by this result'})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['noiseless','noisy']);ap.add_argument('--release',type=Path,required=True);args=ap.parse_args();rel=release(args.stage,args.release)
 try:globals()['run_'+args.stage]()
 except Exception as e:save(P/(args.stage+'-failure-fit-private.json'),{'exception':type(e).__name__,'message':str(e),'release_sha256':rel,'protocol_sha256':sha(P/'protocol-v5.json'),'stage':args.stage});raise
