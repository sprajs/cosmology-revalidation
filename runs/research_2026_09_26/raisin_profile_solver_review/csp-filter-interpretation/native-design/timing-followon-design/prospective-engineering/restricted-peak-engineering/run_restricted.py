"""New fixed-domain estimator; no native calls at import or preparation.
Each execution stage requires a separately frozen root release.
"""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,hashlib,subprocess,time,argparse,importlib.util,re,collections
import numpy as np
from astropy.io import fits
Q=Path(__file__).resolve().parent;P=Q.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova');SHORT=R/'phase2/pte'
A=R/'runs/research_2026_09_26/astra_design/raisin_timing_assets/instrumentation_2021'
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
parser=load('native_parser',A/'parse_native.py')
null_checker=load('full_null_checker',P/'noiseless-checker-recovery/full_null_check.py')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def release(stage,path):
 d=json.loads(path.read_text());check(d['protocol_sha256']==sha(Q/'execution-protocol.json'),'release protocol mismatch');check(d['freeze_sha256']==sha(Q/'execution-freeze.json'),'release freeze mismatch');check(stage in d['stages'],'stage not released')
 for f,h in json.loads((Q/'execution-freeze.json').read_text())['files'].items():check(sha(R/f)==h,'input hash mismatch '+f)
 return sha(path)
def native(label,work,mode,dshift=None):
 activity=Q/'native-activity.json';rows=json.loads(activity.read_text()) if activity.exists() else []
 spent=sum(x['wall_seconds'] for x in rows);allowed=min(40.,120.-spent);check(allowed>0,'new estimator120s budget exhausted')
 env=os.environ.copy()
 for k in ['CSP_POSTINIT_DSHIFT','PROSP_LEDGER','PROSP_FIT_SUPPORT','PROSP_HARD_PEAK_DOMAIN']:env.pop(k,None)
 snana=SHORT/'restricted-peak-engineering/build';binary=snana/'bin/snlc_fit.exe'
 env.update(SNANA_DIR=str(snana),SNDATA_ROOT=str(SHORT/'fit-private-lookup/SNDATA_ROOT'),LD_LIBRARY_PATH=str(R/'phase2/official/build/sysroot/usr/lib'),PROSP_FIT_SUPPORT='1')
 for k in ['SNANA_DIR','SNDATA_ROOT']:check(len(env[k])<120,k+' too long')
 if mode is not None:env['PROSP_HARD_PEAK_DOMAIN']=str(mode)
 if dshift is not None:env['CSP_POSTINIT_DSHIFT']=str(dshift)
 t=time.monotonic()
 with (work/'native.log').open('x') as f:
  try:rc=subprocess.run([str(binary),'fit.nml'],cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=allowed).returncode
  except subprocess.TimeoutExpired:rc='timeout'
 row=dict(label=label,returncode=rc,wall_seconds=time.monotonic()-t,timeout_seconds=allowed,binary_sha256=sha(binary),work=str(work.relative_to(R)),SNANA_DIR=env['SNANA_DIR'],SNDATA_ROOT=env['SNDATA_ROOT'],hard_domain_mode=mode,dshift=dshift)
 rows.append(row);save(activity,rows);save(work/'execution.json',row);check(rc==0,'native failed '+label)
def prepare_fit(name,dataset,bands='JH',iterations=12,peak_step=0,peak_start=0):
 w=P/'fits-restricted'/name;w.mkdir(parents=True,exist_ok=False);tpl=(P/'inputs-v2/fit-template-v3.nml').read_text();(w/'fit.nml').write_text(tpl.format(data_path=dataset,iterations=iterations,peak_step=peak_step,peak_initializer=('INIVAL_PEAKMJD = '+format(57707.80078125+peak_start,'.17g')) if peak_step else '! fixed peak uses native HEAD PEAKMJD',bands=bands));save(w/'input-freeze.json',{'NML_sha256':sha(w/'fit.nml'),'dataset':dataset,'protocol_sha256':sha(Q/'execution-protocol.json')});return w

def run_fit(name,dataset,bands='JH',iterations=12,peak_step=0,peak_start=0,dshift=None,mode=1):
 w=prepare_fit(name,dataset,bands,iterations,peak_step,peak_start);native(name,w,mode,dshift);return w
def files(branch):
 d=(P/'readme-adapted'/branch/'PTE') if branch in ['ledger','noiseless'] else (P/'generation-v5'/branch/'output/PTE');heads=list(d.glob('*HEAD.FITS*'));phots=list(d.glob('*PHOT.FITS*'));check(len(heads)==len(phots)==1,'one native FITS file pair expected');return d,heads[0],phots[0]

def read(branch):
 d,h,p=files(branch)
 with fits.open(h) as f:head=f[1].data.copy();hc=f[1].columns
 with fits.open(p) as f:phot=f[1].data.copy();pc=f[1].columns
 return head,phot,hc,pc

def base_block_review(w,expected_n,expected_cids,support=True):
 blocks=parser.parse(w/'native.log');by={};detail=[]
 table=[];names=None
 for line in (w/'fit.FITRES.TEXT').read_text().splitlines():
  v=line.split()
  if v and v[0]=='VARNAMES:':names=v[1:]
  elif v and v[0]=='SN:':check(names is not None and len(v)-1==len(names),'FITRES schema');table.append(dict(zip(names,v[1:])))
 check({x['CID'] for x in table}==set(expected_cids) and len(table)==len(expected_cids),'FITRES complete membership');check(all(int(x['ERRFLAG_FIT'])==0 and int(x['CUTFLAG_SNANA'])==3 for x in table),'native error/cut flag')
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
  remaining=[(str(x['BAND']).strip(),float(x['MJD']),float(x['FLUXCAL']),float(x['FLUXCALERR'])) for x in wanted]
  for x in observed:
   key=(x['band'],x['MJD'],x['dataF'],x['data_error']);check(key in remaining,'fitter native FITS measurement mismatch '+repr(key));remaining.remove(key)
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
  # same native data/order in every callback; FITS exposure MJD remains R8.
  identity=lambda z:[(r['source_epoch'],r['band'],r['MJD'],r['dataF'],r['data_error']) for r in z['rows']]
  check(all(identity(x)==identity(b) for x in bb),'state accepted rows changed');L=np.linalg.cholesky(prev['C']);d=np.linalg.solve(L,b['C']-prev['C']);d=np.linalg.solve(L,d.T).T;metric=float(np.linalg.norm(d,2));check(metric<=.001,'last iterations whitened C')
  x=b['array'];f=x[:,2];y=x[:,4];a=float(f@b['W']@y/(f@b['W']@f));check(a>0,'nonpositive amplitude');da=-2.5*np.log10(a);check(abs(da)<=.001,'fixedC amplitude stationarity');detail.append({'CID':cid,'D':ob['D'],'peak':ob['peak_absolute'],'dataQ':ob['totalQ']-ob['priorQ']-ob['sigmaQ'],'priorQ':ob['priorQ'],'sigmaQ':ob['sigmaQ'],'fixedC_delta_D':da,'last_C_whitened_norm':metric,'max_standardized_model_residual':float(np.max(np.abs((x[:,4]-x[:,2])/x[:,5]))),'mean_calls':int(sup[cid][0]) if support else None})
 return by,detail

def domain_review(w,by,cids):
 mode=json.loads((w/'execution.json').read_text())['hard_domain_mode']
 lines=[x.split() for x in (w/'native.log').read_text().splitlines() if x.startswith(('PROSP_DOMAIN','PROSP_PEAK_BOUNDS'))]
 if mode!=1:
  check(not lines,'disabled mode emitted domain records');return {'enabled':False}
 schema=json.loads((Q/'build-protocol.json').read_text())['schema'];bounds={};summaries={}
 for v in lines:
  check(v[0] in ['PROSP_PEAK_BOUNDS','PROSP_DOMAIN_SUMMARY'],'domain abort/reject or unknown record')
  check(len(v)-1==len(schema[v[0]]),'domain schema mismatch');d=dict(zip(schema[v[0]],v[1:]));c=d['CID']
  check(c in cids,'unexpected domain CID')
  if v[0]=='PROSP_PEAK_BOUNDS':bounds.setdefault(c,[]).append(d)
  else:check(c not in summaries,'duplicate domain summary');summaries[c]=d
 check(set(bounds)==set(summaries)==set(cids),'missing domain records')
 head,phot,_,_=read('noiseless' if w.name.startswith('noiseless') else 'ledger')
 nit=int(re.search(r'NFIT_ITERATION\s*=\s*(\d+)',(w/'fit.nml').read_text()).group(1));out=[]
 for c in cids:
  hh=next(x for x in head if str(x['SNID']).strip()==c);tt=phot[int(hh['PTROBS_MIN'])-1:int(hh['PTROBS_MAX'])]['MJD'];z=float(hh['REDSHIFT_HELIO'])
  lo=float(np.max(tt-(1+z)*70));hi=float(np.min(tt+(1+z)*20));check(len(tt)==117,'domain metadata count')
  bb=bounds[c];ss=summaries[c]
  check(len(bb)==nit and {int(x['iteration']) for x in bb}==set(range(1,nit+1)),'bounds iteration coverage')
  check({b['ITER'] for b in by[c]}==set(range(1,nit+1)),'callback iteration coverage')
  for x in bb:
   check(int(x['n_all_metadata'])==117 and float(x['zHEL'])==z,'bound metadata changed')
   check([float(x['phase_lower']),float(x['phase_upper']),float(x['KCOR_lower']),float(x['KCOR_upper'])]==[-20,70,-20,85],'table phase bounds changed')
   check(float(x['raw_absolute_lower'])==float(x['safe_absolute_lower'])==lo and float(x['raw_absolute_upper'])==float(x['safe_absolute_upper'])==hi,'bound arithmetic/inward adjustment mismatch')
   off=float(x['MJDOFF']);check(off==0,'frozen engineering MJDOFF changed');check(float(x['parameter_lower'])==lo-off and float(x['parameter_upper'])==hi-off,'MJDOFF bound coordinates')
   pk=float(x['initial_absolute_peak']);check(float(x['initial_distance_lower_days'])==pk-lo and float(x['initial_distance_upper_days'])==hi-pk,'initial endpoint distances')
   check(lo<=pk<=hi,'initial peak outside domain')
  check(int(ss['configured'])==1 and int(ss['bounds_assertion_count'])==nit,'domain setup incomplete')
  check(int(ss['physical_mean_calls'])>0 and int(ss['FCN_physical_calls'])>0,'missing guarded physical calls')
  check(float(ss['min_phase'])>=-20 and float(ss['max_phase'])<=70 and float(ss['minimum_phase_edge_distance_rest_days'])>=0,'guarded phase failure')
  check(float(ss['min_FCN_absolute_peak'])>=lo and float(ss['max_FCN_absolute_peak'])<=hi and float(ss['minimum_peak_edge_distance_observer_days'])>=0,'guarded peak failure')
  support=next(x.split() for x in (w/'native.log').read_text().splitlines() if x.startswith('PROSP_SUPPORT '+c+' '));check(int(ss['physical_mean_calls'])==int(support[2]),'guard versus audit mean-call coverage')
  final=by[c][-1]['objective']['peak_absolute'];check(lo<=final<=hi,'final peak outside domain')
  out.append(dict(CID=c,absolute_bounds=[lo,hi],final_peak=final,final_lower_margin_days=final-lo,final_upper_margin_days=hi-final,final_exactly_on_boundary=final in [lo,hi],search_minimum_edge_distance_days=float(ss['minimum_peak_edge_distance_observer_days']),mean_calls=int(ss['physical_mean_calls']),FCN_physical_calls=int(ss['FCN_physical_calls'])))
 return {'enabled':True,'details':out,'interpretation':'Search endpoint contact is separate from final active-bound status.'}

def block_review(w,expected_n,expected_cids,support=True):
 by,detail=base_block_review(w,expected_n,expected_cids,support)
 domain=domain_review(w,by,expected_cids)
 save(Q/'gates'/(w.name+'-gate.json'),{'gate_pass':True,'details':detail,'domain':domain});return by,detail


def compare(w1,w2,n,cids,tolD=.001,tolT=.01):
 b1,d1=block_review(w1,n,cids,support=w1.name!='noiseless_joint_original');b2,d2=block_review(w2,n,cids)
 for c in cids:
  a,b=b1[c][-1],b2[c][-1];check(abs(a['objective']['D']-b['objective']['D'])<=tolD,'paired numerical D');check(abs(a['objective']['peak_absolute']-b['objective']['peak_absolute'])<=tolT,'paired numerical peak');check(np.array_equal(a['array'][:,[0,4,5]],b['array'][:,[0,4,5]]),'paired measurement identity')
 return b1,b2

def adapt_head(kind,peaks):
 out=P/'derived-restricted'/kind;out.mkdir(parents=True,exist_ok=False);src,h,p=files('ledger');dst=out/'PTE';dst.mkdir()
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
 save(out/'adapter-ledger.json',{'source_HEAD_sha256':sha(h),'source_PHOT_sha256':sha(p),'output_HEAD_sha256':sha(target),'peaks':peaks,'changed_column_only':'PEAKMJD','source_PHOT_same_file':(dst/p.name).resolve()==p.resolve()});return '../../derived-restricted/'+kind

def run_disabled():
 check(json.loads((P/'generation-gate.json').read_text())['gate_pass'],'generation gate')
 cid=[json.loads((P/'generation-gate.json').read_text())['branches']['noiseless'][0]['CID']]
 reference=P/'fits-readme/noiseless_joint';results=[]
 for name,mode in [('noiseless_disabled_absent',None),('noiseless_disabled_zero',0)]:
  w=run_fit(name,'../../readme-adapted/noiseless','grizJH',peak_step=2,mode=mode);block_review(w,117,cid)
  r=null_checker.compare(reference,w,cid,117,True);save(Q/(name+'-full-identity.json'),r);results.append(r)
 save(Q/'disabled-gate.json',{'gate_pass':True,'comparisons':results})

def run_noiseless():
 check(json.loads((Q/'disabled-gate.json').read_text())['gate_pass'],'disabled identity gate')
 cid=[json.loads((P/'generation-gate.json').read_text())['branches']['noiseless'][0]['CID']]
 joint=run_fit('noiseless_joint_active','../../readme-adapted/noiseless','grizJH',peak_step=2);_,d1=block_review(joint,117,cid)
 nir=run_fit('noiseless_NIR_active','../../readme-adapted/noiseless');_,d2=block_review(nir,6,cid)
 truth=json.loads((P/'generation-gate.json').read_text())['branches']['noiseless'][0]['DLMU_true']
 for d in d1+d2:
  check(abs(d['D']-truth)<=.001,'noiseless D recovery');check(abs(d['peak']-57707.80078125)<=.01,'noiseless peak recovery');check(d['max_standardized_model_residual']<=.02,'noiseless native mean mismatch')
 save(Q/'noiseless-gate.json',{'gate_pass':True,'truth_D':truth,'joint':d1,'NIR':d2,'disabled_instrumentation_exact':True})


def run_noisy():
 check(json.loads((Q/'noiseless-gate.json').read_text())['gate_pass'],'noiseless gate');cids=[x['CID'] for x in json.loads((P/'generation-gate.json').read_text())['branches']['ledger']];data='../../readme-adapted/ledger';joint=run_fit('joint12',data,'grizJH',peak_step=2);jb,jd=block_review(joint,117,cids)
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
 null=run_fit('NIR_null12',nulldata);compare(P/'fits-restricted/NIR_true12',null,6,cids,0,0)
 exact=null_checker.compare(P/'fits-restricted/NIR_true12',null,cids,6,True);save(Q/'full-null-identity.json',exact)
 # Exact domain diagnostics are also invariant for a null HEAD adapter.
 for tag in ['PROSP_PEAK_BOUNDS','PROSP_DOMAIN_SUMMARY']:
  a=[x for x in (P/'fits-restricted/NIR_true12/native.log').read_text().splitlines() if x.startswith(tag)]
  b=[x for x in (null/'native.log').read_text().splitlines() if x.startswith(tag)]
  check(a==b,'null domain diagnostics differ')
 rows=[]
 for c in cids:
  truth=next(x['DLMU_true'] for x in json.loads((P/'generation-gate.json').read_text())['branches']['ledger'] if x['CID']==c);d0=finals['true'][c]['D'];d1=finals['estimated'][c]['D'];rows.append({'CID':c,'peak_true':57707.80078125,'peak_fitted_native_R4':float(np.float32(peaks[c])),'D_truth':truth,'D_NIR_truepeak':d0,'D_NIR_fittedpeak':d1,'paired_delta_D':d1-d0,'error_truepeak':d0-truth,'error_fittedpeak':d1-truth})
 save(Q/'engineering-result.json',{'gate_pass':True,'scope':'8 engineering draws, conditional native model/noise/weights and fixed cadence; no population or survey bias estimate','rows':rows,'next64':'not frozen or authorized by this result'})
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['disabled','noiseless','noisy']);ap.add_argument('--release',type=Path,required=True);args=ap.parse_args();rel=release(args.stage,args.release)
 try:globals()['run_'+args.stage]()
 except Exception as e:save(Q/(args.stage+'-failure.json'),{'exception':type(e).__name__,'message':str(e),'release_sha256':rel,'protocol_sha256':sha(Q/'execution-protocol.json'),'stage':args.stage});raise
