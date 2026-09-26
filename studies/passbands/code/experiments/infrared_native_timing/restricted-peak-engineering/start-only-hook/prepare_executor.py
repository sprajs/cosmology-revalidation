from pathlib import Path
import ast
H=Path(__file__).resolve().parent;Q=H.parent;P=Q.parent
old=(P/'noiseless-checker-recovery/full_null_check.py').read_text()
old=old.replace('def compare(left,right,cids,n,support_both=True):','def compare(left,right,cids,n,support_both=True,iterations=12):').replace('set(range(1,13))','set(range(1,iterations+1))')
(H/'full_identity.py').write_text(old)
r=(Q/'run_restricted.py').read_text();tree=ast.parse(r);fn={x.name:ast.get_source_segment(r,x) for x in tree.body if isinstance(x,ast.FunctionDef)}
prefix='''"""Additive start-only executor. Import/preparation never invokes native code.
Each identity/starts/nir stage requires a separately frozen root release.
"""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,hashlib,subprocess,time,argparse,importlib.util,re
import numpy as np
from astropy.io import fits
H=Path(__file__).resolve().parent;Q=H.parent;P=Q.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova');SHORT=R/'phase2/pte'
def load(name,path):
 s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
previous=load('frozen_restricted',Q/'run_restricted.py')
identity=load('full_identity',H/'full_identity.py')
parser=previous.parser
read=previous.read;files=previous.files
CIDS=[str(i) for i in range(1,9)]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\\n')
def check(ok,msg):
 if not ok:raise RuntimeError(msg)
def release(stage,path):
 d=json.loads(path.read_text());check(d['protocol_sha256']==sha(H/'execution-protocol.json'),'release protocol mismatch');check(d['freeze_sha256']==sha(H/'execution-freeze.json'),'release freeze mismatch');check(stage in d['stages'],'stage not released')
 for f,h in json.loads((H/'execution-freeze.json').read_text())['files'].items():check(sha(R/f)==h,'input hash mismatch '+f)
 return sha(path)
def native(label,w,peak_shift=None,dshift=None):
 carried=json.loads((H/'carried-native-activity.json').read_text());activity=H/'native-activity.json'
 rows=json.loads(activity.read_text()) if activity.exists() else carried.copy()
 check(rows[:len(carried)]==carried,'carried activity changed')
 spent=sum(x['wall_seconds'] for x in rows);allowed=min(40.,120.-spent);check(allowed>0,'120s total budget exhausted')
 env=os.environ.copy()
 for key in ['CSP_POSTINIT_DSHIFT','PROSP_LEDGER','PROSP_FIT_SUPPORT','PROSP_HARD_PEAK_DOMAIN','PROSP_MINUIT_PEAK_SHIFT']:env.pop(key,None)
 snana=SHORT/'restricted-peak-engineering/start-only-hook/build';binary=snana/'bin/snlc_fit.exe'
 env.update(SNANA_DIR=str(snana),SNDATA_ROOT=str(SHORT/'fit-private-lookup/SNDATA_ROOT'),LD_LIBRARY_PATH=str(R/'phase2/official/build/sysroot/usr/lib'),PROSP_FIT_SUPPORT='1',PROSP_HARD_PEAK_DOMAIN='1')
 for k in ['SNANA_DIR','SNDATA_ROOT']:check(len(env[k])<120,k+' too long')
 if peak_shift is not None:env['PROSP_MINUIT_PEAK_SHIFT']=str(peak_shift)
 if dshift is not None:env['CSP_POSTINIT_DSHIFT']=str(dshift)
 t=time.monotonic()
 with (w/'native.log').open('x') as out:
  try:rc=subprocess.run([str(binary),'fit.nml'],cwd=w,env=env,stdout=out,stderr=subprocess.STDOUT,timeout=allowed).returncode
  except subprocess.TimeoutExpired:rc='timeout'
 row=dict(label=label,returncode=rc,wall_seconds=time.monotonic()-t,timeout_seconds=allowed,binary_sha256=sha(binary),work=str(w.relative_to(R)),SNANA_DIR=env['SNANA_DIR'],SNDATA_ROOT=env['SNDATA_ROOT'],hard_domain_mode=1,minuit_peak_shift=peak_shift,dshift=dshift)
 rows.append(row);save(activity,rows);save(w/'execution.json',row);check(rc==0,'native failed '+label)
def prepare_fit(name,ds,bands='JH',iterations=12,peak_step=0):
 w=P/'fits-start-only'/name;w.mkdir(parents=True,exist_ok=False)
 tpl=(P/'inputs-v2/fit-template-v3.nml').read_text()
 (w/'fit.nml').write_text(tpl.format(data_path=ds,iterations=iterations,peak_step=peak_step,peak_initializer='INIVAL_PEAKMJD = 57707.80078125' if peak_step else '! fixed peak uses native HEAD PEAKMJD',bands=bands))
 save(w/'input-freeze.json',{'NML_sha256':sha(w/'fit.nml'),'protocol_sha256':sha(H/'execution-protocol.json'),'dataset':ds});return w

def run_fit(name,ds,bands='JH',iterations=12,peak_step=0,peak_shift=None,dshift=None):
 w=prepare_fit(name,ds,bands,iterations,peak_step)
 if bands=='grizJH':
  baseline=P/'fits-restricted'/('joint9' if iterations==9 else 'joint12')
  check((w/'fit.nml').read_bytes()==(baseline/'fit.nml').read_bytes(),'nominal NML identity')
 native(name,w,peak_shift,dshift);return w

def records(path,tag):return [x.split() for x in path.read_text().splitlines() if x.startswith(tag+' ')]
def shift_review(w,by,cids,domain):
 execution=json.loads((w/'execution.json').read_text());shift=execution.get('minuit_peak_shift')
 a=records(w/'native.log','PROSP_MNPARM_PEAK');b=records(w/'native.log','PROSP_MNPOUT_PEAK')
 aborts=records(w/'native.log','PROSP_MNPARM_ABORT');check(not aborts,'MINUIT hook abort')
 if shift in [None,0]:check(not a and not b,'default hook emitted start records');return {'active':False}
 check(shift in [-2,2],'undeclared peak shift')
 nit=int(re.search(r'NFIT_ITERATION\\s*=\\s*(\\d+)',(w/'fit.nml').read_text()).group(1));check(nit==12,'start branch iteration count')
 schema=json.loads((H/'schema.json').read_text());parsed={}
 for tag,rows in [('PROSP_MNPARM_PEAK',a),('PROSP_MNPOUT_PEAK',b)]:
  out={}
  for row in rows:
   check(len(row)-1==len(schema[tag]),'start record schema');d=dict(zip(schema[tag],row[1:]));key=(d['CID'],int(d['iteration']));check(key not in out,'duplicate start record');out[key]=d
  check(set(out)=={(c,it) for c in cids for it in range(1,nit+1)},'start record coverage');parsed[tag]=out
 original=P/'fits-restricted/joint12'
 # Full source entry identity at iteration1, before the start-only hook.
 own_entries={x[1]:x for x in records(w/'native.log','CSP_ENTRY:') if int(x[2])==1}
 base_entries={x[1]:x for x in records(original/'native.log','CSP_ENTRY:') if int(x[2])==1}
 check(own_entries==base_entries,'first common initialization/prior state changed')
 own_bounds=[x for x in records(w/'native.log','PROSP_PEAK_BOUNDS') if int(x[2])==1]
 base_bounds=[x for x in records(original/'native.log','PROSP_PEAK_BOUNDS') if int(x[2])==1]
 check(own_bounds==base_bounds,'first pre-grid metadata/initialization changed')
 details=[]
 for c in cids:
  dd=next(x for x in domain['details'] if x['CID']==c);lo,hi=dd['absolute_bounds']
  for it in range(1,nit+1):
   aa=parsed['PROSP_MNPARM_PEAK'][c,it];bb=parsed['PROSP_MNPOUT_PEAK'][c,it]
   cb=next(x for x in by[c] if x['ITER']==it);entry=cb['entry']['peak_entry_absolute']
   # Objective schema position12 is the actual common INIVAL prior center.
   obj=next(x for x in records(w/'native.log','CSP_OBJECTIVE:') if x[1]==c and int(x[2])==it)
   prior=float(obj[12]);source=float(aa['native_source_parameter']);stored=float(bb['stored_MINUIT_parameter']);proposed=float(aa['proposed_MNPARM_parameter']);expected=shift if it==1 else 0
   check(source==entry==prior==float(aa['unchanged_prior_center_parameter'])==float(bb['unchanged_prior_center_parameter'])==float(bb['native_source_parameter']),'shared prior source changed')
   check(float(aa['requested_shift_days'])==shift and float(aa['applied_shift_days'])==expected,'applied shift metadata')
   check(proposed==source+expected and proposed==float(bb['proposed_MNPARM_parameter']),'proposed shift arithmetic')
   check(abs(stored-source-expected)<.004 and abs(stored-proposed)<.004,'actual stored start separation')
   if it>1:check(stored==source,'later optimizer entry changed')
   check(bb['parameter_name']=='PKMJD' and int(bb['internal_index'])>0,'MINUIT stored parameter identity')
   check(float(aa['lower_bound'])==float(bb['input_lower_bound'])==float(bb['stored_lower_bound'])==lo,'lower bound changed')
   check(float(aa['upper_bound'])==float(bb['input_upper_bound'])==float(bb['stored_upper_bound'])==hi,'upper bound changed')
   check(float(aa['distance_lower_days'])==proposed-lo and float(aa['distance_upper_days'])==hi-proposed,'entry endpoint margins')
   if it==1:details.append({'CID':c,'source_and_prior':source,'actual_MINUIT_start':stored,'actual_shift_days':stored-source,'requested_shift_days':shift})
 return {'active':True,'first_source_CSP_ENTRY_and_pregrid_bounds_exact':True,'details':details,'scope':'First objective function/state preserved; W may depend on coordinates. Later C/prior histories may differ.'}

def block_review(w,n,cids):
 by,detail=previous.base_block_review(w,n,cids,True);domain=previous.domain_review(w,by,cids);starts=shift_review(w,by,cids,domain)
 save(H/'gates'/(w.name+'-gate.json'),{'gate_pass':True,'source_workdir':str(w.relative_to(R)),'existing_reference':w.parent==P/'fits-restricted','details':detail,'domain':domain,'starts':starts});return by,detail

def compare(w1,w2,n,cids,tolD=.001,tolT=.01):
 a,_=block_review(w1,n,cids);b,_=block_review(w2,n,cids)
 for c in cids:
  x,y=a[c][-1],b[c][-1]
  check(abs(x['objective']['D']-y['objective']['D'])<=tolD,'numerical distance agreement')
  check(abs(x['objective']['peak_absolute']-y['objective']['peak_absolute'])<=tolT,'numerical peak agreement')
  check(np.array_equal(x['array'][:,[0,4,5]],y['array'][:,[0,4,5]]),'paired measurement identity')
 return a,b
'''
adapt=fn['adapt_head'].replace("P/'derived-restricted'/kind","P/'derived-start-only'/kind").replace("'../../derived-restricted/'","'../../derived-start-only/'")
stages='''
def run_identity():
 check(json.loads((Q/'noiseless-gate.json').read_text())['gate_pass'],'existing active noiseless gate')
 reports=[]
 for it in [12,9]:
  original=P/'fits-restricted'/('joint'+str(it))
  for mode,label in [(None,'absent'),(0,'zero')]:
   w=run_fit('identity_joint'+str(it)+'_'+label,'../../readme-adapted/ledger','grizJH',it,2,peak_shift=mode);block_review(w,117,CIDS)
   rep=identity.compare(original,w,CIDS,117,True,it)
   for tag in ['PROSP_PEAK_BOUNDS','PROSP_DOMAIN_SUMMARY']:check(records(original/'native.log',tag)==records(w/'native.log',tag),'disabled domain identity')
   save(H/(w.name+'-identity.json'),rep);reports.append(rep)
 compare(P/'fits-restricted/joint12',P/'fits-restricted/joint9',117,CIDS)
 save(H/'identity-gate.json',{'gate_pass':True,'reports':reports,'existing_joint12_and_joint9_reusable':True})

def run_starts():
 check(json.loads((H/'identity-gate.json').read_text())['gate_pass'],'disabled identity gate')
 reports=[]
 for label,shift in [('minus',-2),('plus',2)]:
  w=run_fit('joint_start_'+label,'../../readme-adapted/ledger','grizJH',12,2,peak_shift=shift)
  compare(P/'fits-restricted/joint12',w,117,CIDS);reports.append(str(w.relative_to(R)))
 save(H/'starts-gate.json',{'gate_pass':True,'jobs':reports,'nominal_joint12_reused':True,'offset_scope':'first MINUIT entry only; common first prior/state unchanged'})

def run_nir():
 check(json.loads((H/'identity-gate.json').read_text())['gate_pass'] and json.loads((H/'starts-gate.json').read_text())['gate_pass'],'identity and start gates')
 joint=P/'fits-restricted/joint12';jb,_=block_review(joint,117,CIDS)
 peaks={c:jb[c][-1]['objective']['peak_absolute'] for c in CIDS};ndata=adapt_head('estimated',peaks);nulldata=adapt_head('null',{c:57707.80078125 for c in CIDS})
 finals={}
 for arm,ds in [('true','../../readme-adapted/ledger'),('estimated',ndata)]:
  w=run_fit('NIR_'+arm+'12',ds);_,dd=block_review(w,6,CIDS);finals[arm]={x['CID']:x for x in dd}
  for suffix,it,shift in [('9',9,None),('minus',12,-.2),('plus',12,.2)]:
   v=run_fit('NIR_'+arm+suffix,ds,iterations=it,dshift=shift);a,b=compare(w,v,6,CIDS)
   if shift is not None:
    for c in CIDS:check(abs((b[c][0]['entry']['D_entry']-a[c][0]['entry']['D_entry'])-shift)<1e-4,'D start not applied')
 null=run_fit('NIR_null12',nulldata);compare(P/'fits-start-only/NIR_true12',null,6,CIDS,0,0)
 exact=identity.compare(P/'fits-start-only/NIR_true12',null,CIDS,6,True,12);save(H/'full-null-identity.json',exact)
 for tag in ['PROSP_PEAK_BOUNDS','PROSP_DOMAIN_SUMMARY']:check(records(P/'fits-start-only/NIR_true12/native.log',tag)==records(null/'native.log',tag),'null domain identity')
 rows=[]
 for c in CIDS:
  truth=next(x['DLMU_true'] for x in json.loads((P/'generation-gate.json').read_text())['branches']['ledger'] if x['CID']==c);d0=finals['true'][c]['D'];d1=finals['estimated'][c]['D']
  rows.append({'CID':c,'peak_true':57707.80078125,'peak_fitted_native_R4':float(np.float32(peaks[c])),'D_truth':truth,'D_NIR_truepeak':d0,'D_NIR_fittedpeak':d1,'paired_delta_D':d1-d0,'error_truepeak':d0-truth,'error_fittedpeak':d1-truth})
 save(H/'engineering-result.json',{'gate_pass':True,'scope':'Same existing8 engineering photons, conditional native template/cadence/noise/weights; no population/survey correction','rows':rows,'next64':'not frozen or authorized'})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['identity','starts','nir']);ap.add_argument('--release',type=Path,required=True);args=ap.parse_args();rel=release(args.stage,args.release)
 try:globals()['run_'+args.stage]()
 except Exception as e:save(H/(args.stage+'-failure.json'),{'exception':type(e).__name__,'message':str(e),'release_sha256':rel,'protocol_sha256':sha(H/'execution-protocol.json'),'stage':args.stage});raise
'''
(H/'run_start_only.py').write_text(prefix+'\n\n'+adapt+'\n\n'+stages);ast.parse((H/'run_start_only.py').read_text())
