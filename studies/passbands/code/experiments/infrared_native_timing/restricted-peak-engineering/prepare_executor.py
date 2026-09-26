from pathlib import Path
import ast
Q=Path(__file__).resolve().parent;P=Q.parent
old=(P/'run_engineering_v5_checker_resume.py').read_text();tree=ast.parse(old)
functions={n.name:ast.get_source_segment(old,n) for n in tree.body if isinstance(n,ast.FunctionDef)}
prefix='''"""New fixed-domain estimator; no native calls at import or preparation.
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
def save(p,d):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\\n')
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
'''
# Reuse audited arithmetic and native R8/R4 measurement semantics, not its execution or outputs.
helpers='\n\n'.join(functions[n] for n in ['files','read'])+'\n\n'
b=functions['block_review'].replace('def block_review(', 'def base_block_review(')
b=b.replace("save(P/'noiseless-checker-recovery/gates'/(w.name+'-gate.json'),{'gate_pass':True,'details':detail});return by,detail",'return by,detail')
assert "noiseless-checker-recovery/gates" not in b
helpers+=b+'\n\n'
wrapper='''def domain_review(w,by,cids):
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
 nit=int(re.search(r'NFIT_ITERATION\\s*=\\s*(\\d+)',(w/'fit.nml').read_text()).group(1));out=[]
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
'''
helpers+=wrapper+'\n\n'+functions['compare']+'\n\n'
a=functions['adapt_head'].replace("P/'derived'/kind","P/'derived-restricted'/kind").replace("'../../derived/'","'../../derived-restricted/'")
helpers+=a+'\n\n'
stages='''def run_disabled():
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
'''
noisy=functions['run_noisy'].replace("P/'noiseless-gate.json'","Q/'noiseless-gate.json'").replace("P/'fits-resume/NIR_true12'","P/'fits-restricted/NIR_true12'").replace("P/'engineering-result.json'","Q/'engineering-result.json'")
oldnull="null=run_fit('NIR_null12',nulldata);compare(P/'fits-restricted/NIR_true12',null,6,cids,0,0)"
newnull=oldnull+"\n exact=null_checker.compare(P/'fits-restricted/NIR_true12',null,cids,6,True);save(Q/'full-null-identity.json',exact)\n # Exact domain diagnostics are also invariant for a null HEAD adapter.\n for tag in ['PROSP_PEAK_BOUNDS','PROSP_DOMAIN_SUMMARY']:\n  a=[x for x in (P/'fits-restricted/NIR_true12/native.log').read_text().splitlines() if x.startswith(tag)]\n  b=[x for x in (null/'native.log').read_text().splitlines() if x.startswith(tag)]\n  check(a==b,'null domain diagnostics differ')"
assert oldnull in noisy;noisy=noisy.replace(oldnull,newnull)
end='''
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['disabled','noiseless','noisy']);ap.add_argument('--release',type=Path,required=True);args=ap.parse_args();rel=release(args.stage,args.release)
 try:globals()['run_'+args.stage]()
 except Exception as e:save(Q/(args.stage+'-failure.json'),{'exception':type(e).__name__,'message':str(e),'release_sha256':rel,'protocol_sha256':sha(Q/'execution-protocol.json'),'stage':args.stage});raise
'''
(Q/'run_restricted.py').write_text(prefix+helpers+stages+'\n\n'+noisy+end)
ast.parse((Q/'run_restricted.py').read_text())
