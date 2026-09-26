"""Additive joint-only convergence diagnostics; no NIR/science-score path."""
from pathlib import Path
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,hashlib,argparse,importlib.util,sys,gc,subprocess,time
D=Path(__file__).resolve().parent;O=D.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova');sys.path.insert(0,str(O))
import checker as c
spec=importlib.util.spec_from_file_location('frozen_original_runner',O/'run.py');old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
CIDS=[str(x) for x in range(1,65)]

def release(path):
 x=json.loads(path.read_text());c.check(x['protocol_sha256']==sha(D/'protocol.json') and x['freeze_sha256']==sha(D/'freeze.json') and x['stage']=='remaining_joint_diagnostics','diagnostic release mismatch')
 for manifest in [O/'freeze.json',D/'freeze.json']:
  for name,h in json.loads(manifest.read_text())['files'].items():c.check(sha(R/name)==h,'frozen file mismatch '+name)
 c.check(sha(O/'fits-native-activity.json')==json.loads((D/'protocol.json').read_text())['starting_activity_sha256'],'native activity changed before diagnostic release')
 return sha(path)

def review(w):
 by,result=c.review_saved(w,CIDS);small=old.compact(by,result['details']);del by;gc.collect();old.resources()
 c.save(D/(w.name+'-check.json'),result);return small,result

def native(name,w,pshift,carried):
 ledger=O/'fits-native-activity.json';rows=json.loads(ledger.read_text());remaining=carried+120.-sum(x['wall_seconds'] for x in rows);limit=min(90.,remaining);c.check(limit>0,'additional120s diagnostic cap exhausted');old.resources()
 env=os.environ.copy()
 for k in ['CSP_POSTINIT_DSHIFT','PROSP_LEDGER','PROSP_FIT_SUPPORT','PROSP_HARD_PEAK_DOMAIN','PROSP_MINUIT_PEAK_SHIFT']:env.pop(k,None)
 snana=R/'phase2/pte/restricted-peak-engineering/start-only-hook/build';binary=snana/'bin/snlc_fit.exe';sndata=R/'phase2/pt64/fit-sndata'
 env.update(SNANA_DIR=str(snana),SNDATA_ROOT=str(sndata),PROSP_FIT_SUPPORT='1',PROSP_HARD_PEAK_DOMAIN='1',LD_LIBRARY_PATH=str(R/'phase2/official/build/sysroot/usr/lib'))
 if pshift is not None:env['PROSP_MINUIT_PEAK_SHIFT']=str(pshift)
 c.check(len(str(snana))<120 and len(str(sndata))<120,'native env buffers');t=time.monotonic()
 with (w/'native.log').open('x') as f:
  try:rc=subprocess.run([str(binary),'fit.nml'],cwd=w,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=limit).returncode
  except subprocess.TimeoutExpired:rc='timeout'
 row=dict(label=name,returncode=rc,wall_seconds=time.monotonic()-t,timeout_seconds=limit,binary_sha256=sha(binary),kind='fits',hard_domain_mode=1,minuit_peak_shift=pshift,dshift=None,work=str(w.relative_to(R)),SNANA_DIR=str(snana),SNDATA_ROOT=str(sndata),purpose='convergence diagnostics only')
 rows.append(row);c.save(ledger,rows);c.save(w/'execution.json',row);c.check(rc==0,'native failure '+name);old.resources()

def main(release_path):
 rel=release(release_path);carried=sum(x['wall_seconds'] for x in json.loads((O/'fits-native-activity.json').read_text()));baseline,bcheck=review(O/'fits/joint12');out={};nominal=(O/'fits/joint12/fit.nml').read_text()
 for name,it,shift in [('diagnostic_joint9',9,None),('diagnostic_joint_minus',12,-2),('diagnostic_joint_plus',12,2)]:
  w=O/'fits'/name;w.mkdir(exist_ok=False);s=nominal.replace('NFIT_ITERATION = 12','NFIT_ITERATION = 9') if it==9 else nominal;c.check(s.count('NFIT_ITERATION = '+str(it))==1,'iteration replacement');(w/'fit.nml').write_text(s)
  c.save(w/'input-freeze.json',{'NML_sha256':sha(w/'fit.nml'),'diagnostic_protocol_sha256':sha(D/'protocol.json'),'nominal_NML_exact_except_iteration9':True})
  native(name,w,shift,carried);small,result=review(w);old.compare(baseline,small);out[name]={'native_return_warnings':len(result['native']['abnormal_minimize_returns']),'original_cross_run_D_peak_thresholds_pass':True}
 c.save(D/'diagnostic-result.json',{'completed_prescribed_joint_diagnostics':True,'original_scientific_timing_stage_released':False,'NIR_executed':False,'convergence_certified':False,'baseline_warnings':len(bcheck['native']['abnormal_minimize_returns']),'jobs':out,'release_sha256':rel,'required_next_decision':'Review all occurrence/status/state results and amplitude-plus-peak fixedC local-minimum evidence; no automatic NIR continuation.'})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--release',type=Path,required=True);a=ap.parse_args()
 try:main(a.release)
 except Exception as e:c.save(D/'diagnostic-failure.json',{'exception':type(e).__name__,'message':str(e),'action':'Stop, retain every64 draw/occurrence and all original failures; no NIR/reseed/refill/science effect.'});raise
