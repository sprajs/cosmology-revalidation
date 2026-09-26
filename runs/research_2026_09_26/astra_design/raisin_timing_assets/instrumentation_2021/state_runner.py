from pathlib import Path
import json,hashlib,shutil,os,time,subprocess,argparse
P=Path(__file__).resolve().parent;R=Path('/home/szymon/Documents/ChatGPT/supernova');B=P.parent/'baseline_2021'
CONDITIONS={'absent':None,'zero':'0','plus':'0.2','minus':'-0.2'}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def prepare():
 assert not (P/'state-freeze.json').exists()
 for name in CONDITIONS:
  w=P/'fits'/name;w.mkdir(parents=True,exist_ok=False);src=B/'fits/baseline'
  shutil.copytree(src/'data',w/'data');shutil.copy2(src/'fit.nml',w/'fit.nml')
  (w/'kcor.fits').symlink_to((src/'kcor.fits').resolve());shutil.copytree(src/'snoopy.B18',w/'snoopy.B18',symlinks=True)
 paths=[P/'state-protocol.json',P/'state_runner.py',P/'state_gate.py',P/'parse_native.py',P/'build/bin/snlc_fit.exe',P/'build-verification.json',P/'instrumented-snlc_fit.car',P/'ported-source.patch',B/'baseline-gate.json',B/'resource-amendment.json']
 paths+=list((P/'fits').rglob('*.DAT'))+list((P/'fits').rglob('*.nml'))
 save(P/'state-freeze.json',{'files':{str(p.relative_to(R)):sha(p) for p in paths},'protocol_sha256':sha(P/'state-protocol.json'),'conditions':CONDITIONS})
 print(sha(P/'state-freeze.json'))
def run(name):
 assert name in CONDITIONS
 protocol=json.loads((P/'state-protocol.json').read_text());freeze=json.loads((P/'state-freeze.json').read_text());release=json.loads((P/'execution-release.json').read_text())
 assert release['protocol_sha256']==sha(P/'state-protocol.json') and release['freeze_sha256']==sha(P/'state-freeze.json')
 for f,h in freeze['files'].items():assert sha(R/f)==h,f
 w=P/'fits'/name;assert not (w/'fit.log').exists();resource=protocol['resources']
 prior=resource['earlier_native_seconds'];spent=sum(json.loads(p.read_text())['wall_seconds'] for p in (P/'fits').glob('*/execution.json'))
 allowed=min(resource['per_process_seconds'],resource['this_stage_seconds']-spent,resource['original_total_pilot_seconds']-prior-spent);assert allowed>0
 env=os.environ.copy();env.pop('CSP_POSTINIT_DSHIFT',None)
 if CONDITIONS[name] is not None:env['CSP_POSTINIT_DSHIFT']=CONDITIONS[name]
 env.update(SNANA_DIR=str(P/'build'),SNDATA_ROOT=str(R/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(R/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
 start=time.monotonic()
 with (w/'fit.log').open('x') as f:
  try:code=subprocess.run([str(P/'build/bin/snlc_fit.exe'),'fit.nml'],cwd=w,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=allowed).returncode
  except subprocess.TimeoutExpired:code='timeout'
 save(w/'execution.json',{'condition':name,'postinit_distance_shift':CONDITIONS[name],'returncode':code,'wall_seconds':time.monotonic()-start,'timeout_seconds':allowed,'earlier_native_seconds':prior,'stage_seconds_before':spent,'protocol_sha256':sha(P/'state-protocol.json'),'freeze_sha256':sha(P/'state-freeze.json'),'release_sha256':sha(P/'execution-release.json'),'binary_sha256':sha(P/'build/bin/snlc_fit.exe')});assert code==0,code
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('action',choices=['prepare']+list(CONDITIONS));x=a.parse_args();prepare() if x.action=='prepare' else run(x.action)
