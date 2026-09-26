from pathlib import Path
import json,hashlib,shutil,os,subprocess,time,re
import numpy as np
from parse_native import parse
ROOT=Path.cwd();O=Path(__file__).resolve().parent;B=O/'SNANA-v11_04k-output';P=O/'probes-v3';S=ROOT/'runs/research_2026_09_26/csp_native_filter_response/cohort-preparation/pilot';sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert not P.exists();P.mkdir();(P/'data').symlink_to(S/'data',target_is_directory=True)
base=(S/'fits/nominal/fit.nml').read_text();ref=parse(O/'pilot-v2/fits/nominal/fit.log');last={b['CID']:b for b in ref};arms=[]
for name,shift in [('nominal',None),('start_minus',-.2),('start_plus',.2)]:arms.append(dict(name=name,nml=base,shift=shift,type='fitted'))
for cid,b in last.items():
 for label,offset in [('center',0),('minus',-.15),('plus',.15)]:
  dseed=float(np.float32(b['objective']['D']+offset));nml=re.sub(r"(?m)^(\s*SNCCID_LIST\s*=).*$",lambda m:m[1]+f" '{cid}'",base)
  nml=nml.replace('&FITINP',f'&FITINP\n LFIXPAR_ALL = T\n INIVAL_DLMAG = {dseed:.17g}\n INISTP_DLMAG = 0.0')
  arms.append(dict(name=f'{cid}_{label}',nml=nml,shift=None,type='fixed',CID=cid,requested_offset=offset,D_seed_float32=dseed,reference_D=b['objective']['D']))
inputs=[Path(__file__),O/'parse_native.py',B/'bin/snlc_fit.exe',B/'src/snlc_fit.car',O/'start-hook-v3-protocol.json',O/'pilot-v2/analysis-v2/result.json']
for a in arms:
 d=P/'fits'/a['name'];d.mkdir(parents=True);(d/'fit.nml').write_text(a.pop('nml'));shutil.copyfile(S/'fits/nominal/vpec.list',d/'vpec.list');inputs += [d/'fit.nml',d/'vpec.list']
proto={'scope':'Nominal-only numerical gates; no filter relabel. Default hook-off equivalence first; then first-iteration native post-init shifts and exact fixed native amplitude probes.','arms':arms,'hashes':{str(p.relative_to(ROOT)):sha(p) for p in inputs},'per_invocation_timeout_seconds':30,'total_invocation_cap_seconds':120,'gates':{'default_science_and_LCPLOT':'byte-identical to output-v2 and original','actual_first_entry_shifts':'native entry differences equal ±.2 within1e−12','final_D_multistart':.001,'fixed_native_mean_scaling_relative':2e-8,'all_rows':'exact band,MJD,y,error accepted identity per matched iteration','Qclosure_absolute':1e-8,'support':'all callback phases within source grid−20..70, shape inside.7..1.3; final fitwindow−15..45; rest throughput endpoints inspected separately; no clamp'}}
(P/'protocol.json').write_text(json.dumps(proto,indent=2)+'\n');print('frozen',sha(P/'protocol.json'),flush=True)
env=os.environ.copy();env.update(SNANA_DIR='/tmp/snana-csp-audit-8c5f0d',SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1');env.pop('CSP_POSTINIT_DSHIFT',None)
used=0;results=[]
for a in arms:
 d=P/'fits'/a['name'];e=env.copy()
 if a['shift'] is not None:e['CSP_POSTINIT_DSHIFT']=str(a['shift'])
 t=time.monotonic()
 with (d/'fit.log').open('x') as f:r=subprocess.run([str(B/'bin/snlc_fit.exe'),'fit.nml'],cwd=d,env=e,stdout=f,stderr=subprocess.STDOUT,timeout=min(30,120-used))
 elapsed=time.monotonic()-t;used+=elapsed;x=dict(name=a['name'],returncode=r.returncode,seconds=elapsed,log_sha256=sha(d/'fit.log'));results.append(x);(d/'execution.json').write_text(json.dumps(x,indent=2)+'\n');(P/'execution-results.json').write_text(json.dumps({'arms':results,'total_seconds':used},indent=2)+'\n');print(x,flush=True)
 if r.returncode:raise RuntimeError('Native failure preserved; later arms not run')
 if a['name']=='nominal':
  rd=O/'pilot-v2/fits/nominal';science=lambda p:[x for x in p.read_text().splitlines() if x.startswith(('SN:','VARNAMES:'))]
  assert science(d/'fit.FITRES.TEXT')==science(rd/'fit.FITRES.TEXT')
  assert (d/'fit.LCPLOT.TEXT').read_bytes()==(rd/'fit.LCPLOT.TEXT').read_bytes()
