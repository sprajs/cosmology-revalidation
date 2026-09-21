#!/usr/bin/env python3
"""Run one immutable forward input or record an already completed manual run."""
from pathlib import Path
import argparse,hashlib,json,os,subprocess,datetime
R=Path(__file__).resolve().parents[3];O=R/'phase2/literature/simulations'
p=argparse.ArgumentParser();p.add_argument('model');p.add_argument('--tag',default='pilot02');p.add_argument('--record-existing',action='store_true');a=p.parse_args()
version=f'PH2_{a.tag}_{a.model}';inp=O/'inputs'/(version+'.input');output=O/'outputs'/version;log=O/'logs'/(version+'-attempt01.log');manifest=O/'manifests'/(version+'.json')
exe=R/'phase2/official/build/SNANA-2fe0f56/bin/snlc_sim.exe'
def sha(f):
 h=hashlib.sha256()
 with f.open('rb') as stream:
  for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
 return h.hexdigest()
def obj(f):return {'path':str(f.relative_to(R)),'sha256':sha(f),'bytes':f.stat().st_size}
if manifest.exists():raise RuntimeError('Refusing duplicate run manifest')
if not inp.exists():raise FileNotFoundError(inp)
env=os.environ.copy();env.update(SNANA_DIR=str(exe.parents[1]),SNDATA_ROOT=str(R/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(R/'phase2/official/build/sysroot/usr/lib'))
record={'version':version,'record_time_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'input':obj(inp),'executable':obj(exe),'source_commit':'2fe0f564361a873860661ff61080b4db9c607edf','source_tree_audit':'phase2/official/diagnostics/historical-tree-audit.json','environment':{k:env[k] for k in ['SNANA_DIR','SNDATA_ROOT','LD_LIBRARY_PATH']},'command':[str(exe),str(inp)],'acquisition_archive':'phase2/official/inputs/SNDATA_ROOT_2024-07-03.tar.gz','constant_alpha_compatibility':'phase2/literature/simulations/inputs/constant-alpha-compatibility.json'}
if not a.record_existing:
 if output.exists() or log.exists():raise RuntimeError('Refusing to overwrite existing output/log')
 with log.open('w') as stream:record['exit_code']=subprocess.run(record['command'],env=env,cwd=R,stdout=stream,stderr=subprocess.STDOUT).returncode
else:record['exit_code']=None;record['adopted_manual_run']=True
record['log']=obj(log);record['outputs']=[obj(f) for f in sorted(output.glob('*')) if f.is_file()]
record['success_marker']= 'SUMMARY' in log.read_text() and any(x['path'].endswith('.README') for x in record['outputs'])
manifest.parent.mkdir(exist_ok=True);manifest.write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'manifest':str(manifest),'exit_code':record['exit_code'],'output_files':len(record['outputs'])}))
