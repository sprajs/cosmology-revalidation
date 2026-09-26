from pathlib import Path
import json,hashlib,os,tempfile,shutil
R=Path('/home/szymon/Documents/ChatGPT/supernova');P=R/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering'
A=P/'restart-v5-output-parent';A.mkdir(exist_ok=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
frozen=json.loads((P/'freeze-v5.json').read_text())
for f,h in frozen['files'].items():assert sha(R/f)==h,f
log=P/'generation-v5/original/native.log';text=log.read_text()
assert 'Begin Generating Lightcurves' not in text
assert 'mkdir: cannot create directory' in text and 'No such file or directory' in text
before={}
for branch in ['original','ledger','noiseless']:
 w=P/'generation-v5'/branch
 before[branch]=[str(x.relative_to(P)) for x in sorted(w.rglob('*'))]
 assert not (w/'output').exists(),branch
 assert not list(w.glob('*.FITS*')) and not list(w.rglob('*.DUMP')),branch
activity=json.loads((P/'native-activity-v5.json').read_text());assert len(activity)==3 and activity[-1]['returncode']==-11
assert activity[-1]['work'].endswith('generation-v5/original')
shutil.copyfile(P/'native-activity-v5.json',A/'native-activity-before-restart.json')
shutil.copyfile(P/'root-release-generation-v5.json',A/'root-release-used.json')
archive=A/'failed-attempt';archive.mkdir()
moves=[]
for old,dest in [(log,archive/'native.log'),(P/'generation-v5/original/execution.json',archive/'execution.json'),(P/'generation-failure-v5.json',archive/'generation-failure-v5.json')]:
 h=sha(old);old.rename(dest);assert sha(dest)==h
 moves.append({'old':str(old.relative_to(R)),'archived':str(dest.relative_to(R)),'sha256':h,'bytes':dest.stat().st_size})
checks=[]
for branch in ['original','ledger','noiseless']:
 w=P/'generation-v5'/branch;out=w/'output';out.mkdir()
 with tempfile.NamedTemporaryFile(prefix='.writable-preflight-',dir=out,delete=True) as q:q.write(b'writable\n');q.flush()
 inp=w/'sim.input';values={}
 for line in inp.read_text().splitlines():
  if ':' in line and not line.startswith('#'):
   k,v=line.split(':',1);values[k]=v.strip()
 assert Path(values['PATH_SNDATA_SIM']).resolve()==out.resolve()
 source={}
 for k in ['SIMLIB_FILE','GENMODEL','KCOR_FILE','FLUXERRMODEL_FILE']:
  q=(w/values[k]).resolve();assert q.exists() and os.access(q,os.R_OK),(branch,k)
  source[k]={'resolved':str(q),'readable':True,'sha256':sha(q) if q.is_file() else {p.name:sha(p) for p in sorted(q.iterdir()) if p.is_file()}}
 binary=P/'repair-v5'/('original-build' if branch=='original' else 'ledger-build')/'bin/snlc_sim.exe'
 assert os.access(binary,os.X_OK)
 assert not (out/'PTE').exists() and not (w/'native.log').exists() and not (w/'execution.json').exists()
 checks.append({'branch':branch,'input_sha256':sha(inp),'output_parent':str(out.relative_to(R)),'output_parent_writable':True,'native_version_subdirectory_absent':True,'binary_sha256':sha(binary),'source_inputs':source})
for f,h in frozen['files'].items():assert sha(R/f)==h,f
result={'status':'Administrative restart prepared only; no native execution','failed_attempt':'v5 output parent was absent; native init_simFiles uses mkdir without -p, before the main event loop. The log never reaches Begin Generating Lightcurves and no output directory, FITS or DUMP event exists. Source ordering establishes no GENFLUX_DRIVER invocation in this attempt.','source_order':{'source':'runs/research_2026_09_26/astra_design/raisin_timing_assets/snana_v11_04d/source/src/snlc_sim.c','init_simFiles_before_event_banner':[204,221],'mkdir_without_p':[27392,27395],'GENFLUX_DRIVER_in_event_loop':308},'frozen_input_files_reverified':len(frozen['files']),'input_binary_seed_gate_changes':False,'before_inventory':before,'archived_paths':moves,'preflight':checks,'native_seconds_carried':sum(x['wall_seconds'] for x in activity),'native_seconds_remaining':120-sum(x['wall_seconds'] for x in activity),'native_activity_not_reset_or_moved':True,'restart':'Invoke the unchanged run_engineering_v5.py generation with the same release only when root chooses to resume; the original executor can open fresh x-mode logs after the disclosed archive moves. No wrapper or recompilation is needed.'}
save(A/'preflight.json',result)
shutil.copyfile('/tmp/prepare_v5_directory_restart.py',A/'prepare_restart.py')
paths=[p for p in A.rglob('*') if p.is_file()]
save(A/'manifest.json',{'files':{str(p.relative_to(R)):sha(p) for p in paths},'unchanged_protocol_sha256':sha(P/'protocol-v5.json'),'unchanged_freeze_sha256':sha(P/'freeze-v5.json'),'unchanged_executor_sha256':sha(P/'run_engineering_v5.py')})
print(json.dumps({'preflight_sha256':sha(A/'preflight.json'),'manifest_sha256':sha(A/'manifest.json'),'native_seconds_carried':result['native_seconds_carried'],'remaining':result['native_seconds_remaining'],'moves':moves},indent=2))
