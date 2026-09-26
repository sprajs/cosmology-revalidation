from pathlib import Path
import json,hashlib,shutil,ast,difflib,subprocess
R=Path('/home/szymon/Documents/ChatGPT/supernova');P=R/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering';O=P/'fit-readme-adapter';O.mkdir(exist_ok=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
# Relevant exact crash metadata only; no core extraction or native rerun.
z=subprocess.run(['coredumpctl','--no-pager','info','490235'],capture_output=True,text=True);(O/'coredump-info.txt').write_text(z.stdout+z.stderr)
assert z.returncode==0 and '2026-09-26 08:36:29' in z.stdout and 'check_file_docana' in z.stdout
rr=subprocess.run(['journalctl','-k','--since','2026-09-26 08:36:00','--until','2026-09-26 08:37:00','--no-pager'],capture_output=True,text=True)
(O/'kernel-interval.txt').write_text(rr.stdout+rr.stderr)
assert not any(k in rr.stdout.lower() for k in ['out of memory','oom-kill','killed process'])
source=R/'runs/research_2026_09_26/astra_design/raisin_timing_assets/snana_v11_04d/source/src';c=(source/'sntools.c').read_text();f=(source/'snana.car').read_text()
assert 'if ( FIRSTLINE && nline > 5 ) { break; }' in c
records=[]
for branch in ['noiseless','ledger']:
 src=P/'generation-v5'/branch/'output/PTE';dst=P/'readme-adapted'/branch/'PTE';dst.mkdir(parents=True,exist_ok=False)
 for q in sorted(src.iterdir()):
  if q.name=='PTE.README':
   old=q.read_bytes();assert old.splitlines()[0].strip()==b'DOCUMENTATION:'
   newline=old.index(b'\n')+1;new=old[:newline]+b'\n'+old[newline:];(dst/q.name).write_bytes(new)
   assert new.replace(b'\n\n',b'\n',1)==old and old.split()==new.split()
   before=[len(x) for x in b'\n'.join(old.splitlines()[:6]).split()];after=[len(x) for x in b'\n'.join(new.splitlines()[:6]).split()]
   assert max(before)>=60 and max(after)<60
   assert b'DOCUMENTATION:' in b'\n'.join(new.splitlines()[:6]).split()
   records.append({'branch':branch,'source_readme':str(q.relative_to(R)),'source_sha256':sha(q),'adapted_readme':str((dst/q.name).relative_to(R)),'adapted_sha256':sha(dst/q.name),'inserted_bytes':1,'only_change':'one blank line after DOCUMENTATION:','all_nonwhitespace_tokens_identical':True,'max_token_length_first_six_before':max(before),'max_token_length_first_six_after':max(after)})
  elif q.is_file():(dst/q.name).symlink_to(q.resolve())
 for name in ['PTE.LIST','PTE_HEAD.FITS','PTE_PHOT.FITS']:
  assert (dst/name).resolve()==(src/name).resolve()
 # C search has private path, private/PTE, empty private registry, default privateSIM.
 parent=dst.parent;privateSIM=P/'fit-private-lookup/SNDATA_ROOT/SIM'
 candidates=[parent/'PTE.LIST',dst/'PTE.LIST',privateSIM/'PTE/PTE.LIST'];assert [x for x in candidates if x.is_file()]==[dst/'PTE.LIST']
 # Parse DOCANA sub-block using available YAML parser; post-DOCANA legacy text is unchanged.
 import yaml
 a=old.split(b'DOCUMENTATION_END:')[0];b=new.split(b'DOCUMENTATION_END:')[0];assert yaml.safe_load(a)==yaml.safe_load(b)
 (O/(branch+'-readme.patch')).write_text(''.join(difflib.unified_diff(old.decode().splitlines(True),new.decode().splitlines(True),fromfile=branch+'/original.README',tofile=branch+'/adapter.README')))
s=(P/'run_engineering_v5_fit_private.py').read_text();n=s.replace("'freeze-fit-private.json'","'freeze-fit-readme.json'").replace("P/'fits-private'","P/'fits-readme'").replace("P/'fits-private/NIR_true12'","P/'fits-readme/NIR_true12'").replace("'-failure-fit-private.json'","'-failure-fit-readme.json'")
n=n.replace("d=P/'generation-v5'/branch/'output/PTE';heads=", "d=(P/'readme-adapted'/branch/'PTE') if branch in ['ledger','noiseless'] else (P/'generation-v5'/branch/'output/PTE');heads=")
n=n.replace("data='../../generation-v5/noiseless/output'","data='../../readme-adapted/noiseless'").replace("data='../../generation-v5/ledger/output'","data='../../readme-adapted/ledger'")
assert n!=s;ast.parse(n);(P/'run_engineering_v5_fit_readme.py').write_text(n)
(O/'executor.patch').write_text(''.join(difflib.unified_diff(s.splitlines(True),n.splitlines(True),fromfile='run_engineering_v5_fit_private.py',tofile='run_engineering_v5_fit_readme.py')))
shutil.copyfile('/tmp/prepare_readme_adapter.py',O/'prepare.py')
activity=json.loads((P/'native-activity-v5.json').read_text());shutil.copyfile(P/'native-activity-v5.json',O/'activity-at-amendment.json')
for q in [P/'noiseless-failure-fit-private.json']+list((P/'fits-private/noiseless_joint_original').glob('*')):
 if q.is_file():shutil.copyfile(q,O/('failed-'+q.name))
proto={'status':'Metadata-only adapter frozen before retry; no native fitting or generation in preparation','original_scientific_protocol_sha256':sha(P/'protocol-v5.json'),'cause':{'exact_pid':490235,'function':'check_file_docana','source':'sntools.c:9697-9710 uses key[60]; get_PARSE_WORD copies each token from store_PARSE_WORDS FIRSTLINE. FIRSTLINE reads six physical lines, because nline is incremented and parsed before testing nline>5 at2063.','trigger':'Generated README line6 has SNDATA_ROOT path length75. Copying this token exceeds key[60], consistent with exact process stack-canary abort. No data/model fit began.','limits':'Core metadata and native function backtrace plus source/trigger correspondence; no new core extraction or source-level local-variable claims.'},'adapter':records,'consumer_audit':{'getinfo_photometry':'snana.car:3548 calls CHECK_FILE_DOCANA; the file supplies documentary validation, not fitted data.','other_native_snana_consumers':'DUMP_README cats the full file to stdout (call in read-data path commented at5430); reformat close concatenates original README. Neither interprets numeric photometry from this file. sntools_dataformat_fits only obtains its name in getInfo_PHOTOMETRY_VERSION. No fitter README numerical consumer found.','yaml':'Parsed DOCUMENTATION block before DOCUMENTATION_END is exactly equal; all whitespace-separated tokens in complete file are exactly equal. Original path provenance is retained, not abbreviated or invented.'},'changes':['Symlink all generated data/list assets into readme-adapted/{noiseless,ledger}/PTE and copy only README with exactly one extra newline after first line.','Use same-depth fits-readme/ work folders; versioned executor points to these metadata adapters and inherits private empty SIM registry.','Future estimated/null HEAD adapters copy the safe documentary README while retaining original PHOT file; only HEAD PEAKMJD intervention remains the frozen scientific change.'],'unchanged':'Both historical fitter binaries and instrumented support binary, all means/model/noise/covariance settings, HEAD/PHOT/LIST bytes, template, masks, controls, thresholds, scientific protocol, all previous failures and activity ledger. No shared file is mutated.','budget':{'native_total_seconds':120,'carried_seconds':sum(x['wall_seconds'] for x in activity)},'release':'A new root release binds unchanged protocol-v5.json and freeze-fit-readme.json; only noiseless/noisy stages. This amendment does not authorize execution itself.'}
save(O/'amendment.json',proto)
frozen=json.loads((P/'freeze-fit-private.json').read_text());frozen['status']='Metadata-only README adapter before fit retry';frozen['authoritative_executor']='run_engineering_v5_fit_readme.py';frozen['external_execution_release_required']['freeze_sha256']='SHA256 of freeze-fit-readme.json'
paths=[P/'run_engineering_v5_fit_readme.py']+[x for x in O.rglob('*') if x.is_file()]+[x for x in (P/'readme-adapted').rglob('*') if x.is_file()]
for q in paths:frozen['files'][str(q.relative_to(R))]=sha(q)
save(P/'freeze-fit-readme.json',frozen)
for name,h in frozen['files'].items():assert sha(R/name)==h,name
print(json.dumps({'amendment_sha256':sha(O/'amendment.json'),'executor_sha256':sha(P/'run_engineering_v5_fit_readme.py'),'freeze_sha256':sha(P/'freeze-fit-readme.json'),'files':len(frozen['files']),'carried_native_seconds':proto['budget']['carried_seconds'],'readme_records':records},indent=2))
