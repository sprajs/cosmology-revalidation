from pathlib import Path
import json,hashlib,shutil,difflib,ast,os,tempfile
R=Path('/home/szymon/Documents/ChatGPT/supernova');P=R/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering';S=R/'phase2/official/inputs/SNDATA_ROOT'
O=P/'fit-private-lookup';O.mkdir(exist_ok=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
root=O/'SNDATA_ROOT';root.mkdir();(root/'SIM').mkdir();registry=root/'SIM/PATH_SNDATA_SIM.LIST';registry.write_text('')
refs=[]
for q in sorted(S.iterdir()):
 if q.name=='SIM':continue
 target=root/q.name;target.symlink_to(q.resolve(),target_is_directory=q.is_dir());refs.append({'path':str(target.relative_to(R)),'target':str(q.resolve()),'same_inode_target':target.resolve()==q.resolve(),'file_sha256':sha(q) if q.is_file() else None})
shutil.copyfile(S/'SIM/PATH_SNDATA_SIM.LIST',O/'shared-registry-observed.txt')
shared_hash=sha(S/'SIM/PATH_SNDATA_SIM.LIST')
# Match the C search order, rather than canonicalizing duplicate spellings away.
def candidates(datadir,registry_paths,sndata):
 return [datadir/'PTE.LIST',datadir/'PTE/PTE.LIST']+[x/'PTE/PTE.LIST' for x in registry_paths]+[sndata/'SIM/PTE/PTE.LIST']
preflight=[]
for branch in ['noiseless','ledger']:
 ds=P/'generation-v5'/branch/'output';cs=candidates(ds,[],root);found=[q for q in cs if q.is_file()];assert found==[ds/'PTE/PTE.LIST']
 names=(found[0]).read_text().split();assert names
 for name in names:assert (found[0].parent/name).exists(),name
 preflight.append({'dataset':str(ds.relative_to(R)),'candidate_paths':[str(q.relative_to(R)) for q in cs],'unique_found':str(found[0].relative_to(R)),'listed_files':{x:sha(found[0].parent/x) for x in names},'list_sha256':sha(found[0])})
# Test future adapter lookup structure with copied LIST only; no HEAD adaptation or photons.
probe=O/'adapter-layout-preflight';probe.mkdir()
for kind in ['estimated','null']:
 ds=probe/kind;d=ds/'PTE';d.mkdir(parents=True);src=P/'generation-v5/ledger/output/PTE'
 shutil.copyfile(src/'PTE.LIST',d/'PTE.LIST')
 for name in (d/'PTE.LIST').read_text().split():(d/name).symlink_to((src/name).resolve())
 cs=candidates(ds,[],root);found=[q for q in cs if q.is_file()];assert found==[d/'PTE.LIST'];preflight.append({'future_adapter_layout':kind,'separate_test_folder_only':True,'unique_list':str(found[0].relative_to(R)),'no_scientific_HEAD_adapter_created':True})
short=R/'phase2/pte/fit-private-lookup/SNDATA_ROOT';assert len(str(short))<120
s=(P/'run_engineering_v5.py').read_text();n=s.replace("'freeze-v5.json'","'freeze-fit-private.json'")
a="env.update(SNANA_DIR=str(snana),SNDATA_ROOT=str(R/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(R/'phase2/official/build/sysroot/usr/lib'))"
b="env.update(SNANA_DIR=str(snana),SNDATA_ROOT=str(SHORT/'fit-private-lookup/SNDATA_ROOT') if binary.name=='snlc_fit.exe' else str(R/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(R/'phase2/official/build/sysroot/usr/lib'))"
assert a in n;n=n.replace(a,b).replace("P/'fits'","P/'fits-private'").replace("P/'fits/NIR_true12'","P/'fits-private/NIR_true12'").replace("'-failure-v5.json'","'-failure-fit-private.json'").replace("choices=['generation','noiseless','noisy']","choices=['noiseless','noisy']")
ast.parse(n);(P/'run_engineering_v5_fit_private.py').write_text(n)
(O/'executor.patch').write_text(''.join(difflib.unified_diff(s.splitlines(True),n.splitlines(True),fromfile='run_engineering_v5.py',tofile='run_engineering_v5_fit_private.py')))
assert not (P/'fits-private').exists()
activity=json.loads((P/'native-activity-v5.json').read_text());shutil.copyfile(P/'native-activity-v5.json',O/'activity-at-amendment.json')
for name in ['noiseless-failure-v5.json']:
 shutil.copyfile(P/name,O/name)
for name in ['native.log','execution.json','fit.nml','input-freeze.json']:
 q=P/'fits/noiseless_joint_original'/name
 if q.exists():shutil.copyfile(q,O/('failed-'+name))
proto={'status':'Administrative lookup amendment after pre-data-fit failure; no native fits executed in preparation','scientific_protocol_sha256':sha(P/'protocol-v5.json'),'failure':'Global simulation-path registry is searched even with a private dataset path; six existing PTE LIST candidates were counted. The failed fit did not load/fit data. Original failed fit folder and failure record remain untouched.','source':{'file':'runs/research_2026_09_26/astra_design/raisin_timing_assets/snana_v11_04d/source/src/sntools.c','private_plus_global_search_lines':[5291,5375],'count_existing_LIST_paths_lines':[5390,5420],'registry_write_lines':[5471,5519],'simulator_registers_explicit_path':'snlc_sim.c:5071-5074'},'changes':['Fit processes only use private SNDATA_ROOT with identical non-SIM reference targets and an empty private SIM/PATH_SNDATA_SIM.LIST.','Fit work folders move from fits/ to fits-private/ at the same relative depth.','Versioned executor reads freeze-fit-private.json and uses its own failure filename; only noiseless/noisy stages accepted.'], 'unchanged':'Native fitter source/binaries, NML template, generated HEAD/PHOT, row order, numerical/physical gates, seed, scientific protocol and activity ledger. No shared-registry mutation and no generation rerun.','reference_links':refs,'private_registry_sha256':sha(registry),'shared_registry_snapshot_sha256':shared_hash,'preflight':preflight,'runtime_SNDATA_ROOT':str(short),'runtime_path_length':len(str(short)),'budget':{'carried_native_seconds':sum(x['wall_seconds'] for x in activity),'cap_seconds':120,'activity_ledger':'native-activity-v5.json'},'release':'A new root release binds original protocol-v5.json SHA and new freeze-fit-private.json SHA; stages noiseless/noisy only. Source-aligned path isolation is administrative, not a changed estimator.'}
save(O/'amendment.json',proto)
shutil.copyfile('/tmp/prepare_private_fit_lookup.py',O/'prepare.py')
f=json.loads((P/'freeze-v5.json').read_text());f['status']='Additive fit-only lookup isolation before native retry';f['authoritative_executor']='run_engineering_v5_fit_private.py'
paths=[P/'run_engineering_v5_fit_private.py']+[q for q in O.rglob('*') if q.is_file() and not q.is_symlink()]
for q in paths:f['files'][str(q.relative_to(R))]=sha(q)
f['external_execution_release_required']['freeze_sha256']='SHA256 of freeze-fit-private.json';f['external_execution_release_required']['stages']=['noiseless','noisy']
save(P/'freeze-fit-private.json',f)
assert sha(S/'SIM/PATH_SNDATA_SIM.LIST')==shared_hash
for name,h in f['files'].items():assert sha(R/name)==h,name
print(json.dumps({'amendment_sha256':sha(O/'amendment.json'),'executor_sha256':sha(P/'run_engineering_v5_fit_private.py'),'freeze_sha256':sha(P/'freeze-fit-private.json'),'files':len(f['files']),'native_seconds_carried':proto['budget']['carried_native_seconds'],'SNDATA_ROOT_length':len(str(short))},indent=2))
