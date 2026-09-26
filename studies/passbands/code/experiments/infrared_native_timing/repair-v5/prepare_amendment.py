from pathlib import Path
import hashlib,json,difflib,shutil,ast
R=Path('/home/szymon/Documents/ChatGPT/supernova')
P=R/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering'
V=P/'repair-v5'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def patch(a,b,fa,fb):return ''.join(difflib.unified_diff(a.splitlines(True),b.splitlines(True),fromfile=fa,tofile=fb))
build=json.loads((V/'build-result.json').read_text());assert len(build)==2 and all(x['returncode']==0 for x in build)
# A byte-exact inverse ensures that the instrumentation was not revised with the repair.
checks=[]
for name in ['original','ledger']:
 old=P/'instrumentation'/('snlc_sim.original.c' if name=='original' else 'snlc_sim.instrumented.c')
 new=V/(name+'.c');s=new.read_text()
 s=s.replace(',SNRMAX_FILT[MXCUTWIN_SNRMAX+1][MXFILTINDX]',',SNRMAX_FILT[MXCUTWIN_SNRMAX][MXFILTINDX]')
 s=s.replace('MEM=NEP; if (NEP < MXFILTINDX) { MEM=MXFILTINDX; }\n  MEM++; // one-based TLIST/INDEX_SORT require indices 1..MEM-1','MEM=NEP; if (NEP < MXFILTINDX) { MEM=MXFILTINDX+1; }')
 assert s==old.read_text(),name
 assert sha(new)==sha(V/(name+'-build/src/snlc_sim.c'))
 checks.append({'variant':name,'undo_exact_two_bounds_repairs_recovers_previous_source':True,'build_source_exact':True,'binary_sha256':sha(V/(name+'-build/bin/snlc_sim.exe'))})
save(V/'source-identity-check.json',checks)
# Only PATH_SNDATA_SIM changes between the preserved v4 and new v5 input files.
inputs=[]
for branch in ['original','ledger','noiseless']:
 old=P/'generation-v4'/branch/'sim.input';new=P/'generation-v5'/branch/'sim.input';new.parent.mkdir(parents=True,exist_ok=False)
 s=old.read_text();n=s.replace('generation-v4/','generation-v5/');assert n.count('generation-v5/')==1
 new.write_text(n);inputs.append({'branch':branch,'old_sha256':sha(old),'new_sha256':sha(new),'only_output_directory_changed':True})
 (V/(branch+'-input.patch')).write_text(patch(s,n,str(old.relative_to(R)),str(new.relative_to(R))))
old=P/'run_engineering_v4.py';s=old.read_text();n=s.replace('protocol-v4.json','protocol-v5.json').replace('freeze-v4.json','freeze-v5.json').replace('native-activity-v4.json','native-activity-v5.json').replace('generation-v4','generation-v5').replace('-failure-v4.json','-failure-v5.json')
a="  w=P/'generation-v5'/branch;binary=SHORT/('instrumentation/snlc_sim.original.exe' if branch=='original' else 'build-v11_04d/bin/snlc_sim.exe')\n  native('generate_'+branch,w,binary,SHORT/'build-v11_04d',['sim.input'],ledger=branch!='original')"
b="  w=P/'generation-v5'/branch;snana=SHORT/'repair-v5'/('original-build' if branch=='original' else 'ledger-build');binary=snana/'bin/snlc_sim.exe'\n  native('generate_'+branch,w,binary,snana,['sim.input'],ledger=branch!='original')"
assert a in n;n=n.replace(a,b);ast.parse(n);(P/'run_engineering_v5.py').write_text(n)
(V/'executor.patch').write_text(patch(s,n,str(old.relative_to(R)),str((P/'run_engineering_v5.py').relative_to(R))))
prior=json.loads((P/'native-activity-v4.json').read_text());assert len(prior)==2 and [x['returncode'] for x in prior]==[11,-6]
spent=sum(x['wall_seconds'] for x in prior);shutil.copyfile(P/'native-activity-v4.json',P/'native-activity-v5.json')
amend={'status':'Frozen after v4 crash diagnosis and two private builds; before any repaired-generator execution',
 'previous_protocol_sha256':sha(P/'protocol-v4.json'),'previous_freeze_sha256':sha(P/'freeze-v4.json'),
 'source_identity':'SNANA v11_04d 10ec91297e4482d593cb5d3d055b10d4aa915071 plus exactly two gen_cutwin array-extent repairs; not unpatched historical generator reproduction',
 'repair_build_protocol_sha256':sha(V/'build-protocol.json'),'repair_build_result_sha256':sha(V/'build-result.json'),
 'repairs':['Allocate SNRMAX_FILT[MXCUTWIN_SNRMAX+1][MXFILTINDX] for the existing inclusive 0..MXCUTWIN_SNRMAX convention.','Allocate TLIST and INDEX_SORT with max(NEP,MXFILTINDX)+1 elements for existing one-based 1..NEP/1..NFILT accesses.'],
 'scope':'No mean, noise, RNG, seed, cadence, truth, cut statistic, optimizer, scientific gate or selection setting changes. The repaired comparison and ledger generators have identical two repairs; output-only instrumentation is byte-identical after undoing these repairs.',
 'unexercised_defect':'Global INPUTS.CUTWIN_SNRMAX family has its own maximum-user-cut index defect. Default NCUTWIN_SNRMAX=0 and absence of any CUTWIN_SNRMAX input remain frozen. No global repair or user-cut experiment is performed.',
 'failure_provenance':{'v3':'IDLOCK setup failed before a native light-curve mean/noise event; zero events written.','v4':'One internal event reached GENFLUX and gen_cutwin. Stack canary aborted after mean/noise generation; NGENLC_TOT=1, NGENFLUX_DRIVER=1, NGENLC_WRITE=0 in the exact core. No FITS event or science fit outcome was written or selected. This was not a pre-photon failure.','replay':'Same frozen seed replays the intended eight engineering realizations. The failed internal first realization is not a ninth independent draw; no reseeding, new draws or outcome-selected membership.'},
 'inputs':inputs,'binary':{x['variant']:x['binary_sha256'] for x in build},
 'resources':{'total_native_cap_seconds':120,'prior_native_seconds':spent,'remaining_native_seconds':120-spent,'prior_native_ledger_sha256':sha(P/'native-activity-v4.json'),'total_generator_build_cap_seconds':180,'cumulative_generator_build_seconds':build[-1]['total_generator_build_seconds']},
 'execution':'Requires a new root release binding protocol-v5.json and freeze-v5.json. No repaired-generator execution by this preparation.'}
save(V/'amendment.json',amend)
proto=json.loads((P/'protocol-v4.json').read_text())
proto['status']='Additive two-bounds-repair amendment after post-noise v4 stack crash; new root release required'
proto['source_identity']=amend['source_identity']
proto['generation']['extra_controls']='One replay of the same eight attempts with the identically bounds-repaired comparison generator for exact output-hook identity, plus one SMEARFLAG_FLUX0 deterministic mean control. Planned completed native attempts17 including replay; nine unique engineered records, eight stochastic. The v4 aborted internal first event is disclosed and replayed, not counted as an extra independent draw.'
proto['generation']['v4_failure']=amend['failure_provenance']['v4']
proto['generation']['v5_replay']=amend['failure_provenance']['replay']
proto['generation']['source_repairs']=amend['repairs']
proto['generation']['unpatched_historical_generator']='Unavailable on this compiler because native gen_cutwin invokes out-of-bounds writes. The prospective corrected-source comparison is not exact unpatched historical reproduction.'
proto['generation']['NCUTWIN_SNRMAX']=0
proto['gates']['instrumentation']='All native generator HEAD/PHOT table columns exact between two identically bounds-repaired variants differing only by the preserved output hooks; original vs support-fitter FITRES science rows and every saved mean/W/objective exact.'
proto['resources']['builds_separate_seconds']['generator_total']=build[-1]['total_generator_build_seconds']
proto['resources']['earlier_native_seconds']=spent;proto['resources']['remaining_native_seconds']=120-spent
proto['resources']['activity']='native-activity-v5.json begins with both unchanged v3/v4 failure records; overall120s cap unchanged.'
proto['authoritative_inputs']['generation']='generation-v5/*/sim.input';proto['authoritative_inputs']['executor']='run_engineering_v5.py';proto['authoritative_inputs']['source_amendment']='repair-v5/amendment.json'
proto['authoritative_inputs']['generator_comparison']='repair-v5/original-build/bin/snlc_sim.exe';proto['authoritative_inputs']['generator_ledger']='repair-v5/ledger-build/bin/snlc_sim.exe'
save(P/'protocol-v5.json',proto)
# Preserve the existing authority chain, and add current versioned active paths and source/build evidence.
f=json.loads((P/'freeze-v4.json').read_text());files=dict(f['files'])
paths=[P/'protocol-v5.json',P/'run_engineering_v5.py',P/'protocol-v4.json',P/'freeze-v4.json',P/'native-activity-v4.json',P/'crash-v4-review/result.json',P/'crash-v4-review/gdb.txt',P/'crash-v4-review/core-generation-counts.txt']
for branch in ['original','ledger','noiseless']:paths.append(P/'generation-v5'/branch/'sim.input')
for name in ['build-protocol.json','build-result.json','amendment.json','source-identity-check.json','original-bounds.patch','ledger-bounds.patch','ledger-minus-original.patch','original.c','ledger.c','original-build.sh','ledger-build.sh','executor.patch']:paths.append(V/name)
for variant in ['original','ledger']:paths.extend([V/(variant+'-build/src/snlc_sim.c'),V/(variant+'-build/bin/snlc_sim.exe')])
for p in paths:files[str(p.relative_to(R))]=sha(p)
f['status']='Frozen after disclosed internal first-event crash; before repaired v5 replay'
f['protocol_sha256']=sha(P/'protocol-v5.json');f['files']=files;f['authoritative_executor']='run_engineering_v5.py'
f['external_execution_release_required']={'protocol_sha256':f['protocol_sha256'],'freeze_sha256':'SHA256 of freeze-v5.json','stages':['generation','noiseless','noisy']}
save(P/'freeze-v5.json',f)
print(json.dumps({'protocol_sha256':sha(P/'protocol-v5.json'),'executor_sha256':sha(P/'run_engineering_v5.py'),'freeze_sha256':sha(P/'freeze-v5.json'),'amendment_sha256':sha(V/'amendment.json'),'files':len(files),'prior_native_seconds':spent},indent=2))
