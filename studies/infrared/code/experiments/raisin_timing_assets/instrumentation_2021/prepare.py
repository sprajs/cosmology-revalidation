from pathlib import Path
import shutil,subprocess,json,hashlib,difflib,datetime
R=Path('/home/szymon/Documents/ChatGPT/supernova');P=Path(__file__).resolve().parent;V=P.parent/'snana_v11_04d';I=R/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/instrumentation'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
assert not (P/'build').exists()
shutil.copytree(V/'source',P/'build')
logs=[]
for f in ['output-only-v2.patch','start-hook-v3.patch']:
 shutil.copy2(I/f,P/f)
 r=subprocess.run(['patch','--fuzz=0','-p1','--input',str(P/f)],cwd=P/'build',capture_output=True,text=True)
 logs.append({'patch':f,'sha256':sha(P/f),'exit':r.returncode,'stdout':r.stdout,'stderr':r.stderr})
 save(P/'patch-application.json',logs)
 assert r.returncode==0,logs
old=(V/'source/src/snlc_fit.car').read_text();new=(P/'build/src/snlc_fit.car').read_text()
(P/'ported-source.patch').write_text(''.join(difflib.unified_diff(old.splitlines(keepends=True),new.splitlines(keepends=True),fromfile='v11_04d/src/snlc_fit.car',tofile='instrumented_v11_04d/src/snlc_fit.car')))
shutil.copy2(P/'build/src/snlc_fit.car',P/'instrumented-snlc_fit.car')
shutil.copy2(I/'parse_native.py',P/'parse_native.py')
shell=(V/'build.sh').read_text().replace(str(V/'build'),str(P/'build')).replace('raisin_timing_assets/snana_v11_04d/build','raisin_timing_assets/instrumentation_2021/build')
(P/'build.sh').write_text(shell)
shutil.copy2(V/'run_build.py',P/'run_build.py')
changed=[str(f.relative_to(V/'source')) for f in (V/'source').rglob('*') if f.is_file() and sha(f)!=sha(P/'build'/f.relative_to(V/'source'))]
assert changed==['src/snlc_fit.car'],changed
assert len(str(P/'build'))<160
save(P/'build-protocol.json',{'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Private v11_04d output-only state port and absent/default-zero postinitialization distance-shift hook; no native fits authorized by this build protocol.','source_commit':'10ec91297e4482d593cb5d3d055b10d4aa915071','allowed_scientific_state_change':'Only explicitly nonempty CSP_POSTINIT_DSHIFT diagnostic environment modifies first-iteration D initialization after FITPAR_PREP. Absent variable changes no scientific state.','port_method':'Both previously independently reviewed v11_04k patches apply to v11_04d with zero fuzz; offsets recorded. No model function is called again by diagnostics.','full_state':'Accepted rows, exact model means/data/errors, active W or diagonal weights, total/prior/sigma objectives, entering and previous coordinates, peak prior center, fixed state and MJDOFF.','build_timeout_seconds':180,'no_shared_binary_changes':True,'prebuild_changed_files':changed,'inputs':{str(p.relative_to(R)):sha(p) for p in [V/'source/src/snlc_fit.car',I/'output-only-v2.patch',I/'start-hook-v3.patch',I/'parse_native.py',P/'ported-source.patch',P/'instrumented-snlc_fit.car',P/'build.sh',P/'run_build.py',P/'prepare.py']},'future_required_gates':['Default absent hook science FITRES rows and LCPLOT bytes exactly match passing baseline2021.','Explicit zero hook same results.','Only actual post-initialization starts ±0.2 mag; verify logged entry before minimization.','Full-state objective reconstruction and covariance positive definiteness, masks and coordinates checked before any timing interpretation.']})
print(json.dumps({'patches':logs,'build_protocol_sha256':sha(P/'build-protocol.json'),'SNANA_DIR_length':len(str(P/'build'))},indent=2))
