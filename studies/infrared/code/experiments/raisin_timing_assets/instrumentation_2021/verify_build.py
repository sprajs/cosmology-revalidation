from pathlib import Path
import json,hashlib,difflib,subprocess,os
P=Path(__file__).resolve().parent;V=P.parent/'snana_v11_04d';R=Path('/home/szymon/Documents/ChatGPT/supernova')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
a=json.loads((V/'acquisition.json').read_text());changes=[]
for r in a['files']:
 p=V/'source'/r['path'];assert sha(p)==r['sha256'];b=P/'build'/r['path']
 if p.read_bytes()!=b.read_bytes():changes.append(r['path'])
assert set(changes)=={'src/snlc_fit.car','src/sntools_output.h'},changes
assert sha(P/'build/src/snlc_fit.car')==sha(P/'instrumented-snlc_fit.car')
assert (P/'build/src/sntools_output.h').read_bytes()==(V/'build/src/sntools_output.h').read_bytes()
assert (V/'build/bin/snlc_fit.exe').is_file() and sha(V/'build/bin/snlc_fit.exe')=='d42c96c18984434e4e9ddc8bfd80f7cc3a18c4b64ede0467f98116719dcc5b51'
b=P/'build/bin/snlc_fit.exe';env=os.environ.copy();env['LD_LIBRARY_PATH']=str(R/'phase2/official/build/sysroot/usr/lib');ldd=subprocess.check_output(['ldd',str(b)],text=True,env=env);assert 'not found' not in ldd
result={'pass':True,'source_files_verified':len(a['files']),'modified_source_files':changes,'source_patch_sha256':sha(P/'ported-source.patch'),'binary_sha256':sha(b),'binary_bytes':b.stat().st_size,'ldd':ldd,'source_preserved':True,'original_binary_preserved':True,'native_equivalence':'Not yet executed; mandatory next gate.','build_seconds':json.loads((P/'build-result.json').read_text())['elapsed_seconds']}
(P/'build-verification.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='ldd'},indent=2))
