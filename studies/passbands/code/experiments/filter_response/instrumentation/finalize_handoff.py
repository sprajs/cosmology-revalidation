from pathlib import Path
import json,hashlib
import numpy as np
from parse_native import parse,ROW
O=Path(__file__).resolve().parent;ROOT=Path.cwd();sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for run in ['probes-v3','changed-pilot-probes-v3']:
 P=O/run;arms=json.loads((P/('derived-fixed-probes.json' if run.startswith('changed') else 'protocol.json')).read_text())['arms'];arrays={};index=[]
 for a in arms:
  for i,b in enumerate(parse(P/'fits'/a['name']/'fit.log')):
   key=a['name']+f'__b{i:03d}';arrays[key+'_rows']=b['array'];arrays[key+'_W']=b['W'];arrays[key+'_C']=b['C'];arrays[key+'_bands']=np.array([r['band'] for r in b['rows']]);arrays[key+'_source_epoch']=np.array([r['source_epoch'] for r in b['rows']]);index.append({'key':key,'CID':b['CID'],'ITER':b['ITER'],'objective':b['objective'],'entry':b['entry']})
 np.savez_compressed(P/'native-arrays.npz',**arrays);(P/'array-index.json').write_text(json.dumps({'rows_schema':ROW,'blocks':index},indent=2)+'\n')
B=O/'SNANA-v11_04k-output';old=ROOT/'phase2/official/build/SNANA-v11_04k';differences=[];checked=0
for p in (old/'src').iterdir():
 if not p.is_file() or p.suffix not in ['.car','.c','.h','.F']:continue
 checked+=1;q=B/'src'/p.name
 if sha(p)!=sha(q):differences.append(p.name)
assert differences==['snlc_fit.car'],differences
x={'binary_sha256':sha(B/'bin/snlc_fit.exe'),'original_binary_sha256':sha(old/'bin/snlc_fit.exe'),'source_sha256':sha(B/'src/snlc_fit.car'),'source_files_checked':checked,'source_files_different':differences,'v1_failed_compile_log_duration_seconds_from_birth_to_last_write':4.493056182,'v2_compile_log_duration_seconds_from_birth_to_last_write':15.191194598,'v1_v2_timing_caveat':'File timestamps bound log-writing intervals, not instrumented subprocess timers; neither hit its timeout. Originallogs preserved.','v3_monotonic_compile_seconds':json.loads((O/'build-v3-result.json').read_text())['seconds'],'build_timeout_limits_seconds':[180,120,120],'short_alias':'/tmp/snana-csp-audit-8c5f0d','alias_resolves_to':str(Path('/tmp/snana-csp-audit-8c5f0d').resolve()),'native_seconds':json.loads((O/'resource-ledger.json').read_text())['all_active_native_seconds'],'default_environment':'Explicitly remove CSP_POSTINIT_DSHIFT. Only deliberately requested starts set ±0.2.','scope':'Two pilots only; native mean/matrix and starts verified. No full42 changed response here.'}
(O/'build-and-handoff.json').write_text(json.dumps(x,indent=2)+'\n')
files=[p for p in O.rglob('*') if p.is_file() and 'SNANA-v11_04k-output' not in p.parts and p.name!='handoff-manifest.json']+[B/'src/snlc_fit.car',B/'bin/snlc_fit.exe',ROOT/'docs/research-2026-09-26/csp-filter-interpretation-review.md']
(O/'handoff-manifest.json').write_text(json.dumps({str(p.relative_to(ROOT)):sha(p) for p in files},indent=2)+'\n');print(json.dumps(x,indent=2));print('manifest',sha(O/'handoff-manifest.json'))
