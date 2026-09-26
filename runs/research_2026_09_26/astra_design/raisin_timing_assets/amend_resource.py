from pathlib import Path
import json,hashlib
R=Path.cwd();P=Path(__file__).parent;O=P/'baseline_2021';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
oldroot=R/'runs/research_2026_09_26/raisin_simulation_timing_pilot';prior=[]
for p in sorted(oldroot.rglob('execution.json')):
 d=json.loads(p.read_text());prior.append({'path':str(p.relative_to(R)),'sha256':sha(p),'wall_seconds':float(d['wall_seconds'])})
a={'scope':'Resource-only amendment before any2021baseline fit; v1inputs,protocol,runner,freeze preserved','per_native_process_seconds':60,'combined_2021_baseline_pair_seconds':120,'original_total_pilot_seconds':600,'earlier_native_ledger':prior,'earlier_native_wall_seconds':sum(r['wall_seconds'] for r in prior),'carry_rule':'Before each branch compute spent2021pair time from immutable completed branch execution.json; perprocess timeout=min(60,120-pairspent,600-earlierspent-pairspent). Failbeforelaunch whennonpositive. No resetbetweenbaselineandcopy.','science_changes':None,'SNANA_DIR_length':len(str((P/'snana_v11_04d/build').resolve()))};assert a['SNANA_DIR_length']<160
ap=O/'resource-amendment.json';assert not ap.exists();ap.write_text(json.dumps(a,indent=2)+'\n');ah=sha(ap)
s=(O/'baseline_runner.py').read_text().replace('baseline-freeze.json','baseline-freeze-v2.json').replace('baseline_runner.py','baseline_runner_v2.py').replace('baseline_gate.py','baseline_gate_v2.py')
s=s.replace("    assert sha(DESIGN) == PROTOCOL_SHA",f"    assert sha(DESIGN) == PROTOCOL_SHA\n    assert sha(OUT / 'resource-amendment.json') == '{ah}'")
s=s.replace("    assert release['protocol_sha256'] == PROTOCOL_SHA", "    assert release['protocol_sha256'] == PROTOCOL_SHA\n    assert release['freeze_sha256'] == sha(OUT / 'baseline-freeze-v2.json')")
s=s.replace("    begin = time.monotonic()", "    resource = json.loads((OUT / 'resource-amendment.json').read_text())\n    for entry in resource['earlier_native_ledger']:\n        assert sha(ROOT / entry['path']) == entry['sha256']\n    pair_spent = sum(json.loads(p.read_text())['wall_seconds'] for p in (OUT / 'fits').glob('*/execution.json'))\n    prior_spent = resource['earlier_native_wall_seconds']\n    allowed = min(resource['per_native_process_seconds'], resource['combined_2021_baseline_pair_seconds'] - pair_spent, resource['original_total_pilot_seconds'] - prior_spent - pair_spent)\n    assert allowed > 0, ('native resource cap exhausted', prior_spent, pair_spent)\n    begin = time.monotonic()")
s=s.replace('timeout=600','timeout=allowed')
s=s.replace("'wall_seconds': time.monotonic() - begin,", "'wall_seconds': time.monotonic() - begin, 'timeout_seconds': allowed, 'earlier_native_seconds': prior_spent, 'pair_seconds_before_branch': pair_spent, 'resource_amendment_sha256': sha(OUT / 'resource-amendment.json'),")
(O/'baseline_runner_v2.py').write_text(s)
g=(O/'baseline_gate.py').read_text().replace('baseline-freeze.json','baseline-freeze-v2.json').replace('baseline_runner.py','baseline_runner_v2.py');(O/'baseline_gate_v2.py').write_text(g)
f=json.loads((O/'baseline-freeze.json').read_text());f['status']='v2resource-only freeze before native fits';f['preserved_v1_freeze_sha256']=sha(O/'baseline-freeze.json');f['resource_amendment_sha256']=ah
for n in ['baseline_runner.py','baseline_gate.py']:f['files'].pop(str((O/n).relative_to(R)))
for n in ['baseline_runner_v2.py','baseline_gate_v2.py','resource-amendment.json']:f['files'][str((O/n).relative_to(R))]=sha(O/n)
fp=O/'baseline-freeze-v2.json';fp.write_text(json.dumps(f,indent=2)+'\n')
for n,h in f['files'].items():assert sha(R/n)==h
print(json.dumps({'resource_amendment_sha256':ah,'freeze_v2_sha256':sha(fp),'prior_native_ledger_count':len(prior),'prior_native_seconds':a['earlier_native_wall_seconds'],'SNANA_DIR_length':a['SNANA_DIR_length']},indent=2))
