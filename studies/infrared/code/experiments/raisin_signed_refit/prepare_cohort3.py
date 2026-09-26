from pathlib import Path
import sys,json,subprocess
O=Path('runs/research_2026_09_26/astra_design/raisin_signed_refit').resolve();sys.path.insert(0,str(O));import paired_refit as p
C=O/'fixed-c-profile/cohort10'
s=(O/'cohort_profile_engine_v2.py').read_text().replace('from check_results import table','from check_results import table\nfrom oracle_prefix_gate import check as prefix_gate')
s=s.replace("oracle=Oracle(G/'stream',pars=pars);assert oracle.n==len(kb)","oracle=Oracle(G/'stream',pars=pars);assert oracle.n==len(kb)\n    try:\n      gate=prefix_gate(G/'stream/fit.log',nom,CID,pars);p.save(G/'stream-prefix-gate.json',gate)\n    except Exception:\n      oracle.close();raise")
(O/'cohort_profile_engine_v3.py').write_text(s)
# Gate completed benchmarks retrospectively without modifying their saved manifests.
from oracle_prefix_gate import check
import numpy as np
for cid,rel in [('DES16C1cim','DES16C1cim'),('DES16C3cmy','DES16C3cmy/v2')]:
 d=C/rel;nom=np.load(d/'fixed_reference/native.npz');r=check(d/'stream/fit.log',nom,cid,nom['parameters']);r['before_first_query']=False;r['timing']='Retrospective prefix validation of completed benchmark; outputs unchanged';p.save(C/(cid+'-retrospective-prefix-gate.json'),r)
proto=json.loads((C/'protocol.json').read_text());am=dict(status='Frozen before remaining8 profiles; benchmarks already inspected',change='Exact object-specific native reference parameters, epoch order, flux/errors, mean and precision checked from initial native export before first stream query. No scientific setting changed.',engine_sha256=p.sha(O/'cohort_profile_engine_v3.py'),prefix_gate_sha256=p.sha(O/'oracle_prefix_gate.py'),controller_sha256=p.sha(O/'cohort_profile_remaining_v3.py'),prior_amendment_sha256=p.sha(C/'oracle-coordinate-amendment.json'),resource='Same900 total active seconds including failed attempt, one worker',benchmark='Both completed benchmark initial exports retrospectively match exactly; no rerun')
p.save(C/'prefix-gate-amendment.json',am)
idx=json.loads((C/'active-attempts.json').read_text());p.save(C/'active-attempts-v2-preserved.json',idx)
for cid in proto['cohort'][2:]:
 d=C/cid/'v3';d.mkdir();pr=json.loads((C/cid/'v2/protocol.json').read_text());pr['engine_sha256']=p.sha(O/'cohort_profile_engine_v3.py');pr['hashes'].pop(str((O/'cohort_profile_engine_v2.py').relative_to(p.ROOT)))
 for f in [O/'cohort_profile_engine_v3.py',O/'oracle_prefix_gate.py',C/'prefix-gate-amendment.json']:pr['hashes'][str(f.relative_to(p.ROOT))]=p.sha(f)
 p.save(d/'protocol.json',pr);idx['attempts'][cid]=str(Path(cid)/'v3')
idx['prefix_gate_amendment_sha256']=p.sha(C/'prefix-gate-amendment.json');p.save(C/'active-attempts.json',idx)
(C/'prefix-gate.patch').write_text(subprocess.run(['diff','-u',str(O/'cohort_profile_engine_v2.py'),str(O/'cohort_profile_engine_v3.py')],capture_output=True,text=True).stdout)
print(json.dumps(am))
