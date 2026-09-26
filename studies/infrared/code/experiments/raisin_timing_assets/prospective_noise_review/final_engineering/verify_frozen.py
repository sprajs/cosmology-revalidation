"""Independent metadata/hash gate only. Never imports or invokes the native executor."""
from pathlib import Path
import re,json,hashlib,ast,csv
from collections import Counter
R=Path.cwd();P=R/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering';O=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
protocol=json.loads((P/'protocol.json').read_text());freeze=json.loads((P/'freeze.json').read_text())
checked=[]
for name,h in freeze['files'].items():
 p=R/name;assert sha(p)==h,(name,'hash mismatch');checked.append(name)
runner=P/'run_engineering.py';ast.parse(runner.read_text())
src=R/'runs/research_2026_09_26/raisin_sign_source/sim/simlibs/DES_RAISIN.simlib';s=src.read_text();block=re.search(r'^LIBID:\s+11\s.*?^END_LIBID:\s*11.*?$',s,re.M|re.S).group()
old=[l for l in block.splitlines() if l.startswith('S:')];cad=P/'inputs-v2/cadence.simlib';new=[l for l in cad.read_text().splitlines() if l.startswith('S:')];hdr=cad.read_text()
z=float(re.search(r'REDSHIFT:\s+(\S+)',hdr).group(1));t0=float(re.search(r'PEAKMJD:\s+(\S+)',hdr).group(1))
expected=[l for l in old if -15<=(float(l.split()[1])-t0)/(1+z)<=45]
assert new==expected and len(new)==117;assert Counter(l.split()[3] for l in new)==Counter(g=27,r=28,i=28,z=28,J=3,H=3)
assert re.search(r'FIELD:\s+(\S+)',hdr).group(1)=='DES16E1dcx';assert int(re.search(r'NOBS:\s+(\S+)',hdr).group(1))==117
# Source-selected signed exposure metadata retained byte for byte; no photometry values exist yet.
r={'status':'PASS_frozen_hashes_and_independent_cadence_metadata','native_invocations':0,'executor_imports':0,'manifest_entries':len(checked),'source_block_rows':len(old),'selected_rows':len(new),'selected_band_counts':dict(Counter(l.split()[3] for l in new)),'source_exposure_lines_order_and_bytes_exact':True,'protocol_sha256':sha(P/'protocol.json'),'freeze_sha256':sha(P/'freeze.json'),'runner_sha256':sha(runner),'reviewer_sha256':sha(Path(__file__)),'checked_files':checked}
(O/'frozen-check.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k!='checked_files'},indent=2))
