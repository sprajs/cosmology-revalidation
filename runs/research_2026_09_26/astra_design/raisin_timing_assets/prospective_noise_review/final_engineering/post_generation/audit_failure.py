from pathlib import Path
import hashlib,json,re
R=Path.cwd();O=Path(__file__).resolve().parent;P=R/'runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering';S=R/'runs/research_2026_09_26/astra_design/raisin_timing_assets/snana_v11_04d/source/src/snlc_sim.c'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
fail=json.loads((P/'generation-failure.json').read_text());activity=json.loads((P/'native-activity.json').read_text());log=P/'generation-v3/original/native.log';inp=P/'generation-v3/original/sim.input';cad=P/'inputs-v2/cadence.simlib';output=P/'generation-v3/original/output/PTE'
assert len(activity)==1 and activity[0]['returncode']==11
assert 'NTRY=100000 exceeds bound of MXREAD_SIMLIB=100000' in log.read_text()
assert re.search(r'^SIMLIB_IDLOCK:\s+1\s*$',inp.read_text(),re.M)
ids=[int(x) for x in re.findall(r'^LIBID:\s+(\d+)',cad.read_text(),re.M)];assert ids==[11]
paths=sorted(p for p in output.iterdir() if p.is_file());sizes={p.name:p.stat().st_size for p in paths};assert sizes['PTE_HEAD.FITS']==sizes['PTE_PHOT.FITS']==0
text=(output/'PTE.DUMP').read_text();event_rows=[l for l in text.splitlines() if l.startswith(('SN:','ROW:'))];assert not event_rows
assert not (P/'generation-gate.json').exists() and not (P/'generation-v3/ledger/native.log').exists() and not (P/'generation-v3/noiseless/native.log').exists()
lines=S.read_text().splitlines();evidence=[]
for a,b in [(355,369),(14997,15012),(15046,15070),(15858,15876)]:evidence.append({'first':a,'last':b,'lines':lines[a-1:b]})
inputs=[P/'protocol.json',P/'freeze.json',P/'run_engineering.py',P/'generation-failure.json',P/'native-activity.json',P/'generation-v3/original/execution.json',log,inp,cad,S]+paths
r={'status':'STOP_independent_generation_failure_confirmed','native_invocations_by_reviewer':0,'returncode':11,'native_wall_seconds_consumed':activity[0]['wall_seconds'],'native_budget_remaining_seconds':120-sum(x['wall_seconds'] for x in activity),'configured_IDLOCK':1,'available_LIBIDs':ids,'output_sizes_bytes':sizes,'written_DUMP_event_rows':0,'instrumentation_comparison_available':False,'noise_algebra_audit_available':False,'subsequent_stages_ran':False,'source_mechanism':'findStart assigns literal INPUTS.IDLOCK to GENLC at end; keep_SIMLIB_HEADER rejects LIBID11 against lock1. Mainloop first-accepted special case is reached only after successful generation and cannot repair this rejection.','interpretation':'Source-control preparation failure; no emitted epoch photometry, timing/noise result, or historical-simulation inference. Static preflight review missed this source interaction.','source_excerpts':evidence,'input_hashes':{str(p.relative_to(R)):sha(p) for p in inputs}}
(O/'failure-review.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:v for k,v in r.items() if k not in ['input_hashes','source_excerpts']},indent=2))
