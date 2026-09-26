from pathlib import Path
import hashlib,json
R=Path.cwd();P=Path(__file__).parent
report=R/'docs/research-2026-09-26/raisin-simulation-replay-coherence.md'
inputs=[
'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres/nir.FITRES.gz',
'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres/lcplot-feasibility/FITOPT000.LCPLOT.gz',
'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres/acquisition.json',
'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres/lcplot-feasibility/acquisition.json',
'runs/research_2026_09_26/raisin_sign_source/sim/simlibs/DES_RAISIN.simlib',
'runs/research_2026_09_26/raisin_simulation_execution_review/author/output/fit_nir/DES_RAISIN_NIR_SIM/SUBMIT.INFO',
'runs/research_2026_09_26/raisin_simulation_execution_review/author/SIMLOGS_DES_SPEC/SIMGEN_MASTER_DESSPEC.INPUT',
'runs/research_2026_09_26/raisin_simulation_assets/timing/REFAC_DES_RAISIN_NIR.nml',
'runs/research_2026_09_26/raisin_simulation_timing_pilot/baseline-gate-authorgrid.json',
'runs/research_2026_09_26/raisin_simulation_timing_pilot/baseline-freeze-authorgrid.json',
'runs/research_2026_09_26/raisin_simulation_timing_pilot/fits_authorgrid/baseline/fit.log',
'runs/research_2026_09_26/raisin_simulation_timing_pilot/fits_authorgrid/baseline/fit.FITRES.TEXT',
'runs/research_2026_09_26/raisin_simulation_timing_pilot/fits_authorgrid/baseline/fit.LCPLOT.TEXT',
'sources/repos/RickKessler__SNANA@v11_04k/src/snlc_fit.car',
'sources/repos/RickKessler__SNANA@v11_04k/src/sntools_output_text.c']
files=[R/x for x in inputs]+[report]+[p for p in P.iterdir() if p.is_file() and p.name!='manifest.json']
d={str(p.relative_to(R)):{'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size} for p in sorted(set(files))}
result={'status':'complete independently verified read-only archive integrity audit; event identity not recoverable from reviewed products','files':d}
(P/'manifest.json').write_text(json.dumps(result,indent=2)+'\n')
for name,v in d.items():assert hashlib.sha256((R/name).read_bytes()).hexdigest()==v['sha256']
print(json.dumps({'files':len(d),'manifest_sha256':hashlib.sha256((P/'manifest.json').read_bytes()).hexdigest(),'report_sha256':hashlib.sha256(report.read_bytes()).hexdigest()},indent=2))
