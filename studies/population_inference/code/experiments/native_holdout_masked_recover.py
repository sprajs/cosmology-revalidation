"""Preserve failed file-resolution logs; complete fixed-mask derivative-only check."""
from pathlib import Path
import os,sys,re,json,hashlib,subprocess
from concurrent.futures import ThreadPoolExecutor
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
ROOT=Path(__file__).resolve().parents[3];P=Path(__file__).resolve().parent;O=P/'hold_jm';E=ROOT/'runs/research_2026_09_26/simulation_holdout_control';EXE=ROOT/'phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe';PY=ROOT/'phase2/env-official/bin/python';EXPORT=ROOT/'scripts/research_2026_09_26/export_flux_objectives.py';sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();selection={'P21':[19675,14533],'G10':[980,6801]};steps={1:.001,2:.0001,3:.01}
for arm in selection:
 d=O/'inputs'/('HE_'+arm);src=d/f'HE_{arm}_PHOT.FITS';dst=d/f'PH2_pilot02_{arm}_PHOT.FITS';os.link(src,dst);assert sha(src)==sha(dst)
(O/'file-resolution-recovery.json').write_text(json.dumps({'reason':'Copied HEAD primary header retains original PHOTFILE basename; add byte-identical hardlink in isolated derivative-input directory. No original data or scientific input changes.','source_sha256':sha(__file__)},indent=2)+'\n')
jobs=sorted(p.parent for p in O.glob('*/*/p*/fit.nml'))+sorted(p.parent for p in O.glob('*/baseline/fit.nml'))
env=os.environ.copy();env.update(SNANA_DIR=str(ROOT/'phase2/official/build/SNANA-audit-v3'),SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
def run(d):
 if (d/'fit.log').exists():(d/'fit.log').rename(d/'fit-missingPHOT.log')
 with (d/'fit.log').open('x') as f:r=subprocess.run([str(EXE),str(d/'fit.nml')],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 assert r.returncode==0,d
 r=subprocess.run([str(PY),str(EXPORT),'--log',str(d/'fit.log'),'--output',str(d/'objectives'),'--expected-cids',str(d/'cids.txt')],cwd=ROOT,env=env,capture_output=True,text=True);(d/'export.log').write_text(r.stdout+r.stderr);assert r.returncode==0,(d,r.stderr)
 return str(d.relative_to(O))
with ThreadPoolExecutor(max_workers=2) as pool:
 for d in pool.map(run,jobs):print('finished',d,flush=True)
s=(P/'native_holdout_masked_gate.py').read_text();exec(compile(s[s.index('# Comparison is separate'):],str(P/'native_holdout_masked_gate.py'),'exec'),globals())
