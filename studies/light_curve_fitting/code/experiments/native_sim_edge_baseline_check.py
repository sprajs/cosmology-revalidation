"""Check that fixed-coordinate native mean equals free-fit exported baseline at mean edges."""
from pathlib import Path
import os,re,json,subprocess,hashlib
import numpy as np
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parents[3];P=Path(__file__).resolve().parent;O=P/'sim_edge_J';E=ROOT/'runs/research_2026_09_26/simulation_residual_control/full256';EXE=ROOT/'phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe';PY=ROOT/'phase2/env-official/bin/python';EXP=ROOT/'scripts/research_2026_09_26/export_flux_objectives.py'
sel={'P21':[19612],'G10':[4058,7349,24318]};jobs=[]
for arm,ids in sel.items():
 for law in ['approx_minus99','exact_99']:
  d=O/arm/law/'base';d.mkdir();lines=['VARNAMES: CID PKMJD x0 x1 c']
  for cid in ids:
   x0,x1,c,t0=np.load(E/arm/law/f'objectives/objective_{cid}.npz')['parameters_x0_x1_c_t0'];lines.append(f'SN: {cid} {t0:.17g} {x0:.17g} {x1:.17g} {c:.17g}')
  seed=d/'seed.FITRES';seed.write_text('\n'.join(lines)+'\n');s=(O/arm/law/'p1_plus/fit.nml').read_text();s=re.sub(r"(?m)^\s*SNCID_LIST_FILE\s*=.*$",f" SNCID_LIST_FILE = '{seed}'",s);s=re.sub(r"(?m)^\s*TEXTFILE_PREFIX\s*=.*$",f" TEXTFILE_PREFIX = '{d/'fit'}'",s);(d/'fit.nml').write_text(s);(d/'cids.txt').write_text('\n'.join(map(str,ids))+'\n');jobs.append((d,arm,law,ids))
env=os.environ.copy();env.update(SNANA_DIR=str(ROOT/'phase2/official/build/SNANA-audit-v3'),SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
def run(job):
 d,arm,law,ids=job
 with (d/'fit.log').open('x') as f:r=subprocess.run([str(EXE),str(d/'fit.nml')],cwd=d,env=env,stdout=f,stderr=subprocess.STDOUT)
 assert r.returncode==0
 r=subprocess.run([str(PY),str(EXP),'--log',str(d/'fit.log'),'--output',str(d/'objectives'),'--expected-cids',str(d/'cids.txt')],cwd=ROOT,env=env,capture_output=True,text=True);(d/'export.log').write_text(r.stdout+r.stderr);assert r.returncode==0
 rows=[]
 for cid in ids:
  a=np.load(E/arm/law/f'objectives/objective_{cid}.npz');b=np.load(d/f'objectives/objective_{cid}.npz')
  for k in ['MJD','band','data_flux','data_fluxerr','parameters_x0_x1_c_t0','model_flux']:assert np.array_equal(a[k],b[k]),(arm,law,cid,k)
  rows.append(dict(arm=arm,law=law,CID=cid,all_identity_exact=True))
 return rows
with ThreadPoolExecutor(max_workers=2) as pool:rows=sum(pool.map(run,jobs),[])
r={'status':'PASS: exact native fixed-coordinate/free-fit baseline identity','checks':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(O/'baseline-identity.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
