from pathlib import Path
import os,subprocess,time,json
ROOT=Path(__file__).resolve().parents[3];BASE=Path(__file__).resolve().parent/'shared43'
env=os.environ.copy();env.update(SNANA_DIR=str(ROOT/'phase2/official/build/SNANA-audit-v3'),SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
for k in range(1,10):
 p=BASE/f'model{k:03d}';t=time.monotonic()
 assert not (p/'fit.log').exists()
 with (p/'fit.log').open('w') as log:subprocess.run([str(ROOT/'phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe'),str(p/'fit.nml')],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'scripts/research_2026_09_26/export_flux_objectives.py'),'--log',str(p/'fit.log'),'--output',str(p/'objectives'),'--expected-cids',str(BASE/'cids.txt')],cwd=ROOT,env=env,check=True)
 print('model',k,'wall_seconds',time.monotonic()-t,flush=True)
