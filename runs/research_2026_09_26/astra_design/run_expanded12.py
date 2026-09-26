from pathlib import Path
import os,subprocess,time
ROOT=Path(__file__).resolve().parents[3];P=Path(__file__).resolve().parent/'expanded12';env=os.environ.copy();env.update(SNANA_DIR=str(ROOT/'phase2/official/build/SNANA-audit-v3'),SNDATA_ROOT=str(ROOT/'phase2/official/inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(ROOT/'phase2/official/build/sysroot/usr/lib'),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
for cohort in ['discovery','validation']:
 for name in ['CALSPEC','MWEBV','COLORLAW']:
  p=P/cohort/name;t=time.monotonic();assert not (p/'fit.log').exists()
  with (p/'fit.log').open('w') as log:subprocess.run([str(ROOT/'phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe'),str(p/'fit.nml')],env=env,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
  subprocess.run([str(ROOT/'.venv/bin/python'),str(ROOT/'scripts/research_2026_09_26/export_flux_objectives.py'),'--log',str(p/'fit.log'),'--output',str(p/'objectives'),'--expected-cids',str(P/cohort/'cids.txt')],env=env,cwd=ROOT,check=True)
  print(cohort,name,'wall_seconds',time.monotonic()-t,flush=True)
