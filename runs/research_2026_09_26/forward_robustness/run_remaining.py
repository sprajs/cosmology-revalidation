import json,os,subprocess,time
from pathlib import Path
R=Path(__file__).resolve().parents[3];out=R/'runs/research_2026_09_26/forward_robustness'
variants=['recovered-mask','strict-cuts','bandwidth075','bandwidth150','reference-qplus05','reference-qminus1']
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
for name in variants:
 t=time.time()
 with (out/f'{name}.stdout.txt').open('w') as stdout,(out/f'{name}.stderr.txt').open('w') as stderr:
  proc=subprocess.run([str(R/'.venv/bin/python'),str(R/'scripts/research_2026_09_26/forward_robustness.py'),'--variant',name],cwd=R,env=env,stdout=stdout,stderr=stderr)
 summary=out/name/'summary.json'
 cohort=json.loads(summary.read_text())['cohort'] if summary.exists() else None
 print(json.dumps({'variant':name,'exit_code':proc.returncode,'elapsed_seconds':round(time.time()-t,1),'cohort':cohort}),flush=True)
