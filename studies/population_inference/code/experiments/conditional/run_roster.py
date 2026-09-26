"""Serial registered conditional primary roster; print diagnostics, never scores."""
import json,os,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'runs/research_2026_09_26/conditional'
RUNNER=ROOT/'scripts/phase2/hierarchy/conditional.py'
PYTHON=ROOT/'phase2/hierarchy/.venv/bin/python'
models=['none','stretch','colour','tripp','host','host_colour','broken_colour','evolution','flexible']
env=os.environ.copy();env.update(OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',XLA_FLAGS='--xla_cpu_multi_thread_eigen=false intra_op_parallelism_threads=1')
for noise in ('gaussian','student4'):
    for model in models:
        name=f'main_{noise}_{model}'
        argv=[str(PYTHON),str(RUNNER),'--data','phase2/hierarchy/data/conditioned-multistart-best','--name',name,'--model',model,'--noise',noise,'--warmup','600','--draws','1000','--chains','4','--seed','2026092160','--output-root','runs/research_2026_09_26/conditional']
        t=time.time()
        with (OUT/f'{name}.stdout.txt').open('w') as stdout,(OUT/f'{name}.stderr.txt').open('w') as stderr:
            proc=subprocess.run(argv,cwd=ROOT,env=env,stdout=stdout,stderr=stderr)
        folder=OUT/name;gate=folder/'convergence.json'
        diag=json.loads(gate.read_text()) if gate.exists() else {}
        state={'name':name,'exit_code':proc.returncode,'elapsed_seconds':round(time.time()-t,1),'gate':diag.get('valid_for_scoring'),'max_r_hat':diag.get('max_r_hat'),'min_n_eff':diag.get('min_n_eff'),'divergences':diag.get('divergences'),'max_steps':diag.get('max_steps'),'mean_acceptance':diag.get('mean_acceptance'),'reasons':diag.get('reasons')}
        print(json.dumps(state),flush=True)
