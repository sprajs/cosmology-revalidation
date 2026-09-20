"""Replay saved fit configurations with preserved preliminary summaries."""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2]

def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=3);p.add_argument('--archive-first',action='store_true')
    args=p.parse_args();jobs=[]
    for f in sorted((ROOT/'runs/cosmology').glob('*/configuration.json')):
        c=json.loads(f.read_text())
        if c.get('action')!='fit':continue
        if args.archive_first:
            dest=ROOT/'runs/cosmology/preliminary'/c['name'];dest.mkdir(parents=True,exist_ok=True)
            for name in ['configuration.json','summary.json','manifest.json','q_history.csv']:
                source=f.parent/name
                if source.exists() and not (dest/name).exists():shutil.copy2(source,dest/name)
        cmd=[sys.executable,str(ROOT/'scripts/cosmology/run.py'),'fit']
        for key in ['name','model','data','amplitude','correction','column','zcolumn','scale','tau','steps','burn','walkers','seed','magnitude_table']:
            if c.get(key) is not None:cmd.extend(['--'+key.replace('_','-'),str(c[key])])
        jobs.append((c['name'],cmd))
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
    def run(job):
        name,cmd=job;log=ROOT/'runs/cosmology'/name/'replay.log'
        with log.open('w') as out:
            proc=subprocess.run(cmd,cwd=ROOT,env=env,stdout=out,stderr=subprocess.STDOUT)
        return name,proc.returncode
    failed=[]
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for f in as_completed([pool.submit(run,j) for j in jobs]):
            name,code=f.result();print(json.dumps({'experiment':name,'exit_code':code}),flush=True)
            if code:failed.append(name)
    if failed:raise SystemExit('Failed: '+', '.join(failed))

if __name__=='__main__':main()
