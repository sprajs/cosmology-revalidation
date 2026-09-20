#!/usr/bin/env python3
"""Archive code/plan snapshots by digest and inventory final local evidence."""
from pathlib import Path
import json,hashlib,datetime
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'runs/age_signal';DEST=OUT/'provenance';DEST.mkdir(exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
    manifests=[json.loads(p.read_text()) for p in OUT.glob('*-manifest.json')]
    needed_plan={j['plan_sha256'] for j in manifests};needed_code={j['code_sha256'] for j in manifests}
    p=(ROOT/'docs/experiments/age_signal-plan.md').read_bytes();found_plan=set()
    for end in range(len(p)):
        if p[end:end+1]!=b'\n':continue
        candidate=p[:end+1];h=sha(candidate)
        if h in needed_plan: (DEST/f'{h}.plan.md').write_bytes(candidate);found_plan.add(h)
    code=(ROOT/'scripts/age_signal/analyse.py').read_text(); candidates=[code]
    pilot=code.replace('miniter=20000,maxiter=100000','miniter=5000,maxiter=30000').replace("'miniter':20000,'maxiter':100000","'miniter':5000,'maxiter':30000").replace("for pol in ['G11_first','R19_first']:","for pol in ['G11_first','R19_first','all_rows']:").replace("if pol!='all_rows': variants +=","if pol=='R19_first': variants +=").replace("v=np.asarray(values).reshape(4,-1); n2=v.shape[1]//2; v=np.concatenate([v[:,:n2],v[:,-n2:]],axis=0); n=v.shape[1]","v=np.asarray(values).reshape(4,-1); n=v.shape[1]")
    candidates.append(pilot)
    old=pilot.replace('nchains=4,parallelize=True,seed=20260920+i','nchains=4,parallelize=False,seed=20260920+i').replace('lm.run_mcmc(miniter=5000,maxiter=30000,silent=False)','lm.run_mcmc(miniter=5000,maxiter=30000,silent=True)').replace("'rhat_retained':{k:rhat_four(chain[k]) for k in ['alpha','beta','sigsqr']}","'rhat':{str(k):float(v) for k,v in zip(['alpha','beta','log_sigsqr','ximean','xivar','atanh_corr'],lm._get_Rhat())}")
    start=old.index('def rhat_four(values):');end=old.index('def linmix_run():');old=old[:start]+old[end:];candidates.append(old)
    found_code=set()
    for c in candidates:
        h=sha(c.encode())
        if h in needed_code:(DEST/f'{h}.analyse.py').write_text(c);found_code.add(h)
    for file in (ROOT/'scripts/age_signal').glob('*.py'):(DEST/f'{sha(file.read_bytes())}.{file.name}').write_bytes(file.read_bytes())
    (DEST/f'{sha(p)}.plan.md').write_bytes(p)
    outputs=[q for q in OUT.rglob('*') if q.is_file() and q.suffix!='.npz' and 'provenance' not in q.parts and q.name!='final-inventory.json']
    result={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'plan_snapshots_resolved':sorted(found_plan),'plan_snapshots_unresolved':sorted(needed_plan-found_plan),'analyse_snapshots_resolved':sorted(found_code),'analyse_snapshots_unresolved':sorted(needed_code-found_code),'files':{str(f.relative_to(ROOT)):sha(f.read_bytes()) for f in outputs},'note':'Large numerical chains are preserved locally. Historical script candidates were archived only when their bytes match the recorded SHA-256 exactly.'}
    (OUT/'final-inventory.json').write_text(json.dumps(result,indent=2));print({k:v for k,v in result.items() if k!='files'})
if __name__=='__main__':main()
