"""Reconstruct the bounded historical capacity-regression fixtures.

Requires the first twelve generated training shards and fixed evaluation pool.
Runs four baseline cases and the twenty nominal pilot cases with the ORIGINAL
BBC executable; the caller replays the six failures and four baselines with the
separate repaired executable. It does not generate new photons or change seeds.
"""
import json
from concurrent.futures import ThreadPoolExecutor
import subprocess
import numpy as np
from common import RESULTS, WORK, sha, native_env
import scaled_bbc as bbc
import scaled_bootstrap as bootstrap
from scaled_text_tables import write_clones


def restore():
    root = bbc.SCALE / 'capacity-fixtures'
    assert not root.exists(), 'Refuse to overwrite restored regression fixtures'
    root.mkdir(parents=True)
    old_results = bbc.RESULTS
    (root/'scaled-campaign.json').write_text((RESULTS/'scaled-campaign.json').read_text())
    try:
        bbc.RESULTS=root
        bbc.merge(12)
    finally:
        bbc.RESULTS=old_results
    rec=json.loads((root/'scaled-inputs-k012.json').read_text())
    tables=rec['tables']
    nominal=bbc.fitres(tables['eval_nominal']['path'])
    age=bbc.fitres(tables['eval_age']['path'])
    ids=nominal.index.intersection(age.index)
    ev={arm:bbc.write_table(f.loc[ids],root/f'eval-{arm}-common.FITRES')
        for arm,f in [('nominal',nominal),('age',age)]}
    specifications=[]
    for name,arm,target in [('nominal','nominal','nominal_literal'),
                            ('age_nominal','age','nominal_literal'),
                            ('age_literal','age','age_literal'),
                            ('age_retained','age','age_retained_age')]:
        specifications.append(('baseline-'+name,root/'original'/('baseline-'+name),
                               (ev[arm]['path'],tables[target]['path']),True))
    old_results=bootstrap.RESULTS
    try:
        bootstrap.RESULTS=root
        frames,attempts=bootstrap.prepare(12,True)
    finally:
        bootstrap.RESULTS=old_results
    for i in range(20):
        folder=root/'original'/f'pilot-r{i:03d}'
        rng=np.random.default_rng(972000+i)
        drawn={pool:rng.integers(n,size=n) for pool,n in attempts.items()}
        files={}
        for pool in ['train','eval']:
            f=bootstrap.draw_table(frames[pool+'_nominal'],drawn[pool],
                                   10000000 if pool=='train' else 20000000)
            key='nominal_literal' if pool=='train' else 'eval_nominal'
            files[pool]=write_clones(f,folder/(pool+'.FITRES'),frames['_caches'][key])
        specifications.append((f'pilot-r{i:03d}',folder,(files['eval']['path'],files['train']['path']),False))
    original_exe=WORK/'SNANA-legacy-eff/bin/SALT2mu.exe'
    def run(spec):
        label,folder,(data,train),baseline=spec
        folder.mkdir(parents=True,exist_ok=True)
        text,source,changes=bbc.source_config(data,train,True)
        inp=folder/'bbc.input'
        inp.write_text(text)
        with (folder/'bbc.log').open('w') as stream:
            p=subprocess.run([str(original_exe),str(inp)],cwd=folder,env=native_env(),
                             stdout=stream,stderr=subprocess.STDOUT)
        log=(folder/'bbc.log').read_text()
        good=p.returncode==0 and 'Done.' in log[-1000:] and 'FATAL ERROR ABORT' not in log
        r={'name':label,'graceful':good,'executable_sha256':sha(original_exe),
           'input_sha256':sha(inp),'log_sha256':sha(folder/'bbc.log'),
           'seed':None if baseline else 972000+int(label[-3:]),
           'scope':'Original unrepaired BBC; four baseline and twenty nominal pilot fixtures only.'}
        (folder/'original-run.json').write_text(json.dumps(r,indent=2)+'\n')
        if baseline: assert good,label
        return spec,good
    with ThreadPoolExecutor(max_workers=3) as pool: completed=list(pool.map(run,specifications))
    failures=[spec for spec,good in completed if not good]
    assert [x[0] for x in failures]==[f'pilot-r{i:03d}' for i in [5,7,8,16,17,19]]
    (root/'restoration.json').write_text(json.dumps({
        'code_sha256':sha(__file__),'original_executable':str(original_exe),
        'original_executable_sha256':sha(original_exe),
        'original_native_runs':24,'bootstrap_seed_base':972000,'pilot_replicates':20,
        'failures':[x[0] for x in failures],
        'scope':'Reconstructed numerical fixtures; timestamps and provenance hashes may differ from archived run. No photons or physical assumptions changed.'},indent=2)+'\n')
    return [spec for spec,good in completed if spec[3] or not good]
