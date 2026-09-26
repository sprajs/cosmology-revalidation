#!/usr/bin/env python3
"""Prespecified numerical basin diagnostic on conditioned-mask anomalies.

Three shape starts x three time offsets, plus reference initialization. Results
are retained individually; selection never uses closeness to published values.
"""
from pathlib import Path
import concurrent.futures,hashlib,json,os,re,subprocess
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'phase2/official'
q=pd.read_csv(OUT/'results/conditioned_comparison.csv',dtype={'CID':str})
q=q[(q.IDSURVEY==10)&((q.delta_mB.abs()>.01)|(q.delta_Tripp_fixed_coefficients.abs()>.01))].copy()
base=(OUT/'inputs/snana_des_conditioned0.nml').read_text()
base=re.sub(r'(?m)^\s*MXLC_PLOT.*$','',base).replace("FITRES(text:host) LCPLOT(text:col)","FITRES(text:host)")
tasks=[]
for x in [-2,0,2]:
    for dt in [-4,0,4]:tasks.append((f'x{x:+d}_t{dt:+d}'.replace('+','p').replace('-','m'),x,dt))
tasks.append(('reference',None,None))
manifest={'selection':'DES with conditioned |delta_mB|>.01 or |delta_Tripp|>.01; all anomalies retained, no exclusion.','CID':q.CID.tolist(),'seed_grid':'x1=-2,0,2; t0=defaultconditionedfit_t0 +(-4,0,4)days; x0,c from defaultconditionedfit; plus one publishedreference seed.','same_mask':True,'selection_criterion':'Report all objective values; candidate min reported SNANA FITCHI2, never closest reference. Different iteratively frozen covariance matrices mean min-chi2 across branches is not a fully normalized common likelihood. Quantify branch ambiguity before selecting final inference arm.','tasks':[]}
for label,x,dt in tasks:
    vals=[]
    for row in q.itertuples():
        if x is None:vals.append((row.CID,row.PKMJD_published,row.x0_published,row.x1_published,row.c_published))
        else:vals.append((row.CID,row.t0_double+dt,row.x0_double,x,row.c_double))
    seed=OUT/f'inputs/conditioned_seed_{label}.FITRES'
    text='VARNAMES: CID PKMJD x0 x1 c\n'+''.join('SN: '+cid+' '+' '.join(f'{v:.17g}' for v in par)+'\n' for cid,*par in vals)
    seed.write_text(text)
    nml=OUT/f'inputs/snana_conditioned_seed_{label}.nml'
    s=base.replace('OPT_SNCID_LIST = 1','OPT_SNCID_LIST = 3')
    s=re.sub(r"SNCID_LIST_FILE = '[^']+'",f"SNCID_LIST_FILE = '{seed}'",s)
    s=s.replace('snana_des_conditioned0',f'snana_conditioned_seed_{label}')
    nml.write_text(s)
    manifest['tasks'].append({'label':label,'seed':str(seed.relative_to(ROOT)),'seed_sha256':hashlib.sha256(seed.read_bytes()).hexdigest(),'config':str(nml.relative_to(ROOT)),'config_sha256':hashlib.sha256(nml.read_bytes()).hexdigest()})
p=OUT/'inputs/conditioned-multistart-contract.json';p.write_text(json.dumps(manifest,indent=2)+'\n')
env=os.environ.copy();env.update(SNANA_DIR=str(OUT/'build/SNANA-current'),SNDATA_ROOT=str(OUT/'inputs/SNDATA_ROOT'),LD_LIBRARY_PATH=str(OUT/'build/sysroot/usr/lib'))
def run(t):
    label=t['label'];log=OUT/f'diagnostics/snana-conditioned-seed-{label}.log'
    assert not log.exists(),f'Preserve existing log: {log}'
    with log.open('w') as f:r=subprocess.run([str(OUT/'build/SNANA-current/bin/snlc_fit.exe'),str(ROOT/t['config'])],stdout=f,stderr=subprocess.STDOUT,env=env)
    return {'label':label,'returncode':r.returncode,'graceful':'ENDING PROGRAM GRACEFULLY.' in log.read_text()}
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:results=list(pool.map(run,manifest['tasks']))
(OUT/'diagnostics/conditioned-multistart-status.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))
