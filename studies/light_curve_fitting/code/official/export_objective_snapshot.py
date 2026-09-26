#!/usr/bin/env python3
"""Export read-only final objective dumps and verify their self-consistency."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'phase2/official';DEST=OUT/'portable_pilot'
log=OUT/'diagnostics/snana-audit-pilot.log'
assert 'ENDING PROGRAM GRACEFULLY.' in log.read_text()
records={};cov={};contexts={}
for line in log.read_text().splitlines():
    q=line.split()
    if not q:continue
    if q[0]=='PHASE2_FLUX:':
        cid,i,band=q[1],int(q[2]),q[3]
        records.setdefault(cid,{})[i]=[band]+[float(x) for x in q[4:]]
    if q[0]=='PHASE2_OBJECTIVE:':contexts[q[1]]={'n':int(q[2]),'values':[float(x) for x in q[3:]]}
    if q[0]=='PHASE2_COVINV:':cov.setdefault(q[1],{})[int(q[2])]=[float(x) for x in q[3:]]
rows=[];checks={}
for cid,ctx in contexts.items():
    n=ctx['n'];ordered=[records[cid][i] for i in range(1,n+1)]
    arr=np.array([x[1:] for x in ordered]);ci=np.array([cov[cid][i] for i in range(1,n+1)])
    assert ci.shape==(n,n)
    residual=arr[:,4]-arr[:,2]
    direct=float(residual@ci@residual)
    v=ctx['values'];chi2,prior=v[:2]
    c=np.linalg.inv(ci)
    np.savez_compressed(DEST/f'objective_{cid}.npz',inverse_frozen_flux_covariance=ci,frozen_flux_covariance=c,parameters_x0_x1_c_t0=np.array(v[2:6]),chi2=np.array(chi2),prior_chi2=np.array(prior),initial_search_peakmjd=np.array(v[6]),MJD=arr[:,0],rest_phase=arr[:,1],model_flux=arr[:,2],model_magerr=arr[:,3],data_flux=arr[:,4],data_fluxerr=arr[:,5],zHEL=arr[:,6],MWEBV=arr[:,7],band=np.array([x[0] for x in ordered]))
    checks[cid]={'n':n,'chi2_dump':chi2,'prior_chi2':prior,'reconstructed_data_chi2':direct,'objective_residual':direct+prior-chi2,'min_cov_eigenvalue':float(np.linalg.eigvalsh(c)[0]),'inverse_asymmetry':float(np.max(np.abs(ci-ci.T)))}
    for i,(band,*values) in enumerate(ordered,1):rows.append([cid,i,band]+values)
cols=['CID','fit_row_one_based','BAND','MJD','rest_phase','model_flux','model_magerr','data_flux','data_fluxerr','zHEL','MWEBV']
pd.DataFrame(rows,columns=cols).to_csv(DEST/'exact_epoch_predictions.csv',index=False,float_format='%.17g')
from audit_fits import read_fit
a=read_fit(OUT/'results/snana_audit_pilot.FITRES.TEXT');b=read_fit(OUT/'results/snana_mask32_pilot.FITRES.TEXT')
cols=['CID','x0','x1','c','PKMJD','FITCHI2','NDOF'];a=a[cols].sort_values('CID').reset_index(drop=True);b=b[cols].sort_values('CID').reset_index(drop=True)
assert a.equals(b),'Output-only build changed fit result'
assert max(abs(x['objective_residual']) for x in checks.values())<1e-8
record={'description':'Exact per-accepted-observation REAL8 forward predictions and final frozen inverse flux covariance from last MINUIT function call in iteration3. Covariance was constructed between iterations, not recomputed at final parameters. Model_magerr evaluated at final flux parameters but x1_for_error comes from previous iteration.','unchanged_parameters_and_chi2_verified':True,'checks':checks,'hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [log,OUT/'build/SNANA-audit/bin/snlc_fit.exe',OUT/'build/SNANA-audit/src/snlc_fit.F90',OUT/'inputs/snana_audit_pilot.nml']}}
(DEST/'exact_objective_snapshot.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
