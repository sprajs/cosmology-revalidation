#!/usr/bin/env python3
from pathlib import Path
import json,hashlib
import pandas as pd
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'phase2/official';DEST=OUT/'portable_pilot'
log=OUT/'diagnostics/snana-audit-pilot-v3.log';assert 'ENDING PROGRAM GRACEFULLY.' in log.read_text()
records={}
for line in log.read_text().splitlines():
    if line.startswith('PHASE2_COVSET:'):
        q=line.split();records[(q[1],int(q[2]))]=[float(x) for x in q[3:]]
cols=['CID','fit_row_one_based','frozen_MWXT_FLUXERR','frozen_COVMAG_ERR','frozen_MODELFLUX','DATAFLUX_ERR','FUDGEFLUX_ERR']
d=pd.DataFrame([[cid,i]+v for (cid,i),v in records.items()],columns=cols)
d.to_csv(DEST/'exact_covariance_components.csv',index=False,float_format='%.17g')
from audit_fits import read_fit
a=read_fit(OUT/'results/snana_audit_pilot_v3.FITRES.TEXT');b=read_fit(OUT/'results/snana_mask32_pilot.FITRES.TEXT')
cols=['CID','x0','x1','c','PKMJD','FITCHI2','NDOF'];assert a[cols].sort_values('CID').reset_index(drop=True).equals(b[cols].sort_values('CID').reset_index(drop=True))
record={'n_rows':len(d),'n_objects':d.CID.nunique(),'description':'FITINI_COV finaliteration exact covariance ingredients before matrix inversion. MW covariance is outer product of frozen_MWXT_FLUXERR; diagonal data variance is DATAFLUX_ERR squared plus FUDGEFLUX_ERR squared. Model covariance is remainder after subtracting these from supplied frozen covariance. All variables preserve actual source precision; these common-block inputs are generallyfloat32 except COVMAG_ERR.','unchanged_fit_verified':True,'hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [log,OUT/'build/SNANA-audit-v3/bin/snlc_fit.exe',OUT/'build/SNANA-audit-v3/src/snlc_fit.F90',OUT/'inputs/snana-audit-v3-output-only.patch']}}
(DEST/'exact_covariance_components.json').write_text(json.dumps(record,indent=2)+'\n')
p=DEST/'contract.json';c=json.loads(p.read_text());c['file_sha256']={str(p.relative_to(DEST)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(DEST.rglob('*')) if p.is_file() and p.name!='contract.json'};p.write_text(json.dumps(c,indent=2)+'\n');print(json.dumps(record,indent=2))
