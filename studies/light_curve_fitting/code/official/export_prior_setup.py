#!/usr/bin/env python3
"""Append true final-iteration initialization, bounds and covariance source parameters."""
from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[3];OUT=ROOT/'phase2/official';DEST=OUT/'portable_pilot'
log=OUT/'diagnostics/snana-audit-pilot-v2.log'
assert 'ENDING PROGRAM GRACEFULLY.' in log.read_text()
data={}
for line in log.read_text().splitlines():
    if not line.startswith('PHASE2_SETUP:'):continue
    q=line.split();v=np.array([float(x) for x in q[2:]]);assert len(v)==16
    data[q[1]]={'final_iteration_INIVAL_x0_x1_c_t0':v[:4].tolist(),'previous_iteration_FITVAL_x0_x1_c_t0':v[4:8].tolist(),'bounds_x0_x1_c_t0':v[8:].reshape(4,2).tolist(),'actual_t0_prior_centre':float(v[3]),'x1_for_model_error':float(v[5]),'note':'Final-iteration covariance frozen using preceding fit; actual INIVAL t0 prior centre differs from original flux-search peak.'}
from audit_fits import read_fit
a=read_fit(OUT/'results/snana_audit_pilot_v2.FITRES.TEXT');b=read_fit(OUT/'results/snana_mask32_pilot.FITRES.TEXT')
cols=['CID','x0','x1','c','PKMJD','FITCHI2','NDOF']
assert a[cols].sort_values('CID').reset_index(drop=True).equals(b[cols].sort_values('CID').reset_index(drop=True))
record={'parameter_order':['x0','x1','c','t0'],'unchanged_standard_fit_verified':True,'objects':data,'source_hashes':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [log,OUT/'build/SNANA-audit-v2/bin/snlc_fit.exe',OUT/'build/SNANA-audit-v2/src/snlc_fit.F90',OUT/'inputs/snana-audit-v2-output-only.patch']}}
(DEST/'exact_prior_setup.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps({'objects':len(data),'output':str(DEST/'exact_prior_setup.json')}))
# Update the fixture contract without regenerating original files.
p=DEST/'contract.json';c=json.loads(p.read_text());c['exact_objective_extension']={'prediction_file':'exact_epoch_predictions.csv','per_object_covariance':'objective_<CID>.npz','setup':'exact_prior_setup.json','precision':'REAL8 predictions and exact accepted MJDs; internal data flux and frozen model-flux covariance ingredients can be float32 as in SNANA.','validation':'r^T C^-1 r + prior equals SNANA objective to <6e-14; both output-only audit builds reproduce unchanged fit parameters/chi2.'}
c['known_limits']=[x for x in c['known_limits'] if not x.startswith('LCPLOT model grid')];c['file_sha256']={str(p.relative_to(DEST)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(DEST.rglob('*')) if p.is_file() and p.name!='contract.json'}
p.write_text(json.dumps(c,indent=2)+'\n')
