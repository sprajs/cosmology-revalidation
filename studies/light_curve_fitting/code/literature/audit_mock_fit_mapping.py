#!/usr/bin/env python3
"""Quantify noisy fitted-SALT versus physical source+dust proxy differences.
These selected/noisy residuals do not identify an exact noiseless mapping.
"""
from pathlib import Path
import json,hashlib,io
import numpy as np,pandas as pd
R=Path(__file__).resolve().parents[3];p=R/'phase2/official/results/snana_mock0001.FITRES.TEXT'
lines=p.read_text().splitlines();cols=next(x.split()[1:] for x in lines if x.startswith('VARNAMES:'))
d=pd.read_csv(io.StringIO('\n'.join(x[4:] for x in lines if x.startswith('SN:'))),sep=r'\s+',names=cols)
m=(d.x1.abs()<3)&(d.c.abs()<.3)&(d.x1ERR<1)&(d.PKMJDERR<2)&(d.FITPROB>.001)&(d.zHD>.025)
d['E_proxy']=d.SIM_AV/d.SIM_RV
d['delta_c_proxy']=d.c-d.SIM_c-d.E_proxy
d['delta_mB_proxy']=d.mB-d.SIM_mB-(d.SIM_RV+1)*d.E_proxy
d['delta_mB_proxy_offset_accounted']=d.delta_mB_proxy+0.12
d['delta_x1']=d.x1-d.SIM_x1

def stats(a):
 a=np.asarray(a);a=a[np.isfinite(a)]
 return {'n':len(a),'mean':float(a.mean()),'std':float(a.std()),'quantiles':dict(zip(['p01','p05','p16','p50','p84','p95','p99'],map(float,np.quantile(a,[.01,.05,.16,.5,.84,.95,.99]))))}
out={'source':str(p.relative_to(R)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'total_fitter_accepted':len(d),'quality_subset':int(m.sum()),'quality_cut':'|x1|<3, |c|<0.3, x1ERR<1, PKMJDERR<2, FITPROB>0.001, zHD>0.025; not full BBC membership','proxy_definitions':{'delta_c_proxy':'fitted_c - SIM_c - SIM_AV/SIM_RV','delta_mB_proxy':'fitted_mB - SIM_mB - (SIM_RV+1)*SIM_AV/SIM_RV','delta_mB_proxy_offset_accounted':'fitted_mB - SIM_mB - (SIM_RV+1)*SIM_AV/SIM_RV +0.12 (known job GENMAG_OFF_GLOBAL=-0.12)','delta_x1':'fitted_x1-SIM_x1'},'all_fit':{c:stats(d[c]) for c in ['delta_c_proxy','delta_mB_proxy','delta_mB_proxy_offset_accounted','delta_x1']},'quality':{c:stats(d.loc[m,c]) for c in ['delta_c_proxy','delta_mB_proxy','delta_mB_proxy_offset_accounted','delta_x1']},'quality_low_noise':{c:stats(d.loc[m&(d.cERR<.03)&(d.mBERR<.04),c]) for c in ['delta_c_proxy','delta_mB_proxy','delta_mB_proxy_offset_accounted','delta_x1']},'limits':['Physical dust proxy is deliberately an unvalidated approximation, not fitted latent truth.','The raw magnitude proxy omits the known -0.12 job luminosity offset; offset-accounted values are reported separately.','Differences include fitting noise, selection, spectral model projection, lensing/peculiar/redshift effects and possibly code-version difference.','This is not proof of biased data, cosmological evolution, or an isolated dust projection bias.','Calibration of noiseless effective SALT summaries needs deterministic high-SNR forward light curves or a conditional response model, not substituting these proxies as truth.']}
(R/'phase2/literature/mock-fit-mapping-audit.json').write_text(json.dumps(out,indent=2)+'\n')
d.loc[:,['CID','zHD','SIM_HOSTLIB_LOGMASS_TRUE','SIM_RV','E_proxy','c','cERR','mBERR','delta_c_proxy','delta_mB_proxy','delta_mB_proxy_offset_accounted','delta_x1']].assign(quality=m).to_csv(R/'phase2/literature/mock-fit-mapping.csv.gz',index=False,compression={'method':'gzip','mtime':0})
print(json.dumps(out,indent=2))
