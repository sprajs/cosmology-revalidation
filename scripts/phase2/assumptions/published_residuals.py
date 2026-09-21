#!/usr/bin/env python3
from pathlib import Path
import json,numpy as np,pandas as pd
R=Path(__file__).resolve().parents[3];O=R/'phase2/assumptions';P=pd.read_csv(R/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES5YR_SALT3_LCFIT.LCPLOT.gz',sep=r'\s+',comment='#',dtype={'SNID':str});rows=[]
for (cid,band),d in P.groupby(['SNID','BAND']):
 m=d[d.DATA_MODEL==0].sort_values('MJD').drop_duplicates('MJD');v=d[d.DATA_MODEL!=0].copy();v=v[v.MJD.between(m.MJD.min(),m.MJD.max())]
 if not len(v):continue
 v['interpolated_model_flux']=np.interp(v.MJD,m.MJD,m.FLUXCAL);v['measurement_pull']=(v.FLUXCAL-v.interpolated_model_flux)/v.FLUXCALERR;v['signed_sqrt_diagonal_chi2']=np.sign(v.measurement_pull)*np.sqrt(np.maximum(v.CHI2,0));rows.append(v)
D=pd.concat(rows,ignore_index=True);D.to_csv(O/'published_residual_rows.csv.gz',index=False,compression={'method':'gzip','mtime':0});out=[]
for (dm,band),d in D.groupby(['DATA_MODEL','BAND']):
 for col in ['measurement_pull','signed_sqrt_diagonal_chi2']:
  x=d[col].to_numpy();out.append({'DATA_MODEL':int(dm),'BAND':band,'statistic':col,'n':len(x),'objects':d.SNID.nunique(),'mean':float(np.mean(x)),'rms':float(np.sqrt(np.mean(x*x))),'mad_sigma':float(1.4826*np.median(abs(x-np.median(x)))),'abs_gt3_fraction':float(np.mean(abs(x)>3)),'abs_gt5_fraction':float(np.mean(abs(x)>5))})
(O/'published_residual_statistics.json').write_text(json.dumps({'warning':'In-sample and clipped. Model curves are sampled every2days and linearly interpolated; this approximation is not an independent SALT refit. Published CHI2 is diagonal total-error contribution, not full covariance likelihood.','input_data_rows':int((P.DATA_MODEL!=0).sum()),'interpolation_supported_rows':len(D),'statistics':out},indent=2)+'\n')
print(json.dumps(out,indent=2))
