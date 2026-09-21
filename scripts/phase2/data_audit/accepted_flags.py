#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np,h5py,json,hashlib
R=Path(__file__).resolve().parents[3];O=R/'phase2/data_audit';p=R/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES5YR_SALT3_LCFIT.LCPLOT.gz'
pub=pd.read_csv(p,sep=r'\s+',comment='#',dtype={'SNID':str});pub=pub[pub.DATA_MODEL!=0]
rows=[]
with h5py.File(O/'photometry_audit.h5') as h:
 g=h['original/DES'];ids=np.char.strip(g['head/SNID'][:].astype('U'));ix={cid:j for j,cid in enumerate(ids)}
 for cid,l in pub.groupby('SNID',sort=False):
  j=ix[cid];a=int(g['head/PTROBS_MIN'][j])-1;b=int(g['head/PTROBS_MAX'][j]);mjd=g['phot/MJD'][a:b];mjd=np.array([float(f'{x:.3f}') for x in mjd.astype(np.float32)]);band=np.char.strip(g['phot/BAND'][a:b].astype('U'));flag=g['phot/PHOTFLAG'][a:b];flux=g['phot/FLUXCAL'][a:b]
  look={}
  for k,key in enumerate(zip(mjd,band)):look.setdefault(key,[]).append(k)
  for r in l.itertuples():
   kk=look.get((r.MJD,r.BAND),[])
   if len(kk)>1:kk=[k for k in kk if np.isclose(flux[k],r.FLUXCAL,rtol=5e-5,atol=1e-4)]
   vals=np.unique(flag[kk]);rows.append({'CID':cid,'MJD':r.MJD,'BAND':r.BAND,'DATA_MODEL':r.DATA_MODEL,'raw_match_count':len(kk),'unique_flag_count':len(vals),'PHOTFLAG':int(vals[0]) if len(vals)==1 else None,'possible_flags':','.join(str(int(v)) for v in vals)})
d=pd.DataFrame(rows);d.to_csv(O/'published_lcplot_flags.csv.gz',index=False,compression={'method':'gzip','mtime':0});summary={}
for dm in [1,-1]:
 t=d[d.DATA_MODEL==dm];v=t[t.unique_flag_count==1].PHOTFLAG.astype(int).to_numpy();summary[str(dm)]={'rows':len(t),'unmatched':int((t.raw_match_count==0).sum()),'ambiguous_distinct_flags':int((t.unique_flag_count>1).sum()),'bit_counts':{str(1<<k):int(((v&(1<<k))!=0).sum()) for k in range(16)},'flags':{str(int(k)):int(c) for k,c in zip(*np.unique(v,return_counts=True))}}
(O/'published_lcplot_flags.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
