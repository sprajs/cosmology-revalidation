#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np,h5py,json
R=Path(__file__).resolve().parents[3];O=R/'phase2/data_audit';d=pd.read_csv(O/'duplicate_epochs.csv.gz',dtype={'snid':str});d=d[d.dataset=='original/DES'];acc=pd.read_csv(O/'duplicate_exposure_published_acceptance.csv.gz',dtype={'snid':str});records=[]
with h5py.File(O/'photometry_audit.h5') as h:
 p=h['original/DES/phot']
 for key,g in d.groupby(['snid','IMGNUM','MJD','band']):
  idx=np.sort(g.phot_row.to_numpy(int)-1);differ=[];unique={}
  for c in p:
   a=p[c][idx];n=len(np.unique(a));unique[c]=n
   if n>1:differ.append(c)
  cid,img,mjd,band=key;accept=acc[(acc.snid==cid)&(acc.IMGNUM==img)&(acc.MJD==mjd)&(acc.band==band)]
  records.append({'CID':cid,'IMGNUM':int(img),'MJD':mjd,'BAND':band,'raw_rows':len(idx),'source_rows':','.join(str(x+1) for x in idx),'differing_columns':','.join(differ),'identical_every_phot_column':not differ,'identical_flux':unique['FLUXCAL']==1,'identical_detector_identity':all(unique[x]==1 for x in ['CCDNUM','FIELD','XPIX','YPIX']),'fluxerr_min':float(p['FLUXCALERR'][idx].min()),'fluxerr_max':float(p['FLUXCALERR'][idx].max()),'published_accepted_rows':int(accept.published_accepted_rows.sum()),'representative_first_row':int(idx[0]+1),'representative_maxerr_row':int(idx[np.argmax(p['FLUXCALERR'][idx])]+1)})
r=pd.DataFrame(records);r.to_csv(O/'duplicate_identity.csv.gz',index=False,compression={'method':'gzip','mtime':0});s={'groups':len(r),'all_fields_identical_groups':int(r.identical_every_phot_column.sum()),'identical_flux_and_detector_groups':int((r.identical_flux&r.identical_detector_identity).sum()),'different_flux_groups':int((~r.identical_flux).sum()),'different_detector_groups':int((~r.identical_detector_identity).sum()),'differing_field_patterns':r.differing_columns.value_counts(dropna=False).to_dict(),'groups_with_published_accepted':int((r.published_accepted_rows>0).sum()),'distinct_published_affected_objects':int(r.loc[r.published_accepted_rows>0,'CID'].nunique()),'max_fractional_error_spread':float((r.fluxerr_max/r.fluxerr_min-1).max()),'no_mutation':'Source rows have not been removed or altered; refit arms remain gated.'};(O/'duplicate_identity.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2))
