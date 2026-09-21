#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,numpy as np,h5py,json
R=Path(__file__).resolve().parents[3];O=R/'phase2/data_audit';d=pd.read_csv(O/'duplicate_epochs.csv.gz',dtype={'snid':str});d=d[d.dataset=='original/DES'].sort_values('phot_row');idx=d.phot_row.to_numpy(int)-1
with h5py.File(O/'photometry_audit.h5') as h:
 p=h['original/DES/phot'];f=pd.DataFrame({c:p[c][idx] for c in p});f.insert(0,'CID',d.snid.to_numpy());hashes=pd.util.hash_pandas_object(f,index=False).to_numpy();f['row_hash']=hashes;f['phot_row']=idx+1
 # Verify complete equality within hash groups; hashes only accelerate grouping.
 records=[]
 for _,g in f.groupby('row_hash'):
  if len(g)<2:continue
  cols=[c for c in g if c not in ['row_hash','phot_row']];assert (g[cols]==g[cols].iloc[0]).all().all();rep=int(g.phot_row.min())
  for r in g.itertuples():records.append({'CID':r.CID,'phot_row':int(r.phot_row),'representative_phot_row':rep,'redundant_exact_row':int(r.phot_row)!=rep,'group_size':len(g),'MJD':float(r.MJD),'IMGNUM':int(r.IMGNUM),'BAND':r.BAND.decode().strip(),'PHOTFLAG':int(r.PHOTFLAG)})
r=pd.DataFrame(records);r.to_csv(O/'exact_duplicate_rows.csv.gz',index=False,compression={'method':'gzip','mtime':0})
a=pd.read_csv(O/'duplicate_exposure_published_acceptance.csv.gz',dtype={'snid':str});accepted={(x.snid,x.MJD,x.band) for x in a[a.published_accepted_rows>0].itertuples()};r['in_published_accepted_raw_group']=[(x.CID,x.MJD,x.BAND) in accepted for x in r.itertuples()]
s={'full_row_exact_duplicate_groups':int(r.representative_phot_row.nunique()),'exact_rows_total':len(r),'redundant_exact_rows':int(r.redundant_exact_row.sum()),'affected_objects':int(r.CID.nunique()),'exact_duplicate_groups_with_public_accepted_raw_group':int(r.loc[r.in_published_accepted_raw_group,'representative_phot_row'].nunique()),'redundant_exact_rows_in_public_accepted_raw_groups':int((r.redundant_exact_row&r.in_published_accepted_raw_group).sum()),'interpretation':'Broad four-row exposure-key groups usually comprise two distinct image-metadata records repeated twice. Exact-all-column duplicates can be challenged separately without merging distinct metadata records; scientific effects still require controlled refit.'};(O/'exact_duplicates.json').write_text(json.dumps(s,indent=2)+'\n');print(json.dumps(s,indent=2))
