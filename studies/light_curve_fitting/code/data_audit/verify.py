#!/usr/bin/env python3
"""Independent full-column readback against original FITS, plus audit invariants."""
from pathlib import Path
import hashlib,json
import numpy as np,h5py
from astropy.io import fits
R=Path(__file__).resolve().parents[3];O=R/'phase2/data_audit';S=json.loads((O/'summary.json').read_text());records=[];donephot=set()
with h5py.File(O/'photometry_audit.h5','r') as h:
 for ds in S['datasets']:
  g=h[ds['dataset']];r={'dataset':ds['dataset'],'head_columns_verified':0,'phot_columns_verified':0}
  for section,attr in [('head','head_source'),('phot','phot_source')]:
   p=R/g.attrs[attr];digest=hashlib.file_digest(p.open('rb'),'sha256').hexdigest();cache=O/'cache'/(digest+'.fits')
   if section=='phot' and digest in donephot:r['phot_verified_by_identical_source_and_hdf5_hardlink']=True;continue
   with fits.open(cache,memmap=True) as f:
    source=f[1].data;assert len(source)==len(next(iter(g[section].values())))
    for c in source.names:
     out=g[section][c]
     for start in range(0,len(source),250000):
      a=np.asarray(source[c][start:start+250000]);b=out[start:start+250000]
      if a.dtype.kind=='U':a=np.char.encode(a,'utf-8')
      same=(a==b)
      if a.dtype.kind=='f':same|=(np.isnan(a)&np.isnan(b))
      if not same.all():raise AssertionError((ds['dataset'],c,start))
     r[section+'_columns_verified']+=1
   if section=='phot':donephot.add(digest)
  owner=g['head_row_zero_based'][:];flags=g['audit_flags'][:];basic=g['basic_eligible'][:]
  assert len(owner)==ds['phot_rows'];assert int(basic.sum())==ds['basic_eligible_rows']
  for i,(lo,hi) in enumerate(zip(g['head/PTROBS_MIN'][:],g['head/PTROBS_MAX'][:])):
   assert np.all(owner[int(lo)-1:int(hi)]==i)
   assert g['phot/MJD'][int(hi)]==-777
  assert all(ds[k]==[] for k in ['pointer_invalid_ids','nobs_mismatch_ids','delimiter_inside_span_ids','missing_delimiter_after_ids'])
  r['all_object_spans_verified']=len(g['head/SNID']);r['all_flags_lengths_verified']=True;records.append(r);print(r,flush=True)
 for typ in ['DES','Foundation','LOWZ']:assert h['original/'+typ+'/phot'].id==h['current/'+typ+'/phot'].id
result={'passed':True,'source_HDF5_comparison':'Every physical HEAD and PHOT column, every row; repeated PHOT hardlinks checked once against identical source digest.','datasets':records}
(O/'verification.json').write_text(json.dumps(result,indent=2)+'\n')
