#!/usr/bin/env python3
from pathlib import Path
import json,hashlib,datetime
import numpy as np,pandas as pd,h5py
from astropy.io import fits
R=Path(__file__).resolve().parents[3];O=R/'phase2/assumptions'
def entry(p):return {'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()}
checks={}
old=json.load(open(R/'phase2/data_audit/outputs_manifest.json'))
files=old.get('files',old.get('outputs',[])) if isinstance(old,dict) else old
fail=[]
for f in files:
 p=R/f['path']
 if entry(p)['sha256']!=f['sha256']:fail.append(str(p))
checks['first_wave_output_hash_mismatches']=fail;checks['first_wave_outputs_verified']=len(files)
S=R/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES';C=O/'conditional_mask';bit=1<<29
with fits.open(S/'DES-SN5YR_DES_PHOT.FITS.gz',character_as_bytes=True) as src,fits.open(C/'DES-SN5YR_DES_CONDITIONAL_PHOT.FITS',character_as_bytes=True) as dst:
 for name in src[1].data.names:
  a=np.asarray(src[1].data[name]);b=np.asarray(dst[1].data[name]);b=(b&~bit) if name=='PHOTFLAG' else b
  assert np.array_equal(a,b),name
 checks['all_original_PHOT_fields_preserved_except_added_mask_bit']=True
 keep=(dst[1].data['PHOTFLAG']&bit)==0;assign=pd.read_csv(C/'assignments.csv.gz');assert len(assign)==sum(keep)==62702;assert np.array_equal(np.sort(assign.phot_row_one_based.to_numpy()-1),np.flatnonzero(keep));checks['conditional_mask_acceptance_count']=int(sum(keep))
with fits.open(S/'DES-SN5YR_DES_HEAD.FITS.gz',character_as_bytes=True) as src,fits.open(C/'DES-SN5YR_DES_CONDITIONAL_HEAD.FITS',character_as_bytes=True) as dst:
 for name in src[1].data.names:assert np.array_equal(src[1].data[name],dst[1].data[name]),name
 checks['all_original_HEAD_physical_columns_preserved']=True
D=pd.read_csv(O/'off_signal_rows.csv.gz',dtype={'CID':str});assert np.all(abs(D.phase_days)>200);assert np.all(np.isfinite(D.pull));checks['off_signal_extracted_rows']=len(D)
coverage=[]
for ds in ['original/DES','diffimg/DES']:
 for win in [365,200]:
  for mask in ['basic','quality']:
   q=D[(D.dataset==ds)&(abs(D.phase_days)>win)&D[mask]];coverage.append({'dataset':ds,'abs_phase_greater_than_days':win,'mask':mask,'n_rows':len(q),'n_objects':q.CID.nunique(),'available':bool(len(q))})
(O/'test_coverage.json').write_text(json.dumps(coverage,indent=2)+'\n')
T=pd.read_csv(O/'baseline_heldout_rows.csv.gz',dtype={'CID':str});assert not T.duplicated(['CID','phot_row']).any();assert np.all(abs(T.phase_days)>365);assert np.all(T.train_n>=5);checks['heldout_unique_records']=len(T)
P=json.load(open(O/'pixel_noise_results.json'));assert len(P)==3 and all(p['exact_science_crop_verified'] for p in P);checks['independently_acquired_science_weight_matched_exposures']=[p['exposure'] for p in P]
# Finalize conditional manifest after its run log stopped growing.
cm=json.load(open(C/'manifest.json'));cm['files']=[entry(R/f['path']) for f in cm['files'] if (R/f['path']).is_file()];(C/'manifest.json').write_text(json.dumps(cm,indent=2)+'\n')
checks['passed']=not fail;(O/'verification.json').write_text(json.dumps(checks,indent=2)+'\n')
inputs=[R/'phase2/data_audit/photometry_audit.h5',R/'phase2/data_audit/original_metadata.csv.gz',R/'phase2/data_audit/outputs_manifest.json',R/'phase2/data_audit/duplicate_identity.json',R/'phase2/data_audit/exact_duplicates.json',R/'phase2/data_audit/lcplot_audit_mask32.json',R/'phase2/data_audit/published_lcplot_flags.json',R/'data/des-diffimg/DES-SN5YR_DIFFIMG.README',R/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES5YR_SALT3_LCFIT.LCPLOT.gz',S/'DES-SN5YR_DES.README']
inputs=[p for p in inputs if p.exists()]+[R/'sources/repos/RickKessler__SNANA/src/snlc_fit.F90',R/'sources/repos/des-science__DES-SN5YR@1.3/7_PIPPIN_FILES/base_files/lcfit/lcfit_desSMP_5yr.nml',R/'papers/text/2406.05046v1.txt',R/'phase2/official/fitter-assumptions.json',R/'phase2/official/build/SNANA-2fe0f56/src/snlc_sim.c',R/'phase2/official/inputs/SNDATA_ROOT/simlib/DES/DES-SN5YR_DES_FLUXERRMODEL_SIM.DAT']
(O/'inputs_manifest.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':[entry(p) for p in inputs]},indent=2)+'\n')
outputs=sorted([p for p in O.rglob('*') if p.is_file() and p.name!='outputs_manifest.json' and not p.name.endswith('.log')]+list((R/'scripts/phase2/assumptions').glob('*.py'))+[R/'docs/phase2/snana-assumptions.md'])
(O/'outputs_manifest.json').write_text(json.dumps({'schema':'measurement-assumptions-output-manifest-v1','utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'note':'Includes immutable downloaded science/calibrated image planes and all machine outputs; logs excluded because redirected finalizer log can change during hashing. Download URLs/statuses in nested image_access manifests.','files':[entry(p) for p in outputs]},indent=2)+'\n')
print(json.dumps(checks,indent=2));print('Manifested',len(outputs),'files')
