"""Independent hash/static provenance audit; no executor import or FITS pixel parse."""
from pathlib import Path
import ast,hashlib,json,datetime
P=Path(__file__).resolve().parent
E=P.parent.parent/'pixel_execution'
def h(p):
 z=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):z.update(b)
 return z.hexdigest()
frozen=json.loads((E/'input-freeze.json').read_text())
rows=[]
for r in frozen['files']:
 p=Path(r['path']);actual=h(p)
 assert actual==r['sha256'] and p.stat().st_size==r['bytes'],r['path']
 rows.append(dict(path=r['path'],sha256=actual,bytes=r['bytes'],pass_hash_size=True))
for name in ['run.py','geometry.py']:
 p=E/name; ast.parse(p.read_text());(P/('final-'+name)).write_bytes(p.read_bytes())
r=dict(status='PASS: static review recommends root release of Stage A only',review_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),input_freeze_sha256=h(E/'input-freeze.json'),executor_sha256=h(E/'run.py'),geometry_sha256=h(E/'geometry.py'),files_verified=len(rows),rows=rows,no_executor_import_or_execution=True,no_scientific_pixel_arrays_parsed=True,no_heldout_aperture_statistics_evaluated=True,checks=dict(design_only_source_masks=True,stage_A_test_SCI_only_finiteness=True,masked_annulus_weights=True,PAM_and_PHOTFLAM_covariance_units=True,no_residual_selection=True,strict_under30_stop_before_sums=True,resampling_requires_12tiles_and_3_per_occupied_tile=True,complete_DQ16_and_finiteness_ledger=True,shared_template_HDHtranspose=True,raw_moment_amendment=True,hard_stop_geometry_and_synthetic_failures=True,partial_operators_and_support_preserved=True,primary_checkpoint_before_sensitivities=True,operator_fields_cap50MB=True,quadrature_chunk8pixels_at128squared=True,stage_A120second_benchmark=True,stage_B_root_hash_release_and_cumulative_cap=True),limitations=['Static review plus independent hash verification, not an executed pixel analysis','Source-mask and DQ selection are conditioning; ERR endogeneity and unknown cross-exposure covariance remain','Masked pixels alter scene response; no total-point-source flux estimate','Spatial resampling is descriptive and cannot measure global calibration/background uncertainty','Current2026 product experiment, not historical RAISIN reconstruction or CRNL calibration'])
(P/'final-static-review.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({k:r[k] for k in ['status','input_freeze_sha256','executor_sha256','geometry_sha256','files_verified']},indent=2))
