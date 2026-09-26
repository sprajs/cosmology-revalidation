from pathlib import Path
import json,hashlib
R=Path.cwd();o=Path(__file__).resolve().parent
H=R/'runs/research_2026_09_26/raisin_hst_pixel_feasibility';s=H/'calwf3_variance_source';p=H/'dark_ramp_raw_pilot'
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
gitblob=lambda f:hashlib.sha1(b'blob '+str(f.stat().st_size).encode()+b'\0'+f.read_bytes()).hexdigest()
checks=[]
for folder in ['raw-processing-source','raw-processing-source-helpers']:
 for r in json.loads((s/folder/'manifest.json').read_text())['files']:
  f=R/r['path'];assert sha(f)==r['sha256'] and gitblob(f)==r['git_blob_sha1'];checks.append(str(f.relative_to(R)))
for r in json.loads((s/'source-3.7.3-manifest.json').read_text()):
 f=s/'source-3.7.3'/Path(r['path']).name
 assert sha(f)==r['sha256'] and gitblob(f)==r['git_blob'];checks.append(str(f.relative_to(R)))
# Use saved metadata only; no FITS arrays are opened.
d=json.loads((p/'header-comparison.json').read_text());pairs=[]
keys=['sampnum','samptime','deltatim','samp_pixvalue','time_pixvalue','time_naxis','samp_naxis','science_geometry']
for r in d['pairs']:
 a=r['dark']['groups'][:8];b=r['science']['groups']
 assert len(a)==len(b)==8
 for x,y in zip(a,b):
  for k in keys:assert x[k]==y[k],(r['visit'],k,x[k],y[k])
 assert [x['extver'] for x in a]==list(range(16,8,-1))
 assert [x['extver'] for x in b]==list(range(8,0,-1))
 pairs.append({'visit':r['visit'],'dark':r['dark_root'],'science':r['science_root'],'times':[x['samptime'] for x in a],'recorded_prefix_fields_exact':True})
paths=[s/f/'manifest.json' for f in ['raw-processing-source','raw-processing-source-helpers']]+[s/'source-3.7.3-manifest.json',p/'header-comparison.json',p/'header-inspection.json',p/'root-prefix-review/result.json']+[s/'source-3.7.3'/f for f in ['refdata.c','imageio.c']]
res={'status':'PASS_SOURCE_HASHES_AND_SAVED_PREFIX_METADATA','source_files_sha_and_git_blob_verified':len(checks),'checked_source_files':checks,'pairs':pairs,'no_network':True,'no_pixel_arrays_read':True,'no_build_or_native_execution':True,'qualification':'Recorded sampling identity is not detector-state or calibrated-operator identity. Unknown template quality remains unresolved.','inputs':{str(f.relative_to(R)):sha(f) for f in paths}}
(o/'metadata-source-result.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps({'status':res['status'],'sourcefiles':len(checks),'paircount':len(pairs)}))
