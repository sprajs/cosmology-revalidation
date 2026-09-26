"""Independent saved CALWF3 replay audit; no native execution, root-code imports or edits."""
from pathlib import Path
import hashlib,json,collections,warnings
import numpy as np
from astropy.io import fits
O=Path(__file__).resolve().parent;P=O.parent;H=P.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()
proto=json.loads((P/'replay-protocol.json').read_text());acq=json.loads((P/'acquisition-protocol.json').read_text());results=json.loads((P/'acquisition-result.json').read_text());inv=json.loads((H/'calwf3_native_feasibility/inventory.json').read_text())
assert sha(H/'calwf3_native_feasibility/inventory.json')==acq['inventory_sha256']
treepath=H/'calwf3_variance_source/tree-3.7.3.json';assert sha(treepath)==acq['tree_sha256'];tree=json.loads(treepath.read_text());assert not tree.get('truncated',False)
blobs={x['path']:x for x in tree['tree'] if x['type']=='blob'};sourcefiles=[p for p in (P/'source').rglob('*') if p.is_file()]
assert len(sourcefiles)==len(blobs)==results['regular_source_files']==812
source_records=[]
for p in sourcefiles:
 assert not p.is_symlink();name=str(p.relative_to(P/'source'));d=p.read_bytes();actual=hashlib.sha1(b'blob '+str(len(d)).encode()+b'\0'+d).hexdigest();assert actual==blobs[name]['sha'],name
 source_records.append({'path':name,'git_blob':actual,'sha256':hashlib.sha256(d).hexdigest()})
for p,h in proto['input_hashes'].items():assert sha(Path(p))==h,p
for row in results['rows']:
 p=P/'files'/row['name'];assert p.stat().st_size==row['bytes'] and sha(p)==row['sha256']
refs=[]
for name,entry in inv['refs'].items():
 p=P/'iref'/name;assert p.is_symlink();target=p.resolve();assert str(target) in proto['input_hashes'];assert sha(target)==proto['input_hashes'][str(target)]
 refs.append({'name':name,'target':str(target),'sha256':sha(target)})
assert len(refs)==14
build=json.loads((P/'build-fortran-result.json').read_text());assert build['pass'] and sha(Path(build['binary']))==build['binary_sha256']==proto['input_hashes'][build['binary']]
for j in build['jobs']:assert j['returncode']==0
failure=json.loads((P/'build-failure.json').read_text());assert failure['completed'][0]['returncode']!=0
bp=json.loads((P/'build-fortran-protocol.json').read_text());assert sha(Path(bp['compiler_path']))==bp['compiler_sha256']
f951=R/'phase2/official/build/sysroot/usr/lib/gcc/x86_64-pc-linux-gnu/16/f951';assert sha(f951)==bp['f951_sha256']
all_results=[];warnings_seen=[];inputs=set(map(Path,proto['input_hashes']))|{treepath,P/'acquisition-protocol.json',P/'acquisition-result.json',P/'build-fortran-protocol.json',P/'build-fortran-result.json',P/'build-failure.json',P/'replay-protocol.json',P/'replay-result.json'}
for root in proto['roots']:
 ex=json.loads((P/(root+'-execution.json')).read_text());assert ex['returncode']==0
 w=P/'work'/root;raw=Path(inv['objects'][root]['raw']['path']);archive=Path(inv['objects'][root]['flt']['path']);output=w/(root+'_flt.fits')
 assert sha(raw)==sha(w/raw.name)==inv['objects'][root]['raw']['sha256'];assert sha(archive)==inv['objects'][root]['flt']['sha256']
 for name,row in ex['outputs'].items():assert (w/name).stat().st_size==row['bytes'] and sha(w/name)==row['sha256']
 text=(w/'native.log').read_text();assert 'Version 3.7.3 (Jan-07-2026)' in text
 components=[];headercheck={};headerdiff=[]
 with warnings.catch_warnings(record=True) as ww:
  warnings.simplefilter('always')
  with fits.open(output,memmap=False) as a, fits.open(archive,memmap=False) as b:
   keys=['CAL_VER',*inv['objects'][root]['raw_refs'],*inv['objects'][root]['archive_flt_switches']]
   for k in keys:headercheck[k]={'replay':a[0].header.get(k),'archive':b[0].header.get(k),'exact':a[0].header.get(k)==b[0].header.get(k)}
   for k in sorted(set(a[0].header)|set(b[0].header)):
    if k not in ['HISTORY','COMMENT',''] and a[0].header.get(k)!=b[0].header.get(k):headerdiff.append({'key':k,'replay':a[0].header.get(k),'archive':b[0].header.get(k)})
   for name in ['SCI','ERR','DQ','SAMP','TIME']:
    x=a[name,1].data;y=b[name,1].data;assert x.shape==y.shape==(1014,1014) and x.dtype==y.dtype
    unit=a[name,1].header.get('BUNIT')==b[name,1].header.get('BUNIT');assert unit
    bits=x.tobytes()==y.tobytes();rec={'extension':name,'pixels':int(x.size),'dtype':str(x.dtype),'unit':a[name,1].header.get('BUNIT'),'units_exact':unit,'bits_exact':bits,'array_hash_replay':hashlib.sha256(x.tobytes()).hexdigest(),'array_hash_archive':hashlib.sha256(y.tobytes()).hexdigest(),'unequal_numeric_pixels':int(np.count_nonzero(x!=y))}
    if name!='DQ':assert bits
    else:
     xx=x.astype(np.uint32);yy=y.astype(np.uint32);xor=xx^yy;uu,cc=np.unique(xor,return_counts=True);m=xor!=0
     rec['xor_histogram']={str(int(k)):int(v) for k,v in zip(uu,cc)};rec['archive_4096_count']=int(np.count_nonzero(yy&4096));rec['native_4096_count']=int(np.count_nonzero(xx&4096));rec['every_other_bit_exact']=bool(np.array_equal(xx&np.uint32(0xffff^4096),yy&np.uint32(0xffff^4096)));rec['differences_only_archive_adds4096']=bool(np.all(yy[m]==(xx[m]|4096)) and np.all((xx[m]&4096)==0))
     assert rec['every_other_bit_exact'] and rec['differences_only_archive_adds4096'] and rec['native_4096_count']==0
     assert rec['unequal_numeric_pixels']=={'icxoi1bcq':3897,'icxoi4hgq':4262}[root]
     coords=np.argwhere(m);np.save(O/(root+'-dq4096-coordinates.npy'),coords)
    components.append(rec)
  warnings_seen+=list(map(lambda w:str(w.message),ww))
 assert all(x['exact'] for x in headercheck.values())
 all_results.append({'root':root,'native_returncode':ex['returncode'],'raw_copy_exact':True,'all_declared_reference_calibration_switch_headers_equal':True,'header_checks':headercheck,'other_primary_header_differences':headerdiff,'components':components,'full_product_exact':False,'SCI_ERR_SAMP_TIME_exact':True})
 inputs|={P/(root+'-execution.json'),raw,archive,*[p for p in w.iterdir() if p.is_file()]}
reported=json.loads((P/'replay-result.json').read_text());assert reported['primary_exact'] is False
for x in reported['comparisons']:
 own=next(r for r in all_results if r['root']==x['root'])
 for m in x['metrics']:
  mine=next(r for r in own['components'] if r['extension']==m['HDU']);assert mine['bits_exact']==m['bit_exact'] and mine['unequal_numeric_pixels']==m['unequal_values']
res={'independent_review_pass':True,'native_calls':0,'source_commit':acq['source_commit'],'source_regular_Git_blobs_verified':812,'pinned_input_hashes_verified':len(proto['input_hashes']),'reference_files_verified':len(refs),'reference_files':refs,'source_binary_build_identity_verified':True,'initial_Fortran_build_failure_preserved':True,'CALWF3_version':'3.7.3','whole_product_gate_remains_false':True,'root_result_arithmetic_reproduced':True,'objects':all_results,'warnings':sorted(set(warnings_seen)),'scope':'Exact current-archive CALWF3 SCI/ERR/SAMP/TIME components for these2 RAWs and fixed refs, not whole archive/DQ/Drizzle or historical author reduction.'}
(O/'result.json').write_text(json.dumps(res,indent=2,allow_nan=False)+'\n');(O/'source-blobs.json').write_text(json.dumps(source_records,indent=2)+'\n')
inputs|={O/'review.py',O/'result.json',O/'source-blobs.json'};inputs|=set(O.glob('*-coordinates.npy'))
(O/'manifest.json').write_text(json.dumps({'files':{str(p.relative_to(R)):sha(p) for p in sorted(inputs)}},indent=2)+'\n')
print(json.dumps({k:v for k,v in res.items() if k not in ['objects','reference_files','warnings']},indent=2))
