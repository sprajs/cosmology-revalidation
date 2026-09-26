"""Independent prefix byte/header review and authorized second science null.
Uses parent's prefix function exactly; never executes or scores the dark pair.
"""
from pathlib import Path
import hashlib,json,os,time,subprocess,importlib.util
import numpy as np
from astropy.io import fits
O=Path(__file__).resolve().parent;H=O.parents[1];R=H.parents[2];OLD=H/'calwf3_native_replay';script=R/'scripts/research_2026_09_26/prepare_calwf3_dark_prefix.py';run=R/'scripts/research_2026_09_26/run_calwf3_dark_prefix.py'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(n,x):(O/n).write_text(json.dumps(x,indent=2)+'\n')
proto=json.loads((OLD/'replay-protocol.json').read_text());constructor=H/'calwf3_dark_prefix/constructor-result.json';records=json.loads(constructor.read_text());binary=OLD/'build-fortran/pkg/wfc3/calwf3.e';root='icxoi4hgq';source=H/'dark_ramp_raw_pilot/files'/f'{root}_raw.fits';reference=OLD/'work'/root/f'{root}_flt.fits'
inputs=[script,run,constructor,OLD/'replay-protocol.json',binary,source,reference,Path(__file__)]+[Path(p) for p in proto['input_hashes']]
assert not (O/'protocol.json').exists()
save('protocol.json',{'status':'Frozen before second native null; parent explicitly authorized one template science replay,60seconds; no dark SCI scores.','source_constructor':str(script),'native_binary':str(binary),'inputs':{str(p):sha(p) for p in inputs},'native_root':root,'environment':proto['environment'],'native_timeout_seconds':60,'output_bytes_cap':350000000,'component_gate':'shape,dtype,unit,data-byte exact SCI ERR DQ SAMP TIME versus previous native template; archive DQ not reference','constructor_gate':'Independently compare all120 retained original-data blocks and full header-card differences in existing3 prepared files; then import unchanged sameprefixfunction for8read template wholefileidentity.'})
for p,h in proto['input_hashes'].items():assert sha(p)==h,p
assert sha(binary)=='230dc5e2c760217e7924f0c11153278dd43e10c22fa73828cad8aff0a9c97d00'
checks=[]
for record in records:
 s=Path(record['source']);t=Path(record['target']);assert sha(s)==record['source_sha256'];assert sha(t)==record['target_sha256'];sb=s.read_bytes();tb=t.read_bytes();different=[];blocks=0
 with fits.open(s,do_not_scale_image_data=True,memmap=True) as a,fits.open(t,do_not_scale_image_data=True,memmap=True) as b:
  n=int(a[0].header['NSAMP']);assert n in [8,16];assert len(b)==41
  # Parse exact card images independently; permit only intended primary keyword changes.
  keys=set(a[0].header.keys())|set(b[0].header.keys());primarychanges=[]
  for key in keys:
   if key in ('','COMMENT','HISTORY'):assert a[0].header.get(key)==b[0].header.get(key);continue
   if a[0].header.get(key)!=b[0].header.get(key):primarychanges.append(key)
  assert set(primarychanges)==(set() if n==8 else {'NSAMP','NEXTEND','EXPTIME','EXPEND'})
  for newver in range(1,9):
   oldver=n-8+newver
   for name in ('SCI','ERR','DQ','SAMP','TIME'):
    old=a[name,oldver];new=b[name,newver];oi=old.fileinfo();ni=new.fileinfo();assert sb[oi['datLoc']:oi['datLoc']+oi['datSpan']]==tb[ni['datLoc']:ni['datLoc']+ni['datSpan']];blocks+=1
    oldcards=list(old.header.cards);newcards=list(new.header.cards);assert len(oldcards)==len(newcards)
    for ac,bc in zip(oldcards,newcards):
     assert ac.keyword==bc.keyword
     if ac.keyword=='EXTVER':assert bc.value==newver
     else:assert ac.image==bc.image,(name,oldver,ac.keyword)
  if n==8:assert sb==tb
  else:
   assert b[0].header['NSAMP']==8 and b[0].header['NEXTEND']==40;assert b[0].header['EXPTIME']==a['SCI',9].header['SAMPTIME'];assert b[0].header['EXPEND']==a[0].header['EXPSTART']+b[0].header['EXPTIME']/86400
 checks.append({'source':str(s),'target':str(t),'data_blocks_exact':blocks,'changed_primary_keywords':sorted(primarychanges),'wholefileidentity_if8':n!=8 or sb==tb})
save('constructor-independent-review.json',{'pass':True,'files':checks,'total_data_blocks':sum(x['data_blocks_exact'] for x in checks)})
# Exact parent's constructor function, no source edits.
spec=importlib.util.spec_from_file_location('parent_prefix_constructor',script);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
work=O/'work'/root;work.mkdir(parents=True,exist_ok=False);copy=work/f'{root}_raw.fits';copy_result=module.prefix(source,copy);assert copy_result['null_whole_file_exact'] and not copy_result['edits'];save('template-constructor-result.json',copy_result)
env=os.environ.copy();env.update(proto['environment']);start=time.monotonic();error=None
try:
 with (work/'native.log').open('x') as f:proc=subprocess.run([str(binary),copy.name],cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=60)
 rc=proc.returncode
except Exception as e:rc=None;error=repr(e)
runtime=time.monotonic()-start;save('execution.json',{'root':root,'returncode':rc,'error':error,'seconds':runtime,'binary_sha256':sha(binary),'raw_sha256':sha(copy),'environment':proto['environment']});assert rc==0,error
native=work/f'{root}_flt.fits';components=[]
with fits.open(native,memmap=False) as a,fits.open(reference,memmap=False) as b:
 for name in ('SCI','ERR','DQ','SAMP','TIME'):
  x,y=a[name,1].data,b[name,1].data;same=x.shape==y.shape and x.dtype==y.dtype and x.tobytes()==y.tobytes();unit=a[name,1].header.get('BUNIT')==b[name,1].header.get('BUNIT');components.append({'name':name,'shape':list(x.shape),'dtype':str(x.dtype),'bit_exact':same,'unit_equal':unit,'unit':a[name,1].header.get('BUNIT'),'unequal_values':int(np.count_nonzero(x!=y))})
result={'pass':all(c['bit_exact'] and c['unit_equal'] for c in components),'native_seconds':runtime,'whole_RAW_byte_identity':sha(source)==sha(copy),'components':components,'native_sha256':sha(native),'reference_sha256':sha(reference),'comparison':'previous native template, no archive DQ test','dark_execution_or_scoring':False,'constructor_review_pass':True,'output_bytes':sum(p.stat().st_size for p in O.rglob('*') if p.is_file())}
save('result.json',result);assert result['pass'];assert result['output_bytes']<=350000000
save('manifest.json',{str(p):sha(p) for p in O.rglob('*') if p.is_file() and p.name!='manifest.json'});print(json.dumps(result,indent=2))
