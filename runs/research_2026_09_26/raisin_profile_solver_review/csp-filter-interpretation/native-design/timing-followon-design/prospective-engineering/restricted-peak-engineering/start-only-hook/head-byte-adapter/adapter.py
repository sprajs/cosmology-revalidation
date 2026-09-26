"""Exact byte-preserving FITS HEAD peak adapter; no native code or serialization."""
from pathlib import Path
import hashlib,json,struct
import numpy as np
from astropy.io import fits

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def require(ok,msg):
 if not ok:raise RuntimeError(msg)

def patch_head(source,target,peaks,require_null=False):
 source,target=Path(source),Path(target)
 require(not target.exists(),'target already exists')
 raw=source.read_bytes()
 with fits.open(source,memmap=False,character_as_bytes=True) as hd:
  require(len(hd)==2 and isinstance(hd[1],fits.BinTableHDU),'expected single BINTABLE')
  require(all('CHECKSUM' not in h.header and 'DATASUM' not in h.header for h in hd),'checksum extension requires distinct implementation')
  t=hd[1];a=t.data;names=a.dtype.names
  require(t.columns['PEAKMJD'].format=='E','PEAKMJD must be 1E float32')
  require(t.columns['PEAKMJD'].bscale is None and t.columns['PEAKMJD'].bzero is None,'scaled PEAKMJD unsupported')
  dtype,offset=a.dtype.fields['PEAKMJD'][:2]
  require(dtype==np.dtype('>f4') and dtype.itemsize==4,'raw peak dtype')
  stride=int(t.header['NAXIS1']);count=int(t.header['NAXIS2']);loc=int(t.fileinfo()['datLoc'])
  require(stride==a.dtype.itemsize and count==len(a) and loc+stride*count<=len(raw),'row layout')
  cids=[bytes(x).decode('ascii').strip() for x in a['SNID']]
  require(len(cids)==len(set(cids)) and set(cids)==set(peaks),'CID coverage')
  oldvalues={c:float(a['PEAKMJD'][i]) for i,c in enumerate(cids)}
  before={n:a[n].copy() for n in names}
 out=bytearray(raw);allowed=set();rows=[]
 for i,c in enumerate(cids):
  val=float(np.float32(peaks[c]));require(np.isfinite(val),'nonfinite peak')
  start=loc+i*stride+offset;allowed.update(range(start,start+4))
  require(struct.unpack('>f',raw[start:start+4])[0]==oldvalues[c],'decoded peak/layout mismatch')
  out[start:start+4]=struct.pack('>f',val)
  rows.append({'CID':c,'row_index':i,'byte_start':start,'byte_stop':start+4,'old_peak':oldvalues[c],'requested_peak_R8':float(peaks[c]),'stored_peak_R4':val})
 diff=[i for i,(x,y) in enumerate(zip(raw,out)) if x!=y]
 require(len(out)==len(raw) and set(diff)<=allowed,'nonpeak raw byte mutation')
 if require_null:require(bytes(out)==raw,'null file is not bit-identical')
 with target.open('xb') as f:f.write(out)
 with fits.open(target,memmap=False,character_as_bytes=True) as hd:
  after=hd[1].data
  require(after.dtype==a.dtype,'HEAD dtype changed')
  for n in names:
   if n!='PEAKMJD':require(np.array_equal(before[n],after[n]),'nonpeak column changed '+n)
  require(all(float(after['PEAKMJD'][i])==row['stored_peak_R4'] for i,row in enumerate(rows)),'peak readback changed')
 return {'source':str(source),'target':str(target),'source_sha256':sha(source),'target_sha256':sha(target),'file_bytes':len(raw),'data_offset':loc,'row_stride':stride,'peak_field_offset':offset,'rows':rows,'allowed_bytes':sorted(allowed),'changed_bytes':diff,'nonpeak_bytes_exact':True,'all_nonpeak_columns_exact':True,'null_bit_exact':bytes(out)==raw}
