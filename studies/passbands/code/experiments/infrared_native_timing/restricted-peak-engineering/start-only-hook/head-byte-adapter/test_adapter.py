"""Read-only input diagnosis plus isolated exact-byte fixtures. No native fits."""
from pathlib import Path
import json,hashlib
import numpy as np
from astropy.io import fits
from adapter import patch_head,sha
O=Path(__file__).resolve().parent;H=O.parent;P=H.parent.parent
src=P/'readme-adapted/ledger/PTE/PTE_HEAD.FITS';failed=P/'derived-start-only/null/PTE/PTE_HEAD.FITS'
raw=src.read_bytes();bad=failed.read_bytes();assert len(raw)==len(bad)
diff=np.flatnonzero(np.frombuffer(raw,'u1')!=np.frombuffer(bad,'u1'));assert len(diff)==712
with fits.open(src,memmap=False,character_as_bytes=True) as h, fits.open(failed,memmap=False,character_as_bytes=True) as k:
 a,b=h[1].data,k[1].data;loc=h[1].fileinfo()['datLoc'];stride=h[1].header['NAXIS1'];strings=[];stringbytes=set();checks={}
 for n in a.dtype.names:
  d,o=a.dtype.fields[n][:2]
  if d.kind=='S':
   strings.append(n)
   for i in range(len(a)):stringbytes.update(range(loc+i*stride+o,loc+i*stride+o+d.itemsize))
  checks[n]=bool(np.array_equal(a[n],b[n]))
 assert all(checks[n] for n in checks if n not in strings)
 assert set(diff)<=stringbytes and all(raw[i]==32 and bad[i]==0 for i in diff)
 nullpeaks={bytes(x['SNID']).decode().strip():float(x['PEAKMJD']) for x in a}
# Use the already frozen nominal joint final peaks, without scoring them.
peaks={}
for line in (P/'fits-restricted/joint12/native.log').read_text().splitlines():
 v=line.split()
 if v and v[0]=='CSP_OBJECTIVE:' and v[2]=='12':peaks[v[1]]=float(v[11])
assert set(peaks)==set(nullpeaks)
fixtures=O/'fixtures';fixtures.mkdir(exist_ok=False)
n=patch_head(src,fixtures/'null_HEAD.FITS',nullpeaks,True)
e=patch_head(src,fixtures/'estimated_HEAD.FITS',peaks)
assert n['source_sha256']==n['target_sha256'] and len(e['allowed_bytes'])==32
# Same strict ordinary-Astropy column comparison as the frozen identity checker.
with fits.open(src) as a,fits.open(fixtures/'null_HEAD.FITS') as b:
 assert np.array_equal(a[1].data,b[1].data)
 for col in a[1].data.dtype.names:assert np.array_equal(a[1].data[col],b[1].data[col])
report={'gate_pass':True,'native_calls':0,'diagnosis':{'changed_bytes':len(diff),'all_changes_space_to_NUL':True,'string_columns':strings,'column_equal':checks,'numeric_columns_equal':True,'source_sha256':sha(src),'failed_null_sha256':sha(failed)},'null':n,'estimated':e}
(O/'test-result.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'gate_pass':True,'string_columns':strings,'all_numeric_columns_exact':True,'bad_padding_bytes':len(diff),'null_hash_exact':True,'estimated_changed_bytes':len(e['changed_bytes']),'allowed_peak_bytes':len(e['allowed_bytes'])},indent=2))
