"""Independent metadata/source audit only; no likelihoods, photometric values or native calls."""
from pathlib import Path
from collections import Counter
import hashlib,json
import numpy as np
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[6]
P=ROOT/'phase2/pte'
OUT=Path(__file__).resolve().parent

def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
rows=[line.split() for line in (P/'inputs-v2/cadence.simlib').read_text().splitlines() if line.startswith('S:')]
cad=sorted((float(x[1]),x[3]) for x in rows)
with fits.open(P/'inputs-v2/snoopy.B18/SNooPy_B18.fits',memmap=True) as h:
 model_phase=np.asarray(h['TREST-GRID'].data['TREST'],float)
with fits.open(P/'inputs-v2/kcor.fits',memmap=True) as h:
 kphase=np.asarray(h['KCOR'].data['Trest'],float)
 kz=np.asarray(h['KCOR'].data['Redshift'],float)
a=max(model_phase.min(),kphase.min()); b=min(model_phase.max(),kphase.max())
ledger=[]
with fits.open(P/'readme-adapted/ledger/PTE/PTE_HEAD.FITS',memmap=True) as hh, fits.open(P/'readme-adapted/ledger/PTE/PTE_PHOT.FITS',memmap=True) as pp:
 H=hh[1].data; M=pp[1].data['MJD']; B=pp[1].data['BAND']
 # No FLUXCAL, FLUXCALERR, model mean, distance or fitted peak columns are read.
 for h in H:
  lo=int(h['PTROBS_MIN'])-1;hi=int(h['PTROBS_MAX']); t=np.asarray(M[lo:hi],float)
  bands=[str(x).strip() for x in B[lo:hi]]
  z=float(h['REDSHIFT_HELIO']); factor=1+z
  L=float(np.max(t-factor*b));U=float(np.min(t-factor*a))
  center=float(h['PEAKMJD']) # predeclared engineering input, not estimated peak.
  assert sorted(zip(t,bands))==cad
  assert len(t)==117 and L<U and kz.min()<=z<=kz.max()
  inward=[float(np.nextafter(L,np.inf)),float(np.nextafter(U,-np.inf))]
  # Source FCN subtraction order is (MJD-MJDOFF)-offsetPeak; take actual MJDOFF
  # from structural CSP_OBJECTIVE fields, never read Q or fitted D/peak.
  off=None
  for line in (P/'fits-resume/joint12/native.log').read_text().splitlines():
   x=line.split()
   if x and x[0]=='CSP_OBJECTIVE:' and x[1]==str(h['SNID']).strip():
    off=float(x[-1]);break
  assert off is not None
  test=[]
  for endpoint_kind,peak in [("raw_lower",L),("raw_upper",U),("inward_lower",inward[0]),("inward_upper",inward[1])]:
   native=(t-off-(peak-off))/factor
   mask=np.float32(np.float32(t-peak)/np.float32(np.float32(1)+np.float32(z)))
   test.append({'endpoint_kind':endpoint_kind,'peak':peak,'R8_phase_min':float(native.min()),'R8_phase_max':float(native.max()),'R4_mask_phase_min':float(mask.min()),'R4_mask_phase_max':float(mask.max()),'all_R8_supported':bool(np.all((native>=a)&(native<=b))),'all_R4_supported':bool(np.all((mask>=a)&(mask<=b)))})
  grid=[center+start+d for start in [-2,0,2] for d in [-4,-2,0,2,4]]
  ledger.append({'CID':str(h['SNID']).strip(),'n':len(t),'band_counts':dict(Counter(bands)),'zHEL':z,'MJDmin':float(t.min()),'MJDmax':float(t.max()),'MJDOFF':off,'exact_interval':[L,U],'inward_R8_endpoints':inward,'width':U-L,'extremal_epochs':{'lower_bound':{'MJD':float(t[t.argmax()]),'band':bands[t.argmax()]},'upper_bound':{'MJD':float(t[t.argmin()]),'band':bands[t.argmin()]}},'relative_declared_center':[L-center,U-center],'declared_initial_grid_supported':all(L<x<U for x in grid),'endpoint_arithmetic':test})
support=[];first=[]
for line in (P/'fits-resume/joint12/native.log').read_text().splitlines():
 x=line.split()
 if x and x[0]=='PROSP_SUPPORT':support.append({'CID':x[1],'mean_calls':int(x[2]),'unsupported_calls':int(x[3]),'nonfinite_calls':int(x[4]),'phase_min':float(x[5]),'phase_max':float(x[6])})
 if x and x[0]=='PROSP_SUPPORT_FIRST_BAD':first.append({'CID':x[1],'iteration':int(x[2]),'filter_index':int(x[3]),'epoch_index':int(x[4]),'Trest':float(x[5]),'Tobs':float(x[6]),'z':float(x[7])})
paths=[P/x for x in ['protocol-v5.json','inputs-v2/cadence.simlib','inputs-v2/snoopy.B18/SNooPy_B18.fits','inputs-v2/kcor.fits','readme-adapted/ledger/PTE/PTE_HEAD.FITS','readme-adapted/ledger/PTE/PTE_PHOT.FITS','fits-resume/joint12/native.log','fits-resume/joint12/fit.nml','noisy-failure-checker-resume.json','fit-support/build/src/snlc_fit.car','fit-support/build/src/snana.car','fit-support/build/src/genmag_snoopy.c','fit-support/schema.json']]
r={'scope':'Metadata/source only. No new fits, source modification, builds, flux scores or distance outcomes.','model_phase':[float(model_phase.min()),float(model_phase.max())],'KCOR_phase':[float(kphase.min()),float(kphase.max())],'KCOR_z':[float(kz.min()),float(kz.max())],'common_phase':[float(a),float(b)],'same117_cadence_for_all8':True,'cohort':ledger,'original_failed_branch_support':support,'first_unsupported_call':first,'inputs':[{'path':str(p.relative_to(ROOT)),'resolved_path':str(p.resolve().relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)} for p in paths]}
(OUT/'metadata-result.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps({'n':len(ledger),'interval':ledger[0]['exact_interval'],'width':ledger[0]['width'],'endpoint_arithmetic':ledger[0]['endpoint_arithmetic'],'original_unsupported_counts':{x['CID']:x['unsupported_calls'] for x in support}},indent=2))
