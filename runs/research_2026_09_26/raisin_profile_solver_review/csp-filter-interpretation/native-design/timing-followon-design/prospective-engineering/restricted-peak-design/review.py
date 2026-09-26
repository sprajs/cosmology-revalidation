from pathlib import Path
import json,hashlib,importlib.util,collections
import numpy as np
from astropy.io import fits
O=Path(__file__).resolve().parent;P=O.parent;R=Path('/home/szymon/Documents/ChatGPT/supernova')
spec=importlib.util.spec_from_file_location('records',P/'noiseless-checker-recovery/full_null_check.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
_,blocks=m.records(P/'fits-resume/joint12/native.log')
source=P/'generation-v5/ledger/output/PTE'
with fits.open(source/'PTE_HEAD.FITS') as f:heads=f[1].data.copy()
with fits.open(source/'PTE_PHOT.FITS') as f:phot=f[1].data.copy()
head=heads[0];raw=phot[int(head['PTROBS_MIN'])-1:int(head['PTROBS_MAX'])];assert len(raw)==117
z=float(head['REDSHIFT_HELIO']);mjd=np.array(raw['MJD'],float);tlo=max(-20.,-20.);thi=min(70.,85.)
# Exact all-row intersection, with inward floating-point guard validated in native operation order.
lo=float(np.max(mjd-(1+z)*thi));hi=float(np.min(mjd-(1+z)*tlo));rawlo,rawhi=lo,hi
while np.max((mjd-lo)/(1+z))>thi:lo=float(np.nextafter(lo,np.inf))
while np.min((mjd-hi)/(1+z))<tlo:hi=float(np.nextafter(hi,-np.inf))
assert lo<hi
for peak in [lo,hi]:assert np.min((mjd-peak)/(1+z))>=tlo and np.max((mjd-peak)/(1+z))<=thi
truth=float(head['PEAKMJD']);starts=[truth+x for x in [-2,0,2]];initialgrid=[x+y for x in starts for y in [-4,-2,0,2,4]];assert min(initialgrid)>=lo and max(initialgrid)<=hi
support={}
for line in (P/'fits-resume/joint12/native.log').read_text().splitlines():
 t=line.split()
 if t and t[0]=='PROSP_SUPPORT':assert t[1] not in support;support[t[1]]=t
final=[]
for cid in [str(x) for x in range(1,9)]:
 bb=[b for b in blocks if b['cid']==cid];b=bb[-1];o=[float(x) for x in b['objective'][5:]];a=np.array([[float(x) for x in row[6:]] for row in b['rows']]);sv=support[cid]
 assert b['n']==117 and o[4]==1 and o[5]==0
 final.append({'CID':cid,'all_mean_calls':int(sv[2]),'unsupported_calls':int(sv[3]),'all_call_phase_bounds':[float(sv[5]),float(sv[6])],'final_callback_phase_bounds':[float(a[:,1].min()),float(a[:,1].max())],'final_peak_in_allrow_table_domain':lo<=o[6]<=hi,'final_peak_within_original4day_gate':abs(o[6]-truth)<=4,'final_rows117':True,'final_Q_arithmetic_error':b['qerr']})
res={'status':'Original unrestricted eight-draw branch FAILS; separate restricted estimator proposal only','native_execution':False,'intersection_basis':'Exposure MJD metadata, fixed zHEL and tabulated model/KCOR domains; no fitted-peak or flux-dependent bound','frame':'Trest=(MJD-peak_absolute)/(1+zHEL), MJDOFF0 in this frozen input','zHEL':z,'model_phase_table':[-20,70],'KCOR_phase_table':[-20,85],'combined_phase_table':[tlo,thi],'all_physical_epochs':len(mjd),'earliest_MJD':float(mjd.min()),'latest_MJD':float(mjd.max()),'raw_peak_bounds_absolute':[rawlo,rawhi],'inward_safe_peak_bounds_absolute':[lo,hi],'inward_rounding_days':[lo-rawlo,rawhi-hi],'peak_bounds_relative_to_frozen_truth_only_for_reporting':[lo-truth,hi-truth],'planned_initial_grid_and_all_three_starts_inside':True,'initial_grid_peak_range':[min(initialgrid),max(initialgrid)],'original_failed_draws':final,'not_physical_validation':'These are tabulation bounds. They do not validate training extrapolation or provide a population timing prior. Valid final states do not undo invalid intermediate means.'}
(O/'result.json').write_text(json.dumps(res,indent=2)+'\n')
paths=[Path(__file__),P/'noiseless-checker-recovery/full_null_check.py',source/'PTE_HEAD.FITS',source/'PTE_PHOT.FITS',P/'fits-resume/joint12/native.log',O/'result.json'];(O/'manifest.json').write_text(json.dumps({'files':{str(q):hashlib.sha256(q.read_bytes()).hexdigest() for q in paths}},indent=2)+'\n')
print(json.dumps(res,indent=2))
