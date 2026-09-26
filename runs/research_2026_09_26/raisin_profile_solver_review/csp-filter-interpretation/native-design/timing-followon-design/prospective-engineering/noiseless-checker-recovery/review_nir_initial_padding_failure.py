from pathlib import Path
import json,hashlib,importlib.util,collections
import numpy as np
from astropy.io import fits
O=Path(__file__).resolve().parent;P=O.parent;W=P/'fits-resume/noiseless_NIR'
spec=importlib.util.spec_from_file_location('independent_records',O/'full_null_check.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
_,blocks=m.records(W/'native.log');assert len(blocks)==12 and [x['iter'] for x in blocks]==list(range(1,13))
d=P/'readme-adapted/noiseless/PTE'
with fits.open(d/'PTE_HEAD.FITS') as f:head=f[1].data.copy()
with fits.open(d/'PTE_PHOT.FITS') as f:phot=f[1].data.copy()
h=head[0];rows=phot[int(h['PTROBS_MIN'])-1:int(h['PTROBS_MAX'])];rows=rows[np.isin(rows['BAND'],['J','H'])];assert len(rows)==6
truthline=next(x.split() for x in (P/'generation-v5/noiseless/native.log').read_text().splitlines() if x.startswith('PROSP_EVENT '));truth=float(truthline[8]);peak=float(h['PEAKMJD'])
expected=collections.Counter((str(x['BAND']).strip(),float(x['MJD']),float(x['FLUXCAL']),float(x['FLUXCALERR'])) for x in rows)
allrows=[];Cs=[];Ws=[];objects=[]
for b in blocks:
 assert b['cid']=='1' and b['n']==6
 rr=b['rows'];a=np.array([[float(v) for v in x[6:]] for x in rr]);wv=np.array([[float(v) for v in x[4:]] for x in b['w']]);weight=wv if b['cov'] else np.diag(wv[:,0]);C=np.linalg.inv(weight)
 observed=collections.Counter((x[5],float(x[6]),float(x[10]),float(x[11])) for x in rr);assert observed==expected
 ob=[float(x) for x in b['objective'][5:]];assert ob[4]==1 and ob[5]==0 and ob[6]==peak
 assert a[:,1].min()>=-20 and a[:,1].max()<=70 and np.linalg.eigvalsh(C).min()>0
 Cs.append(C);Ws.append(weight);allrows.append(a);objects.append(ob)
support=[x.split() for x in (W/'native.log').read_text().splitlines() if x.startswith('PROSP_SUPPORT')];assert len(support)==1 and len(support[0])==11 and support[0][0]=='PROSP_SUPPORT';assert support[0][1]=='1' and int(support[0][2])>0 and int(support[0][3])==int(support[0][4])==0
lines=(W/'fit.FITRES.TEXT').read_text().splitlines();names=next(x.split()[1:] for x in lines if x.startswith('VARNAMES:'));vals=[dict(zip(names,x.split()[1:])) for x in lines if x.startswith('SN:')];assert len(vals)==1 and vals[0]['CID']=='1' and int(vals[0]['CUTFLAG_SNANA'])==3 and int(vals[0]['ERRFLAG_FIT'])==0
f=allrows[-1][:,2];y=allrows[-1][:,4];a=float(f@Ws[-1]@y/(f@Ws[-1]@f));assert a>0;dd=-2.5*np.log10(a)
L=np.linalg.cholesky(Cs[-2]);x=np.linalg.solve(L,Cs[-1]-Cs[-2]);x=np.linalg.solve(L,x.T).T;cnorm=float(np.linalg.norm(x,2))
D=objects[-1][3];pk=objects[-1][6];maxr=float(np.max(np.abs((y-f)/allrows[-1][:,5])))
assert abs(D-truth)<=.001 and abs(pk-peak)<=.01 and maxr<=.02 and abs(dd)<=.001 and cnorm<=.001
assert abs(objects[-1][3]-objects[-2][3])<=.001 and abs(objects[-1][6]-objects[-2][6])<=.01
result={'pass':True,'scope':'Independent pending noiseless NIR gate; no new native calls','callback_count':len(blocks),'rows_per_callback':6,'exact_native_FITS_MJD_flux_error_membership':True,'every_callback_peak_exactly_header':True,'truth_D':truth,'final_D':D,'D_error_mag':D-truth,'peak_error_days':pk-peak,'max_abs_standardized_mean_residual':maxr,'max_abs_Q_reconstruction':max(b['qerr'] for b in blocks),'final_fixedC_amplitude_delta_D_mag':float(dd),'last_two_whitened_C_norm':cnorm,'mean_call_count':int(support[0][2]),'mean_call_phase_bounds':[float(support[0][5]),float(support[0][6])],'CUTFLAG_SNANA':3,'ERRFLAG_FIT':0,'full_noiseless_interpretation':'With existing independently passed joint/control result this closes the declared engineering noiseless gate; it is not a physical validation of empirical model or population bias.'}
(O/'independent-noiseless-NIR.json').write_text(json.dumps(result,indent=2)+'\n')
paths=[Path(__file__),O/'full_null_check.py',W/'native.log',W/'fit.FITRES.TEXT',d/'PTE_HEAD.FITS',d/'PTE_PHOT.FITS',P/'generation-v5/noiseless/native.log',O/'independent-noiseless-NIR.json'];(O/'independent-noiseless-NIR-manifest.json').write_text(json.dumps({'files':{str(x):hashlib.sha256(x.read_bytes()).hexdigest() for x in paths}},indent=2)+'\n')
print(json.dumps(result,indent=2))
