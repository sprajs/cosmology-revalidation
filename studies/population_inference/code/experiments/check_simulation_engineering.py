"""Independent pre-score engineering review; no simulated residual scores."""
from pathlib import Path
import sys,hashlib,json
from collections import Counter
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[3];sys.path.insert(0,str(ROOT))
from scripts.salt_dust_audit import flux_response as fr
from scripts.salt_dust_audit.snana_extinction import SnanaExtinction
P=Path(__file__).resolve().parent;D=ROOT/'runs/research_2026_09_26/simulation_residual_control/engineering8'
ext=SnanaExtinction(ROOT/'runs/salt_dust_audit/snana_extinction/libsnana_extinction.so')
class OldLaw(fr.VariableF99):
 def propagate(self,wave,flux,phase=None):
  ebv,rv=self._parameters
  return flux*10**(-.4*ext(wave,rv,ebv,option=-99))
model,bands,paths,zp=fr.build_model();old=fr.sncosmo.Model(source=model.source,effects=[OldLaw(),fr.VariableF99()],effect_names=['mw','host'],effect_frames=['obs','rest'])
offs={str(r['Filter Name'])[-1]:float(r['Primary Mag']) for r in zp};rows=[];inputs=[]
for arm in ['P21','G10']:
 base=ROOT/f'phase2/literature/simulations/outputs/PH2_pilot02_{arm}'
 h=fits.getdata(base/f'PH2_pilot02_{arm}_HEAD.FITS',1);ph=fits.getdata(base/f'PH2_pilot02_{arm}_PHOT.FITS',1)
 heads={int(r['SNID']):r for r in h};ids=list(map(int,(P/f'simulation_design/{arm}-engineering8-cids.txt').read_text().split()))
 for law in ['approx_minus99','exact_99']:
  d=D/arm/law;setup=json.loads((d/'prepared.json').read_text());assert set(setup['selected_ids'])==set(ids)
  nml=(d/'fit.nml').read_text();assert 'LFIXPAR_ALL' not in nml and 'OPT_SNCID_LIST = 2' in nml
  choice=-99 if law=='approx_minus99' else 99;assert f'OPT_MWCOLORLAW = {choice}' in nml
  manifest=json.loads((d/'objectives/manifest.json').read_text());assert set(map(int,manifest['checks']))==set(ids)
  m=old if choice==-99 else model
  for cid in ids:
   f=d/f'objectives/objective_{cid}.npz';inputs.append(f);a=np.load(f);hh=heads[cid];raw=ph[int(hh['PTROBS_MIN'])-1:int(hh['PTROBS_MAX'])]
   raw_tuples=Counter(zip(raw['MJD'].astype(float),np.char.strip(raw['BAND'].astype(str)),raw['FLUXCAL'].astype(float),raw['FLUXCALERR'].astype(float)))
   observed=Counter(zip(a['MJD'],a['band'],a['data_flux'],a['data_fluxerr']));assert not (observed-raw_tuples)
   C=a['frozen_flux_covariance'];ci=a['inverse_frozen_flux_covariance'];L=np.linalg.cholesky(C)
   res=a['data_flux']-a['model_flux'];closure=float(solve_triangular(L,res,lower=True)@solve_triangular(L,res,lower=True)+a['prior_chi2']-a['chi2']);assert abs(closure)<1e-7
   x0,x1,c,t0=a['parameters_x0_x1_c_t0'];z=float(a['zHEL'][0]);ebv=float(a['MWEBV'][0]);t=a['MJD'];labels=a['band'];bp=np.array([bands[b] for b in labels],dtype=object);conv=np.array([10**(-.4*(.27+offs[b])) for b in labels])
   def flux(th):
    m.set(z=z,t0=t0+th[3],x0=x0*np.exp(-fr.K*th[0]),x1=x1+th[1],c=c+th[2],mwebv=ebv,mwrv=3.1,hostebv=0,hostrv=3.1)
    return conv*m.bandflux(bp,t,zp=27.5,zpsys='ab')
   zero=np.zeros(4);native=flux(zero);steps=np.array([1e-4,1e-3,1e-4,.01]);jac=[];half=[]
   for j,st in enumerate(steps):
    v=np.eye(4)[j]*st;jac.append((flux(v)-flux(-v))/(2*st));half.append((flux(v/2)-flux(-v/2))/st)
   J=np.column_stack(jac);Jh=np.column_stack(half);derivative=float(np.linalg.norm(J-Jh)/np.linalg.norm(Jh));J[:,0]=-fr.K*a['model_flux']
   Jw=solve_triangular(L,J,lower=True);u,s,vt=np.linalg.svd(Jw,full_matrices=False);assert (s>s[0]*1e-10).sum()==4
   diff=solve_triangular(L,native-a['model_flux'],lower=True);perp=diff-u@(u.T@diff)
   rows.append(dict(arm=arm,law=law,CID=cid,epochs=len(t),objective_error=closure,c_inverse_identity_max=float(np.max(abs(C@ci-np.eye(len(t))))),J_rank=4,J_singular_ratio=float(s[-1]/s[0]),halfstep_relative_difference=derivative,mean_difference_projected_norm=float(np.linalg.norm(perp)),mean_difference_quoted_error_max=float(np.max(abs((native-a['model_flux'])/a['data_fluxerr'])))))
f=pd.DataFrame(rows);f.to_csv(P/'simulation-engineering-independent.csv',index=False)
result={'status':'PASS: source/free-fit/raw-epoch/native-law/unit/covariance/tangent engineering gate; no observer residual scores computed','object_law_checks':len(rows),'epochs':int(f.epochs.sum()),'max_objective_error':float(f.objective_error.abs().max()),'max_C_inverse_error':float(f.c_inverse_identity_max.max()),'min_tangent_singular_ratio':float(f.J_singular_ratio.min()),'max_derivative_halfstep_difference':float(f.halfstep_relative_difference.max()),'max_native_mean_difference_projected_norm':float(f.mean_difference_projected_norm.max()),'max_native_mean_difference_quoted_sigma':float(f.mean_difference_quoted_error_max.max()),'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs+[Path(__file__)]}}
(P/'simulation-engineering-independent.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='sha256'},indent=2))
