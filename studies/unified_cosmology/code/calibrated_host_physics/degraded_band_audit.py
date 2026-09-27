"""One declared all-host covariance-preserving spectral degradation audit.

No age fits, new FSPS calls or edits of the original response audit.
"""
import os
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,hashlib,sys,importlib.metadata
import numpy as np
from scipy import sparse
from scipy.linalg import solve_triangular
import extinction
ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).resolve().parent
OUT=ROOT/'.work/unified-cosmology/full-spectrum-design/degraded';OUT.mkdir(parents=True,exist_ok=True)
INPUT=ROOT/'.work/unified-cosmology/calibrated-hosts/spectra'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
parent=ROOT/'studies/unified_cosmology/results/calibrated_host_physics/full-spectrum-audit.json'
parent_result=json.loads(parent.read_text());design_path=HERE/'degraded-band-audit-design.json';design=json.loads(design_path.read_text())
assert sha(parent)==design['parent_result_sha256']
assert sha(HERE/'full-spectrum-audit-design.json')==design['parent_design_sha256']
for kind in ['input_sha256','source_sha256','output_sha256','origin_sha256']:
 for name,digest in parent_result[kind].items():assert sha(ROOT/name)==digest,name
base_design=json.loads((HERE/'full-spectrum-audit-design.json').read_text())
rows=json.loads((ROOT/'.work/unified-cosmology/full-spectrum-design/input-inventory.json').read_text())['rows']
libpath=ROOT/'.work/unified-cosmology/calibrated-host-physics/ssp-library.npz';lib=np.load(libpath)
lw=lib['wavelength_A'];stellar=lib['spectra_Lsun_per_A_per_formed_Msun'].reshape(-1,len(lw));all_lines=lib['emission_line_wavelength_A']
windows=np.asarray(base_design['windows_vacuum_rest_A'],float);reference=base_design['reference_window_index']
dust=[(0.,-.7)]+[(tau,power) for tau in [.5,1.,2.] for power in [-.7,-1.]]

def blur(wave,sigma_kms):
 sigma=sigma_kms/299792.458;step=float(np.median(np.diff(wave)));n=len(wave)
 radius=int(np.ceil(5*sigma*wave[-1]/step))+1;rr=[];cc=[];vv=[]
 for off in range(-radius,radius+1):
  ii=np.arange(max(0,-off),min(n,n-off));jj=ii+off;x=np.log(wave[ii]/wave[jj])/sigma;keep=abs(x)<=5
  ii,jj,x=ii[keep],jj[keep],x[keep];rr.extend(ii);cc.extend(jj);vv.extend(step/wave[ii]/(np.sqrt(2*np.pi)*sigma)*np.exp(-.5*x*x))
 return sparse.csr_matrix((vv,(rr,cc)),shape=(n,n))

def audit(record,extra):
 meta=json.loads((INPUT/(record['targetid']+'.json')).read_text());a=np.load(ROOT/meta['array_file']);z=meta['desi_z'];parts=[]
 for cam,lo,hi in [('B',-np.inf,5780.),('R',5780.,7570.),('Z',7570.,np.inf)]:
  wave=a[cam+'_WAVELENGTH'].astype(float);rest=wave/(1+z);f=a[cam+'_FLUX'].astype(float);iv=a[cam+'_IVAR'].astype(float)
  valid=(a[cam+'_MASK']==0)&np.isfinite(f)&np.isfinite(iv)&(iv>0)
  for line in all_lines[(all_lines>rest[0]*.997)&(all_lines<rest[-1]*1.003)]:valid &= abs(np.log(rest/line))*299792.458>400
  G=blur(wave,base_design['kernel_sigma_km_s']);H=blur(wave,extra) if extra else sparse.eye(len(wave),format='csr');K=H@G
  support=K.copy();support.data[:]=1.;usable=np.asarray(support@(~valid).astype(float)).ravel()==0
  edge=5*(base_design['kernel_sigma_km_s']+extra)/299792.458
  usable &= (wave>wave[0]*np.exp(edge))&(wave<wave[-1]*np.exp(-edge))
  cuts=(wave>=lo)&(wave<hi)
  geom=np.array([np.maximum(0,np.minimum(rest+.4/(1+z),h)-np.maximum(rest-.4/(1+z),l))*cuts for l,h in windows])
  d=a[cam+'_RESOLUTION'].astype(float);k=d.shape[0]//2;R=sparse.dia_matrix((d,np.arange(k,-k-1,-1)),shape=(len(wave),len(wave))).tocsr()
  spec=np.array([np.interp(rest,lw,s) for s in stellar]).T
  fore=10**(-.4*extinction.fitzpatrick99(wave,3.1*meta['fibermap']['EBV']*.86,3.1))/(1+z)
  parts.append(dict(wave=wave,rest=rest,flux=np.where(valid,f,0),var=np.divide(1,iv,out=np.zeros_like(iv),where=valid),G=G,H=H,K=K,R=R,f=spec,fore=fore,widths=geom*usable,geom=geom))
 retained=sum(p['widths'].sum(axis=1) for p in parts);coverage=sum(p['geom'].sum(axis=1) for p in parts)/np.diff(windows,axis=1).ravel()
 indices=np.flatnonzero((retained>0)&(coverage>=.95));ref=int(np.flatnonzero(indices==reference)[0]);n=len(indices)
 y=np.zeros(n);C=np.zeros((n,n));independent=np.zeros_like(C);left=np.zeros((n,len(stellar)*len(dust)));right=left.copy();stress=[]
 for p in parts:
  W=p['widths'][indices]/retained[indices,None];WK=np.asarray(p['K'].T@W.T).T;WH=np.asarray(p['H'].T@W.T).T
  WKR=np.asarray(p['R'].T@WK.T).T;WHR=np.asarray(p['R'].T@WH.T).T
  y+=WK@p['flux'];C+=(WK*p['var'])@WK.T;independent+=np.einsum('ik,k,jk->ij',WK,p['var'],WK,optimize=False)
  for j,(tau,power) in enumerate(dust):
   D=p['fore']*np.exp(-tau*(p['rest']/5500)**power);A=WKR*D;B=np.asarray(p['G'].T@(WHR*D).T).T;sl=slice(j*len(stellar),(j+1)*len(stellar))
   left[:,sl]+=A@p['f'];right[:,sl]+=B@p['f'];stress.append((A-B)/D[None,:]/float(np.median(np.diff(p['wave']))))
 L=np.linalg.cholesky(C);sd=np.sqrt(np.diag(C));assert y[ref]>0 and np.all(right[ref]>0)
 errors=(left-right)*(y[ref]/right[ref])[None,:];white=solve_triangular(L,errors,lower=True)
 single=float(np.max(abs(errors)/sd[:,None]));joint=float(np.max(np.linalg.norm(white,axis=0)))
 sharp=max(float(np.max(np.linalg.norm(solve_triangular(L,x*y[ref],lower=True),axis=0))) for x in stress)
 corr=C/np.outer(sd,sd);np.fill_diagonal(corr,0);check=float(np.max(abs(C-independent))/np.max(abs(C)));assert check<1e-11
 features={}
 for name,ids in [('Hdelta_masked_contrast',[2,3,4]),('Hgamma_G4300_masked_contrast',[5,6,7])]:
  if not all(i in indices for i in ids):features[name]=None;continue
  v=np.zeros(n);v[[int(np.flatnonzero(indices==i)[0]) for i in ids]]=[-.5,1,-.5]
  noise=float(np.sqrt(v@C@v)/y[ref]);shapes=(v@right[:,:len(stellar)])/right[ref,:len(stellar)]
  features[name]={'formal_noise_relative_to_reference':noise,'resolved_unreddened_SSP_contrast_min_max':[float(shapes.min()),float(shapes.max())],'resolved_SSP_contrast_range_over_noise':float(np.ptp(shapes)/noise)}
 return dict(targetid=record['targetid'],extra_sigma_kms=extra,windows_retained=indices.tolist(),retained_rest_width_A=retained.tolist(),signed_broad_flux=y.tolist(),full_broad_covariance=C.tolist(),max_band_error_sigma=single,max_joint_error=joint,max_1A_line_joint_error=sharp,max_abs_cross_correlation=float(np.max(abs(corr))),independent_covariance_relative_error=check,finite_family_passed=single<=.1 and joint<=.1,line_stress_passed=sharp<=.1,all_operator_gates_passed=single<=.1 and joint<=.1 and sharp<=.1,features=features)

records=[]
for row in rows:
 old=audit(row,0.);new=audit(row,design['extra_kernel_sigma_km_s']);records.append(dict(targetid=row['targetid'],original=old,degraded=new))
 print(row['targetid'],'new',new['max_joint_error'],'line',new['max_1A_line_joint_error'],flush=True)
ledger=OUT/'records.json';ledger.write_text(json.dumps(records,indent=2,allow_nan=False)+'\n')
for key in ['input_sha256','source_sha256','output_sha256','origin_sha256']:
 for name,digest in parent_result[key].items():assert sha(ROOT/name)==digest,name
features={}
for name in ['Hdelta_masked_contrast','Hgamma_G4300_masked_contrast']:
 valid=[x for x in records if x['original']['features'][name] and x['degraded']['features'][name]]
 features[name]={'hosts':len(valid),'formal_noise_degraded_over_original_min_median_max':np.percentile([x['degraded']['features'][name]['formal_noise_relative_to_reference']/x['original']['features'][name]['formal_noise_relative_to_reference'] for x in valid],[0,50,100]).tolist(),'template_contrast_range_over_noise_degraded_over_original_min_median_max':np.percentile([x['degraded']['features'][name]['resolved_SSP_contrast_range_over_noise']/x['original']['features'][name]['resolved_SSP_contrast_range_over_noise'] for x in valid],[0,50,100]).tolist()}
previous={x['targetid']:x for x in json.loads((ROOT/'.work/unified-cosmology/full-spectrum-design/response-remedy-ledger.json').read_text())}
baseline_delta=max(max(abs(x['original']['max_joint_error']-previous[x['targetid']]['maximum_column_joint_whitened_error']),abs(x['original']['max_1A_line_joint_error']-previous[x['targetid']]['unit_observed_EW_line_max_joint_error'])) for x in records)
assert baseline_delta<1e-12
validation=dict(original_operator_reproduction_max_absolute_error=baseline_delta,independent_covariance_max_relative_error=max(x['degraded']['independent_covariance_relative_error'] for x in records),minimum_full_covariance_eigenvalue=min(float(np.linalg.eigvalsh(x['degraded']['full_broad_covariance'])[0]) for x in records),retained_summed_band_width_fraction_min_median_max=np.percentile([np.sum(x['degraded']['retained_rest_width_A'])/np.sum(x['original']['retained_rest_width_A']) for x in records],[0,50,100]).tolist(),hosts_losing_at_least_one_window=sum(len(x['degraded']['windows_retained'])<len(x['original']['windows_retained']) for x in records))
summary=dict(status='global_repair_failed_not_adopted_no_age_fit',validation=validation,source_sha256={str(Path(__file__).relative_to(ROOT)):sha(__file__)},design_sha256=sha(design_path),parent_result_sha256=sha(parent),input_identity='Exact input and origin hashes in bound parent result; verified before and after calculation.',output_sha256={str(ledger.relative_to(ROOT)):sha(ledger)},versions={k:importlib.metadata.version(k) for k in ['numpy','scipy','extinction']},hosts=len(records),extra_sigma_kms=design['extra_kernel_sigma_km_s'],original_all_operator_pass=sum(x['original']['all_operator_gates_passed'] for x in records),degraded_all_operator_pass=sum(x['degraded']['all_operator_gates_passed'] for x in records),degraded_finite_family_pass=sum(x['degraded']['finite_family_passed'] for x in records),degraded_line_stress_pass=sum(x['degraded']['line_stress_passed'] for x in records),degraded_max_band_error_sigma=max(x['degraded']['max_band_error_sigma'] for x in records),degraded_max_joint_error=max(x['degraded']['max_joint_error'] for x in records),degraded_max_1A_line_joint_error=max(x['degraded']['max_1A_line_joint_error'] for x in records),degraded_failures=[{'targetid':x['targetid'],'max_1A_line_joint_error':x['degraded']['max_1A_line_joint_error']} for x in records if not x['degraded']['all_operator_gates_passed']],absorption_contrast_sensitivity=features,interpretation='Finite-family and fixed1A unresolved-line numerical ordering checks only. No proof of true library kernel, abundance/calibration/coadd adequacy, stellar-age identification or arbitrary sharp-feature coverage. Additional smoothing is a measurement transformation, not a physical dispersion prior.')
public=ROOT/'studies/unified_cosmology/results/calibrated_host_physics/degraded-band-audit.json';public.write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n');print(json.dumps(summary,indent=2),flush=True)
