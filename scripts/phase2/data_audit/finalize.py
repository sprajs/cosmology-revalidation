#!/usr/bin/env python3
"""Finalize the completed row audit, with schema-aware release comparison and hashes."""
from pathlib import Path
import datetime,json,hashlib
import h5py,numpy as np,pandas as pd
import audit as a

def main():
 a.summaries=json.loads((a.OUT/'progress.json').read_text())
 if len(a.summaries)!=7: raise RuntimeError('All seven dataset audits must complete before finalization')
 mt,md=a.metadat();cross={};npz={}
 for name,d in mt.items():
  err=d[['x0ERR','x1ERR','cERR']].to_numpy(float);cov=np.zeros((len(d),3,3));cov[:,range(3),range(3)]=err**2
  for i,j,col in [(0,1,'COV_x1_x0'),(0,2,'COV_c_x0'),(1,2,'COV_x1_c')]:cov[:,i,j]=cov[:,j,i]=d[col].to_numpy(float)
  jac=np.ones((len(d),3));jac[:,0]=-2.5/(np.log(10)*d.x0)
  cm=cov*jac[:,:,None]*jac[:,None,:]
  corr=cov/(err[:,:,None]*err[:,None,:]);valid=np.linalg.eigvalsh(corr)[:,0]>=-1e-8
  for key,val in {'cid':d.CID.to_numpy(str),'mean_mB_x1_c':d[['mB','x1','c']].to_numpy(float),'cov_x0_x1_c':cov,'cov_mB_x1_c':cm,'naive_covariance_psd':valid}.items():npz[name+'_'+key]=val
 np.savez_compressed(a.OUT/'fitted_parameter_covariances.npz',**npz)
 with h5py.File(a.OUT/'photometry_audit.h5','r') as h:
  for s in a.summaries:
   g=h[s['dataset']]
   for key in ['head_source','phot_source']:a.register(a.ROOT/g.attrs[key],'SNANA_FITS')
  for typ in ['DES','Foundation','LOWZ']:
   g=h['original/'+typ];k=h['current/'+typ];cross[typ]={t+'_bytes_identical':a.sha(a.ROOT/g.attrs[t.lower()+'_source'])==a.sha(a.ROOT/k.attrs[t.lower()+'_source']) for t in ['HEAD','PHOT']}
   ah,bh=g['head'],k['head'];diff={};cross[typ]['head_columns_removed']=sorted(set(ah)-set(bh));cross[typ]['head_columns_added']=sorted(set(bh)-set(ah))
   for c in sorted(set(ah)&set(bh)):
    x,y=ah[c][:],bh[c][:];neq=x!=y
    if x.dtype.kind=='f':neq &= ~(np.isnan(x)&np.isnan(y))
    if neq.any():diff[c]={'changed_rows':int(neq.sum())}
   cross[typ]['head_field_changes']=diff
 a.jsonout('summary.json',{'schema_version':'des-photometry-audit-v1','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'datasets':a.summaries,'fitted_metadata':md,'cross_release_identity':cross,'flag_bits':a.FLAG_BITS,'source_level':'calibrated flux; no detector pixels'})
 a.jsonout('contract.json',{'schema_version':'des-calibrated-photometry-audit-v1','sources':'Original FITS remain immutable source. HDF5 is column-preserving conversion, not new measurements. Hard-linked byte-identical PHOT tables are not independent data.','groups':['original/DES','original/Foundation','original/LOWZ','current/DES','current/Foundation','current/LOWZ','diffimg/DES'],'column_arrays':'<group>/head/<original-column>, <group>/phot/<original-column>','indices':'FITS PTROBS_MIN/MAX 1-based inclusive: python slice [MIN-1:MAX]. HDF5 head_row_zero_based maps every phot row to zero-based HEAD owner; -1 is unowned. External phot_row and head_row are 1-based.','masks':'basic_eligible does not remove negative flux or low SNR. audit_flags retains all data, including delimiters, and describes failures without silently excluding. Duplicate MJD/band flags require exposure-level inspection.','flags':a.FLAG_BITS,'SNR':'signed FLUXCAL/FLUXCALERR','flux_scale':'mag=27.5-2.5log10(FLUXCAL) for positive flux; use signed flux likelihood for all values.','depth':'27.5-2.5log10(5*FLUXCALERR) is an error-derived depth proxy; no independent completeness calibration.','quality_units':'SNANA CUTWIN_PSF uses PSF_SIG1*PIXSIZE*2.355 in FWHM arcsec; ZPNPE=ZEROPT+2.5log10(GAIN), with gain<.01 replaced by .001.','partial_cuts':'Only DES proxies: .05<=zHEL<=1.2,0<=zHELERR<=.01,0<=MWEBV<=.25,>=5 broad-phase epochs,>=2 bands with SNRmax>5, some phase<=5 and>=5. Published PKMJD for metadata matches, header PEAKMJD otherwise. NOT complete published selection; low-z proxy columns are inapplicable.','covariance_npz':'means ordered [mB,x1,c]; unrepaired covariance from public errors and crosscolumns transformed with dmB/dx0=-2.5/(ln10*x0). naive_covariance_psd flags mathematical compatibility only; MINOS and parabolic errors can be mixed. Invalid/rounded matrices must not be silently interpreted as precise Gaussian likelihoods.','fitted_metadata':'Original and current metadata are inferred fit/BBC products; they are not raw photometry and are not independent releases. No full all-transient SALT fit covariance was found in these release metadata files.'})
 for base in [a.ORIG,a.CUR]:
  for sub in ['0_DATA','3_CLASSIFICATION','7_PIPPIN_FILES']:
   for p in sorted((base/sub).rglob('*')):
    if p.is_file() and p.suffix.lower() in ['.md','.readme','.nml','.input','.yml','.ignore','.list']:a.register(p,'source_documentation_or_configuration')
 for p in [a.ROOT/'data/des-diffimg/DES-SN5YR_DIFFIMG.README',a.ROOT/'catalog/des-diffimg-zenodo.json',a.ROOT/'sources/repos/RickKessler__SNANA/src/snana.F90',a.ROOT/'sources/repos/RickKessler__SNANA/src/snlc_fit.F90'] :a.register(p,'source_documentation_or_provenance')
 for name in ['des-science__DES-SN5YR','des-science__DES-SN5YR@1.3','RickKessler__SNANA']:
  p=a.ROOT/'catalog/repositories'/(name+'.commit.json')
  if p.exists():a.register(p,'pinned_repository_commit')
 for p in [a.OUT/'preregistration.json',a.OUT/'audit-executed-v1.py',Path(__file__),Path(a.__file__)]:a.register(p,'analysis_or_preregistration')
 for p in [a.ORIG/'0_DATA/DES5YR_SALT3_LCFIT.LCPLOT.gz',a.ROOT/'phase2/official/results/snana_pilot12_plots.LCPLOT.TEXT',a.ROOT/'phase2/official/results/snana_mask32_checkplots.LCPLOT.TEXT']:
  if p.exists():a.register(p,'supplementary_published_or_refit_LCPLOT')
 a.jsonout('inputs_manifest.json',{'files':list(a.INPUTS.values())})
 # Hashes are finalized separately after validation/reporting/log completion.
 print('Finalized:',len(a.summaries),'datasets',flush=True)
if __name__=='__main__':main()
