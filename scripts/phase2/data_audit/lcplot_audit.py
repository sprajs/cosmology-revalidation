#!/usr/bin/env python3
"""Audit released accepted-point masks against calibrated FITS and refit pilot."""
from pathlib import Path
import pandas as pd,numpy as np,h5py,json,hashlib,sys
R=Path(__file__).resolve().parents[3];O=R/'phase2/data_audit'
p=R/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES5YR_SALT3_LCFIT.LCPLOT.gz'
pub=pd.read_csv(p,sep=r'\s+',comment='#',dtype={'SNID':str});pub=pub[pub.DATA_MODEL!=0].copy()
suffix='_mask32' if '--mask32' in sys.argv else ''
newpath=R/('phase2/official/results/snana_mask32_checkplots.LCPLOT.TEXT' if suffix else 'phase2/official/results/snana_pilot12_plots.LCPLOT.TEXT')
new=pd.read_csv(newpath,sep=r'\s+',names=['SNID','MJD','PHASE','FLUXCAL','FLUXCALERR','DATA_MODEL','BAND','CHI2','IGNOREME'],dtype={'SNID':str});new=new[new.DATA_MODEL!=0]
dup=pd.read_csv(O/'duplicate_epochs.csv.gz',dtype={'snid':str});dup=dup[dup.dataset=='original/DES'].copy()
# LCPLOT serializes SNANA single-precision MJD to F*.3; raw FITS MJD has greater precision.
def lc_mjd(x):return np.array([float(f'{v:.3f}') for v in np.asarray(x,dtype=np.float32)])
dup['lc_mjd']=lc_mjd(dup.MJD);q=dup.groupby(['snid','MJD','band']).agg(raw_rows=('phot_row','size'),IMGNUM=('IMGNUM','first'),lc_mjd=('lc_mjd','first')).reset_index();q['published_rows']=0;q['published_accepted_rows']=0;q['published_excluded_rows']=0
lookup=pub.groupby(['SNID','MJD','BAND']).DATA_MODEL.apply(list).to_dict()
for i,r in q.iterrows():
 vals=lookup.get((r.snid,r.lc_mjd,r.band),[]);q.loc[i,'published_rows']=len(vals);q.loc[i,'published_accepted_rows']=vals.count(1);q.loc[i,'published_excluded_rows']=vals.count(-1)
q.to_csv(O/'duplicate_exposure_published_acceptance.csv.gz',index=False,compression={'method':'gzip','mtime':0})
results={}
with h5py.File(O/'photometry_audit.h5') as h:
 g=h['original/DES'];ids=np.char.strip(g['head/SNID'][:].astype('U'))
 for cid in ['1442085','1896213','1289306','1256426']:
  j=np.flatnonzero(ids==cid)[0];a=int(g['head/PTROBS_MIN'][j])-1;b=int(g['head/PTROBS_MAX'][j]);cols=['MJD','BAND','FLUXCAL','FLUXCALERR','PHOTFLAG','PSF_SIG1','ZEROPT','GAIN','IMGNUM'];raw=pd.DataFrame({c:g['phot/'+c][a:b] for c in cols});raw.BAND=raw.BAND.str.decode('utf-8').str.strip();raw['lc_mjd']=lc_mjd(raw.MJD)
  p=pub[pub.SNID==cid];n=new[new.SNID==cid];rows=[];unmatched=[]
  for r in p.itertuples():
   match=raw[(raw.lc_mjd==r.MJD)&(raw.BAND==r.BAND)]
   if len(match)!=1: unmatched.append({'MJD':r.MJD,'band':r.BAND,'count':len(match)});continue
   k=match.index[0];v=match.iloc[0];nm=n[(abs(n.MJD-v.MJD)<.00251)&(n.BAND==r.BAND)]
   if len(nm)>1:nm=nm.iloc[np.argsort(abs(nm.FLUXCAL-v.FLUXCAL))[:1]]
   row={'cid':cid,'raw_phot_row':int(a+k+1),'raw_mjd':float(v.MJD),'published_mjd':r.MJD,'band':r.BAND,'raw_flux':float(v.FLUXCAL),'published_flux':r.FLUXCAL,'raw_err':float(v.FLUXCALERR),'published_err':r.FLUXCALERR,'published_used':r.DATA_MODEL,'published_chi2':r.CHI2,'published_phase_observer_days':r.PHASE,'raw_photflag':int(v.PHOTFLAG),'raw_psf_sigma_pixels':float(v.PSF_SIG1),'raw_zp':float(v.ZEROPT),'raw_gain':float(v.GAIN)}
   if len(nm)==1:
    vnew=nm.iloc[0];row.update(refit_flux=float(vnew.FLUXCAL),refit_err=float(vnew.FLUXCALERR),refit_used=int(vnew.DATA_MODEL),refit_chi2=float(vnew.CHI2),refit_phase_observer_days=float(vnew.PHASE))
   rows.append(row)
  d=pd.DataFrame(rows);d.to_csv(O/f'lcplot_raw_{cid}{suffix}.csv',index=False);results[cid]={'raw_rows':len(raw),'published_data_rows':len(p),'published_accepted':int((p.DATA_MODEL==1).sum()),'refit_pilot_accepted':int((n.DATA_MODEL==1).sum()) if len(n) else None,'unmatched_public_rows':unmatched,'max_abs_flux_relative_difference':float((abs(d.published_flux-d.raw_flux)/np.maximum(abs(d.raw_flux),1)).max()),'max_abs_err_relative_difference':float((abs(d.published_err-d.raw_err)/d.raw_err).max())}
  extra=[]
  for r in n[n.DATA_MODEL==1].itertuples():
   rm=raw[(abs(raw.MJD-r.MJD)<.00251)&(raw.BAND==r.BAND)]
   if len(rm)==1:
    v=rm.iloc[0];pm=p[(p.MJD==v.lc_mjd)&(p.BAND==r.BAND)]
    if len(pm)==0:extra.append({'raw_mjd':float(v.MJD),'band':r.BAND,'refit_used':int(r.DATA_MODEL),'refit_chi2':float(r.CHI2),'refit_flux':float(r.FLUXCAL),'refit_err':float(r.FLUXCALERR),'raw_photflag':int(v.PHOTFLAG),'raw_psf_sigma':float(v.PSF_SIG1)})
  results[cid]['refit_accepted_absent_from_public_lcplot']=extra
  if 'refit_used' in d:results[cid]['acceptance_disagreements']=d[(d.published_used!=d.refit_used)|d.refit_used.isna()][['raw_mjd','band','published_used','refit_used','published_chi2','refit_chi2','raw_flux','raw_err']].replace({np.nan:None}).to_dict('records')
summary={'note':'Public LCPLOT MJD requires float32 then printed3 decimals; matching unrounded double FITS to3-decimal public values with0.0005-day tolerance would falsely miss most rows. No flux corrections were reapplied.','targeted':results,'duplicate_raw_groups':len(q),'duplicate_groups_on_public_lcplot':int((q.published_rows>0).sum()),'duplicate_groups_with_accepted_rows':int((q.published_accepted_rows>0).sum()),'duplicate_groups_multiple_accepted_rows':int((q.published_accepted_rows>1).sum()),'raw_rows_in_groups_with_any_accepted':int(q.loc[q.published_accepted_rows>0,'raw_rows'].sum()),'public_accepted_rows_in_raw_duplicate_groups':int(q.published_accepted_rows.sum()),'public_duplicate_snid_mjd_band_rows_all_data':int(pub.duplicated(['SNID','MJD','BAND'],keep=False).sum()),'files':[{'path':str(p.relative_to(R)),'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest()} for p in [R/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES5YR_SALT3_LCFIT.LCPLOT.gz',newpath,Path(__file__)]]}
(O/('lcplot_audit'+suffix+'.json')).write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
