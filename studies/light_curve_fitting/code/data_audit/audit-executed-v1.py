#!/usr/bin/env python3
"""All-row audit of public calibrated DES photometry; no cosmology fitting."""
from pathlib import Path
import csv, datetime, gzip, hashlib, json, shutil, sys
import numpy as np
import pandas as pd
import h5py
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'phase2/data_audit'; OUT.mkdir(parents=True,exist_ok=True)
CACHE=OUT/'cache'; CACHE.mkdir(exist_ok=True)
ORIG=ROOT/'sources/repos/des-science__DES-SN5YR@1.3'
CUR=ROOT/'sources/repos/des-science__DES-SN5YR'
SENTINELS=[-9,-99,-999,-9999,-777,99,999]
FLAG_BITS={'delimiter':1,'unowned':2,'invalid_mjd':4,'nonfinite_flux':8,'invalid_fluxerr':16,'invalid_band':32,'duplicate_object_mjd_band':64,'documented_diffimg_error_flag':128,'des_psf_cut_failure':256,'des_electron_zeropoint_cut_failure':512}
INPUTS={}; profiles=[]; summaries=[]; all_objects=[]; duplicates=[]; joins={}

def sha(p):
 return hashlib.file_digest(p.open('rb'),'sha256').hexdigest()
def register(p, role):
 p=Path(p); key=str(p.relative_to(ROOT));
 if key not in INPUTS: INPUTS[key]={'path':key,'bytes':p.stat().st_size,'sha256':sha(p),'role':role}
 return INPUTS[key]
def jsonout(name,d): (OUT/name).write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def native(a):
 a=np.asarray(a)
 if a.dtype.kind=='U': return np.char.encode(a,'utf-8')
 return a.astype(a.dtype.newbyteorder('='),copy=False)
def profile(a, dataset, table, col):
 a=np.asarray(a); r={'dataset':dataset,'table':table,'column':col,'n':len(a),'dtype':str(a.dtype)}
 if a.dtype.kind in 'fiu':
  f=np.isfinite(a); r.update(nonfinite=int((~f).sum()),minimum=float(a[f].min()) if f.any() else None,maximum=float(a[f].max()) if f.any() else None)
  for v in SENTINELS: r['equal_'+str(v)]=int((a==v).sum())
 else:
  s=np.char.strip(a.astype('S'));r['blank']=int((s==b'').sum())
 profiles.append(r)
def dataframe(d):
 return pd.DataFrame({c:np.char.strip(d[c].astype('U')) if d[c].dtype.kind in 'SU' else native(d[c]) for c in d.names})
def readtab(p):
 register(p,'fitted_parameter_or_classifier_table')
 first=p.open().readline()
 if first.startswith('VARNAMES:') or first.startswith('#'):
  lines=[x for x in p.read_text().splitlines() if x.startswith(('VARNAMES:','SN:'))]
  cols=lines[0].split()[1:]; vals=[x.split()[1:] for x in lines[1:]];t=pd.DataFrame(vals,columns=cols)
  for c in cols:
   if c not in ['CID','FIELD','BANDLIST']: t[c]=pd.to_numeric(t[c],errors='coerce')
  return t
 return pd.read_csv(p,dtype={'CID':str,'SNID':str})

def metadat():
 result={}; mt={}
 for name,base in [('original',ORIG),('current',CUR)]:
  p=base/'4_DISTANCES_COVMAT'/('DES-SN5YR_HD+MetaData.csv' if name=='original' else 'DES-Dovekie_Metadata.csv')
  d=readtab(p);d['CID']=d.CID.astype(str);mt[name]=d
  d.to_csv(OUT/f'{name}_metadata.csv.gz',index=False,compression={'method':'gzip','mtime':0})
  for c in d: profile(d[c].to_numpy(),name,'metadata',c)
  # Correlation scaling avoids treating small x0 variances as roundoff-zero eigenvalues.
  err=d[['x0ERR','x1ERR','cERR']].to_numpy(float)
  cov=np.zeros((len(d),3,3));cov[:,range(3),range(3)]=err**2
  for i,j,k in [(0,1,'COV_x1_x0'),(0,2,'COV_c_x0'),(1,2,'COV_x1_c')]: cov[:,i,j]=cov[:,j,i]=d[k].to_numpy(float)
  corr=cov/(err[:,:,None]*err[:,None,:]); eig=np.linalg.eigvalsh(corr)
  result[name]={'n':len(d),'duplicate_cid':int(d.CID.duplicated().sum()),'nonpositive_errors':int((err<=0).any(axis=1).sum()),'nonpositive_x0':int((d.x0<=0).sum()),'non_psd_correlation_count_tol_1e_8':int((eig[:,0]<-1e-8).sum()),'min_correlation_eigenvalue':float(eig[:,0].min()),'zero_COV_x1_x0':int((d.COV_x1_x0==0).sum()),'zero_COV_c_x0':int((d.COV_c_x0==0).sum()),'zero_COV_x1_c':int((d.COV_x1_c==0).sum())}
  pd.DataFrame({'CID':d.CID,'min_corr_eigenvalue':eig[:,0],'max_abs_correlation':np.abs(corr).max(axis=(1,2))}).to_csv(OUT/f'{name}_covariance_audit.csv.gz',index=False,compression={'method':'gzip','mtime':0})
  for p in sorted((base/'3_CLASSIFICATION').glob('*.csv')):
   c=readtab(p);key='CID' if 'CID' in c else 'SNID';c[key]=c[key].astype(str)
   for col in c: profile(c[col].to_numpy(),name,p.name,col)
   result[name][p.name]={'n':len(c),'duplicate_id':int(c[key].duplicated().sum()),'metadata_ids_matched':int(d.CID.isin(c[key]).sum()),'unmatched_classifier_ids':int((~c[key].isin(d.CID)).sum()),'probability_outside_0_1_nonmissing':{col:int(((c[col]<0)|(c[col]>1)).sum()) for col in c if col.startswith('PROB')}}
 return mt,result

def get_uncompressed(p):
 rec=register(p,'SNANA_FITS');target=CACHE/(rec['sha256']+'.fits')
 if not target.exists():
  tmp=target.with_suffix('.part')
  with gzip.open(p,'rb') as f,tmp.open('wb') as g: shutil.copyfileobj(f,g,8*1024*1024)
  tmp.rename(target)
 rec['gzip_crc_read_to_end']=True
 return target

def process(h5,name,hp,pp,mt,photcache):
 dataset=name; print('AUDIT',name,flush=True)
 hf=fits.open(get_uncompressed(hp),memmap=True);pf=fits.open(get_uncompressed(pp),memmap=True)
 hf.verify('exception');pf.verify('exception');h=hf[1].data;p=pf[1].data
 nh,n=len(h),len(p); hdr=dict(hf[0].header);ids=np.char.strip(h['SNID'].astype('U'))
 grp=h5.create_group(name);grp.attrs['head_source']=str(hp.relative_to(ROOT));grp.attrs['phot_source']=str(pp.relative_to(ROOT));grp.attrs['measurement_level']='released calibrated flux; not detector pixels'
 hg=grp.create_group('head')
 for c in h.names:
  arr=native(h[c]);hg.create_dataset(c,data=arr,compression='gzip',compression_opts=1,shuffle=True);profile(arr,name,'HEAD',c)
 grp.attrs['primary_header_json']=json.dumps({k:hdr[k] for k in ['SURVEY','VERSION','FILTERS','DATATYPE','SNANA_VERSION','MWEBV_APPLYFLAG','PHOTFLAG_DETECT'] if k in hdr})
 photsha=register(pp,'SNANA_FITS')['sha256']
 if photsha in photcache: grp['phot']=h5[photcache[photsha]]
 else:
  pg=grp.create_group('phot')
  for c in p.names:
   arr=native(p[c]);pg.create_dataset(c,data=arr,compression='gzip',compression_opts=1,shuffle=True,chunks=True)
  photcache[photsha]=name+'/phot'
 # All PHOT fields profiled, excluding documented delimiter rows for distributions only.
 mjd=np.asarray(p['MJD']); delim=mjd==-777
 for c in p.names: profile(p[c][~delim],name,'PHOT_nondelimiter',c)
 owner=np.full(n,-1,dtype=np.int32);overlap=0;badpointer=[];nobsmismatch=[];sentinelinside=[];missingafter=[]
 lo=np.asarray(h['PTROBS_MIN'],dtype=np.int64)-1;hi=np.asarray(h['PTROBS_MAX'],dtype=np.int64)
 for j,(a,b) in enumerate(zip(lo,hi)):
  if a<0 or b>n or b<a:badpointer.append(str(ids[j]));continue
  overlap+=int((owner[a:b]>=0).sum());owner[a:b]=j
  if b-a!=h['NOBS'][j]:nobsmismatch.append(str(ids[j]))
  if delim[a:b].any():sentinelinside.append(str(ids[j]))
  if b>=n or not delim[b]:missingafter.append(str(ids[j]))
 band=np.char.strip(p['BAND'].astype('S20'));field=np.char.strip(p['FIELD'].astype('S12'))
 flux=np.asarray(p['FLUXCAL']);err=np.asarray(p['FLUXCALERR']);flags=np.zeros(n,dtype=np.uint16)
 flags[delim]|=1; flags[owner<0]|=2;flags[(~np.isfinite(mjd))|(mjd<=0)]|=4
 flags[~np.isfinite(flux)]|=8;flags[(~np.isfinite(err))|(err<=0)]|=16;flags[band==b'']|=32
 # 1/2 warnings and 4096 detection retained; only documented DIFFIMG error bits mask.
 if name=='diffimg/DES': flags[(p['PHOTFLAG'] & (8|16|32|64|128|256|512|1024|2048))!=0]|=128
 basic=(flags&63)==0
 if name.endswith('/DES'):
  pixsize=np.where(owner>=0,np.asarray(h['PIXSIZE'])[np.maximum(owner,0)],np.nan)
  psf_fwhm_arcsec=p['PSF_SIG1']*(pixsize*2.355)
  flags[(psf_fwhm_arcsec<.5)|(psf_fwhm_arcsec>2.75)|~np.isfinite(psf_fwhm_arcsec)]|=256
  with np.errstate(divide='ignore',invalid='ignore'): zpnpe=p['ZEROPT']+2.5*np.log10(np.where(p['GAIN']<.01,.001,p['GAIN']))
  flags[(zpnpe<30.5)|(zpnpe>100)|~np.isfinite(zpnpe)]|=512
 else:zpnpe=np.full(n,np.nan)
 snr=np.divide(flux,err,out=np.full(n,np.nan,dtype=float),where=np.isfinite(err)&(err>0))
 rows=[]; ndupimg=0
 version=name.split('/')[0];meta=mt.get(version,mt['original']);ms=set(meta.CID);index=meta.set_index('CID')
 for j,(a,b) in enumerate(zip(lo,hi)):
  rr={'dataset':name,'head_row':j+1,'snid':ids[j],'source_head':str(hp.relative_to(ROOT)),'source_phot':str(pp.relative_to(ROOT)), 'in_original_hd':bool(ids[j] in set(mt['original'].CID)), 'in_current_hd':bool(ids[j] in set(mt['current'].CID))}
  if not (0<=a<=b<=n):rows.append(rr);continue
  rec=np.rec.fromarrays([mjd[a:b],band[a:b]],names='MJD,BAND');_,inverse,counts=np.unique(rec,return_inverse=True,return_counts=True);dup=counts[inverse]>1
  flags[a:b][dup]|=64
  if dup.any():
   for k in np.flatnonzero(dup): duplicates.append({'dataset':name,'snid':str(ids[j]),'phot_row':int(a+k+1),'MJD':float(mjd[a+k]),'band':band[a+k].decode(),'IMGNUM':int(p['IMGNUM'][a+k]),'flux':float(flux[a+k]),'fluxerr':float(err[a+k])})
  goodimage=p['IMGNUM'][a:b]>0
  ii=np.rec.fromarrays([p['IMGNUM'][a:b][goodimage],band[a:b][goodimage]],names='IMG,BAND');ndupimg+=len(ii)-len(np.unique(ii))
  good=basic[a:b];quality=good&((flags[a:b]&(128|256|512))==0)
  rr.update(nobs=int(b-a),basic_eligible_nobs=int(good.sum()),quality_eligible_nobs=int(quality.sum()),duplicate_mjd_band_rows=int(dup.sum()),n_negative_flux=int((flux[a:b][good]<0).sum()),mjd_min=float(mjd[a:b][good].min()) if good.any() else None,mjd_max=float(mjd[a:b][good].max()) if good.any() else None)
  mx=[]
  for flt in np.unique(band[a:b]):
   if not flt:continue
   bm=quality&(band[a:b]==flt);sn=snr[a:b][bm]; f=flt.decode()
   rr['n_'+f]=int(bm.sum());rr['snrmax_'+f]=float(np.max(sn)) if len(sn) else None
   if len(sn):mx.append(np.max(sn));rr['depth5_median_'+f]=float(np.median(27.5-2.5*np.log10(5*err[a:b][bm])))
  rr['n_bands_snr5']=int(np.sum(np.array(mx)>5))
  z=float(h['REDSHIFT_HELIO'][j]);ze=float(h['REDSHIFT_HELIO_ERR'][j]);ebv=float(h['MWEBV'][j]);pk=float(h['PEAKMJD'][j])
  rr['peakmjd_source']='header_estimate'
  if ids[j] in index.index:pk=float(index.loc[ids[j],'PKMJD']);rr['peakmjd_source']='published_fitted'
  rest=(mjd[a:b]-pk)/(1+z) if z>-1 else np.full(b-a,np.nan);fitwin=quality&(rest>=-15)&(rest<=45);broad=quality&(rest>=-20)&(rest<=60)
  rr.update(nobs_fitted_phase_proxy=int(fitwin.sum()),cut_z_des=bool(.05<=z<=1.2 and 0<=ze<=.01),cut_mwebv_des=bool(0<=ebv<=.25),cut_nepoch_proxy=bool(broad.sum()>=5),cut_snr2bands_proxy=bool(rr['n_bands_snr5']>=2),cut_phase_proxy=bool(broad.any() and rest[broad].min()<=5 and rest[broad].max()>=5))
  rr['partial_des_cut_pass']=all(rr[k] for k in ['cut_z_des','cut_mwebv_des','cut_nepoch_proxy','cut_snr2bands_proxy','cut_phase_proxy'])
  rows.append(rr)
 obj=dataframe(h);obj.insert(0,'head_row',np.arange(nh)+1);obj.insert(0,'dataset',name)
 rr=pd.DataFrame(rows);obj=obj.merge(rr,on=['dataset','head_row'],validate='one_to_one');all_objects.append(obj)
 grp.create_dataset('audit_flags',data=flags,compression='gzip',compression_opts=1,shuffle=True)
 grp.create_dataset('head_row_zero_based',data=owner,compression='gzip',compression_opts=1,shuffle=True)
 grp.create_dataset('basic_eligible',data=basic,compression='gzip',compression_opts=1)
 grp.create_dataset('snr',data=snr,compression='gzip',compression_opts=1,shuffle=True)
 grp.attrs['flag_bits_json']=json.dumps(FLAG_BITS)
 active=~delim
 def hist(a):
  v,c=np.unique(a,return_counts=True);return {str(x.decode() if isinstance(x,bytes) else x):int(y) for x,y in zip(v,c)}
 s={'dataset':name,'head_rows':nh,'phot_rows':n,'nondelimiter_rows':int(active.sum()),'delimiter_rows':int(delim.sum()),'header_columns':len(h.names),'phot_columns':len(p.names),'duplicate_snid':int(nh-len(set(ids))),'pointer_invalid_ids':badpointer,'nobs_mismatch_ids':nobsmismatch,'pointer_overlapping_rows':overlap,'delimiter_inside_span_ids':sentinelinside,'missing_delimiter_after_ids':missingafter,'unowned_nondelimiter_rows':int(((owner<0)&active).sum()),'basic_eligible_rows':int(basic.sum()),'audit_flag_counts_nondelimiter':{k:int(((flags&v!=0)&active).sum()) for k,v in FLAG_BITS.items()},'duplicate_image_band_excess_rows':ndupimg,'bands':hist(band[active]),'fields':hist(field[active]),'photflag':hist(p['PHOTFLAG'][active]),'snr_quantiles_basic':dict(zip(['min','p01','p05','p50','p95','p99','max'],[float(x) for x in np.quantile(snr[basic],[0,.01,.05,.5,.95,.99,1])])),'snr_threshold_counts':{str(x):int((snr[basic]>x).sum()) for x in [0,3,5,10,20]},'negative_flux_rows_basic':int((flux[basic]<0).sum()),'original_hd_members':int(rr.in_original_hd.sum()),'current_hd_members':int(rr.in_current_hd.sum()),'partial_des_cut_pass':int(rr.partial_des_cut_pass.sum()),'partial_des_cut_pass_published_members':int((rr.partial_des_cut_pass&rr[('in_original_hd' if version!='current' else 'in_current_hd')]).sum()),'cuts_warning':'DES partial proxies only; fitted peak for published members, header estimate otherwise; no model wavelength eligibility, iterative outlier rejection, exact SNANA cut ordering, BBC validity or classifier eligibility.'}
 summaries.append(s);jsonout('progress.json',summaries)
 hf.close();pf.close();print('DONE',name,nh,n,flush=True)

def main():
 reg=register(OUT/'preregistration.json','prospective_local_register')
 mt,md=metadat(); photcache={}
 with h5py.File(OUT/'photometry_audit.h5','w') as h5:
  h5.attrs['schema_version']='des-calibrated-photometry-audit-v1';h5.attrs['preregistration_sha256']=reg['sha256'];h5.attrs['index_contract']='HEAD row and phot_row use 1-based external indices; head_row_zero_based is Python mapping, -1 unowned. Source pointer MIN inclusive 1-based, MAX inclusive 1-based.'
  for ver,base in [('original',ORIG),('current',CUR)]:
   for typ in ['DES','Foundation','LOWZ']:
    d=base/'0_DATA'/('DES-SN5YR_'+typ);stem='DES-SN5YR_'+typ
    process(h5,ver+'/'+typ,d/(stem+'_HEAD.FITS.gz'),d/(stem+'_PHOT.FITS.gz'),mt,photcache)
  d=ROOT/'data/des-diffimg';process(h5,'diffimg/DES',d/'DES-SN5YR_DIFFIMG_HEAD.FITS.gz',d/'DES-SN5YR_DIFFIMG_PHOT.FITS.gz',mt,photcache)
 pd.concat(all_objects,ignore_index=True).to_csv(OUT/'objects.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 pd.DataFrame(profiles).to_csv(OUT/'field_profiles.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 pd.DataFrame(duplicates).to_csv(OUT/'duplicate_epochs.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 differences={}
 for typ in ['DES','Foundation','LOWZ']:
  stem='DES-SN5YR_'+typ;differences[typ]={}
  for table in ['HEAD','PHOT']:
   a=ORIG/'0_DATA'/stem/(stem+'_'+table+'.FITS.gz');b=CUR/'0_DATA'/stem/(stem+'_'+table+'.FITS.gz');differences[typ][table+'_bytes_identical']=sha(a)==sha(b)
  with fits.open(get_uncompressed(ORIG/'0_DATA'/stem/(stem+'_HEAD.FITS.gz')),memmap=True) as a,fits.open(get_uncompressed(CUR/'0_DATA'/stem/(stem+'_HEAD.FITS.gz')),memmap=True) as b:
   diff={}
   for c in a[1].columns.names:
    x,y=a[1].data[c],b[1].data[c]; neq=x!=y
    if x.dtype.kind=='f':neq &= ~(np.isnan(x)&np.isnan(y))
    if neq.any():diff[c]={'changed_rows':int(neq.sum())}
   differences[typ]['head_field_changes']=diff
 register(Path(__file__),'analysis_script')
 register(ROOT/'sources/repos/RickKessler__SNANA/src/snana.F90','reference_for_PSF_FWHM_and_electron_zeropoint_cut_units')
 for base in [ORIG,CUR]:
  for sub in ['0_DATA','3_CLASSIFICATION','7_PIPPIN_FILES']:
   for p in sorted((base/sub).rglob('*')):
    if p.is_file() and (p.suffix.lower() in ['.md','.readme','.nml','.input','.yml','.ignore','.list']):register(p,'source_documentation_or_pipeline_configuration')
 for p in [ROOT/'data/des-diffimg/DES-SN5YR_DIFFIMG.README',ROOT/'catalog/des-diffimg-zenodo.json'] :register(p,'source_documentation_or_provenance')
 for name in ['des-science__DES-SN5YR','des-science__DES-SN5YR@1.3']:
  register(ROOT/'catalog/repositories'/(name+'.commit.json'),'pinned_repository_commit')
 jsonout('summary.json',{'schema_version':'des-photometry-audit-v1','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'datasets':summaries,'fitted_metadata':md,'cross_release_identity':differences,'flag_bits':FLAG_BITS,'source_level':'calibrated flux; no detector pixels'})
 jsonout('inputs_manifest.json',{'files':list(INPUTS.values())})
 outputs=[]
 for p in sorted(OUT.iterdir()):
  if p.is_file() and p.name!='outputs_manifest.json': outputs.append({'path':str(p.relative_to(ROOT)),'bytes':p.stat().st_size,'sha256':sha(p)})
 jsonout('outputs_manifest.json',{'files':outputs})
 print('COMPLETE',flush=True)
if __name__=='__main__':main()
