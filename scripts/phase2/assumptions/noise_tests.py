#!/usr/bin/env python3
"""Preregistered off-signal tests of calibrated flux, without any SALT refit."""
from pathlib import Path
import hashlib,json
import h5py,numpy as np,pandas as pd
R=Path(__file__).resolve().parents[3];O=R/'phase2/assumptions'
assert (O/'preregistration.json').exists()
M=pd.read_csv(R/'phase2/data_audit/original_metadata.csv.gz',dtype={'CID':str});M=M[M.IDSURVEY==10].set_index('CID')
def strings(x):return np.char.strip(np.asarray(x).astype('U'))
rows=[]
with h5py.File(R/'phase2/data_audit/photometry_audit.h5') as h:
 for ds in ['original/DES','diffimg/DES']:
  g=h[ds];hd=g['head'];cols=['MJD','BAND','FIELD','FLUXCAL','FLUXCALERR','PHOTFLAG','PSF_SIG1','GAIN','ZEROPT','IMGNUM','CCDNUM']
  arrays={c:g['phot/'+c][:] for c in cols};ids=strings(hd['SNID'][:]);mins=hd['PTROBS_MIN'][:]-1;maxs=hd['PTROBS_MAX'][:];pix=hd['PIXSIZE'][:]
  for j,cid in enumerate(ids):
   if cid not in M.index:continue
   a,b=int(mins[j]),int(maxs[j]); t=arrays['MJD'][a:b]; phase=t-M.loc[cid,'PKMJD']; keep=np.flatnonzero(np.abs(phase)>200)
   if not len(keep):continue
   d=pd.DataFrame({c:arrays[c][a:b][keep] for c in cols});d['dataset']=ds;d['CID']=cid;d['phase_days']=phase[keep];d['phot_row']=a+keep+1
   d.BAND=strings(d.BAND);d.FIELD=strings(d.FIELD);d['pull']=d.FLUXCAL/d.FLUXCALERR;d['basic']=np.isfinite(d.MJD)&(d.MJD>0)&np.isfinite(d.FLUXCAL)&np.isfinite(d.FLUXCALERR)&(d.FLUXCALERR>0)
   mask=sum([8,16,32,64,128,256,512,1024,2048]) if ds=='diffimg/DES' else 32
   psf=d.PSF_SIG1*pix[j]*2.355;zp=d.ZEROPT+2.5*np.log10(np.where(d.GAIN<.01,.001,d.GAIN))
   d['quality']=d.basic & ((d.PHOTFLAG&mask)==0)&psf.between(.5,2.75)&(zp>=30.5)
   sb=np.full(len(d),np.nan)
   for band in 'griz':
    key='HOSTGAL_SB_FLUXCAL_'+band
    if key in hd:
     x=float(hd[key][j]);sb[d.BAND==band]=27.5-2.5*np.log10(x) if x>0 else np.nan
   d['host_sb_mag']=sb;d['host_sb_bin']=np.where(~np.isfinite(sb),'missing',np.where(sb<22,'<22',np.where(sb<24,'22..24','>=24')))
   rows.append(d)
  print(ds,'collected',sum(len(d) for d in rows),flush=True)
D=pd.concat(rows,ignore_index=True);D.to_csv(O/'off_signal_rows.csv.gz',index=False,compression={'method':'gzip','mtime':0})
rng=np.random.default_rng(731150);stats=[];objects=[];pairs=[]
def summarize(d,attrs,bootstrap=False):
 x=d.pull.to_numpy();n=len(x)
 if not n:return
 med=np.median(x);r={**attrs,'n_rows':n,'n_objects':d.CID.nunique(),'mean':float(np.mean(x)),'median':float(med),'std':float(np.std(x)),'rms':float(np.sqrt(np.mean(x*x))),'mad_sigma':float(1.4826*np.median(np.abs(x-med))),'abs_gt3_fraction':float(np.mean(np.abs(x)>3)),'abs_gt5_fraction':float(np.mean(np.abs(x)>5)),'negative_fraction':float(np.mean(x<0))}
 if bootstrap:
  z=d.assign(pull2=x*x,tail3=(np.abs(x)>3).astype(float),tail5=(np.abs(x)>5).astype(float)).groupby('CID').agg(n=('pull','size'),s=('pull','sum'),s2=('pull2','sum'),t3=('tail3','sum'),t5=('tail5','sum')).to_numpy();B=np.empty((1000,4))
  for b in range(1000):
   v=z[rng.integers(0,len(z),len(z))].sum(axis=0);B[b]=[v[1]/v[0],np.sqrt(v[2]/v[0]),v[3]/v[0],v[4]/v[0]]
  for k,name in enumerate(['mean','rms','abs_gt3_fraction','abs_gt5_fraction']):r[name+'_cluster_ci95']=np.quantile(B[:,k],[.025,.975]).tolist()
 stats.append(r)
for ds in ['original/DES','diffimg/DES']:
 for window in [365,200]:
  for mask in ['basic','quality']:
   d=D[(D.dataset==ds)&(abs(D.phase_days)>window)&D[mask]].copy();attrs={'dataset':ds,'window_abs_days':window,'mask':mask,'stratum':'all','value':'all'}
   summarize(d,attrs,True)
   for by in ['BAND','FIELD','host_sb_bin']:
    for val,sub in d.groupby(by,sort=True):summarize(sub,{**attrs,'stratum':by,'value':str(val)})
   if window!=365 or mask!='quality':continue
   for (cid,band),sub in d.groupby(['CID','BAND']):
    x=sub.pull.to_numpy();objects.append({'dataset':ds,'CID':cid,'BAND':band,'n':len(x),'mean':np.mean(x),'rms':np.sqrt(np.mean(x*x))})
    # Exact time ties excluded from pair test; no claim that distinct times are independent.
    sub=sub.sort_values(['MJD','phot_row']).drop_duplicates('MJD');x=sub.pull.to_numpy();dt=np.diff(sub.MJD.to_numpy());prod=x[:-1]*x[1:]
    for lo,hi in [(0,1),(1,10),(10,60),(60,1e9)]:
     sel=(dt>lo)&(dt<=hi)
     if sel.any():pairs.append({'dataset':ds,'CID':cid,'BAND':band,'gap_low':lo,'gap_high':hi,'n_pairs':int(sel.sum()),'sum_products':float(prod[sel].sum()),'sum_left2':float((x[:-1][sel]**2).sum()),'sum_right2':float((x[1:][sel]**2).sum())})
pd.DataFrame(stats).to_json(O/'off_signal_statistics.json',orient='records',indent=2)
pd.DataFrame(objects).to_csv(O/'off_signal_object_band.csv.gz',index=False,compression={'method':'gzip','mtime':0})
P=pd.DataFrame(pairs);P.to_csv(O/'off_signal_pairs.csv.gz',index=False,compression={'method':'gzip','mtime':0});ps=[]
for (ds,band,lo,hi),d in P.groupby(['dataset','BAND','gap_low','gap_high']):
 z=d.groupby('CID')[['n_pairs','sum_products','sum_left2','sum_right2']].sum().to_numpy();v=z.sum(axis=0);boot=[]
 for b in range(1000):
  a=z[rng.integers(0,len(z),len(z))].sum(axis=0);boot.append(a[1]/np.sqrt(a[2]*a[3]))
 ps.append({'dataset':ds,'BAND':band,'gap_low':int(lo),'gap_high':int(hi),'n_pairs':int(v[0]),'n_objects':len(z),'normalized_uncentered_product':float(v[1]/np.sqrt(v[2]*v[3])),'cluster_ci95':np.quantile(boot,[.025,.975]).tolist(),'interpretation':'Common baseline offsets, real variability, repeated exposures and noise correlation are not separately identified by this statistic.'})
(O/'off_signal_pair_statistics.json').write_text(json.dumps(ps,indent=2)+'\n')
print(json.dumps([s for s in stats if s['stratum']=='all'],indent=2))
