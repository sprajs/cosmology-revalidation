#!/usr/bin/env python3
"""Derived flag-only mask, conditional on published accepted epochs, no deduplication."""
from pathlib import Path
import json, hashlib, shutil, datetime, gzip
import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import linear_sum_assignment

R=Path(__file__).resolve().parents[3]
O=R/'phase2/assumptions/conditional_mask'; O.mkdir(parents=True,exist_ok=True)
S=R/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA'
D=S/'DES-SN5YR_DES'; BIT=1<<29
spec={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Conditioned-mask baseline diagnosis, not independently selected cosmology sample','bit':BIT,'assignment':'Per CID, float32 MJD printed 3 decimals and band; one-to-one assignment uses relative flux and error differences. Numerical ties resolved by source row order and explicitly recorded. All raw rows preserved, no deduplication.','unaccepted':'Bit 2^29 set on every source PHOT row not assigned to a published DATA_MODEL=1 row; other bits and every other field unchanged.','ambiguity':'Rounded public flux and error can map multiple raw rows; retain matching multiplicity and disclose nonunique assignments.'}
if not (O/'preregistration.json').exists(): (O/'preregistration.json').write_text(json.dumps(spec,indent=2)+'\n')
p=pd.read_csv(S/'DES5YR_SALT3_LCFIT.LCPLOT.gz',sep=r'\s+',comment='#',dtype={'SNID':str})
p=p[p.DATA_MODEL==1].copy()
h=fits.open(D/'DES-SN5YR_DES_HEAD.FITS.gz',character_as_bytes=True); f=fits.open(D/'DES-SN5YR_DES_PHOT.FITS.gz',character_as_bytes=True); hd=h[1].data; ph=f[1].data
ids={v.decode().strip():i for i,v in enumerate(hd['SNID'])}; accepted=np.zeros(len(ph),bool); assignments=[]; ambiguity=[]
for cid,sub in p.groupby('SNID',sort=False):
 j=ids[cid];a=int(hd['PTROBS_MIN'][j])-1;b=int(hd['PTROBS_MAX'][j]); q=ph[a:b]
 mjd=np.asarray([float(f'{x:.3f}') for x in q['MJD'].astype(np.float32)]);band=np.char.strip(q['BAND'].astype('U'))
 for (t,k),group in sub.groupby(['MJD','BAND'],sort=False):
  cand=np.flatnonzero((mjd==t)&(band==k)); assert len(cand)>=len(group),(cid,t,k,'too few source rows')
  pf=group.FLUXCAL.to_numpy(); pe=group.FLUXCALERR.to_numpy(); rf=q['FLUXCAL'][cand];re=q['FLUXCALERR'][cand]
  df=np.abs(pf[:,None]-rf[None,:])/np.maximum(1,np.abs(rf))[None,:];de=np.abs(pe[:,None]-re[None,:])/re[None,:]
  valid=(df<5e-5)&(de<5e-5);cost=df+de;cost[~valid]=1e6
  rr,cc=linear_sum_assignment(cost); assert np.all(valid[rr,cc]),(cid,t,k,'no valid assignment')
  for u,v in zip(rr,cc):
   raw=a+int(cand[v]);accepted[raw]=True
   assignments.append({'CID':cid,'public_index':int(group.index[u]),'phot_row_one_based':raw+1,'raw_flag':int(ph['PHOTFLAG'][raw]),'possible_matches':int(valid[u].sum())})
  if np.any(valid.sum(axis=1)>1):
   ambiguity.append({'CID':cid,'MJD':t,'BAND':k,'public_accepted':len(group),'candidate_source_rows':[a+int(v)+1 for v in cand],'possible_counts':valid.sum(axis=1).tolist(),'chosen_source_rows':[a+int(cand[v])+1 for v in cc],'possible_flags':ph['PHOTFLAG'][a+cand].astype(int).tolist(),'all_candidates_kept':len(cc)==len(cand)})
orig=ph['PHOTFLAG'].copy(); assert not np.any(orig & BIT)
ph['PHOTFLAG'][~accepted] |= BIT
name='DES-SN5YR_DES_CONDITIONAL'
f[0].header['PHOTFILE']=name+'_PHOT.FITS' if 'PHOTFILE' in f[0].header else name+'_PHOT.FITS'
# SNANA HEAD primary header names PHOTFILE; change only that navigation field.
h[0].header['PHOTFILE']=name+'_PHOT.FITS'
for suffix in ['HEAD','PHOT']:
 with gzip.open(D/('DES-SN5YR_DES_'+suffix+'.FITS.gz'),'rb') as src, (O/(name+'_'+suffix+'.FITS')).open('wb') as dst: shutil.copyfileobj(src,dst)
with fits.open(O/(name+'_HEAD.FITS'),mode='update',character_as_bytes=True) as derived_head:
 derived_head[0].header['PHOTFILE']=name+'_PHOT.FITS'
# Direct strided update avoids Astropy rewriting FITS fixed-width spaces as NUL.
mm=np.memmap(O/(name+'_PHOT.FITS'),mode='r+',offset=f[1].fileinfo()['datLoc'],dtype=ph.dtype,shape=(len(ph),))
mm['PHOTFLAG'][:]=ph['PHOTFLAG'];mm.flush();del mm
(O/(name+'.LIST')).write_text(name+'_HEAD.FITS\n')
(O/(name+'.README')).write_text('DERIVED CONDITIONED-MASK DIAGNOSTIC. Original release fluxes and row multiplicities unchanged. PHOTFLAG bit536870912 rejects rows not matched one-to-one to published DATA_MODEL=1 epochs. See preregistration.json and ambiguity.json.\n')
pd.DataFrame(assignments).to_csv(O/'assignments.csv.gz',index=False,compression={'method':'gzip','mtime':0})
(O/'ambiguity.json').write_text(json.dumps(ambiguity,indent=2)+'\n')
f.close();h.close()
with fits.open(O/(name+'_PHOT.FITS'),character_as_bytes=True) as check, fits.open(D/'DES-SN5YR_DES_PHOT.FITS.gz',character_as_bytes=True) as source:
 for col in ph.names:
  if col!='PHOTFLAG':
   same=np.array_equal(check[1].data[col],source[1].data[col],equal_nan=True) if source[1].data[col].dtype.kind=='f' else np.array_equal(check[1].data[col],source[1].data[col])
   assert same,col
 assert np.array_equal(check[1].data['PHOTFLAG'] & ~BIT,orig)
 assert np.count_nonzero((check[1].data['PHOTFLAG'] & BIT)==0)==len(p)
summary={'accepted_rows':len(p),'accepted_objects':p.SNID.nunique(),'source_phot_rows':len(ph),'ambiguous_groups':len(ambiguity),'ambiguous_subset_groups':sum(not v['all_candidates_kept'] for v in ambiguity),'verified_all_nonflag_columns_unchanged':True,'conditional_not_independent_selection':True,'files':[]}
for path in [Path(__file__),S/'DES5YR_SALT3_LCFIT.LCPLOT.gz',D/'DES-SN5YR_DES_HEAD.FITS.gz',D/'DES-SN5YR_DES_PHOT.FITS.gz']+sorted(O.glob('*')):
 if path.is_file() and path.name!='manifest.json': summary['files'].append({'path':str(path.relative_to(R)),'bytes':path.stat().st_size,'sha256':hashlib.file_digest(path.open('rb'),'sha256').hexdigest()})
(O/'manifest.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='files'},indent=2))
