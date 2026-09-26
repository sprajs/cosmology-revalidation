#!/usr/bin/env python3
from pathlib import Path
import json,numpy as np,pandas as pd
R=Path(__file__).resolve().parents[3];O=R/'phase2/assumptions';assert (O/'off_signal_followup_preregistration.json').exists()
D=pd.read_csv(O/'off_signal_rows.csv.gz',dtype={'CID':str});M=pd.read_csv(R/'phase2/data_audit/original_metadata.csv.gz',dtype={'CID':str}).set_index('CID');D['PIa']=D.CID.map(M.PROB_SNNV19);D['PIa_bin']=np.where(D.PIa>.9,'>0.9','<=0.9');D['side']=np.where(D.phase_days<0,'pre','post');D['season']=np.floor((D.MJD-56474)/365.25).astype(int)
Q=D[(D.dataset=='diffimg/DES')&(abs(D.phase_days)>365)&D.quality].copy();out=[];pred=[];rng=np.random.default_rng(231059)
def stats(d,col,label,extra=None):
 x=d[col].to_numpy(dtype=float);med=np.median(x);r={'label':label,**(extra or {}),'n':len(d),'objects':d.CID.nunique(),'mean':float(x.mean()),'rms':float(np.sqrt(np.mean(x*x))),'mad_sigma':float(1.4826*np.median(abs(x-med))),'abs_gt3_fraction':float(np.mean(abs(x)>3)),'abs_gt5_fraction':float(np.mean(abs(x)>5))}
 if not extra:
  z=d.assign(x=x,x2=x*x,t3=abs(x)>3).groupby('CID').agg(n=('x','size'),s=('x','sum'),s2=('x2','sum'),t3=('t3','sum')).to_numpy();v=[]
  for b in range(1000):
   a=z[rng.integers(0,len(z),len(z))].sum(axis=0);v.append([a[1]/a[0],np.sqrt(a[2]/a[0]),a[3]/a[0]])
  for i,n in enumerate(['mean','rms','abs_gt3_fraction']):r[n+'_cluster_ci95']=np.quantile(np.array(v)[:,i],[.025,.975]).tolist()
 out.append(r)
stats(Q,'pull','original')
for c in ['side','PIa_bin','BAND','host_sb_bin']:
 for v,d in Q.groupby(c):stats(d,'pull','original',{c:str(v)})
for (cid,band),d in Q.groupby(['CID','BAND']):
 for parity in [0,1]:
  train=d[d.season%2==parity];test=d[d.season%2!=parity].copy()
  if len(train)<5 or len(test)<2:continue
  w=1/train.FLUXCALERR.to_numpy(dtype=float)**2;y=train.FLUXCAL.to_numpy(dtype=float);base=np.sum(w*y)/sum(w);formal=1/sum(w)
  scores=pd.Series(w*(y-base),index=np.floor(train.MJD).astype(int)).groupby(level=0).sum();n=len(scores);robust=float(np.sum(scores**2)/sum(w)**2*n/(n-1)) if n>1 else formal
  test['train_n']=len(train);test['train_seasons']=train.season.nunique();test['baseline']=base;test['baseline_var_formal']=formal;test['baseline_var_night_sandwich']=max(formal,robust);test['heldout_pull_formal']=(test.FLUXCAL-base)/np.sqrt(test.FLUXCALERR**2+formal);test['heldout_pull_robust']=(test.FLUXCAL-base)/np.sqrt(test.FLUXCALERR**2+max(formal,robust));pred.append(test)
T=pd.concat(pred,ignore_index=True);assert not T.duplicated(['CID','phot_row']).any();T.to_csv(O/'baseline_heldout_rows.csv.gz',index=False,compression={'method':'gzip','mtime':0})
for col in ['pull','heldout_pull_formal','heldout_pull_robust']:
 stats(T,col,'heldout_support_'+col)
 for c in ['side','PIa_bin','BAND','host_sb_bin']:
  for v,d in T.groupby(c):stats(d,col,'heldout_support_'+col,{c:str(v)})
P=[]
for (cid,band),d in T.groupby(['CID','BAND']):
 d=d.sort_values('MJD').drop_duplicates('MJD');x=d.heldout_pull_robust.to_numpy();dt=np.diff(d.MJD)
 for lo,hi in [(0,1),(1,10),(10,60),(60,1e9)]:
  k=(dt>lo)&(dt<=hi)
  if k.any():P.append({'CID':cid,'BAND':band,'gap_low':lo,'gap_high':hi,'n':int(k.sum()),'xy':float(np.sum(x[:-1][k]*x[1:][k])),'x2':float(np.sum(x[:-1][k]**2)),'y2':float(np.sum(x[1:][k]**2))})
pd.DataFrame(P).to_csv(O/'baseline_heldout_pair_moments.csv.gz',index=False,compression={'method':'gzip','mtime':0});pairout=[]
for (b,lo,hi),d in pd.DataFrame(P).groupby(['BAND','gap_low','gap_high']):
 z=d[['n','xy','x2','y2']].to_numpy();v=z.sum(axis=0);boot=[]
 for k in range(1000):
  a=z[rng.integers(0,len(z),len(z))].sum(axis=0);boot.append(a[1]/np.sqrt(a[2]*a[3]))
 pairout.append({'BAND':b,'gap_low':int(lo),'gap_high':int(hi),'n':int(v[0]),'normalized_uncentered_product':float(v[1]/np.sqrt(v[2]*v[3])),'cluster_ci95':np.quantile(boot,[.025,.975]).tolist()})
# Same pixels, separate reductions: exact source keys, multiplicity-sensitive merge disclosed.
S=D[(D.dataset=='original/DES')&D.quality].copy();F=D[(D.dataset=='diffimg/DES')&D.quality].copy();keys=['CID','IMGNUM','MJD','BAND'];S['key_count']=S.groupby(keys).CID.transform('size');F['key_count']=F.groupby(keys).CID.transform('size')
A=S[S.key_count==1].merge(F[F.key_count==1],on=keys,suffixes=('_SMP','_DIFFIMG'));A['flux_difference']=A.FLUXCAL_SMP-A.FLUXCAL_DIFFIMG;A['error_ratio']=A.FLUXCALERR_DIFFIMG/A.FLUXCALERR_SMP;A.to_csv(O/'cross_reduction_off_signal.csv.gz',index=False,compression={'method':'gzip','mtime':0})
cross={'exact_unique_key_pairs':len(A),'objects':A.CID.nunique(),'smp_eligible_rows':len(S),'diffimg_eligible_rows':len(F),'caution':'Exact time match only; unmatched and multiplicities excluded. Same pixels, errors correlated; no quadrature significance.'}
if len(A):
 cross['median_error_ratio_DIFFIMG_over_SMP']=float(np.median(A.error_ratio));cross['median_flux_difference_SMP_minus_DIFFIMG']=float(np.median(A.flux_difference));cross['median_pull_SMP']=float(np.median(A.pull_SMP));cross['median_pull_DIFFIMG']=float(np.median(A.pull_DIFFIMG))
(O/'off_signal_followup.json').write_text(json.dumps({'statistics':out,'pair_statistics':pairout,'cross_reduction':cross,'baseline_note':'Cross-season predicted constant, uncertainties include formal weighted-mean or max(formal, night-cluster sandwich). Season covariance and contaminating flux remain confounded; this is a diagnostic, not a released-photometry correction.'},indent=2)+'\n');print(json.dumps([x for x in out if not any(k in x for k in ['BAND','side','PIa_bin','host_sb_bin'])],indent=2));print(cross)
