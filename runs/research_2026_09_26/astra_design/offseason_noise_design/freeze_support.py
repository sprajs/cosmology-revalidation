"""Independent selected-row metadata verification and design freeze; no flux access."""
from pathlib import Path
import hashlib,json
from collections import defaultdict,deque
import pandas as pd,numpy as np
from astropy.io import fits
P=Path(__file__).resolve().parent;R=P.parents[3]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
D=pd.read_csv(P/'pre30-epoch-metadata.csv.gz',low_memory=False)
D=D[(D.arm=='SMP')&(D.rest_phase<=-60)].copy()
A=pd.read_csv(P/'object-cadence.csv');A=A[(A.arm=='SMP')&(A.min_rest_phase<=-60)]
M=pd.read_csv(R/'phase2/data_audit/original_metadata.csv.gz',usecols=['CID','IDSURVEY','PKMJD','zHEL']);M=M[M.IDSURVEY==10];M.CID=M.CID.astype(int)
M=M.rename(columns={'PKMJD':'published_fitted_peak','zHEL':'published_zHEL'});A=A.merge(M,on='CID',validate='one_to_one');A['clump_minus_published']=A.clump_peak-A.published_fitted_peak
A.to_csv(P/'early-object-peak-agreement.csv',index=False,float_format='%.17g');A=A.set_index('CID')
base=R/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES'
with fits.open(base/'DES-SN5YR_DES_HEAD.FITS.gz') as hf,fits.open(base/'DES-SN5YR_DES_PHOT.FITS.gz') as pf:
    h=hf[1].data;p=pf[1].data;idx={int(x['SNID']):x for x in h};absolute=[]
    for cid,g in D.groupby('CID',sort=False):
        head=idx[int(cid)];lo=int(head['PTROBS_MIN'])-1;hi=int(head['PTROBS_MAX']);match=defaultdict(deque)
        for j in range(lo,hi):match[(float(p['MJD'][j]),str(p['BAND'][j]).strip(),int(p['IMGNUM'][j]))].append(j)
        for row in g.itertuples():
            # CSV round-trip MJD tolerance checked against exact float metadata.
            possible=[k for k in match if abs(k[0]-row.MJD)<1e-8 and k[1]==row.BAND and k[2]==int(row.IMGNUM)]
            assert len(possible)==1 and match[possible[0]]
            j=match[possible[0]].popleft();absolute.append((row.Index,j+1))
            assert abs((float(p['MJD'][j])-A.loc[cid,'clump_peak'])/(1+float(head['REDSHIFT_HELIO']))-row.rest_phase)<1e-8
            assert int(p['PHOTFLAG'][j])==int(row.PHOTFLAG)
D['phot_row_one_based']=pd.Series(dict(absolute));assert D.phot_row_one_based.is_unique
D['primary_quality']=D.q1016&D.image_quality
D['strict_pre100']=D.rest_phase<=-100
D['published_peak_rest_phase']=(D.MJD-D.CID.map(A.published_fitted_peak))/(1+D.CID.map(A.published_zHEL))
D['head_peak_rest_phase']=(D.MJD-D.CID.map(A.head_peak))/(1+D.CID.map(A.zHEL))
D.to_csv(P/'frozen-pre60-rows.csv',index=False,float_format='%.17g')
summary=[]
for cohort in ['hubble1635','validation1020','common1006']:
    for cut in [-60,-100]:
        d=D[(D[cohort].astype(str)=='True')&D.primary_quality&(D.rest_phase<=cut)]
        n=d.groupby(['CID','BAND']).night.nunique();byid=d.groupby('CID').size();ns=d.groupby('CID').night.nunique()
        eligible=set(n[n>=2].index.get_level_values(0));pairs=[]
        for (cid,band),g in d.groupby(['CID','BAND']):
            t=np.sort(g.groupby('night').MJD.mean());dt=np.diff(t)
            for k,x in enumerate(dt):pairs.append(dict(CID=cid,BAND=band,dt=float(x)))
        pp=pd.DataFrame(pairs)
        summary.append(dict(cohort=cohort,cut=cut,n=len(d),objects=int(d.CID.nunique()),fields=int(d.FIELD.nunique()),nights=int(d.night.nunique()),groups=len(n),groups_ge2=int((n>=2).sum()),groups_ge4=int((n>=4).sum()),objects_with_pairs=len(eligible),residual_intercept_dof=int(sum(n-1)),adjacent_pairs=len(pairs),pairs_gt60=int((pp.dt>60).sum()) if len(pp) else 0,largest_object_rows=int(byid.max()),largest_object_fraction=float(byid.max()/len(d)),all_peak_definitions_pre60=int(((d.published_peak_rest_phase<=-60)&(d.head_peak_rest_phase<=-60)).sum()),all_peak_definitions_pre100=int(((d.published_peak_rest_phase<=-100)&(d.head_peak_rest_phase<=-100)).sum())))
result=dict(status='METADATA PASS; no flux/noise score executed',row_identity='323 source rows uniquely mapped by CID/MJD/BAND/IMGNUM and exact flags; no duplicates on primary window',summary=summary,peak_agreement=dict(max_abs_clump_minus_published=float(abs(A.clump_minus_published).max()),max_abs_clump_minus_head=float(abs(A.clump_minus_head).max())),flag_counts={str(int(k)):int(v) for k,v in D.PHOTFLAG.value_counts().items()},source_provenance_gate='Per-epoch free-vs-fixed SMP status and original DiffImg scene-fit peak are not documented in available per-row columns; not inferred from PHOTFLAG bit1 or flux values.')
(P/'frozen-support.json').write_text(json.dumps(result,indent=2)+'\n')
paths=[P/'metadata_support.py',P/'manifest.json',P/'support.json',P/'object-cadence.csv',P/'pre30-epoch-metadata.csv.gz',P/'frozen-pre60-rows.csv',P/'frozen-support.json',P/'early-object-peak-agreement.csv',Path(__file__),R/'phase2/data_audit/published_lcplot_flags.json']
(P/'freeze-manifest.json').write_text(json.dumps({'sha256':{str(p.relative_to(R)):sha(p) for p in paths}},indent=2)+'\n')
print(json.dumps(result,indent=2))
