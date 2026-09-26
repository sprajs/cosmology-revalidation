"""Metadata-only pre-explosion support audit. Never access FLUXCAL or SIM truth.

FITS record storage necessarily contains other columns; this script accesses only
the enumerated cadence, measured-coordinate and quality metadata columns.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd
from astropy.io import fits

ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
HASH={}
def hashed(p):
    p=ROOT/p if not Path(p).is_absolute() else Path(p)
    HASH[str(p.relative_to(ROOT))]=hashlib.sha256(p.read_bytes()).hexdigest()
    return p
def csv(p,cols):return pd.read_csv(hashed(p),usecols=cols)
def snana(p):
    lines=hashed(p).read_text().splitlines()
    cols=next(x.split()[1:] for x in lines if x.startswith('VARNAMES:'))
    selected=['CID','PKMJDINI','zHEL']
    ii=[cols.index(c)+1 for c in selected]
    return pd.DataFrame([[x.split()[i] for i in ii] for x in lines if x.startswith('SN:')],columns=selected).astype({'CID':int,'PKMJDINI':float,'zHEL':float}).set_index('CID')
M=csv('phase2/data_audit/original_metadata.csv.gz',['CID','IDSURVEY','zHEL'])
real_ids=set(M.loc[M.IDSURVEY==10,'CID'].astype(int));assert len(real_ids)==1635
vids=set(csv('runs/research_2026_09_26/astra_design/validation1020/cohort.csv',['CID']).CID);assert len(vids)==1020
rp=csv('runs/research_2026_09_26/common_classifier_residual/inference/real-probabilities.csv',['CID','gt999'])
cids=set(rp.loc[rp.gt999,'CID']);assert len(cids)==1006 and cids<=vids
simprob=pd.concat([csv(p,['arm','CID','gt999']) for p in ['runs/research_2026_09_26/common_classifier_residual/inference/simulation-probabilities.csv','runs/research_2026_09_26/simulation_holdout_control/inference/holdout-probabilities.csv']])
assert not simprob.duplicated(['arm','CID']).any()
PBASE=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES'
datasets=[('SMP',PBASE/'DES-SN5YR_DES','phase2/classification/reconstruction_20260926/clump_run/DES-SN5YR_DES.SNANA.TEXT')]
for arm in ['P21','G10']:
    v='PH2_pilot02_'+arm
    datasets.append((arm,ROOT/'phase2/literature/simulations/outputs'/v/v,f'runs/research_2026_09_26/common_classifier_residual/clump/{arm}/{v}.SNANA.TEXT'))
rows=[];epochs=[];schemas={};stages={}
PCOLS=['MJD','BAND','FIELD','PHOTFLAG','PSF_SIG1','GAIN','ZEROPT','IMGNUM']
def native(x):
    a=np.array(x)
    return a.astype(a.dtype.newbyteorder('='))
for arm,stem,pkpath in datasets:
    pk=snana(pkpath);assert pk.index.is_unique
    ext='.FITS.gz' if arm=='SMP' else '.FITS'
    hp=hashed(Path(str(stem)+'_HEAD'+ext));pp=hashed(Path(str(stem)+'_PHOT'+ext))
    with fits.open(hp,memmap=False) as f:
        hd=f[1].data
        hc=['SNID','NOBS','PTROBS_MIN','PTROBS_MAX','PIXSIZE','REDSHIFT_HELIO','PEAKMJD']
        hd={c:native(hd[c]) for c in hc if c in hd.names}
    with fits.open(pp,memmap=False) as f:
        schemas[arm]=dict(phot_columns=f[1].columns.names,accessed_phot_columns=PCOLS,accessed_head_columns=list(hd),phot_rows=len(f[1].data))
        phot={c:native(f[1].data[c]) for c in PCOLS if c in f[1].columns.names}
    if arm=='SMP': ss={'hubble1635':real_ids,'validation1020':vids,'common1006':cids}
    else:
        ff=csv(f'phase2/hierarchy/forward/{arm}-fitted.csv.gz',['CID','basic_quality_pass'])
        fid=set(ff.CID);quality=set(ff.loc[ff.basic_quality_pass,'CID'])
        classified=set(simprob.loc[(simprob.arm==arm)&simprob.gt999,'CID'])
        ss={'all_written':set(map(int,hd['SNID'])),'archived_fit':fid,'archived_quality':quality,'quality_common_classifier':quality&classified}
    stages[arm]={k:len(v) for k,v in ss.items()}
    for i,id0 in enumerate(hd['SNID']):
        cid=int(id0)
        if not any(cid in v for v in ss.values()):continue
        if cid not in pk.index:raise ValueError(('missing measured peak',arm,cid))
        p=float(pk.loc[cid,'PKMJDINI']);z=float(hd['REDSHIFT_HELIO'][i])
        a=int(hd['PTROBS_MIN'][i])-1;b=int(hd['PTROBS_MAX'][i]);assert b-a==hd['NOBS'][i]
        d=pd.DataFrame({c:phot[c][a:b] for c in phot})
        d['BAND']=d.BAND.astype(str).str.strip();d['FIELD']=d.FIELD.astype(str).str.strip()
        assert (d.MJD>0).all()
        valid=np.isfinite(p) and p>50000 and np.isfinite(z) and z>0
        if not valid:raise ValueError(('invalid measured peak/redshift',arm,cid,p,z))
        phase=(d.MJD-p)/(1+z);rawpeak=float(hd.get('PEAKMJD',np.full(len(hd['SNID']),np.nan))[i])
        rows.append(dict(arm=arm,CID=cid,zHEL=z,clump_peak=p,head_peak=rawpeak,clump_minus_head=p-rawpeak,nobs=len(d),min_rest_phase=phase.min(),max_rest_phase=phase.max(),min_observer_phase=(d.MJD-p).min(),head_relative_min=(d.MJD-rawpeak).min(),**{k:cid in ids for k,ids in ss.items()}))
        # Keep only phase<=-30 for all later support calculations. No residuals.
        d=d[phase<=-30].copy();d['rest_phase']=phase[phase<=-30]
        if len(d)==0:continue
        d['arm']=arm;d['CID']=cid;d['night']=np.floor(d.MJD-.5).astype(int)
        d['q32']=((d.PHOTFLAG.astype(int)&32)==0)
        d['q1016']=((d.PHOTFLAG.astype(int)&1016)==0)
        psf=d.PSF_SIG1*float(hd['PIXSIZE'][i])*2.355
        zp=d.ZEROPT+2.5*np.log10(np.where(d.GAIN<.01,.001,d.GAIN))
        d['image_quality']=psf.between(.5,2.75)&zp.between(30.5,100)
        for k,ids in ss.items():d[k]=cid in ids
        epochs.append(d)
D=pd.concat(epochs,ignore_index=True);R=pd.DataFrame(rows)
R.to_csv(OUT/'object-cadence.csv',index=False,float_format='%.17g')
D.to_csv(OUT/'pre30-epoch-metadata.csv.gz',index=False,compression={'method':'gzip','mtime':0})
stats=[];breakdowns=[]
for arm,cohorts in stages.items():
    for cohort in cohorts:
        source=D[(D.arm==arm)&D[cohort].fillna(False).astype(bool)]
        for cut in [-100,-60,-40,-30]:
            for mask in ['all_metadata','q32_image','q1016_image']:
                d=source[source.rest_phase<=cut]
                if mask!='all_metadata':d=d[d[mask.split('_')[0]]&d.image_quality]
                nights=d[['CID','BAND','night']].drop_duplicates();nn=nights.groupby(['CID','BAND']).size()
                paircounts={str(k):0 for k in ['0..1','1..10','10..60','>60']}
                for _,g in d.groupby(['CID','BAND']):
                    # Mean cadence of a night is metadata only, no flux weights.
                    dt=np.diff(np.sort(g.groupby('night').MJD.mean()))
                    for k,test in [('0..1',(dt>0)&(dt<=1)),('1..10',(dt>1)&(dt<=10)),('10..60',(dt>10)&(dt<=60)),('>60',dt>60)]:paircounts[k]+=int(test.sum())
                stats.append(dict(arm=arm,cohort=cohort,rest_phase_max=cut,mask=mask,rows=len(d),objects=d.CID.nunique(),fields=d.FIELD.nunique(),nights=d.night.nunique(),object_band_nights=len(nights),object_band_groups=len(nn),object_band_ge2nights=int((nn>=2).sum()),object_band_ge4nights=int((nn>=4).sum()),adjacent_distinct_night_pairs=paircounts))
                if mask=='all_metadata':
                    for (band,field),g in d.groupby(['BAND','FIELD']):breakdowns.append(dict(arm=arm,cohort=cohort,rest_phase_max=cut,band=band,field=field,rows=len(g),objects=g.CID.nunique(),nights=g.night.nunique()))
pd.DataFrame(breakdowns).to_csv(OUT/'support-band-field.csv',index=False)
(OUT/'support.json').write_text(json.dumps(dict(stages=stages,schemas=schemas,support=stats,flux_columns_accessed=False,simulation_truth_columns_accessed=False),indent=2)+'\n')
hashed(Path(__file__))
for p in ['papers/text/2406.05046v1.txt','docs/phase2/snana-assumptions.md','sources/repos/des-science__DES-SN5YR@1.3/0_DATA/README.md','sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES/DES-SN5YR_DES.README','phase2/literature/simulations/inputs/PH2_pilot02_P21.input','phase2/literature/simulations/inputs/PH2_pilot02_G10.input']:hashed(p)
(OUT/'manifest.json').write_text(json.dumps(dict(inputs_sha256=HASH,outputs_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name not in ['manifest.json','metadata_support.py']}),indent=2)+'\n')
print(json.dumps([s for s in stats if s['mask']=='all_metadata' and s['cohort'] in ['hubble1635','common1006','quality_common_classifier']],indent=2))
