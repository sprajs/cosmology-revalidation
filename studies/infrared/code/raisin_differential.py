#!/usr/bin/env python3
"""Pinned RAISIN w-distance accounting and frozen descriptive paired contrasts.

Run `freeze` before `audit` and `contrast`. No cosmological parameters are fitted.
The LCDM curve in audit only verifies the already-exported mures column.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import nnls

ROOT = Path(__file__).resolve().parents[2]
REL = ROOT / 'sources/repos/djones1040__RAISIN_DataRelease@a383c4b'
OUT = ROOT / 'runs/research_2026_09_26/raisin_differential'
MAN = ROOT / 'sources/updates/2026-09-26-raisin/acquisition-manifest.json'
BRANCHES = ['nir', 'optical', 'opticalnir']
GROUPS = dict(mwebv=[1], hstcal=[2], pecvel=[3], massstep=[4,5], lcfitter=[6],
              lowzcal=list(range(7,15)), photcal=[2]+list(range(7,23)),
              biascorst=[23,25], biascorav=[24,26], biascor=[23,24,25,26],
              kcor=[27], tmpl=[28], all=list(range(1,29)))

def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(name, obj):
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/name).write_text(json.dumps(obj, indent=2, allow_nan=False)+'\n')
def table(branch, variant=0):
    p=REL/f'distances/w/{branch}_dist/RAISIN_combined_FITOPT{variant:03d}.FITRES'
    return fitres_file(p)
def fitres_file(p):
    lines=p.read_text().splitlines()
    cols=next(x.split()[1:] for x in lines if x.startswith('VARNAMES:'))
    d=pd.DataFrame([x.split()[1:] for x in lines if x.startswith('SN:')],columns=cols)
    assert not d.CID.duplicated().any()
    d=d.set_index('CID')
    for col in d:
        try: d[col]=pd.to_numeric(d[col])
        except ValueError: pass
    return d
def lcparams(branch, group):
    return np.loadtxt(REL/f'distances/w/{branch}_syst/RAISIN_{group}_lcparams_cosmosis.txt')
def offcov(branch, group):
    x=np.loadtxt(REL/f'distances/w/{branch}_syst/RAISIN_{group}.covmat')
    return x[1:].reshape(int(x[0]),int(x[0]))
def syscov(branch, group):
    """Restore diagonal folded into the group-specific lcparams dmb."""
    return offcov(branch,group)+np.diag(lcparams(branch,group)[:,5]**2-lcparams(branch,'stat')[:,5]**2)
def photometry():
    result={}
    # One calibrator outside w (2007sr) has two alternate peak-header files.
    for p in sorted((REL/'photometry/RAISIN').rglob('*')):
        if not p.is_file() or p.suffix.lower()!='.dat': continue
        heads={}; obs=[]; cols=None; duplicates=[]
        for line in p.read_text().splitlines():
            if ':' not in line or line.startswith('#'): continue
            key,val=line.split(':',1); val=val.strip()
            if key=='VARLIST': cols=val.split()
            elif key=='OBS': obs.append(dict(zip(cols,val.split())))
            elif key in heads: duplicates.append(key)
            else: heads[key]=val
        if 'SNID' not in heads: continue
        cid=heads['SNID'].split()[0]
        if cid in result:
            assert cid not in table('nir').index, 'Ambiguous w-sample light curve'
            result[cid].setdefault('alternate_paths_outside_w_sample',[]).append(str(p.relative_to(ROOT)))
            continue
        result[cid]=dict(CID=cid, path=str(p.relative_to(ROOT)), sha256=sha(p),
            survey=heads['SURVEY'].split()[0], RA=float(heads['RA'].split()[0]),
            DEC=float(heads['DEC'].split()[0]), zHEL_phot=float(heads['REDSHIFT_HELIO'].split()[0]),
            peak_phot=float(heads['PEAKMJD'].split()[0]), nobs=len(obs),
            declared_nobs=int(heads['NOBS']), bands=dict(Counter(x['FLT'] for x in obs)),
            negative_flux_count=sum(float(x['FLUXCAL'])<0 for x in obs),
            nonpositive_error_count=sum(float(x['FLUXCALERR'])<=0 for x in obs), duplicate_header_keys=duplicates)
    return result

def freeze():
    if (OUT/'protocol.json').exists(): raise RuntimeError('Protocol exists; do not overwrite it.')
    manifest=json.loads(MAN.read_text()); verify=[]
    for row in manifest['files']:
        p=ROOT/row['path']; b=p.read_bytes()
        ok=sha(p)==row['sha256'] and hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()==row['git_blob']
        verify.append(ok)
    assert all(verify) and len(verify)==490
    phot=photometry(); tabs={b:table(b) for b in BRANCHES}; base=tabs['nir']
    assert all(list(d.index)==list(base.index) for d in tabs.values())
    members=[]
    for cid,row in base.iterrows():
        assert cid in phot
        for br,d in tabs.items():
            assert float(d.loc[cid,'zHD'])==float(row.zHD)
            assert float(d.loc[cid,'HOST_LOGMASS'])==float(row.HOST_LOGMASS)
        members.append(dict(CID=cid, survey=phot[cid]['survey'], RA=phot[cid]['RA'],DEC=phot[cid]['DEC'],
             zHEL=float(row.zHEL), zHD=float(row.zHD), host_logmass=float(row.HOST_LOGMASS),
             stratum='low' if row.zHD<.1 else 'high', high_host_mass=bool(row.HOST_LOGMASS>=10),
             photometry_path=phot[cid]['path']))
    assert len(members)==79 and sum(x['stratum']=='low' for x in members)==42
    assert all(x['zHD']>.2 or x['zHD']<.1 for x in members)
    pd.DataFrame(members).to_csv(OUT/'frozen-membership.csv',index=False)
    protocol=dict(created_utc=datetime.now(timezone.utc).isoformat(), stage='before paired redshift outcomes',
        release_commit=manifest['commit'],manifest_sha256=sha(MAN),verified_files=490,
        membership_sha256=sha(OUT/'frozen-membership.csv'),membership_rule='Exact common physical CIDs of nominal w optical/NIR/combined branches; coordinate-confirmed raw light curves; no outcome exclusion.',
        primary='T = mean_high(DLMAG_optical-DLMAG_nir) - mean_low(DLMAG_optical-DLMAG_nir), unweighted; high zHD>0.2, low zHD<0.1; units mag.',
        cancellation='Each paired difference cancels common distance-redshift relation; high-minus-low cancels any branch-constant intercept. No cosmology is fitted.',
        correction_variants={'released':'DLMAG', 'without_mass':'DLMAG - MASS_CORR',
             'without_bias':'DLMAG + DLMAG_biascor', 'without_bias_or_mass':'DLMAG + DLMAG_biascor - MASS_CORR'},
        correction_scope='These undo exported additive terms only, leaving original light-curve fitting and selection fixed. They are not raw flux analyses, alternative validated corrections, or regenerated selection corrections.',
        secondary='Same estimand for opticalnir minus nir; descriptive only, strongly overlapping information.',
        metadata_splits=['PS1 high vs all low CSP; DES high vs all low CSP (shared low denominator)',
             'host logmass>=10 and <10, each high vs low; threshold fixed by published host-mass imbalance, no threshold search'],
        uncertainty=['Paired empirical SE sqrt(s_high^2/n_high+s_low^2/n_low), descriptive finite-sample scatter, not a full systematic error.',
             'Published statistical marginal errors with unknown per-object cross-branch correlation: rho=0 illustration and sharp [-1,+1] pairwise bounds, conditional on independence across SNe.',
             'Systematic marginal variance bound: |sqrt(a Copt a)-sqrt(a Cnir a)| to sum, without inventing cross covariance.',
             'Only physically documented signed systematic mappings may define a conditional shared-response covariance; separately report any failure of export closure.'],
        no_claims=['No p-value, cosmology fit, independent multiplication of branch/sample evidence, dust-law inference, correction tuning, or inferred cross-branch intrinsic covariance.'])
    write('protocol.json',protocol)
    write('photometry-inventory.json',list(phot.values()))
    write('provenance-verification.json',dict(verified_sha256_and_git_blob=len(verify), manifest_sha256=sha(MAN),
        files_by_top_directory=dict(Counter(Path(x['path']).relative_to(REL.relative_to(ROOT)).parts[0] for x in manifest['files']))))
    print(json.dumps({'frozen':79,'counts':dict(Counter(x['survey'] for x in members)), 'protocol_sha256':sha(OUT/'protocol.json')}))

def overlap(members):
    """1 arcsec coordinate matching, with redshift differences retained for audit."""
    def skysep(ra,dec,ras,decs):
        r1,d1,r2,d2=np.deg2rad([ra,dec,0,0]); r2=np.deg2rad(ras);d2=np.deg2rad(decs)
        return np.rad2deg(2*np.arcsin(np.sqrt(np.clip(np.sin((d1-d2)/2)**2+np.cos(d1)*np.cos(d2)*np.sin((r1-r2)/2)**2,0,1))))*3600
    pp=ROOT/'sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/4_DISTANCES_AND_COVAR/Pantheon+SH0ES.dat'
    p=pd.read_csv(pp,sep=r'\s+')
    dh=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz'
    with fits.open(dh) as h:
        d=pd.DataFrame({k:np.asarray(h[1].data[k]).astype(str if k in ['SNID','IAUC'] else float) for k in ['SNID','IAUC','RA','DEC','REDSHIFT_HELIO']})
    hd=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/4_DISTANCES_COVMAT/DES-SN5YR_HD+MetaData.csv'
    fitids=set(pd.read_csv(hd).CID.astype(str))
    rows=[]
    for row in members.itertuples():
        for name,tab,cidcol,zcol in [('Pantheon+',p,'CID','zHEL'),('DES5YR_photometry',d,'SNID','REDSHIFT_HELIO')]:
            sep=skysep(row.RA,row.DEC,tab.RA.to_numpy(),tab.DEC.to_numpy()); ix=np.flatnonzero(sep<=1)
            for i in ix:
                t=tab.iloc[i]; target=str(t[cidcol]).strip()
                rows.append(dict(CID=row.CID,source=name,matched_ID=target,sep_arcsec=float(sep[i]),
                    delta_zHEL=float(row.zHEL-t[zcol]), in_DES5YR_distance=target in fitids if name.startswith('DES') else False,
                    IDSURVEY=str(t.IDSURVEY) if name=='Pantheon+' else '',IAUC=str(t.IAUC) if name.startswith('DES') else ''))
    matches=pd.DataFrame(rows);matches.to_csv(OUT/'overlap-matches.csv',index=False)
    write('overlap-summary.json',dict(match_radius_arcsec=1, reference_hashes={str(x.relative_to(ROOT)):sha(x) for x in [pp,dh,hd]},
        Pantheon_unique_RAisin_objects=int(matches.loc[matches.source=='Pantheon+','CID'].nunique()),
        Pantheon_matched_rows=int((matches.source=='Pantheon+').sum()),
        DES5YR_photometry_unique=int(matches.loc[matches.source=='DES5YR_photometry','CID'].nunique()),
        DES5YR_distances_unique=int(matches.loc[matches.in_DES5YR_distance,'CID'].nunique()),
        note='Multiple Pantheon+ rows for one physical supernova are not independent objects. DES photometry and distance membership reported separately. Redshift differences do not affect the coordinate match and are exposed in the ledger.'))

def audit():
    protocol=json.loads((OUT/'protocol.json').read_text()); assert sha(OUT/'frozen-membership.csv')==protocol['membership_sha256']
    members=pd.read_csv(OUT/'frozen-membership.csv');overlap(members)
    cosm=FlatLambdaCDM(H0=70,Om0=.3); summary={}; modes=[]
    for br in BRANCHES:
        d=table(br); stat=lcparams(br,'stat'); mu=cosm.distmod(d.zHD.to_numpy()).value
        intercept=(d.DLMAG-d.MASS_CORR-d.mures-mu).to_numpy()
        branch=dict(n=len(d),cosmosis_max_abs_mb_minus_DLMAG_plus19_36=float(np.max(abs(stat[:,4]-d.DLMAG+19.36))),
            cosmosis_max_abs_redshift_minus_zHD=float(np.max(abs(stat[:,1:3]-d.zHD.to_numpy()[:,None]))),
            stat_error_export_max_abs=float(np.max(abs(stat[:,5]-d.DLMAGERR))),
            mures_identity='mures = DLMAG - MASS_CORR - mu_LCDM_H070_Om0.3 - offset',
            mures_offset=float(np.mean(intercept)),mures_identity_max_error=float(np.max(abs(intercept-np.mean(intercept)))),
            nir_fixed_parameter_check={key:sorted(d[key].unique().tolist()) for key in ['STRETCH','AV','PKMJDERR']} if br=='nir' else {},
            correction_only_variant_closure={},covariance={})
        for k in [4,5,23,24,25,26]:
            s=table(br,k);closure=s.DLMAG-d.DLMAG+s.DLMAG_biascor-d.DLMAG_biascor-(s.MASS_CORR-d.MASS_CORR)
            branch['correction_only_variant_closure'][str(k)]=float(np.max(abs(closure)))
        for group,ks in GROUPS.items():
            C=syscov(br,group);ev=np.linalg.eigvalsh(C); off=offcov(br,group)
            pars=lcparams(br,group)
            branch['covariance'][group]=dict(offdiag_max_diagonal=float(abs(np.diag(off)).max()),
                max_mb_difference_from_stat=float(abs(pars[:,4]-stat[:,4]).max()),
                reconstructed_diag_range=[float(np.diag(C).min()),float(np.diag(C).max())],
                reconstructed_eigenvalue_range=[float(ev[0]),float(ev[-1])],symmetry_max_error=float(abs(C-C.T).max()))
            # Project out the intercept on both sides. This removes arbitrary
            # centering, but does not alter the predeclared high-low estimand.
            n=len(d);H=np.eye(n)-np.ones((n,n))/n; target=(H@C@H)
            for flavor in ['DLMAG','HD','without_mass']:
                vs=[]
                for k in ks:
                    s=table(br,k).loc[d.index];v=(s.DLMAG-d.DLMAG).to_numpy()
                    if flavor=='HD': v-=cosm.distmod(s.zHD.to_numpy()).value-mu
                    if flavor=='without_mass':v-=(s.MASS_CORR-d.MASS_CORR).to_numpy()
                    vs.append(H@v)
                X=np.stack([np.outer(v,v).ravel() for v in vs],axis=1)
                if np.linalg.norm(target)==0: weights=np.zeros(len(ks));error=0
                else:
                    weights,_=nnls(X,target.ravel(),maxiter=10000)
                    error=np.linalg.norm(X@weights-target.ravel())/np.linalg.norm(target)
                modes.append(dict(branch=br,group=group,flavor=flavor,relative_covariance_residual=float(error),
                    variant_indices=ks,nonnegative_outer_product_weights=weights.tolist(),
                    qualification='Numerically recovered marginal export weights, not inferred physical priors or guaranteed cross-branch coupling.'))
        summary[br]=branch
    write('algebra-covariance-audit.json',summary);write('systematic-reconstruction.json',modes)
    inventory=[]
    for p in sorted(REL.glob('**/*.fits')):
        with fits.open(p) as hdus:
            inventory.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),hdus=[dict(name=h.name,naxis1=h.header.get('NAXIS1'),naxis2=h.header.get('NAXIS2'),tfields=h.header.get('TFIELDS'),duplicate_column_names=len(h.columns.names)!=len(set(h.columns.names)) if hasattr(h,'columns') else False) for h in hdus]))
    write('fits-inventory.json',inventory)
    supplementary_audit(members)
    print(json.dumps({'audit':'complete','branches':list(summary),'paired_outcomes':'not computed by this audit command'}))

def supplementary_audit(members):
    phot=photometry();nir=table('nir'); nml=[]; filters=[]
    for p in sorted((REL/'lcfitting').glob('*.nml')):
        lines=p.read_text().splitlines(); modes=[]; settings={}
        for line in lines:
            m=re.match(r'\s*-\s*/([^/]+)/\s*(.*)',line)
            if m:modes.append(dict(FITOPT=len(modes)+1,tag=m[1],command=m[2]))
            if line.strip().startswith(('#','!')):continue
            if '=' in line:
                key,val=line.split('=',1)
                if key.strip() in ['FILTLIST_FIT','TREST_REJECT','FITMODEL_NAME','KCOR_FILE','HEADER_OVERRIDE_FILE','INIVAL_SHAPE','INISTP_SHAPE','INIVAL_AV','INISTP_AV','INISTP_PEAKMJD','INIVAL_RV','OPT_COVAR_FLUX']:settings[key.strip()]=val.strip()
        nml.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),modes=modes,settings=settings))
    write('nml-signed-variants.json',nml)
    for p in sorted((REL/'kcor').glob('*.fits')):
        with fits.open(p) as h:
            dat=h['FilterTrans'].data;wave=dat.field(0)
            for name in dat.columns.names[1:]:
                tr=dat[name];assert np.isfinite(tr).all() and np.all(np.diff(wave)>0)
                filters.append(dict(file=p.name,filter=name,grid_min_A=float(wave.min()),grid_max_A=float(wave.max()),
                    response_min=float(tr.min()),response_max=float(tr.max()),integral_A=float(np.trapezoid(tr,wave))))
    pd.DataFrame(filters).to_csv(OUT/'kcor-filter-inventory.csv',index=False)
    diffs={br:np.abs(table(br).PKMJD-nir.PKMJD).to_numpy() for br in ['optical','opticalnir']}
    pd.DataFrame(dict(CID=nir.index,nir_peak=nir.PKMJD,phot_header_peak=[phot[c]['peak_phot']for c in nir.index],
         optical_peak=table('optical').PKMJD,combined_peak=table('opticalnir').PKMJD)).to_csv(OUT/'peak-origin-audit.csv',index=False)
    write('supplementary-audit.json',dict(photometry_physical_objects=len(phot),photometry_file_count=len(list((REL/'photometry').rglob('*.DAT')))+len(list((REL/'photometry').rglob('*.dat'))),
        photometry_by_survey=dict(Counter(x['survey']for x in phot.values())),
        nominal_nir_peak_max_abs_difference_from_header=float(max(abs(phot[c]['peak_phot']-nir.loc[c,'PKMJD'])for c in nir.index)),
        nominal_nir_vs_exported_other_peak_abs_quantiles={br:np.quantile(x,[0,.5,1]).tolist()for br,x in diffs.items()},
        nir_FITCHI2_per_NDOF_quantiles=np.quantile(nir.FITCHI2/nir.NDOF,[0,.5,.9,1]).tolist(),
        nir_DLMAGERR_range=[float(nir.DLMAGERR.min()),float(nir.DLMAGERR.max())],
        mass_divide_variant_check={br:dict(hosts_between10_and10_44=int(((table(br).HOST_LOGMASS>=10)&(table(br).HOST_LOGMASS<10.44)).sum()),
            sign_flips_in_fitopt004=int(((table(br).MASS_CORR*table(br,4).MASS_CORR)<0).sum()),
            baseline_values=table(br).MASS_CORR.unique().tolist(),fitopt004_values=table(br,4).MASS_CORR.unique().tolist())for br in BRANCHES},
        missing_referenced_systematic_kcor=[x for x in ['kcor_CSPDR3_BD17_sys.fits','kcor_DES_NIR_sys.fits','kcor_PS1MD_NIR_sys.fits']if not(REL/'kcor'/x).exists()],
        archive_has_fitopt029=any(REL.glob('distances/w/*_dist/*FITOPT029*')), all_nml_has_29_variants=all(len(x['modes'])==29 for x in nml),
        covariance_group_nonadditivity_rel_frobenius={br:float(np.linalg.norm(sum(syscov(br,g) for g in ['photcal','massstep','lcfitter','biascor','kcor','tmpl','pecvel','mwebv'])-syscov(br,'all'))/np.linalg.norm(syscov(br,'all')))for br in BRANCHES}))
    dh=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz'
    nearest=[]
    with fits.open(dh) as h:
        d=h[1].data;c=SkyCoord(d.RA*u.deg,d.DEC*u.deg)
        for row in members.query("survey=='DES'").itertuples():
            sep=c.separation(SkyCoord(row.RA*u.deg,row.DEC*u.deg)).arcsec;i=np.argmin(sep);cid=str(d.SNID[i]).strip()
            nearest.append(dict(CID=row.CID,DES_CID=cid,sep_arcsec=float(sep[i]),zHEL=float(d.REDSHIFT_HELIO[i]),
                peakMJD=float(d.PEAKMJD[i]),RAISIN_phot_peak=phot[row.CID]['peak_phot'],
                in_validation1020=(ROOT/'runs/research_2026_09_26/astra_design/validation1020/analysis/objects'/f'{cid}.npz').exists()))
    pd.DataFrame(nearest).to_csv(OUT/'des-nearest-metadata.csv',index=False)
    roster=ROOT/'sources/updates/2026-09-26-training-roster/roster.csv'
    if roster.exists():
        t=pd.read_csv(roster,dtype=str).fillna('')
        ras=t.RA.map(lambda x:float(x.split()[0])).to_numpy();decs=t.DEC.map(lambda x:float(x.split()[0])).to_numpy()
        sky=SkyCoord(ras*u.deg,decs*u.deg);matches=[]
        for row in members.itertuples():
            sep=sky.separation(SkyCoord(row.RA*u.deg,row.DEC*u.deg)).arcsec
            for i,train in enumerate(t.itertuples()):
                same_name=row.CID.lower()==train.SNID.lower()
                if sep[i]<=1 or (same_name and sep[i]<=3):
                    matches.append(dict(CID=row.CID,training_SNID=train.SNID,training_survey=train.SURVEY,sep_arcsec=float(sep[i]),
                        match_rule='sky<=1arcsec' if sep[i]<=1 else 'exact_name_and_sky<=3arcsec',
                        input_list=train.input_list,training_file=train.archive_path))
        df=pd.DataFrame(matches);df.to_csv(OUT/'training-input-overlap.csv',index=False)
        write('training-input-overlap-summary.json',dict(roster_sha256=sha(roster),listed_files=len(t),
            RAISIN_physical_objects=int(df.CID.nunique()),strict_1arcsec=int(df.loc[df.sep_arcsec<=1,'CID'].nunique()),
            exact_name_rounding_additions=int(df.loc[df.sep_arcsec>1,'CID'].nunique()),
            counts_by_RAisin_survey=df.merge(members[['CID','survey']],on='CID').groupby('survey').CID.nunique().to_dict(),
            qualification='Public K21 config-listed input proxy, not an execution-linked accepted SALT3.DES5YR training roster; optical calibration changes and CFA-U exclusions documented by DES. No light-curve outcomes used to set membership.'))
    # Primary-paper provenance and explicit-table sign check (three rows copied
    # from Jones2022 Table7). Do not treat the rounded table as a raw FITRES dump.
    table7={'PScB480464':40.188,'PScH540087':40.621,'PScF520188':40.350}
    raw_offsets=[float(nir.loc[c,'DLMAG']+nir.loc[c,'DLMAG_biascor']-nir.loc[c,'MASS_CORR']-v)for c,v in table7.items()]
    write('primary-paper-provenance.json',dict(papers=[dict(url=f'https://arxiv.org/pdf/{tag}',sha256=sha(OUT/f'{tag}.pdf'))for tag in ['2201.07801v2','2402.18624v2']],
        appendix_table7_raw_values=table7,offsets_released_plus_bias_minus_mass_minus_table_raw=raw_offsets,
        note='Offsets agree within published 0.001mag rounding, supporting the released minus-bias sign despite appendix prose saying plus; original raw fit export remains absent.'))

def contrast():
    protocol=json.loads((OUT/'protocol.json').read_text())
    assert sha(OUT/'frozen-membership.csv')==protocol['membership_sha256']
    assert (OUT/'algebra-covariance-audit.json').exists()
    members=pd.read_csv(OUT/'frozen-membership.csv').set_index('CID')
    tabs={b:table(b).loc[members.index] for b in BRANCHES}
    low=(members.stratum=='low').to_numpy(); high=~low
    selectors={'all':(high,low),
       'PS1_high_CSP_low':((members.survey=='PS1MD').to_numpy(),low),
       'DES_high_CSP_low':((members.survey=='DES').to_numpy(),low),
       'host_ge10':(high&(members.host_logmass>=10).to_numpy(),low&(members.host_logmass>=10).to_numpy()),
       'host_lt10':(high&(members.host_logmass<10).to_numpy(),low&(members.host_logmass<10).to_numpy())}
    variants={name:(int('without_bias' in name),-int(name in ['without_mass','without_bias_or_mass'])) for name in protocol['correction_variants']}
    rows=[]; result=[]; sensitivities=[]; conditional=[]
    for branch in ['optical','opticalnir']:
        o,n=tabs[branch],tabs['nir']
        for name,(bs,ms) in variants.items():
            delta=(o.DLMAG-n.DLMAG+bs*(o.DLMAG_biascor-n.DLMAG_biascor)+ms*(o.MASS_CORR-n.MASS_CORR)).to_numpy()
            for i,cid in enumerate(members.index):rows.append(dict(CID=cid,branch_minus_nir=branch,correction_variant=name,delta_mag=float(delta[i])))
            for split,(hi,lo) in selectors.items():
                a=hi/hi.sum()-lo/lo.sum()
                # Each statistical pair has an unknown rho in [-1,1].
                so=o.DLMAGERR.to_numpy();sn=n.DLMAGERR.to_numpy()
                empirical_se=np.sqrt(np.var(delta[hi],ddof=1)/hi.sum()+np.var(delta[lo],ddof=1)/lo.sum())
                item=dict(branch_minus_nir=branch,correction_variant=name,split=split,high_n=int(hi.sum()),low_n=int(lo.sum()),
                   high_mean_mag=float(delta[hi].mean()),low_mean_mag=float(delta[lo].mean()),contrast_mag=float(a@delta),
                   empirical_paired_se_mag=float(empirical_se),statistical_rho0_se_mag=float(np.sqrt(np.sum(a*a*(so*so+sn*sn)))),
                   statistical_per_object_rho_bounds_se_mag=[float(np.sqrt(np.sum(a*a*(so-sn)**2))),float(np.sqrt(np.sum(a*a*(so+sn)**2)))])
                if name=='released':
                    vo=float(a@syscov(branch,'all')@a);vn=float(a@syscov('nir','all')@a)
                    assert vo>=0 and vn>=0
                    item.update(systematic_marginal_se_mag=[np.sqrt(vo),np.sqrt(vn)],
                      systematic_unknown_cross_covariance_se_bounds_mag=[abs(np.sqrt(vo)-np.sqrt(vn)),np.sqrt(vo)+np.sqrt(vn)],
                      total_branchwise_stat_plus_sys_unknown_cross_se_bounds_mag=[abs(np.sqrt(vo+np.sum(a*a*so*so))-np.sqrt(vn+np.sum(a*a*sn*sn))),np.sqrt(vo+np.sum(a*a*so*so))+np.sqrt(vn+np.sum(a*a*sn*sn))])
                result.append(item)
        a=high/high.sum()-low/low.sum()
        for k in range(1,29):
            ov=table(branch,k).loc[members.index];nv=table('nir',k).loc[members.index]
            shift=(ov.DLMAG-o.DLMAG)-(nv.DLMAG-n.DLMAG)
            sensitivities.append(dict(branch_minus_nir=branch,variant=k,contrast_shift_mag=float(a@shift.to_numpy()),
                interpretation='Released signed variation, not a standardized independent sigma unless documented.'))
        # Physical common-input mapping is documented by NML and paper, rather
        # than choosing relative signs from the observed paired distances.
        # Keep the components as sensitivities, not a complete joint covariance.
        for group,ks in {'calibration':[2]+list(range(7,23)),'MW':[1],
                         'population_shape':[23,25],'population_AV':[24,26]}.items():
            shifts=[x['contrast_shift_mag'] for x in sensitivities if x['branch_minus_nir']==branch and x['variant'] in ks]
            conditional.append(dict(branch_minus_nir=branch,component=group,variants=ks,
                independent_unit_variant_rss_mag=float(np.linalg.norm(shifts)),
                status='Conditional shared signed-input response under unit documented variation sizes. Partial sensitivity only; not full paired systematic uncertainty.'))
    pd.DataFrame(rows).to_csv(OUT/'paired-object-ledger.csv',index=False)
    pd.DataFrame(result).to_csv(OUT/'paired-contrasts.csv',index=False)
    pd.DataFrame(sensitivities).to_csv(OUT/'systematic-contrast-shifts.csv',index=False)
    write('results.json',dict(protocol_sha256=sha(OUT/'protocol.json'),
         interpretation='Descriptive conditional distances, not new cosmology or independent dust/correction evidence.',
         contrasts=result,conditional_shared_response_components=conditional,
         covariance_warning='No cross-branch statistical/intrinsic covariance is supplied. Marginal systematic export reconstruction is not globally closed. Report bounds rather than a full paired likelihood.'))
    print(json.dumps([x for x in result if x['split']=='all']))

def freeze_mass_counterfactual():
    name='mass-threshold/counterfactual-protocol.json'
    if (OUT/name).exists():raise RuntimeError('Counterfactual protocol exists; do not overwrite.')
    code=OUT/'mass-threshold/code/raisin_cosmo/cosmo_sys.py'
    write(name,dict(created_utc=datetime.now(timezone.utc).isoformat(),
        label='Conditional correction-propagation repair, not a refit or cosmology measurement',
        source_commit='b888214a5cbae38ac0bf886488ce7734f5b77e87',source_sha256=sha(code),
        source_evidence='cosmo_sys.py 550-566 estimates mixture membership at10.44 for MASS_DIVIDE;621-624 applies fitted step at hard-coded10;659-672 adds centered signed variant outer products with unit weight.',
        frozen_release_manifest_sha256=sha(MAN),primary_protocol_sha256=sha(OUT/'protocol.json'),
        rule='Keep every FITOPT004 fitted half-step amplitude A4, every fitted nuisance/error, all nominal distances and every other variant unchanged. Construct MASS_CORR4_cf=+A4 if HOST_LOGMASS>10.44 else -A4, retaining the author > and <= semantics. Set DLMAG4_cf=DLMAG4_published+MASS_CORR4_cf-MASS_CORR4_published.',
        affected_rule='10<HOST_LOGMASS<=10.44; verify published nominal sign uses>10. No probabilistic host-mass reweighting or new step fit.',
        estimand='Same frozen high37-minus-low42 linear contrast for each branch variant and paired optical-minus-NIR; all outputs are deterministic conditional sensitivities.',
        covariance='Reproduce massdivide (=release massstep) from centered FITOPT004/005 signed response outer products with unit weights. Center each response by inverse-variance mean using that variant DLMAGERR. zHD unchanged. Only after closure, replace old004 outer product by cf004 outer product in supplied massstep and all marginal exports, retaining all other contributions.',
        closure_gate='Offdiagonal differences must be within the %.5e export rounding half-bin plus1e-12; reconstructed diagonal differences within(dmb_group+dmb_stat)*1e-6+1e-12 (the two %.6f error exports). Otherwise do not release a repaired covariance.',
        no_claim='No published source/input overwritten, no recomputation of nominal corrections, no new cosmology likelihood or dust/population inference. New covariance is conditional on retaining all other released contributions and on this source-linked propagation fix.'))
    print(json.dumps({'frozen_counterfactual_protocol_sha256':sha(OUT/name)}))

def mass_counterfactual():
    """Deterministic source-linked sensitivity; no optimization or input writes."""
    proto=OUT/'mass-threshold/counterfactual-protocol.json';q=json.loads(proto.read_text())
    assert q['source_sha256']==sha(OUT/'mass-threshold/code/raisin_cosmo/cosmo_sys.py')
    members=pd.read_csv(OUT/'frozen-membership.csv').set_index('CID');high=(members.stratum=='high').to_numpy();low=~high
    a=high/high.sum()-low/low.sum();H=np.eye(len(a))-np.ones((len(a),len(a)))/len(a)
    summary={};ledger=[];arrays={};cosm=FlatLambdaCDM(H0=70,Om0=.3)
    def closure(pred,br,group):
        export=syscov(br,group);off=offcov(br,group);mask=~np.eye(len(a),dtype=bool)
        # Exact literal print-precision gate; does not use a fitted tolerance.
        tol=np.zeros_like(off)+1e-12
        nz=off!=0;tol[nz]+=.5*10.**(np.floor(np.log10(abs(off[nz])))-5)
        dmb=lcparams(br,group)[:,5];stat=lcparams(br,'stat')[:,5]
        dtol=(dmb+stat)*1e-6+1e-12
        de=np.abs(np.diag(pred-export));oe=np.abs((pred-export)[mask])
        return dict(pass_gate=bool(np.all(oe<=tol[mask])and np.all(de<=dtol)),
            offdiag_max_abs_mag2=float(oe.max()),diagonal_max_abs_mag2=float(de.max()),
            max_offdiag_error_over_rounding_bound=float(np.max(oe/tol[mask])),
            max_diagonal_error_over_rounding_bound=float(np.max(de/dtol)),
            relative_frobenius_error=float(np.linalg.norm(pred-export)/np.linalg.norm(export)))
    for br in BRANCHES:
        d=table(br).loc[members.index];s=table(br,4).loc[members.index];s5=table(br,5).loc[members.index]
        mass=d.HOST_LOGMASS.to_numpy();A0=float(abs(d.MASS_CORR).max());A4=float(abs(s.MASS_CORR).max())
        nominal_sign=np.where(mass>10,1.,-1.);new_sign=np.where(mass>10.44,1.,-1.)
        assert np.max(abs(d.MASS_CORR-A0*nominal_sign))<1e-12
        assert np.max(abs(s.MASS_CORR-A4*nominal_sign))<1e-12
        affected=(mass>10)&(mass<=10.44);assert affected.sum()==16
        expected_mass=A4*new_sign;cf=s.DLMAG.to_numpy()+expected_mass-s.MASS_CORR.to_numpy()
        delta4=(s.DLMAG-d.DLMAG).to_numpy();delta5=(s5.DLMAG-d.DLMAG).to_numpy();delta_cf=cf-d.DLMAG.to_numpy()
        assert np.max(abs(s.zHD-d.zHD))==0 and np.max(abs(s5.zHD-d.zHD))==0
        centered4=delta4-np.average(delta4,weights=1/s.DLMAGERR.to_numpy()**2)
        centered5=delta5-np.average(delta5,weights=1/s5.DLMAGERR.to_numpy()**2)
        centered_cf=delta_cf-np.average(delta_cf,weights=1/s.DLMAGERR.to_numpy()**2)
        C4=np.outer(centered4,centered4);C5=np.outer(centered5,centered5);Ccf=np.outer(centered_cf,centered_cf)
        mass_gate=closure(C4+C5,br,'massstep')
        full_checks={};base_mu=cosm.distmod(d.zHD.to_numpy()).value;all_pred=np.zeros_like(C4)
        for k in range(1,30):
            if k<=28:v=table(br,k).loc[members.index]
            else:v=fitres_file(OUT/f'mass-threshold/code/output/cosmo_fitres_{br}/RAISIN_combined_FITOPT029_new.FITRES').loc[members.index]
            diff=(v.DLMAG-d.DLMAG).to_numpy();global_off=np.average(diff,weights=1/v.DLMAGERR.to_numpy()**2)
            response=diff-global_off-(cosm.distmod(v.zHD.to_numpy()).value-base_mu)
            all_pred+=np.outer(response,response)
            if k in [28,29]:full_checks[str(k)]=closure(all_pred,br,'all')
        correction=Ccf-C4
        arrays[f'{br}_nominal_distance']=d.DLMAG.to_numpy();arrays[f'{br}_fitopt004_published']=s.DLMAG.to_numpy()
        arrays[f'{br}_fitopt004_counterfactual']=cf;arrays[f'{br}_delta4_published']=delta4;arrays[f'{br}_delta4_counterfactual']=delta_cf
        arrays[f'{br}_centered4_published']=centered4;arrays[f'{br}_centered4_counterfactual']=centered_cf
        arrays[f'{br}_covariance_contribution_change']=correction
        arrays[f'{br}_massstep_source_reconstruction']=C4+C5
        if mass_gate['pass_gate']:
            arrays[f'{br}_massstep_published']=syscov(br,'massstep')
            arrays[f'{br}_massstep_counterfactual']=syscov(br,'massstep')+correction
        full_pass=mass_gate['pass_gate']and any(x['pass_gate']for x in full_checks.values())
        if full_pass:
            arrays[f'{br}_all_sys_published']=syscov(br,'all')
            arrays[f'{br}_all_sys_counterfactual']=syscov(br,'all')+correction
        summary[br]=dict(fitted_full_step_nominal=2*A0,fitted_full_step_variant004=2*A4,
            affected_high=int((affected&high).sum()),affected_low=int((affected&low).sum()),
            actual004_highlow_shift_mag=float(a@delta4),counterfactual004_highlow_shift_mag=float(a@delta_cf),
            counterfactual_minus_actual_highlow_mag=float(a@(delta_cf-delta4)),
            mass_covariance_gate=mass_gate,all_covariance_reconstruction=full_checks,
            all_counterfactual_covariance_released=full_pass,
            old_mass_group_contrast_sd_mag=float(np.sqrt(a@syscov(br,'massstep')@a)),
            counterfactual_mass_group_contrast_sd_mag=float(np.sqrt(a@(syscov(br,'massstep')+correction)@a))if mass_gate['pass_gate']else None,
            old_mass_response_rank=int(np.linalg.matrix_rank(np.stack([H@delta4,H@delta5]),tol=1e-10)),
            new_mass_response_rank=int(np.linalg.matrix_rank(np.stack([H@delta_cf,H@delta5]),tol=1e-10)))
        for i,cid in enumerate(members.index):
            if affected[i]:ledger.append(dict(branch=br,CID=cid,survey=members.loc[cid,'survey'],stratum=members.loc[cid,'stratum'],HOST_LOGMASS=mass[i],
                nominal_mass_corr=float(d.MASS_CORR.iloc[i]),fitopt004_actual_mass_corr=float(s.MASS_CORR.iloc[i]),fitopt004_declared_threshold_mass_corr=float(expected_mass[i]),
                actual_delta4_mag=float(delta4[i]),counterfactual_delta4_mag=float(delta_cf[i]),counterfactual_minus_actual_mag=float(delta_cf[i]-delta4[i])))
    arrays['CID']=members.index.to_numpy(str);arrays['a_high_minus_low']=a
    np.savez_compressed(OUT/'mass-threshold/counterfactual-arrays.npz',**arrays)
    pd.DataFrame(ledger).to_csv(OUT/'mass-threshold/affected16-by-branch.csv',index=False)
    pair={}
    for br in ['optical','opticalnir']:
        pair[br]=dict(actual_paired_FITOPT004_change_mag=summary[br]['actual004_highlow_shift_mag']-summary['nir']['actual004_highlow_shift_mag'],
            counterfactual_paired_FITOPT004_change_mag=summary[br]['counterfactual004_highlow_shift_mag']-summary['nir']['counterfactual004_highlow_shift_mag'],
            counterfactual_minus_actual_mag=summary[br]['counterfactual_minus_actual_highlow_mag']-summary['nir']['counterfactual_minus_actual_highlow_mag'])
    write('mass-threshold/counterfactual-result.json',dict(protocol_sha256=sha(proto),branches=summary,paired_sensitivity=pair,
        note='Frozen amplitude/parameters only; no raw refit. Full all-matrix replacement withheld unless exact code-formula closure passes the frozen print-precision gate.'))
    print(json.dumps({'branches':summary,'paired_sensitivity':pair}))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('command',choices=['freeze','audit','contrast','freeze-mass','mass-counterfactual']);a=p.parse_args()
    {'freeze':freeze,'audit':audit,'contrast':contrast,'freeze-mass':freeze_mass_counterfactual,'mass-counterfactual':mass_counterfactual}[a.command]()
