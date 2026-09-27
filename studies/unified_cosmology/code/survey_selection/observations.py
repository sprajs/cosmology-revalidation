#!/usr/bin/env python3
"""Check physical-event, flux-pointer and classification support in public DES data."""
import ast,json
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table
from common import ROOT,WORK,RESULTS,sha
from dovekie import load,marginalized_loglike

def strings(x):return np.char.strip(np.asarray(x).astype(str))
def audit_pair(head_path,phot_path):
    with fits.open(head_path) as hh, fits.open(phot_path) as ph:
        h=hh[1].data;p=ph[1].data
        ids=strings(h['SNID']);lo=h['PTROBS_MIN'].astype(int)-1;hi=h['PTROBS_MAX'].astype(int)
        n=h['NOBS'].astype(int)
        assert len(set(ids))==len(ids) and np.array_equal(hi-lo,n)
        assert lo.min()==0 and hi.max()<=len(p)
        assert np.array_equal(lo[1:],hi[:-1]+1)
        assert np.all(p['MJD'][hi[:-1]]==-777)
        event=np.repeat(np.arange(len(h)),n)
        usable=np.ones(len(p),bool);usable[hi[hi<len(p)]]=False
        assert usable.sum()==n.sum()
        f=np.asarray(p['FLUXCAL'])[usable];e=np.asarray(p['FLUXCALERR'])[usable]
        assert np.isfinite(f).all() and np.isfinite(e).all()
        bad=np.bincount(event,weights=e<=0,minlength=len(h)).astype(int)
        bands=strings(p['BAND'])[usable]
        per_bands=np.zeros(len(h),int)
        for b in np.unique(bands):per_bands+=np.bincount(event,weights=bands==b,minlength=len(h))>0
        r=pd.DataFrame({'CID':ids,'RA':np.asarray(h['RA']).astype(float),'DEC':np.asarray(h['DEC']).astype(float),'SNTYPE':h['SNTYPE'].astype(int),'NOBS':n,'bands':per_bands,'epochs_nonpositive_error':bad,'host_specz':np.asarray(h['HOSTGAL_SPECZ']).astype(float),'peakmjd':np.asarray(h['PEAKMJD']).astype(float)})
        result={'events':len(h),'phot_rows_including_separators':len(p),'science_phot_rows':int(n.sum()),'bands':{str(b):int(np.sum(bands==b)) for b in np.unique(bands)},'event_band_counts':r.bands.value_counts().to_dict(),'events_with_nonpositive_fluxerr':int((bad>0).sum()),'epochs_nonpositive_fluxerr':int((e<=0).sum()),'negative_flux_epochs':int((f<0).sum()),'host_specz_positive':int((r.host_specz>0).sum()),'spectroscopic_type_counts':r.SNTYPE.value_counts().to_dict(),'phot_columns':list(p.names),'head_nobs_pointer_separator_checks':True,'hdu_names':[(x.name,x.__class__.__name__) for x in ph]}
        if 'PHOTFLAG' in p.names:result['PHOTFLAG_counts']={str(k):int(v) for k,v in zip(*np.unique(p['PHOTFLAG'][usable],return_counts=True))}
    return r,result

def main():
    base=WORK/'release/0_DATA/DES-SN5YR_DES';follow=WORK/'followup'
    paths=[base/'DES-SN5YR_DES_HEAD.FITS.gz',base/'DES-SN5YR_DES_PHOT.FITS.gz',follow/'DES-SN5YR_DIFFIMG_HEAD.FITS.gz',follow/'DES-SN5YR_DIFFIMG_PHOT.FITS.gz']
    smp,smpcheck=audit_pair(*paths[:2]);diff,diffcheck=audit_pair(*paths[2:])
    cpath=next((WORK/'release/3_CLASSIFICATION').glob('*.csv'));classif=pd.read_csv(cpath,dtype={'SNID':str}).rename(columns={'SNID':'CID'});assert classif.CID.is_unique
    pc=[c for c in classif if c!='CID'];assert len(pc)==1;classif=classif.rename(columns={pc[0]:'pIa'})
    assert classif.pIa.between(0,1).all()
    ledger=pd.read_csv(WORK/'normalized/dovekie-ledger.csv',dtype={'CID':str});selected=set(ledger[ledger.survey_family.eq('DES')].CID)
    r=diff.merge(smp[['CID','RA','DEC','NOBS','bands']],on='CID',how='outer',suffixes=('_DIFFIMG','_SMP'),validate='one_to_one',indicator=True)
    assert (r._merge=='right_only').sum()==0
    common=r._merge.eq('both');dra=(r.loc[common,'RA_DIFFIMG']-r.loc[common,'RA_SMP'])*np.cos(np.deg2rad(r.loc[common,'DEC_SMP']));ddec=r.loc[common,'DEC_DIFFIMG']-r.loc[common,'DEC_SMP']
    posdelta=np.hypot(dra,ddec)*3600;assert posdelta.max()<.1
    r=r.merge(classif,on='CID',how='left',validate='one_to_one');r['in_Dovekie']=r.CID.isin(selected);r['in_SMP']=common.to_numpy();assert r.in_Dovekie.sum()==1623
    r.to_csv(follow/'candidate-ledger.csv',index=False)
    spec=r.SNTYPE.eq(1);classified=r.pIa.notna();cosmo=r.in_Dovekie
    classifier={'rows':len(classif),'all_IDs_in_SMP':set(classif.CID)<=set(smp.CID),'accepted_Dovekie_missing_probability':int((cosmo & ~classified).sum()),'accepted_Dovekie_pIa_below_0_5':int((cosmo&r.pIa.lt(.5)).sum()),'classified_spec_normal_Ia':int((spec&classified).sum()),'classified_spec_normal_Ia_pIa_above_0_5':int((spec&classified&r.pIa.gt(.5)).sum()),'classified_other_spec_objects':int((~spec&r.SNTYPE.ne(0)&classified).sum()),'classification_is_not_unbiased_population_purity':True}
    # Extract and execute only the mathematical likelihood function from author source.
    lp=WORK/'release/4_DISTANCES_COVMAT/DES-Dovekie-SN_Likelihood.py';tree=ast.parse(lp.read_text());func=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='cov_log_likelihood');namespace={'np':np};exec(compile(ast.Module(body=[func],type_ignores=[]),str(lp),'exec'),namespace)
    d=load();model=d['MU']+np.linspace(-.12,.16,len(d['MU']));one=np.ones(len(model));normalization=.5*np.log((one@d['precision']@one)/(2*np.pi));delta=abs(namespace['cov_log_likelihood'](model,d['MU'],d['precision'])-(marginalized_loglike(model,d)-normalization));assert delta<1e-8
    try:
        bad=Table.read(WORK/'release/4_DISTANCES_COVMAT/DES-Dovekie_HD.csv',format='ascii.csv');bad['zHD'];parser_error=None
    except Exception as ex:parser_error=type(ex).__name__+': '+str(ex)
    assert parser_error is not None
    result={'status':'passed','code_sha256':sha(__file__),'SMP':smpcheck,'DIFFIMG':diffcheck,'physical_ID_join':{'SMP_without_DIFFIMG':int((r._merge=='right_only').sum()),'detected_without_SMP':int((r._merge=='left_only').sum()),'shared_max_coordinate_delta_arcsec':float(posdelta.max()),'selected_in_detected':int(cosmo.sum())},'classifier':classifier,'author_likelihood_algebra_difference':delta,'author_csv_parser_failure':parser_error,'author_default_filenames_exist':[ (lp.parent/n).exists() for n in ['DES-SN5YR_HD.csv','STAT+SYS.txt.gz'] ],'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths+[cpath,lp]},'outputs':{str((follow/'candidate-ledger.csv').relative_to(ROOT)):sha(follow/'candidate-ledger.csv')},'limits':['DIFFIMG is a detected-candidate denominator; authors explicitly disallow its flux for cosmology.','SMP contains calibrated flux and individual errors; no cross-epoch covariance HDU is released here. Calibration cross-object covariance exists in corrected-distance products and cannot automatically be transplanted into rawflux likelihood.','Known spectroscopic types form a selected follow-up sample; classifier acceptance there is not unbiased purity or CC contamination in the untyped sample.','The exact multi-season candidate veto and every original selection label are not reconstructed by identifier membership alone.']}
    (RESULTS/'observed-survey-checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['physical_ID_join','classifier','author_likelihood_algebra_difference','author_csv_parser_failure']},indent=2))
if __name__=='__main__':main()
