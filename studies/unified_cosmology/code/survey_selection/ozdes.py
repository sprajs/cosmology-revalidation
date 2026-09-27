#!/usr/bin/env python3
"""Positional host association to OzDES observed targets, retaining failed redshifts.

Frozen gates: <=1 arcsec, no second target within1arcsec; q=3/4 and
|z_host-z_target|<=.003 for secure spectroscopy. These verify association,
not redshift independence (many DES redshifts originate in OzDES).
"""
import gzip,json,re
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from astropy.io import fits
from common import ROOT,WORK,RESULTS,sha
from acquire import get
BASE='https://cdsarc.cds.unistra.fr/ftp/J/MNRAS/496/19/'

def xyz(ra,dec):
    r=np.deg2rad(ra);d=np.deg2rad(dec)
    return np.column_stack([np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)])

def head(path,cols):
    with fits.open(path) as h:
        return pd.DataFrame({c:np.char.strip(np.array(h[1].data[c]).astype(str)) if c=='SNID' else np.array(h[1].data[c]).astype(float) for c in cols})

def main():
    directory=WORK/'followup';cp=directory/'ozdesdr2.dat.gz';rp=directory/'ReadMe'
    get(BASE+'sp/',directory/'sp-index.html')
    spectra=set(re.findall(r'href="([^"/]+\.fits)"',(directory/'sp-index.html').read_text()))
    rows=[]
    with gzip.open(cp,'rt') as f:
        for line in f:
            c=line.rstrip('\n').split('|');assert len(c) in [9,10],c
            if len(c)==9:c.append('')
            rows.append([v.strip() for v in c[:10]])
    cat=pd.DataFrame(rows,columns=['OzDES_ID','RA','DEC','rmag','type','z','q_z','transient_type','comment','both'])
    for c in ['RA','DEC','rmag','z','q_z']:cat[c]=pd.to_numeric(cat[c],errors='coerce')
    assert len(cat)==38624 and cat.OzDES_ID.is_unique
    cat['host_target']=cat.type.str.contains('SN_host|SN_free_host',regex=True)
    cat['secure_z']=cat.q_z.isin([3,4]) & cat.z.gt(0)
    cat['spectrum_available']=cat.OzDES_ID.add('.fits').isin(spectra)
    cat['spectrum_url']=BASE+'sp/'+cat.OzDES_ID+'.fits'
    cat.to_csv(directory/'ozdes-normalized.csv',index=False)
    hp=WORK/'release/0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_HEAD.FITS.gz'
    h=head(hp,['SNID','RA','DEC','HOSTGAL_RA','HOSTGAL_DEC','HOSTGAL_SPECZ','HOSTGAL_SPECZ_ERR','HOSTGAL_SNSEP','HOSTGAL_DDLR','HOSTGAL_CONFUSION','HOSTGAL_OBJID','HOSTGAL_LOGMASS','HOSTGAL_LOGMASS_ERR','SNTYPE'])
    ledger=pd.read_csv(WORK/'normalized/dovekie-ledger.csv',dtype={'CID':str})
    h['in_Dovekie']=h.SNID.isin(ledger.loc[ledger.survey_family.eq('DES'),'CID'])
    assert h.in_Dovekie.sum()==1623
    valid=h.HOSTGAL_RA.between(0,360) & h.HOSTGAL_DEC.between(-90,90)
    h['valid_host_position']=valid
    tree=cKDTree(xyz(cat.RA.to_numpy(),cat.DEC.to_numpy()))
    dist,ind=tree.query(xyz(h.loc[valid,'HOSTGAL_RA'].to_numpy(),h.loc[valid,'HOSTGAL_DEC'].to_numpy()),k=2)
    a=2*np.arcsin(np.clip(dist/2,0,1))*180/np.pi*3600
    h['nearest_arcsec']=np.nan;h['second_arcsec']=np.nan
    h.loc[valid,'nearest_arcsec']=a[:,0];h.loc[valid,'second_arcsec']=a[:,1]
    matched=cat.iloc[ind[:,0]].copy();matched.index=h.index[valid]
    for c in cat.columns:h.loc[valid,'oz_'+c]=matched[c]
    h['unique_within_1arcsec']=h.nearest_arcsec.le(1) & h.second_arcsec.gt(1)
    h['redshift_delta']=h.HOSTGAL_SPECZ-h.oz_z
    h['redshift_agreement']=h.HOSTGAL_SPECZ.gt(0) & h.redshift_delta.abs().le(.003)
    h['clean_host_spectrum']=h.unique_within_1arcsec & h.oz_host_target.eq(True) & h.oz_secure_z.eq(True) & h.redshift_agreement & h.oz_spectrum_available.eq(True)
    out=directory/'des-ozdes-crosswalk.csv';h.to_csv(out,index=False)
    clean=h[h.clean_host_spectrum & h.in_Dovekie].copy()
    clean.to_csv(directory/'dovekie-ozdes-clean.csv',index=False)
    counts={}
    for name,mask in [('all_SMP',np.ones(len(h),dtype=bool)),('Dovekie',h.in_Dovekie),('outside_Dovekie',~h.in_Dovekie)]:
        d=h.loc[mask]
        counts[name]={'rows':len(d),'valid_host_position':int(d.valid_host_position.sum()),'nearest_within_1arcsec':int(d.nearest_arcsec.le(1).sum()),'unique_within_1arcsec':int(d.unique_within_1arcsec.sum()),'unique_host_targets':int((d.unique_within_1arcsec & d.oz_host_target.eq(True)).sum()),'clean_host_spectrum':int(d.clean_host_spectrum.sum()),'clean_unique_targets':int(d.loc[d.clean_host_spectrum,'oz_OzDES_ID'].nunique()),'ambiguous_two_targets_within_1arcsec':int(d.second_arcsec.le(1).sum()),'secure_positional_matches_redshift_disagree':int((d.unique_within_1arcsec & d.oz_secure_z.eq(True) & d.HOSTGAL_SPECZ.gt(0) & ~d.redshift_agreement).sum())}
    bins=[]
    targets=cat[cat.host_target]
    for lo,hi in zip([0,20,21,22,23,24,25],[20,21,22,23,24,25,40]):
        d=targets[targets.rmag.ge(lo)&targets.rmag.lt(hi)];n=len(d);k=int(d.secure_z.sum())
        # Wilson 95% interval; conditional on observation, not complete survey selection.
        z=1.95996398454;p=k/n if n else None
        if n:
            den=1+z*z/n;mid=(p+z*z/(2*n))/den;rad=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
        bins.append({'rmag_low':lo,'rmag_high':hi,'observed_targets':n,'secure_redshifts':k,'fraction':p,'wilson95':[mid-rad,mid+rad] if n else None})
    record={'code_sha256':sha(__file__),'status':'passed','association_gates':{'radius_arcsec':1,'second_neighbor_outside_radius':True,'secure_q_z':[3,4],'abs_redshift_delta_max':.003,'requires_host_target_tag':True},'catalogue_rows':len(cat),'spectra_listing_entries':len(spectra),'catalogue_quality_counts':cat.q_z.value_counts().to_dict(),'host_target_rows':len(targets),'host_target_quality_counts':targets.q_z.value_counts().to_dict(),'crossmatch':counts,'observed_host_target_success_by_aperture_rmag':bins,'clean_Dovekie_host_redshift_range':[float(clean.oz_z.min()),float(clean.oz_z.max())],'clean_Dovekie_z_above_0_5':int(clean.oz_z.gt(.5).sum()),'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [cp,rp,hp,directory/'sp-index.html',WORK/'normalized/dovekie-ledger.csv']},'outputs':{str(p.relative_to(ROOT)):sha(p) for p in [directory/'ozdes-normalized.csv',out,directory/'dovekie-ozdes-clean.csv']},'limits':['Crossmatch uses host positions, not transient positions; explicit host association uncertainty (DDLR/confusion) remains in ledger.','Matching redshifts are often the same OzDES measurements; agreement checks identity, not independent cosmological validation.','Secure-redshift fraction conditions on observed OzDES SN-host targets; it is not the targeting or detection efficiency.','r magnitude in OzDES is two-arcsec aperture and is not interchangeable with Kron magnitude in survey efficiency tables.','Spectra may have wavelength-dependent response, aperture and calibration limitations; preserve variance/mask and fit those nuisances.']}
    representative=clean.iloc[0]
    spectrum_path=directory/(representative.oz_OzDES_ID+'.fits')
    acquired=get(representative.oz_spectrum_url,spectrum_path)
    with fits.open(spectrum_path) as f:
        assert f[1].name=='VARIANCE' and f[2].name=='BADPIX'
        assert f[0].data.shape==f[1].data.shape==f[2].data.shape
        keep=f[2].data==0
        assert np.isfinite(f[0].data[keep]).all() and np.isfinite(f[1].data[keep]).all()
        assert np.all(f[1].data[keep]>0)
        record['representative_spectrum']={'acquisition':acquired,'OzDES_ID':representative.oz_OzDES_ID,'HDU_count':len(f),'pixels':len(f[0].data),'unmasked_pixels':int(keep.sum()),'flux_unit':f[0].header.get('BUNIT'),'variance_unit':f[1].header.get('BUNIT'),'wavelength_unit':f[0].header.get('CUNIT1'),'wavelength_WCS':{k:f[0].header[k] for k in ['CRPIX1','CRVAL1','CDELT1']},'same_stack_and_variance_and_mask_dimensions':True,'physical_flux_calibration_established':False}
    (RESULTS/'ozdes-crossmatch.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ['crossmatch','host_target_rows','clean_Dovekie_host_redshift_range','clean_Dovekie_z_above_0_5']},indent=2))
if __name__=='__main__':main()
