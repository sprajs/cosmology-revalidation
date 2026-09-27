#!/usr/bin/env python3
"""Corrected public DES global host photometry and release-history audit."""
from pathlib import Path
import sys,json,hashlib,datetime
import numpy as np,pandas as pd
from astropy.io import fits,ascii
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT))
from lib.records import fitres
W=ROOT/'.work/host-transport';O=ROOT/'studies/host_ages/results/host_transport'
def sha(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    p=W/'des/DES-SN5YR_DES_HEAD.FITS.gz';mp=W/'des/DES-Dovekie_Metadata.csv';dp=W/'des-deep-crosswalk.csv'
    with fits.open(p) as f:
        a=f[1].data;cols=['SNID','NAME_IAUC','NAME_TRANSIENT','HOSTGAL_RA','HOSTGAL_DEC','HOSTGAL_OBJID','HOSTGAL_DDLR','HOSTGAL_NMATCH','HOSTGAL_NMATCH2','REDSHIFT_HELIO','HOSTGAL_SPECZ','MWEBV','MWEBV_ERR']+[f'HOSTGAL_{v}_{b}' for v in ['MAG','MAGERR'] for b in 'griz']
        d=pd.DataFrame({k:np.char.strip(np.array(a[k]).astype(str)) if a[k].dtype.kind in 'US' else np.array(a[k]).astype(float) for k in cols})
    m=fitres(mp);selected=m[m.IDSURVEY==10];assert set(selected.index).issubset(set(d.SNID));d['selected_Dovekie']=d.SNID.isin(selected.index)
    d['redshift']=d.SNID.map(m.zHEL).fillna(d.REDSHIFT_HELIO)
    alias_source=W/'des/DES-SN5YR_DIFFIMG_HEAD.FITS.gz'
    with fits.open(alias_source) as f:
        a=f[1].data
        aliases=dict(zip(np.char.strip(np.array(a.SNID).astype(str)),np.char.strip(np.array(a.IAUC).astype(str))))
    d['legacy_alias']=d.SNID.map(aliases)
    d.to_csv(W/'des-smp-hosts.csv',index=False);rows=[]
    for _,r in d.iterrows():
        for b in 'griz':
            mag=float(r['HOSTGAL_MAG_'+b]);e=float(r['HOSTGAL_MAGERR_'+b]);valid=np.isfinite([mag,e]).all() and 0<mag<50 and 0<e<10
            flux=10**(.4*(22.5-mag)) if valid else None
            rows.append(dict(sn_id=r.SNID,host_id=int(r.HOSTGAL_OBJID),host_ra_deg=r.HOSTGAL_RA,host_dec_deg=r.HOSTGAL_DEC,redshift=r.redshift,selected_Dovekie=bool(r.selected_Dovekie),host_dlr=r.HOSTGAL_DDLR,host_nmatch=r.HOSTGAL_NMATCH,aperture='global host catalogue',filter_family='DECam',band=b,native_magnitude=mag if valid else None,native_magnitude_error=e if valid else None,magnitude_system='AB',flux_ab_nanomaggies=flux,flux_error_ab_nanomaggies=flux*np.log(10)/2.5*e if valid else None,measurement_available=bool(valid),native_likelihood='Gaussian catalogue magnitude; flux error is first-order approximation',mw_ebv_sn=r.MWEBV,mw_ebv_sn_error=r.MWEBV_ERR,foreground_state='Released global host magnitude; exact host correction must be checked separately from SN foreground metadata',covariance_available=False))
    pd.DataFrame(rows).to_csv(W/'des-smp-photometry.csv',index=False)
    q=pd.read_csv(dp,dtype={'SNID':str});q=q.merge(d,on='SNID',suffixes=('_old','_smp'),validate='one_to_one');q=q[q.primary_match&q.selected_Dovekie_smp&q.deep_quality&q.host_dlr_lt4]
    comparisons=[]
    for b in 'griz':
        f=q[f'deep_BDF_FLUX_CALIB_{b.upper()}'];old=q[f'HOSTGAL_MAG_{b}_old'];new=q[f'HOSTGAL_MAG_{b}_smp'];ok=(f>0)&(old>0)&(old<50)&(new>0)&(new<50)
        deep=30-2.5*np.log10(f[ok]);delta=deep-new[ok]
        comparisons.append({'band':b,'n':int(ok.sum()),'deep_minus_old_mag_q025_q50_q975':np.quantile(deep-old[ok],[.025,.5,.975]).tolist(),'deep_minus_corrected_mag_q025_q50_q975':np.quantile(delta,[.025,.5,.975]).tolist(),'old_minus_corrected_mag_q025_q50_q975':np.quantile(old[ok]-new[ok],[.025,.5,.975]).tolist(),'median_corrected_magerr':float(q.loc[ok,f'HOSTGAL_MAGERR_{b}_smp'].median())})
    wp=W/'wiseman/tableb1.dat';wr=W/'wiseman/ReadMe';w=ascii.read(wp,readme=str(wr),format='cds').to_pandas();j=d.merge(w,left_on='legacy_alias',right_on='Name');paper=[]
    for b in 'griz':
        delta=j[f'HOSTGAL_MAG_{b}']-j[b+'mag'];ok=np.isfinite(delta)&(j[b+'mag']>0)&(j[b+'mag']<50)&(j[f'HOSTGAL_MAG_{b}']>0)&(j[f'HOSTGAL_MAG_{b}']<50)
        paper.append({'band':b,'matched':int(ok.sum()),'corrected_minus_Wiseman_q025_q50_q975':np.quantile(delta[ok],[.025,.5,.975]).tolist() if ok.any() else None})
    r={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'candidate_rows':len(d),'selected_DES_metadata_rows':len(selected),'corrected_host_mag_source':'Pinned SMP HEAD; its official README logs the December2025 correction. DIFFIMG HEAD and Dovekie metadata retain older discrepant HOST_MAG values; do not treat those as current host photometry. No SN distance effect inferred here.','deep_comparison':comparisons,'Wiseman2020_comparison':paper,'coordinate_check_max_change_deg':float(np.max(abs(q[['HOSTGAL_RA_old','HOSTGAL_DEC_old']].to_numpy()-q[['HOSTGAL_RA_smp','HOSTGAL_DEC_smp']].to_numpy()))),'photometry_scope':'No outcome-based magnitude cut; native magnitudes/errors retained, not original signed forced flux. Missing measurements are not interpreted as upper limits. Model/measurement aperture difference and interband covariance remain unresolved.','output_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [W/'des-smp-hosts.csv',W/'des-smp-photometry.csv']},'input_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [p,mp,dp,wp,wr,alias_source,W/'des/DES-SN5YR_DES.README']},'code_sha256':{str(Path(__file__).relative_to(ROOT)):sha(Path(__file__))}}
    (O/'des-smp-interface.json').write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');print(json.dumps(r,indent=2))
if __name__=='__main__':main()
