#!/usr/bin/env python3
"""Independent calibration, identity, weight and provenance checks."""
from pathlib import Path
import argparse,datetime,hashlib,importlib.metadata,json,py_compile,sys
import numpy as np,pandas as pd
from astropy.io import fits
ROOT=Path(__file__).resolve().parents[4];W=ROOT/'.work/host-transport';O=ROOT/'studies/host_ages/results/host_transport';H=Path(__file__).parent
sys.path.insert(0,str(ROOT));from lib.records import fitres

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--include-physical-review',action='store_true');args=parser.parse_args()
    pinned=json.loads((O/'inputs.json').read_text())['files'];checked={}
    for r in pinned:
        p=ROOT/r['path'];s=sha(p);assert s==r['sha256'] and p.stat().st_size==r['bytes'],p;checked[r['path']]=s
    md5=[]
    for field in ['C3','E2','X3']:
        d=json.loads((W/f'des/deep-schema-{field}.json').read_text())
        for name,r in d['file_info']['files'].items():
            p=W/'des'/name
            with p.open('rb') as stream:s=hashlib.file_digest(stream,'md5').hexdigest()
            assert s==r['md5sum'] and p.stat().st_size==r['size'];md5.append({'path':str(p.relative_to(ROOT)),'publisher_md5':s})
    records={}
    for p in O.glob('*.json'):
        if p.name in ['validation.json','inputs.json','availability.json']:continue
        if p.name=='physical-operator-review.json' and not args.include_physical_review:continue
        r=json.loads(p.read_text());records[p.name]=sha(p)
        for kind in ['input_sha256','code_sha256','output_sha256']:
            m=r.get(kind,{})
            if isinstance(m,dict):
                for path,s in m.items():assert sha(ROOT/path)==s,(p,path)
    if args.include_physical_review:
        review=json.loads((O/'physical-operator-review.json').read_text())
        assert review['passed'] and review['max_abs_colour_error_mag']<review['tolerance_mag']
        assert sha(H/'review_physical_operator.py')==review['review_code_sha256']
        assert sha(ROOT/'studies/host_ages/code/physical_ages/model.py')==review['model_sha256']
    r=json.loads((O/'roman-interface.json').read_text())
    for pathkey,hashkey in [('output_path','output_sha256'),('properties_path','properties_sha256')]:assert sha(ROOT/r[pathkey])==r[hashkey]
    rf=json.loads((O/'filter-interface.json').read_text());assert sha(H/'filters.py')==rf['code_sha256'];assert sha(ROOT/rf['npz_path'])==rf['npz_sha256']
    for b in rf['filters']:
        for kind in ['input','output']:assert sha(ROOT/b[kind+'_path'])==b[kind+'_sha256']
    s=json.loads((O/'des-selection.json').read_text());assert sha(W/'des-eightband-availability.csv')==s['output_sha256']
    p=pd.read_csv(W/'roman-photometry.csv');local=p[(p.aperture=='local')&p.measurement_available];zps={}
    for family,target in [('SDSS',0),('CFHT_MegaCam_SNLS',30)]:
        t=local[(local.filter_family==family)&(local.native_flux>0)&(local.native_magnitude<90)]
        if len(t):
            zp=t.native_magnitude+2.5*np.log10(t.native_flux);err=float(abs(zp-target).max());assert err<1e-5;zps[family]={'n':len(t),'target':target,'max_abs_error_mag':err}
    # Compare retained signed native fluxes directly to the FITS rows.
    with fits.open(W/'roman/snprop.fits.gz') as f:a=f[1].data
    direct_negative=sum(int(np.sum((a['flux_local_'+b]<0)&(a['eflux_local_'+b]>0))) for b in 'ugriz') if 'flux_local_u' in a.names else None
    assert int(local.negative_flux.sum())==direct_negative==48
    q=pd.read_csv(W/'des-deep-crosswalk.csv',dtype={'SNID':str});g=q[q.primary_match&q.selected_Dovekie&q.deep_quality&q.host_dlr_lt4];assert len(g)==331 and g.deep_ID.nunique()==331
    selected_path=W/'des-deep-selected-photometry.csv';sel=pd.read_csv(selected_path,float_precision='round_trip');src=pd.read_csv(W/'des-deep-crosswalk.csv',float_precision='round_trip');src=src[src.primary_match&src.selected_Dovekie&src.deep_quality&src.host_dlr_lt4]
    assert len(sel)==2648 and sel.measurement_available.all()
    for b in ['u','g','r','i','z','j','h','ks']:
        a=sel[sel.band==b];native=src['deep_BDF_FLUX_CALIB_'+b.upper()].to_numpy();assert not np.any(abs(native)==9.999e9)
        assert a.sn_id.astype(str).to_list()==src.SNID.astype(str).to_list();assert (a.flux_ab_nanomaggies.to_numpy()==native/1000).all()
    integrity={'proof':'Only revised availability predicate excludes abs(native value)==9.999e9. None of the331 rows per band meets that condition, so old/revised branches are identical on all2648 selected photometry rows. Native values and row order independently equal selected source crosswalk.','old_full_interface_sha256':'6d08bbc3a9a6502494667c75207f0bcf54c82beed393c7fa3d718366512c2dc9','selected_objects':331,'selected_rows':len(sel),'selected_csv_sha256':sha(selected_path),'selected_csv_path':str(selected_path.relative_to(ROOT))}
    (O/'selected-photometry-integrity.json').write_text(json.dumps(integrity,indent=2)+'\n')

    meta=fitres(W/'des/DES-Dovekie_Metadata.csv');sm=pd.read_csv(W/'des-smp-hosts.csv',dtype={'SNID':str});ids=set(meta.loc[meta.IDSURVEY==10].index.astype(str));assert len(ids)==1623 and ids<=set(sm.SNID)
    angles=np.deg2rad(q[['HOSTGAL_RA','HOSTGAL_DEC','deep_RA','deep_DEC']].to_numpy());ra,dec,ra2,dec2=angles.T
    sep=np.arctan2(np.hypot(np.cos(dec2)*np.sin(ra2-ra),np.cos(dec)*np.sin(dec2)-np.sin(dec)*np.cos(dec2)*np.cos(ra2-ra)),np.sin(dec)*np.sin(dec2)+np.cos(dec)*np.cos(dec2)*np.cos(ra2-ra))*206264.80624709636
    max_sep=float(abs(sep-q.match_sep_arcsec).max());assert max_sep<1e-7
    cal=[]
    for b,k in zip(['U','G','R','I','Z','J','H','KS'],[3.9627,3.186,2.140,1.569,1.196,.705,.441,.308]):
        flux=q['deep_BDF_FLUX_DERED_CALIB_'+b];mag=q['deep_BDF_MAG_DERED_CALIB_'+b];ok=(flux>0)&np.isfinite(mag)&(mag<90)
        err=float(abs(mag[ok]+2.5*np.log10(flux[ok])-30).max());assert err<1e-10
        native=q['deep_BDF_FLUX_CALIB_'+b];pred=native*10**(.4*k*q.deep_EBV_SFD98);available=np.isfinite(native)&np.isfinite(flux)&(abs(native)!=9.999e9)&(abs(flux)!=9.999e9);errfg=float(np.nanmax((abs(pred-flux)/np.maximum(abs(flux),1e-20))[available]));assert errfg<1e-10
        cal.append({'band':b,'mag_identity_error':err,'foreground_relative_error':errfg})
    raw=pd.read_csv(W/'roman-properties.csv');valid=np.isfinite(raw[['uv_local','uv_global','euv_local','euv_global','stellar_mass','redshift']]).all(axis=1)&(raw.euv_local>0)&(raw.euv_global>0)&(raw.uv_local<90)&(raw.uv_global<90)&(raw.stellar_mass>0)
    raw=raw[valid];tests=[]
    for source,target in [('SDSS','SNLS5'),('SNLS5','SDSS')]:
        d=pd.read_csv(W/f'roman-transport-{source}-to-{target}.csv');t=raw[(raw.survey==target)&(raw.redshift>=.2)&(raw.redshift<.4)];cols=['redshift','uv_global','stellar_mass'];w=d.transport_weight.to_numpy();assert abs(w.sum()-1)<1e-12 and (w>0).all()
        error=float(abs((w@d[cols].to_numpy()-t[cols].mean().to_numpy())/d[cols].std(ddof=0).to_numpy()).max());assert error<1e-6
        tests.append({'source':source,'target':target,'independent_standardized_balance_error':error,'ess':float(1/np.sum(w**2))})
    d=pd.read_csv(W/'des-eightband-availability.csv');fieldtests=[]
    for r in s['heldout_fields']:
        t=d[d.field==r['heldout_field']];x=t[['redshift','host_i','host_g_minus_i']].to_numpy();coef=np.array(r['coefficients']);eta=coef[0]+((x-r['training_mean'])/r['training_std'])@coef[1:];p=1/(1+np.exp(-eta));assert np.max(abs(p-t.heldout_availability_probability))<1e-12
        # Direct all-positive/all-negative pair comparison, independent of rank formula.
        y=t.available.to_numpy();diff=p[y,None]-p[~y][None,:];auc=float(np.mean((diff>0)+.5*(diff==0)));assert abs(auc-r['AUC'])<1e-12
        fieldtests.append({'field':r['heldout_field'],'pairwise_AUC':auc,'predicted_probability_max_error':float(np.max(abs(p-t.heldout_availability_probability)))})
    for p in H.glob('*.py'):py_compile.compile(str(p),doraise=True)
    result={'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'passed':True,'physical_operator_review_included':args.include_physical_review,'pinned_sources_verified':len(pinned),'published_parquet_md5':md5,'record_sha256':records,'roman_native_zero_points':zps,'roman_negative_measurements_retained':int(local.negative_flux.sum()),'DES_selected_IDSURVEY10':len(ids),'DES_quality_eightband_hosts':len(g),'atan2_sky_separation_max_error_arcsec':max_sep,'DES_photometric_identities':cal,'independent_entropy_weight_checks':tests,'independent_pairwise_auc_checks':fieldtests,'environment':{p:importlib.metadata.version(p) for p in ['numpy','scipy','pandas','astropy','pyarrow']},'code_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(H.glob('*')) if p.is_file()},'source_registry_sha256':sha(O/'inputs.json')}
    (O/'validation.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
