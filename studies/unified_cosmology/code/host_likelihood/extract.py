"""Read public HDF5 arrays without deserializing author pickles.

Physical quantities use the same native draw index. Author postprocessed NPY
arrays independently resample properties and therefore cannot preserve their
joint covariance. All ages here retain the author's WMAP9 clock and SED prior.
"""
from __future__ import annotations
import ast,io,json,zipfile
from datetime import datetime,timezone
from pathlib import Path
import h5py,numpy as np,pandas as pd
from astropy.cosmology import WMAP9
from scipy.interpolate import PchipInterpolator
from scipy.spatial import cKDTree
from acquire import ROOT,WORK,OUT,sha

def vectors(ra,dec):
    r,d=np.deg2rad(ra),np.deg2rad(dec)
    return np.column_stack([np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)])

def headers(path):
    result={}
    for line in path.read_text().splitlines():
        if line.startswith('OBS:'):break
        if ':' in line and not line.startswith('#'):
            key,value=line.split(':',1);result[key]=value.strip()
    return result

def age_and_sfh(chain,labels,tuniv):
    n=len(chain);tuniv=np.broadcast_to(tuniv,(n,))
    # Exact independent algebra for constant-SFR bins and their formed masses.
    logedges=np.column_stack([np.zeros(n),np.full(n,7.4772),
      np.stack([8+(np.log10(.9*tuniv*1e9)-8)*f for f in np.linspace(0,1,5)],axis=1),np.log10(tuniv*1e9)])
    edges=10**logedges/1e9;dt=np.diff(edges,axis=1)
    ratios=chain[:,[labels.index('logsfr_ratios_'+str(i)) for i in range(1,7)]]
    relative=np.column_stack([np.ones(n),10**(-np.cumsum(ratios,axis=1))])
    fraction=dt*relative;fraction/=fraction.sum(axis=1)[:,None]
    age=np.sum(fraction*(edges[:,:-1]+edges[:,1:])/2,axis=1)
    assert np.allclose(fraction.sum(axis=1),1) and np.all(age>0) and np.all(age<tuniv)
    return age,fraction

def main():
    table=WORK/'frankenblast_spectroscopic_SN_sample.csv'
    hosts_path=WORK/'frankenblast_SN_sample_hosts.csv'
    sample=pd.read_csv(table).set_index('Name');hosts=pd.read_csv(hosts_path).set_index('Name')
    assert sample.index.is_unique and hosts.index.is_unique
    archive=WORK/'spectroscopic_SN_host_sbipp_output.zip'
    zgrid=np.linspace(0,1.6,3201);clock=PchipInterpolator(zgrid,WMAP9.age(zgrid).value)
    rng=np.random.default_rng(927081)
    zz=rng.uniform(0,1.5,100)
    clock_error=float(np.max(np.abs(clock(zz)-WMAP9.age(zz).value)))
    assert clock_error<1e-6
    chain_rows=[];phot_rows=[];selected_files={};invalid=[]
    yse={p.name.split('.snana.dat')[0]:(p,headers(p)) for p in (WORK/'yse/yse_dr1_zenodo').glob('*.snana.dat')}
    export=WORK/'joint-host-draws';export.mkdir(exist_ok=True)
    with zipfile.ZipFile(archive) as bundle:
        entries={Path(name).parent.name:name for name in bundle.namelist() if name.endswith('.h5') and '__MACOSX' not in name}
        for number,(name,entry) in enumerate(sorted(entries.items())):
            payload=bundle.read(entry)
            with h5py.File(io.BytesIO(payload)) as handle:
                if 'sampling/chain' not in handle:
                    invalid.append({'name':name,'reason':'native HDF5 has no posterior chain'});continue
                chain=handle['sampling/chain'][:].astype(float)
                labels=json.loads(handle['sampling'].attrs['theta_labels'])
                obs=handle['obs'];filters=json.loads(obs.attrs['filters'])
                redshift=json.loads(obs.attrs['redshift'])
                flux,error,mask=[obs[key][:] for key in ['maggies','maggies_unc','phot_mask']]
            assert chain.shape[1]==len(labels) and np.all(np.isfinite(chain))
            z=chain[:,labels.index('zred')]
            if np.any((z<0)|(z>1.6)):
                invalid.append({'name':name,'reason':'posterior redshift outside declared clock grid'});continue
            age,fraction=age_and_sfh(chain,labels,float(WMAP9.age(redshift).value) if redshift is not None else clock(z))
            age_drawz,_=age_and_sfh(chain,labels,clock(z))
            dust2=chain[:,labels.index('dust2')];dust1=chain[:,labels.index('dust1_fraction')]
            physical=np.column_stack([age,1.086*dust2,1.086*dust2*(1+dust1),chain[:,labels.index('logmass')],chain[:,labels.index('logzsol')],z])
            cov=np.cov(physical,rowvar=False,ddof=1)
            meta=sample.loc[name] if name in sample.index else None
            host=hosts.loc[name] if name in hosts.index else None
            row=dict(name=name,archive_member=entry,draws=len(chain),parameters=len(labels),
                     spec_class=None if meta is None else meta['Spec. Class'],
                     sn_ra=None if meta is None else float(meta.RA),sn_dec=None if meta is None else float(meta.Dec),
                     source_z=None if meta is None else float(meta.Redshift),obs_z=redshift,
                     host_ra=None if host is None else float(host.RA),host_dec=None if host is None else float(host.Dec),
                     host_probability=None if host is None else float(host.host_prob),
                     age_mean=float(age.mean()),age_sd=float(age.std(ddof=1)),age_median=float(np.median(age)),
                     Av_diffuse_mean=float(physical[:,1].mean()),Av_young_mean=float(physical[:,2].mean()),
                     logmass_formed_mean=float(physical[:,3].mean()),logZ_mean=float(physical[:,4].mean()),
                     corr_age_Av_diffuse=float(cov[0,1]/np.sqrt(cov[0,0]*cov[1,1])),
                     corr_age_logmass=float(cov[0,3]/np.sqrt(cov[0,0]*cov[3,3])),
                     posterior_z_sd=float(z.std(ddof=1)),drawz_minus_fixedz_age_mean=float((age_drawz-age).mean()),
                     band_count=len(filters),valid_bands=int((mask&np.isfinite(flux)&np.isfinite(error)&(error>0)).sum()))
            chain_rows.append(row)
            assert len(filters)==len(flux)==len(error)==len(mask)
            for band,f,e,keep in zip(filters,flux,error,mask):
                phot_rows.append(dict(name=name,filter=band,flux_maggies=float(f),error_maggies=float(e),
                                     native_mask=bool(keep),valid=bool(keep and np.isfinite(f) and np.isfinite(e) and e>0),
                                     MW_state='author_F99_Rv3.1_corrected',uncertainty_state='reported_with_1percent_floor'))
            if name in yse:
                path=export/(name+'.npz')
                np.savez_compressed(path,raw_chain=chain,theta_labels=np.asarray(labels),physical_draws=physical,
                  physical_labels=np.asarray(['mass_weighted_age_Gyr','Av_diffuse_mag','Av_young_mag','log10_mass_formed_Msun','log10_Z_Zsun','redshift']),
                  covariance=cov,mass_fraction=fraction,age_at_draw_redshift_Gyr=age_drawz,
                  flux_maggies=flux,error_maggies=error,phot_mask=mask,filters=np.asarray(filters))
                selected_files[str(path.relative_to(ROOT))]=sha(path)
            if number%1000==0:print('hosts',number,flush=True)
    summary=pd.DataFrame(chain_rows);summary.to_csv(WORK/'host-summary.csv',index=False)
    phot=pd.DataFrame(phot_rows);phot.to_csv(WORK/'host-photometry.csv',index=False)
    normal=summary.spec_class.isin(['SN Ia','SNIa-norm','SN Ia-norm'])
    cross=[]
    for name,(path,meta) in sorted(yse.items()):
        if name not in set(summary.name):continue
        row=summary.set_index('name').loc[name]
        ysera=float(meta['RA'].split()[0]);ysedec=float(meta['DECL'].split()[0])
        sep=float(np.linalg.norm(vectors([ysera],[ysedec])[0]-vectors([row.sn_ra],[row.sn_dec])[0])*206264.806)
        z=float(meta['REDSHIFT_FINAL'].split()[0]);hostz=float(meta['VETTED_HOST_GALAXY_REDSHIFT'].split()[0])
        hostra,hostdec=ast.literal_eval(meta['VETTED_HOST_GALAXY_COORDS(RA, DECL)[deg]'])
        hostsep=float(np.linalg.norm(vectors([hostra],[hostdec])[0]-vectors([row.host_ra],[row.host_dec])[0])*206264.806) if hostra>=0 and abs(hostdec)<=90 else np.nan
        data=dict(name=name,lightcurve=str(path.relative_to(ROOT)),lightcurve_sha256=sha(path),
                  host_draws=str((export/(name+'.npz')).relative_to(ROOT)),sn_separation_arcsec=sep,
                  spec_class=meta['SPEC_CLASS'],normal_Ia=meta['SPEC_CLASS'] in ['SN Ia','SN Ia-norm'],
                  z_hel=z,host_z_hel=hostz,redshift_source=meta['REDSHIFT_FINAL'],host_probability=row.host_probability,
                  mwebv=float(meta['MWEBV'].split()[0]),discovery_mjd=float(meta['SEARCH_PEAKMJD']),
                  host_status=meta['VETTED_HOST_GALAXY_STATUS'],source_z_difference=z-row.source_z)
        data['vetted_host_separation_arcsec']=hostsep
        data['host_centroids_agree']=bool(np.isfinite(hostsep) and hostsep<1)
        data['host_spec_z_available']=hostz>0
        data['primary_candidate']=bool(data['normal_Ia'] and sep<1 and abs(data['source_z_difference'])<.002 and row.host_probability>=.9 and data['host_centroids_agree'] and hostz>0 and .01<hostz<.25)
        cross.append(data)
    cross=pd.DataFrame(cross);cross.to_csv(WORK/'yse-host-crosswalk.csv',index=False)
    selected=summary[normal&(summary.host_probability>=.9)]
    def q(v):return dict(zip(['p05','median','p95'],map(float,np.quantile(v,[.05,.5,.95]))))
    result=dict(completed_utc=datetime.now(timezone.utc).isoformat(),status='recovered_observed_photometry_and_joint_model_conditioned_posteriors',
      code_sha256=sha(__file__),inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in [table,hosts_path,archive]},
      source_sha256={str((WORK/'frankenblast-host'/n).relative_to(ROOT)):sha(WORK/'frankenblast-host'/n) for n in ['fit_host_sed.py','postprocess_sbi.py','sbi_pp.py']},
      outputs_sha256={str((WORK/n).relative_to(ROOT)):sha(WORK/n) for n in ['host-summary.csv','host-photometry.csv','yse-host-crosswalk.csv']},
      selected_draws_sha256=selected_files,hosts=len(summary),archive_hdf5_members=len(entries),invalid=invalid,
      normal_Ia_hosts=int(normal.sum()),normal_Ia_host_probability_ge09=len(selected),
      photometry_rows=len(phot),valid_fluxes=int(phot.valid.sum()),valid_negative_fluxes=int((phot.valid&(phot.flux_maggies<0)).sum()),
      YSE_exact_name_common=len(cross),YSE_normal_Ia=int(cross.normal_Ia.sum()),YSE_primary_candidates=int(cross.primary_candidate.sum()),
      selected_normal_Ia_summary={key:q(selected[key].dropna()) for key in ['source_z','age_mean','age_sd','corr_age_Av_diffuse','corr_age_logmass','posterior_z_sd']},
      clock_interpolation_max_error_Gyr=clock_error,
      semantics={'age':'formed-mass weighted global host age, constant SFR bins, WMAP9 cosmic-time ceiling; not progenitor delay',
                 'dust':'diffuse stellar attenuation and young-star total attenuation; neither is SN line-of-sight extinction',
                 'mass':'total formed mass; surviving stellar mass would require model-dependent mass loss',
                 'posterior':'SBI++ conditioned on empirical training mixture and photometric selection; not prior-free likelihood',
                 'joint_draws':'derived quantities preserve HDF5 row alignment; independent author NPY resampling is bypassed',
                 'photometry':'author foreground-corrected global aperture flux in maggies;1percent error floor included; no released interband covariance',
                 'YSE':'exact identifiers and SN positions checked; actual redshift uncertainty and epoch/sample selection remain additional requirements'})
    (OUT/'host-interface.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({key:result[key] for key in ['hosts','normal_Ia_hosts','YSE_exact_name_common','YSE_normal_Ia','YSE_primary_candidates','valid_negative_fluxes']}))

if __name__=='__main__':main()
