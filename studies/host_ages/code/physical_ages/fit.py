#!/usr/bin/env python3
"""Fit all eligible recovered fluxes; keep failed physical-model checks visible."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import argparse, json, hashlib, time
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from scipy.stats import chi2, norm
from sfdmap2 import sfdmap
from model import Grid, Library, age_bound, feasible_fit, fieller

ROOT=Path(__file__).resolve().parents[4]; WORK=ROOT/'.work/physical-ages'
OUT=ROOT/'studies/host_ages/results/physical_ages'
GRID=None

def initialize(sample='sdss'):
    global GRID
    GRID=Grid(Library(response_path=ROOT/'.work/host-transport/filters/des-deep-responses.npz')) if sample=='des' else Grid()

def summarize_interval(values):
    return [float(min(values[0])),float(max(values[1]))] if values[0] and values[1] else [None,None]

def fit_one(item):
    grid=GRID; z=item['z']; ages,matrices,_,_=grid.at_redshift(z,item['ebv'])
    flux=np.asarray(item['flux']); error=np.asarray(item['error'])
    normalizer=max(float(np.median(abs(flux))),1e-8); flux=flux/normalizer;error=error/normalizer
    dn=grid.observed_dn(z,item['ebv'])
    threshold=chi2.ppf(.95,len(flux)); joint_threshold=chi2.ppf(.975,len(flux))
    rows=[]
    for floor in [0.,.03]:
        sigma=np.hypot(error,.4*np.log(10)*floor*abs(flux))
        best=None; age_values=[[],[]]; dn_values=[[],[]]; joint_values=[[],[]]
        feasible=0; joint_feasible=0; numerical_failures=0
        for j,a in enumerate(matrices):
            coeff,chisq=feasible_fit(a,flux,sigma)
            if best is None or chisq<best['chi2']:
                best={'chi2':chisq,'nuisance':j,'age_Gyr':float(ages@coeff/coeff.sum()) if coeff.sum()>0 else None,
                      'Dn4000':float(dn[j,1]@coeff/(dn[j,0]@coeff)) if coeff.sum()>0 else None}
            if chisq <= threshold:
                feasible+=1
                pair=[age_bound(a,flux,sigma,ages,threshold,maximize) for maximize in [False,True]]
                if all(x is not None for x in pair):
                    for k in range(2):age_values[k].append(pair[k]['age'])
                else:numerical_failures+=1
                if item.get('dn_interval95') is not None:
                    # Coefficients here represent blue-band light, not mass.
                    pair=[age_bound(a/dn[j,0],flux,sigma,dn[j,1]/dn[j,0],threshold,maximize) for maximize in [False,True]]
                    if all(x is not None for x in pair):
                        for k in range(2):dn_values[k].append(pair[k]['age'])
                    else:numerical_failures+=1
            if chisq<=joint_threshold and item.get('dn_interval975') is not None:
                dlo,dhi=item['dn_interval975']
                inequalities=np.vstack([dn[j,1]-dhi*dn[j,0],dlo*dn[j,0]-dn[j,1]])
                diagnostics=[[],[]]
                pair=[age_bound(a,flux,sigma,ages,joint_threshold,maximize,inequalities,diagnostics[k]) for k,maximize in enumerate([False,True])]
                if all(x is not None for x in pair):
                    joint_feasible+=1
                    for k in range(2):joint_values[k].append(pair[k]['age'])
                else:
                    # A failed solver is not an infeasibility certificate,
                    # including when BOTH endpoints fail numerically.
                    proven_infeasible=all(d and d[-1].get('status')=='PrimalInfeasible' for d in diagnostics)
                    if not proven_infeasible:numerical_failures+=1
        alo,ahi=summarize_interval(age_values); dlo,dhi=summarize_interval(dn_values); jlo,jhi=summarize_interval(joint_values)
        if numerical_failures:
            # A missing extremum must never silently shrink the reported union.
            alo=ahi=dlo=dhi=jlo=jhi=None
        row={k:v for k,v in item.items() if k not in ['flux','error','dn_interval95','dn_interval975']}
        row.update(model_discrepancy_mag=floor,age_low_Gyr=alo,age_high_Gyr=ahi,
            joint_age_low_Gyr=jlo,joint_age_high_Gyr=jhi,photometry_Dn_low=dlo,photometry_Dn_high=dhi,
            nuisance_cells_feasible=feasible,joint_nuisance_cells_feasible=joint_feasible,numerical_failures=numerical_failures,
            minimum_chi2=best['chi2'],best_photometry_age_Gyr=best['age_Gyr'],best_photometry_Dn=best['Dn4000'],
            best_nuisance=best['nuisance'],cosmic_age_ceiling_Gyr=float(ages.max()))
        interval=item.get('dn_interval95')
        row['heldout_Dn_intervals_disjoint']=bool(dhi<interval[0] or dlo>interval[1]) if dlo is not None and interval else None
        rows.append(row)
    return rows

def inputs(sample):
    dust=sfdmap.SFDMap(str(WORK/'sfd'),scaling=.86)
    source_files=[]; items=[]; excluded={}
    if sample=='sdss':
        source=ROOT/'.work/galaxy-validation/host-spectrum-age-ledger.csv';source_files.append(source)
        table=pd.read_csv(source,dtype={k:str for k in ['specObjID','bestObjID','photometric_objID','photometric_flags']})
        table=table[table.specObjID.notna()].drop_duplicates('specObjID').sort_values('specObjID')
        for _,row in table.iterrows():
            flux=row[[f'fiberFlux_{b}' for b in 'ugriz']].to_numpy(float)
            ivar=row[[f'fiberFluxIvar_{b}' for b in 'ugriz']].to_numpy(float)
            if not (np.isfinite(flux).all() and np.isfinite(ivar).all() and (ivar>0).all()):
                excluded['invalid_five_band_flux_or_error']=excluded.get('invalid_five_band_flux_or_error',0)+1;continue
            ab=10**(-.4*np.array([-.04,0,0,0,.02]))
            ebv=float(dust.ebv(row.photometric_ra,row.photometric_dec))
            item=dict(id=row.specObjID,source_id=row.source_id,z=float(row.spectral_z),ebv=ebv,
                flux=(flux*ab).tolist(),error=(ivar**-.5*ab).tolist(),
                published_global_age_Gyr=float(row.sed_age) if np.isfinite(row.sed_age) else None,
                photometric_clean=bool(row.photometric_clean==1),SN_to_fibre_arcsec=float(row.SN_to_fibre_arcsec) if np.isfinite(row.SN_to_fibre_arcsec) else None,
                brightness_sample=bool(row.in_401_TITAN_ZTF_brightness_sample),raw_Dn4000=float(row.raw_Dn4000) if np.isfinite(row.raw_Dn4000) else None)
            if row.valid is True or row.valid==True:
                if min(row.blue_coverage_fraction,row.red_coverage_fraction)>=.995:
                    for level in [.95,.975]:
                        key='dn_interval'+('95' if level==.95 else '975')
                        item[key]=fieller(row.raw_red_Fnu_uJy,row.raw_blue_Fnu_uJy,row.raw_red_variance_uJy2,row.raw_blue_variance_uJy2,row.raw_blue_red_covariance_uJy2,norm.ppf((1+level)/2))
            items.append(item)
    elif sample=='des':
        source=ROOT/'.work/host-transport/des-deep-photometry.csv';source_files.append(source)
        source_files.append(ROOT/'.work/host-transport/filters/des-deep-responses.npz')
        table=pd.read_csv(source,dtype={'sn_id':str,'host_id':str})
        table=table[table.selected_Dovekie & table.deep_quality & table.host_dlr_lt4]
        for name,group in table.groupby('host_id',sort=True):
            # One host may contain multiple events; never duplicate its light.
            group=group.drop_duplicates('band').set_index('band').reindex(['u','g','r','i','z','j','h','ks'])
            flux=group.flux_dered_ab_nanomaggies.to_numpy(float);error=group.flux_error_dered_ab_nanomaggies.to_numpy(float)
            if not (np.isfinite(flux).all() and np.isfinite(error).all() and (error>0).all()):
                excluded['invalid_eight_band_flux_or_error']=excluded.get('invalid_eight_band_flux_or_error',0)+1;continue
            row=group.iloc[0];items.append(dict(id=name,source_id=row.sn_id,z=float(row.redshift),ebv=0.,flux=flux.tolist(),error=error.tolist(),foreground='Released DES dereddened AB flux; no second foreground correction'))
    elif sample=='roman':
        source=ROOT/'.work/host-transport/roman-photometry.csv';source_files.append(source)
        table=pd.read_csv(source);table=table[(table.aperture=='local')&(table.filter_family=='SDSS')]
        for name,group in table.groupby('sn_id',sort=True):
            group=group.set_index('band').reindex(list('ugriz'));flux=group.flux_ab_nanomaggies.to_numpy(float);error=group.flux_error_ab_nanomaggies.to_numpy(float)
            if not (np.isfinite(flux).all() and np.isfinite(error).all() and (error>0).all()):
                excluded['invalid_five_band_flux_or_error']=excluded.get('invalid_five_band_flux_or_error',0)+1;continue
            row=group.iloc[0];items.append(dict(id=name,z=float(row.redshift),ebv=float(dust.ebv(row.sn_ra_deg,row.sn_dec_deg)),flux=flux.tolist(),error=error.tolist()))
    return items,source_files,excluded

def quantiles(series):
    x=pd.to_numeric(series,errors='coerce').dropna()
    return dict(zip(['p05','median','p95'],map(float,np.quantile(x,[.05,.5,.95])))) if len(x) else None

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--sample',choices=['sdss','roman','des'],required=True);parser.add_argument('--workers',type=int,default=4);parser.add_argument('--limit',type=int)
    args=parser.parse_args();start=time.perf_counter();items,sources,excluded=inputs(args.sample)
    sources += [WORK/'ssp-library.npz', WORK/'sfd/SFD_dust_4096_ngp.fits', WORK/'sfd/SFD_dust_4096_sgp.fits']
    if args.limit:items=items[:args.limit]
    allrows=[]
    with ProcessPoolExecutor(args.workers,initializer=initialize,initargs=(args.sample,)) as pool:
        for i,rows in enumerate(pool.map(fit_one,items,chunksize=1)):
            allrows.extend(rows)
            if (i+1)%25==0:print(f'{args.sample}: {i+1}/{len(items)} objects',flush=True)
    table=pd.DataFrame(allrows);tag=args.sample+('-pilot' if args.limit else '')
    output=WORK/(tag+'-age-bounds.csv');table.to_csv(output,index=False)
    summaries=[]
    for floor,group in table.groupby('model_discrepancy_mag'):
        good=group.age_low_Gyr.notna();joint=group.joint_age_low_Gyr.notna();tested=group.heldout_Dn_intervals_disjoint.notna()
        summaries.append(dict(model_discrepancy_mag=floor,objects=len(group),photometry_feasible=int(good.sum()),
            age_interval_width_Gyr=quantiles(group.age_high_Gyr-group.age_low_Gyr),lower_age_Gyr=quantiles(group.age_low_Gyr),upper_age_Gyr=quantiles(group.age_high_Gyr),
            joint_feasible=int(joint.sum()),joint_interval_width_Gyr=quantiles(group.joint_age_high_Gyr-group.joint_age_low_Gyr),
            heldout_Dn_tested=int(tested.sum()),heldout_Dn_disjoint=int(group.loc[tested,'heldout_Dn_intervals_disjoint'].astype(bool).sum()),
            numerical_failures=int(group.numerical_failures.sum()),minimum_chi2=quantiles(group.minimum_chi2)))
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),sample=args.sample,objects=len(items),excluded=excluded,workers=args.workers,
        seconds=time.perf_counter()-start,summaries=summaries,assumptions='Finite SSP/dust/metallicity grid; formed-mass age; fixed cosmological clock; diagonal flux errors; central fibres not necessarily SN-local populations. 0.03mag floor is a sensitivity, not a measured model covariance. Index constraints include only >=99.5% coverage in both bands.',
        confidence='Statistical-only photometric95% absolute Gaussian mean-flux ellipsoid; joint95% via97.5% photometry and97.5% Fieller ratio Bonferroni. Exact coverage requires correct fixed Gaussian variances and an included physical model. The0.03mag observed-flux-dependent floor gives assumed compatibility sensitivity, not calibrated95% coverage. Heldout95%+95% disjointness has at least90% joint coverage by Bonferroni (at most10% false incompatibility under the assumptions); do not call it a95% rejection. None is a Bayesian posterior credible interval.',
        source_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        output=str(output.relative_to(ROOT)),output_sha256=hashlib.sha256(output.read_bytes()).hexdigest(),
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('model.py')]},design_sha256=hashlib.sha256(Path(__file__).with_name('design.json').read_bytes()).hexdigest())
    (OUT/(tag+'-age-bounds.json')).write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
