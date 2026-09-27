#!/usr/bin/env python3
"""Conservative host-to-published-galaxy membership inventory for ZTF/TITAN.

Group membership is inherited only from a spatially/redshift matched individual
host galaxy. Being close to a group on the sky never establishes membership.
"""
from pathlib import Path
import datetime,gzip,hashlib,json,sys,urllib.request
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import ascii
from scipy.stats import norm

ROOT=Path(__file__).resolve().parents[4]
HERE=Path(__file__).parent
OUT=ROOT/'studies/host_ages/results/population_transport'
WORK=ROOT/'.work/population-transport'
TULLY=WORK/'tully'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')


def main():
    registry=json.loads((OUT/'group-inputs.json').read_text())
    for row in registry['files']:
        p=ROOT/row['path']
        if not p.exists():
            p.parent.mkdir(parents=True,exist_ok=True)
            p.write_bytes(urllib.request.urlopen(row['url'],timeout=60).read())
        assert sha(p)==row['sha256'],p
        if p.suffix=='.gz':p.with_suffix('').write_bytes(gzip.decompress(p.read_bytes()))
    plan=json.loads((HERE/'group-scenarios.json').read_text())
    source=ROOT/'.work/environment-validation/titan-ztf-join.csv'
    d=pd.read_csv(source)
    d=d[d.selected_titan].copy().reset_index(drop=True)
    assert len(d)==401 and d.ztfname.is_unique
    gal=ascii.read(TULLY/'table5.dat',readme=str(TULLY/'ReadMe'),format='cds').to_pandas()
    groups=ascii.read(TULLY/'table3.dat',readme=str(TULLY/'ReadMe'),format='cds').to_pandas()
    assert len(gal)==43038 and gal.PGC.is_unique and groups.Nest.is_unique
    assert np.isfinite(d[['ra_host','dec_host','redshift']]).all().all()
    hostcoords=SkyCoord(d.ra_host.to_numpy()*u.deg,d.dec_host.to_numpy()*u.deg,frame='icrs')
    galaxycoords=SkyCoord(l=gal.GLON.to_numpy()*u.deg,b=gal.GLAT.to_numpy()*u.deg,frame='galactic').icrs
    nearest,separation,_=hostcoords.match_to_catalog_sky(galaxycoords)
    _,second,_=hostcoords.match_to_catalog_sky(galaxycoords,nthneighbor=2)
    m=gal.iloc[nearest].reset_index(drop=True)
    joined=pd.concat([d,m.add_prefix('tully_')],axis=1)
    joined['host_galaxy_separation_arcsec']=separation.arcsec
    joined['second_galaxy_separation_arcsec']=second.arcsec
    joined['heliocentric_velocity_difference_kms']=299792.458*joined.redshift-joined.tully_HV
    statuses=[]
    for radius in [1,3,5]:
        angular=separation.arcsec<=radius
        unique=second.arcsec>radius
        velocity=np.abs(joined.heliocentric_velocity_difference_kms.to_numpy())<=300
        accepted=angular&unique&velocity
        group=accepted&(joined.tully_Nmb.to_numpy()>=2)&(joined.tully_Nest.to_numpy()>0)
        selected=joined[group]
        count=selected.groupby('tully_Nest').ztfname.nunique()
        repeated=count[count>=2].index
        distinct=selected.groupby('tully_Nest').tully_PGC.nunique()
        cousins=distinct[distinct>=2].index
        statuses.append(dict(radius_arcsec=radius,within_radius=int(angular.sum()),
            ambiguous_more_than_one_galaxy=int((angular&~unique).sum()),
            within_radius_velocity_rejected=int((angular&unique&~velocity).sum()),
            accepted_individual_host_galaxies=int(accepted.sum()),
            group_member_events=int(group.sum()),unique_groups=int(len(count)),
            repeated_groups=int(len(repeated)),events_in_repeated_groups=int(selected.tully_Nest.isin(repeated).sum()),
            distinct_host_cousin_groups=int(len(cousins)),events_in_distinct_host_cousin_groups=int(selected.tully_Nest.isin(cousins).sum())))
        if radius==3:
            joined['accepted_galaxy_match']=accepted;joined['accepted_group_member']=group
            primary=selected.copy()
    repeated=primary.groupby('tully_Nest').ztfname.nunique()
    repeated_ids=repeated[repeated>=2].index
    rows=[]
    for nest,events in primary[primary.tully_Nest.isin(repeated_ids)].groupby('tully_Nest'):
        group=groups[groups.Nest==nest].iloc[0]
        distance=10**((float(group.DM)-25)/5)
        extent=1.5*float(group.R2t)
        small_sigma=5/np.log(10)*extent/distance/np.sqrt(3)
        span=5*np.log10((distance+extent)/(distance-extent)) if distance>extent else None
        # Relative redshift-distance scatter if one used member z as distance;
        # not an independent term to add to group-intercept apparent magnitudes.
        pv_sigma=5/np.log(10)*float(group.sigP)/float(group['<Vcmba>'])
        a=events.mass_weighted_age_50.to_numpy()
        ae=(events.mass_weighted_age_84.to_numpy()-events.mass_weighted_age_16.to_numpy())/2
        rows.append(dict(Nest=int(nest),events=len(events),distinct_PGC=int(events.tully_PGC.nunique()),
            IAU_names=events.iau_name.tolist(),PGC=events.tully_PGC.astype(int).tolist(),
            z_heliocentric=events.redshift.tolist(),global_age_summary_gyr=a.tolist(),
            age_half_68_width_gyr=ae.tolist(),age_range_gyr=float(a.max()-a.min()),
            centred_age_sum_squares_gyr2=float(np.sum((a-a.mean())**2)),
            distinct_host_global_age_design_contrast=bool(events.tully_PGC.nunique()>1),
            group_members=int(group.Nmb),group_R2t_Mpc=float(group.R2t),group_distance_Mpc=distance,
            ratio_R2t_to_distance=float(group.R2t)/distance,
            conditional_uniform_depth_sigma_mag=small_sigma,
            maximum_endpoint_depth_pair_span_mag=span,
            depth_scope='Offsets of distinct group members under an imposed uniform-depth model. Events with the same PGC share the host distance and its depth offset; these are not independent depth errors for siblings.',
            catalogue_velocity_dispersion_kms=float(group.sigP),
            member_redshift_distance_sigma_mag_for_comparison=pv_sigma))
    # Independent great-circle check in ICRS, including transform round-trip.
    r1=np.deg2rad(d.ra_host.to_numpy());d1=np.deg2rad(d.dec_host.to_numpy())
    r2=galaxycoords.ra.radian[nearest];d2=galaxycoords.dec.radian[nearest]
    hav=np.sin((d2-d1)/2)**2+np.cos(d1)*np.cos(d2)*np.sin((r2-r1)/2)**2
    direct=2*np.arcsin(np.sqrt(np.clip(hav,0,1)))*180/np.pi*3600
    separation_error=float(np.max(np.abs(direct-separation.arcsec)))
    assert separation_error<1e-7
    joined.to_csv(WORK/'ztf-titan-tully-crosswalk.csv',index=False)
    group_info=sum(r['centred_age_sum_squares_gyr2']/(.12**2+r['conditional_uniform_depth_sigma_mag']**2) for r in rows if r['distinct_PGC']>1)
    repeats=primary[primary.tully_Nest.isin(repeated_ids)].copy()
    # Global host ages must be shared by siblings; inconsistent posterior centres
    # for the same PGC are not two distinct underlying physical global ages.
    age_per_host=repeats.groupby('tully_PGC').mass_weighted_age_50.transform('mean').to_numpy()
    membership=np.column_stack([(repeats.tully_Nest.to_numpy()==nest).astype(float) for nest in repeated_ids])
    design=np.column_stack([membership,repeats.x1,repeats.c,age_per_host])
    rank=int(np.linalg.matrix_rank(design))
    rank_without_age=int(np.linalg.matrix_rank(design[:,:-1]))
    centred_age=age_per_host-membership@np.linalg.lstsq(membership,age_per_host,rcond=None)[0]
    controls=np.column_stack([membership,repeats.x1,repeats.c])
    age_residual=age_per_host-controls@np.linalg.lstsq(controls,age_per_host,rcond=None)[0]
    leverage=dict(rows=len(repeats),parameters=int(design.shape[1]),rank=rank,rank_without_age=rank_without_age,
        age_identifiable_after_group_intercepts_width_colour=bool(rank>rank_without_age),
        residual_age_sum_squares_after_controls_gyr2=float(age_residual@age_residual),
        same_PGC_sibling_global_ages_pooled_for_design=True)
    optimistic=None if group_info==0 else 1/np.sqrt(group_info)
    result=dict(completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        selected_events=401,method=plan,sensitivity=statuses,repeated_group_details=rows,
        optimistic_group_only_known_age_sigma_mag_per_gyr=optimistic,
        optimistic_two_sided_5percent_power_for_abs_003=None if optimistic is None else float(norm.cdf(-1.96-.03/optimistic)+norm.sf(1.96-.03/optimistic)),
        group_width_colour_age_design=leverage,
        same_host_sibling_age_leverage_excluded=True,
        optimistic_information_scope='Assumes exact posterior-centre ages, independent 0.12 mag intrinsic scatter, declared uniform depth, and group intercepts only; omits width/colour/mass and joint host uncertainty. An optimistic error scale conditional on those assumptions, not a calibrated confidence interval or empirical lower bound.',
        association_limit='The actual ZTF host position is matched to an individual Tully member and its heliocentric velocity. The TITAN host identity remains conditional on shared SN identifier: old TITAN summary has no coordinates/PGC. Both DLR cuts do not prove those two host associations are identical.',
        new_age_slope_fit_performed=False,
        reason='Only one distinct-host cousin pair; group-intercept, width, colour and age design is rank deficient, and group sampling uncertainty cannot be calibrated. Same-PGC sibling summary differences are not physical global-age contrasts.' if len(rows)<5 else 'Enough groups for a separately specified conditional regression; see followup record.',
        physical_group_depth_or_velocity_model_validated=False,
        validation=dict(source_rows=43038,unique_PGC=True,unique_group_Nest=True,
            independent_haversine_separation_max_error_arcsec=separation_error,
            repeated_groups_stable_across_1_3_5_arcsec=len({(x['repeated_groups'],x['events_in_repeated_groups'],x['distinct_host_cousin_groups']) for x in statuses})==1),
        code_sha256={str(p.relative_to(ROOT)):sha(p) for p in [Path(__file__),HERE/'group-scenarios.json']},
        input_sha256={str(source.relative_to(ROOT)):sha(source),**{r['path']:r['sha256'] for r in registry['files']}},
        crosswalk_sha256=sha(WORK/'ztf-titan-tully-crosswalk.csv'))
    write(OUT/'group-feasibility.json',result)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
