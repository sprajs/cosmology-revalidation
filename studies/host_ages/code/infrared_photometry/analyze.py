#!/usr/bin/env python3
"""Conservative spherical catalogue association; no physical SED inference."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
from datetime import datetime,timezone
import json,re,sys
from importlib.metadata import version
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.spatial import cKDTree
from acquire import ROOT,HERE,WORK,OUT,TABLES,sha,sample

def xyz(ra,dec):
    r,d=np.deg2rad(np.asarray(ra,float)),np.deg2rad(np.asarray(dec,float))
    return np.column_stack([np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)])

def separation(a,b):
    # Stable small-angle exact spherical separation from unit-vector chord.
    return np.rad2deg(2*np.arcsin(np.clip(np.linalg.norm(a-b,axis=-1)/2,0,1)))*3600

def numeric(value):
    try:return float(value)
    except (ValueError,TypeError):return np.nan

def flags(value):
    match=re.search(r'(-?\d+)$',str(value))
    return int(match.group(1)) if match else -1

def main():
    hosts=sample();acq=json.loads((OUT/'acquisition.json').read_text())
    catalogs={};candidates=[];allmeta={};inputs={OUT/'acquisition.json':sha(OUT/'acquisition.json')}
    for rec in acq['catalogues']:
        table=rec['catalogue'];p=ROOT/rec['path'];assert sha(p)==rec['sha256']
        q=pd.read_csv(p,dtype=str);q['catalogue']=table;q['source_key']=table+':'+q.cntr;q['field']=rec['field']
        q['ra']=q.ra.astype(float);q['dec']=q.dec.astype(float);assert q.source_key.is_unique
        catalogs[table]=q
        candidates.append(q[['catalogue','source_key','cntr','ra','dec','field']])
        # IRSA SERVS descriptions contain unescaped commas/quotes. The first
        # three metadata fields are simple unquoted identifiers/types/units;
        # preserve everything after their third comma as the literal description.
        lines=(WORK/f'{table}-columns.csv').read_text().splitlines()
        assert lines[0]=='column_name,datatype,unit,description'
        meta=pd.DataFrame([line.split(',',3) for line in lines[1:] if line],columns=lines[0].split(','))
        assert meta.column_name.is_unique
        allmeta[table]=meta.set_index('column_name')
    pd.concat(candidates,ignore_index=True).to_csv(WORK/'candidates.csv',index=False)
    # All native source IDs remain string-valued through ingestion and output.
    associations=[];nearby=[];optical_inputs={}
    for field,tables in TABLES.items():
        op=ROOT/f'.work/host-transport/des/Y3_DEEP_FIELDS_PHOTOM-SN-{field}-0000.parquet'
        optical_inputs[str(op.relative_to(ROOT))]=sha(op)
        opt=pq.read_table(op,columns=['ID','RA','DEC','FLAGS','MASK_FLAGS','KNN_CLASS']).to_pandas()
        opt=opt[(opt.FLAGS==0)&(opt.MASK_FLAGS==0)&(opt.KNN_CLASS==1)].copy()
        ox=xyz(opt.RA,opt.DEC);otree=cKDTree(ox)
        h=hosts[hosts.field==field];hx=xyz(h.deep_RA,h.deep_DEC)
        for table in tables:
            q=catalogs[table];cx=xyz(q.ra,q.dec);tree=cKDTree(cx)
            radius=2. if '_24_' in table else 1.
            for j,host in enumerate(h.itertuples()):
                hit=tree.query_ball_point(hx[j],2*np.sin(np.deg2rad(75/3600)/2))
                distances=separation(cx[hit],hx[j]) if hit else np.array([])
                order=np.argsort(distances);hit=np.asarray(hit,dtype=int)[order];distances=distances[order]
                row=dict(SNID=str(host.SNID),deep_ID=str(host.deep_ID),field=field,catalogue=table,
                    family='SERVS' if table.startswith('servs') else ('SWIRE24' if '_24_' in table else 'SWIRE'),
                    primary_radius_arcsec=radius,candidates_within75=len(hit),source_key='',separation_arcsec=np.nan,
                    primary_unique_ir=False,target_nearest_optical=False,other_optical_within_primary=0,
                    other_optical_within3=0,other_optical_within6=0,nearest_other_optical_arcsec=np.nan)
                for r in [1,2,3,6]:row[f'candidates_within{r}']=int((distances<=r).sum())
                for idx,d in zip(hit,distances):
                    nearby.append(dict(SNID=str(host.SNID),deep_ID=str(host.deep_ID),catalogue=table,
                        source_key=q.iloc[idx].source_key,separation_arcsec=float(d)))
                if len(hit):
                    source=q.iloc[hit[0]];v=cx[hit[0]];row['source_key']=source.source_key
                    row['separation_arcsec']=float(distances[0])
                    row['primary_unique_ir']=bool((distances<=radius).sum()==1)
                    near=otree.query_ball_point(v,2*np.sin(np.deg2rad(75/3600)/2))
                    od=separation(ox[near],v);oi=opt.iloc[near].ID.to_numpy()
                    target=oi==int(host.deep_ID);assert target.any() or distances[0]>74
                    other=od[~target];dtarget=float(od[target][0]) if target.any() else float('inf')
                    row['target_nearest_optical']=bool(not len(other) or dtarget<other.min())
                    row['other_optical_within_primary']=int((other<=radius).sum())
                    row['other_optical_within3']=int((other<=3).sum());row['other_optical_within6']=int((other<=6).sum())
                    row['nearest_other_optical_arcsec']=float(other.min()) if len(other) else np.nan
                row['geometry_pass']=bool(row['primary_unique_ir'] and row['target_nearest_optical'] and row['other_optical_within_primary']==0)
                associations.append(row)
    a=pd.DataFrame(associations);assert len(a)==795
    # Same catalogue source assigned to two selected SN hosts is not independent.
    counts=a[a.primary_unique_ir].groupby('source_key').SNID.nunique()
    a['source_shared_by_selected_hosts']=a.primary_unique_ir & (a.source_key.map(counts).fillna(0).astype(int)>1)
    a['geometry_pass'] &= ~a.source_shared_by_selected_hosts
    a.to_csv(WORK/'associations.csv',index=False);pd.DataFrame(nearby).to_csv(WORK/'host-candidates.csv',index=False)
    phot=[]
    for assoc in a.itertuples():
        table=assoc.catalogue;q=catalogs[table].set_index('source_key');source=q.loc[assoc.source_key] if assoc.source_key else None
        bands=[(3.6,'1'),(4.5,'2')] if table.startswith('servs') else ([(24.,'24')] if '_24_' in table else [(3.6,'36'),(4.5,'45'),(5.8,'58'),(8.,'80'),(24.,'24')])
        for band,key in bands:
            servs=table.startswith('servs');fc=f'flux_aper_2_{key}' if servs else f'flux_ap2_{key}'
            ec=f'fluxerr_aper_2_{key}' if servs else f'uncf_ap2_{key}'
            assert allmeta[table].loc[fc,'unit']=='uJy' and allmeta[table].loc[ec,'unit']=='uJy'
            f=numeric(source[fc]) if source is not None else np.nan;e=numeric(source[ec]) if source is not None else np.nan
            valid=bool(np.isfinite(f) and np.isfinite(e) and e>0 and f!=-99 and e!=-99)
            if servs:
                covcol=f'cov_{key}';covkind='hundred_second_frames'
                cov=numeric(source[covcol]) if source is not None else np.nan
                flagtext=source[f'flags_{key}'] if source is not None else ''
                flag=flags(flagtext)
                mask=flags(source[f'mask_{key}']) if source is not None else -1
                extension=np.nan;detid=''
            else:
                covcol=f'cov_avg_{key}' if '_24_' in table else f'cov_{key}'
                covkind='mean_number_frames' if '_24_' in table else 'coverage_flag'
                cov=numeric(source.get(covcol,np.nan)) if source is not None else np.nan
                flagtext=source[f'flgs_{key}'] if source is not None else ''
                flag=flags(flagtext);mask=np.nan
                extension=numeric(source.get(f'ext_fl_{key}',np.nan)) if source is not None else np.nan
                detid=source.get('detid_24','') if source is not None and band==24 else ''
            conflict=bool(table=='chandra_cat_f05' and band!=24)
            extension_conflict=bool(not servs and extension==-1)
            band_geometry=bool(assoc.geometry_pass);mips_link='not_applicable'
            if band==24 and '_24_' not in table:
                # The bandmerged coordinate is IRAC dominated. Its 24um flux
                # needs the actual 24um centroid and association checked too.
                mt=next(t for t in TABLES[assoc.field] if '_24_' in t)
                ma=a[(a.SNID==assoc.SNID)&(a.catalogue==mt)].iloc[0]
                ms=catalogs[mt].set_index('source_key').loc[ma.source_key] if ma.source_key else None
                same=bool(ms is not None and str(detid)==str(ms.detid_24) and str(detid) not in ['', '-99','nan'])
                band_geometry=bool(band_geometry and ma.geometry_pass and same)
                mips_link='same_detection_verified' if band_geometry else 'separate_centroid_or_identity_not_verified'
            within=bool(np.isfinite(assoc.separation_arcsec) and assoc.separation_arcsec<=assoc.primary_radius_arcsec)
            if not within:state='no_catalogue_association_coverage_unknown'
            elif not valid:state='covered_flux_missing' if cov>0 else ('no_coverage' if cov==0 else 'flux_missing_coverage_unknown')
            elif not band_geometry:state='ambiguous_association'
            else:state='catalogue_measurement_geometry_pass'
            phot.append(dict(SNID=assoc.SNID,deep_ID=assoc.deep_ID,field=assoc.field,catalogue=table,source_key=assoc.source_key,
                band_um=band,native_flux_column=fc,native_error_column=ec,native_flux_uJy=f,native_error_uJy=e,
                flux_valid=valid,within_primary=within,geometry_pass=band_geometry,anchor_geometry_pass=assoc.geometry_pass,
                bandmerged24_centroid_link=mips_link,coverage_native=cov,coverage_column=covcol,coverage_kind=covkind,
                native_flag=flag,native_flag_text=flagtext,native_mask=mask,native_extended_flag=extension,detid_24=detid,
                aperture_metadata_conflict=conflict,extension_semantics_conflict=extension_conflict,
                clean_photometry=bool(valid and flag==0 and (not servs or mask==0) and not(extension>0) and not conflict and not extension_conflict),status=state))
    p=pd.DataFrame(phot);p.to_csv(WORK/'photometry.csv',index=False)
    summary=[]
    for (family,field),g in a.groupby(['family','field']):
        summary.append(dict(family=family,field=field,hosts=len(g),radii_arcsec={str(r):dict(any_match=int((g[f'candidates_within{r}']>0).sum()),unique_ir=int((g[f'candidates_within{r}']==1).sum())) for r in [1,2,3,6]},
            primary_unique_ir=int(g.primary_unique_ir.sum()),geometry_pass=int(g.geometry_pass.sum()),
            passing_with_other_optical_within3=int((g.geometry_pass&(g.other_optical_within3>0)).sum()),
            passing_with_other_optical_within6=int((g.geometry_pass&(g.other_optical_within6>0)).sum()),
            target_not_nearest_within_primary=int((g.primary_unique_ir&~g.target_nearest_optical).sum()),
            optical_competition=int((g.primary_unique_ir&(g.other_optical_within_primary>0)).sum()),
            source_shared_by_selected_hosts=int(g.source_shared_by_selected_hosts.sum())))
    bands=[]
    for (table,band),g in p.groupby(['catalogue','band_um']):
        good=g.geometry_pass & g.flux_valid
        bands.append(dict(catalogue=table,band_um=band,geometry_valid_measurements=int(good.sum()),
            clean_geometry_measurements=int((good&g.clean_photometry).sum()),
            extension_semantics_conflict_valid_associations=int((good&g.extension_semantics_conflict).sum()),
            statuses=g.status.value_counts().to_dict(),aperture_metadata_conflict=bool(g.aperture_metadata_conflict.any())))
    union=[dict(band_um=float(band),unique_hosts_with_valid_geometric_association=int(g.loc[g.geometry_pass&g.flux_valid,'SNID'].nunique()),
                interpretation='Union of qualifying catalogue associations, not averaged or independent repeated measurements.')
           for band,g in p.groupby('band_um')]
    consistency=[];ratios=[]
    for field,tables in TABLES.items():
        swire=next(t for t in tables if not t.startswith('servs') and '_24_' not in t)
        servs=next(t for t in tables if t.startswith('servs'))
        for band in [3.6,4.5]:
            left=p[(p.catalogue==swire)&(p.band_um==band)].set_index('SNID')
            right=p[(p.catalogue==servs)&(p.band_um==band)].set_index('SNID')
            for snid in left.index.intersection(right.index):
                l,r=left.loc[snid],right.loc[snid]
                if not(l.geometry_pass and r.geometry_pass):continue
                lc=catalogs[swire].set_index('source_key').loc[l.source_key]
                rc=catalogs[servs].set_index('source_key').loc[r.source_key]
                d=float(separation(xyz([lc.ra],[lc.dec]),xyz([rc.ra],[rc.dec]))[0])
                both=bool(l.flux_valid and r.flux_valid);positive=bool(both and l.native_flux_uJy>0 and r.native_flux_uJy>0)
                ratios.append(dict(SNID=snid,field=field,band_um=band,swire_source_key=l.source_key,servs_source_key=r.source_key,
                    infrared_centroid_separation_arcsec=d,centroids_agree_within1=d<=1,
                    swire_flux_uJy=l.native_flux_uJy,servs_flux_uJy=r.native_flux_uJy,
                    signed_servs_minus_swire_uJy=r.native_flux_uJy-l.native_flux_uJy if both else np.nan,
                    valid_native_pair=both,positive_flux_pair=positive,
                    servs_over_swire=r.native_flux_uJy/l.native_flux_uJy if positive else np.nan,
                    swire_metadata_conflict=l.aperture_metadata_conflict,
                    both_clean_photometry=bool(l.clean_photometry and r.clean_photometry)))
    cr=pd.DataFrame(ratios);cr.to_csv(WORK/'release-consistency.csv',index=False)
    for (field,band),g in cr.groupby(['field','band_um']):
        good=g.centroids_agree_within1 & g.positive_flux_pair
        values=g.loc[good,'servs_over_swire']
        consistency.append(dict(field=field,band_um=band,geometry_pairs=len(g),centroid_disagreements_gt1=int((~g.centroids_agree_within1).sum()),
            missing_or_invalid_pairs=int((~g.valid_native_pair).sum()),nonpositive_valid_pairs=int((g.valid_native_pair&~g.positive_flux_pair).sum()),
            ratio_pairs=len(values),ratio_p05_median_p95=np.quantile(values,[.05,.5,.95]).tolist() if len(values) else None,
            clean_pairs=int((good & g.both_clean_photometry).sum()),
            interpretation='Release consistency only; detection selection, aperture/extension flags, shared IRAC/SWIRE-tied calibration and intrinsic variability remain. Individual exposure overlap was not audited; no independent-calibration or chi-square uncertainty claim.'))
    outputs={str((WORK/n).relative_to(ROOT)):sha(WORK/n) for n in ['candidates.csv','hosts.csv','associations.csv','host-candidates.csv','photometry.csv','release-consistency.csv']}
    result=dict(created_utc=datetime.now(timezone.utc).isoformat(),status='catalogue_associations_only',hosts=265,
        environment=dict(python=sys.version,packages={name:version(name) for name in ['numpy','pandas','scipy','pyarrow']}),
        candidate_sources=sum(len(q) for q in catalogs.values()),associations=len(a),photometry_rows=len(p),
        matching=summary,bands=bands,unique_host_union_by_band=union,release_consistency=consistency,photometry_statuses=p.status.value_counts().to_dict(),
        code_sha256=sha(Path(__file__)),dependencies_sha256={str((HERE/n).relative_to(ROOT)):sha(HERE/n) for n in ['acquire.py','source-clarification.json','acquisition-amendment.json','analysis-addendum.json']},
        design_sha256=sha(HERE/'design.json'),acquisition_sha256=sha(OUT/'acquisition.json'),
        optical_input_sha256=optical_inputs,output_sha256=outputs,
        limits=['No counterpart is not a numerical upper limit or zero flux.',
                'Native aperture2 is point-source corrected; extended/flagged objects retained with caution.',
                'CDFS SWIRE TAP error-aperture radius labels conflict with flux labels and explanatory common-aperture documentation; retained, not repaired.',
                'SWIRE documentation contradicts the meaning of extended-source flag -1; retain it and exclude it from the strict flags-clear diagnostic. Flag0 remains indeterminate, not certified pointlike.',
                'SERVS is calibrated to SWIRE and uses IRAC; individual exposure overlap between these catalogue records was not audited. SWIRE bandmerged and separate24 detections may duplicate.',
                'IRAC includes stellar continuum and is not dust-only luminosity. No physical age, dust luminosity or age correction inferred.',
                'The recovered catalogue columns do not provide a joint cross-band or cross-survey calibration covariance.',
                'Optical competitor screen includes only unflagged deep-catalogue galaxies; undetected/flagged/other-type IR emitters remain possible.'])
    (OUT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
