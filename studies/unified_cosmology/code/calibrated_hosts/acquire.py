"""Acquire bounded public DESI DR1 matches to observed OzDES host positions."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import datetime, hashlib, io, json, time
import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT = Path(__file__).resolve().parents[4]
CODE = Path(__file__).resolve().parent
WORK = ROOT / '.work/unified-cosmology/calibrated-hosts'
RESULT = ROOT / 'studies/unified_cosmology/results/calibrated_hosts'
TAP = 'https://datalab.noirlab.edu/tap/sync'
COLS = 'targetid,mean_fiber_ra,mean_fiber_dec,z,zerr,zwarn,spectype,survey,program,healpix,zcat_primary,coadd_numexp,tsnr2_lrg,tsnr2_elg,tsnr2_bgs'


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    WORK.mkdir(parents=True, exist_ok=True)
    RESULT.mkdir(parents=True, exist_ok=True)
    queries = WORK / 'queries'; queries.mkdir(exist_ok=True)
    design = json.loads((CODE/'design.json').read_text())
    source = ROOT/design['input']
    rows = pd.read_csv(source, dtype={'SNID':str})
    grouped = rows.groupby('oz_OzDES_ID')
    assert grouped.oz_RA.nunique().max() == grouped.oz_DEC.nunique().max() == 1
    hosts = rows.drop_duplicates('oz_OzDES_ID').sort_values('oz_OzDES_ID').reset_index(drop=True)
    hosts.to_csv(WORK/'input-hosts.csv', index=False)
    radius = design['maximum_radius_arcsec']/3600

    def query_batch(start):
        subset = hosts.iloc[start:start+32]
        terms = [f"q3c_radial_query(mean_fiber_ra,mean_fiber_dec,{r.oz_RA:.12f},{r.oz_DEC:.12f},{radius:.15f})='t'" for r in subset.itertuples()]
        sql = f'SELECT TOP 10000 {COLS} FROM desi_dr1.zpix WHERE ' + ' OR '.join(terms)
        qp=queries/f'{start:04d}.sql'; cp=queries/f'{start:04d}.csv'; mp=queries/f'{start:04d}.json'
        if not (cp.exists() and qp.exists() and qp.read_text()==sql and mp.exists() and cp.read_text().startswith('targetid,')):
            if cp.exists():
                (queries/f'{start:04d}.failed-get.txt').write_bytes(cp.read_bytes())
                if mp.exists(): (queries/f'{start:04d}.failed-get.json').write_bytes(mp.read_bytes())
            qp.write_text(sql)
            t=time.monotonic()
            response=requests.post(TAP, data={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':sql},timeout=180)
            cp.write_bytes(response.content)
            meta={'url':TAP,'method':'POST','http_status':response.status_code,'elapsed_seconds':time.monotonic()-t,'query_sha256':sha(qp),'response_sha256':sha(cp),'first_host':start,'hosts':len(subset)}
            mp.write_text(json.dumps(meta,indent=2)+'\n')
            response.raise_for_status()
        assert cp.read_text().startswith('targetid,'), f'Non-CSV TAP response retained: {cp}'
        df=pd.read_csv(cp,dtype={'targetid':str})
        assert len(df)<10000, 'Server row cap reached; split query before continuing.'
        return df

    batches=list(ThreadPoolExecutor(4).map(query_batch,range(0,len(hosts),32)))
    candidates=pd.concat(batches,ignore_index=True).drop_duplicates(['targetid','survey','program'])
    candidates.to_csv(WORK/'desi-candidates.csv',index=False)
    matches=[]; statuses=[]
    if len(candidates):
        coords=SkyCoord(candidates.mean_fiber_ra.to_numpy()*u.deg,candidates.mean_fiber_dec.to_numpy()*u.deg)
    for host in hosts.itertuples():
        selected=candidates.copy()
        if len(selected):
            sep=SkyCoord(host.oz_RA*u.deg,host.oz_DEC*u.deg).separation(coords).arcsec
            selected=selected.assign(separation_arcsec=sep).loc[sep<=design['maximum_radius_arcsec']].copy()
        distinct=selected.targetid.nunique()
        clean=[]
        for row in selected.to_dict('records'):
            row['oz_OzDES_ID']=host.oz_OzDES_ID
            row['oz_z']=float(host.oz_z)
            row['redshift_delta']=float(row['z']-host.oz_z)
            row['distinct_spatial_targets']=int(distinct)
            primary=str(row['zcat_primary']).lower() in ('true','t','1')
            row['quality_ok']=row['zwarn']==0 and row['spectype']=='GALAXY' and primary
            row['clean_match']=bool(distinct==1 and row['quality_ok'] and abs(row['redshift_delta'])<=design['maximum_absolute_redshift_difference'])
            if row['clean_match']:clean.append(row)
            matches.append(row)
        assert len(clean)<=1, f'Multiple released primary spectra for host {host.oz_OzDES_ID}'
        statuses.append({'oz_OzDES_ID':host.oz_OzDES_ID,'input_SN_count':int((rows.oz_OzDES_ID==host.oz_OzDES_ID).sum()),'spatial_spectrum_rows':len(selected),'distinct_spatial_targets':distinct,'clean_spectrum_rows':len(clean),'status':'clean' if clean else ('no_spatial_match' if not len(selected) else ('ambiguous_targets' if distinct>1 else 'quality_or_redshift_failure'))})
    pd.DataFrame(matches).to_csv(WORK/'host-desi-crosswalk.csv',index=False)
    status=pd.DataFrame(statuses);status.to_csv(WORK/'host-match-status.csv',index=False)
    clean=pd.DataFrame([x for x in matches if x['clean_match']]);clean.to_csv(WORK/'clean-matches.csv',index=False)
    outputs=['input-hosts.csv','desi-candidates.csv','host-desi-crosswalk.csv','host-match-status.csv','clean-matches.csv']
    record={'schema':'calibrated-host-acquisition-v1','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'complete','input_SN_rows':len(rows),'unique_hosts':len(hosts),'query_batches':len(batches),'candidate_spectrum_rows':len(candidates),'spatial_associations':len(matches),'clean_hosts':len(clean),'status_counts':status.status.value_counts().to_dict(),'clean_z_range':[float(clean.z.min()),float(clean.z.max())] if len(clean) else None,'query_semantics':'Native Data Lab Q3C cones, radius1arcsec; all redshifts/qualities returned before local quality cuts. Coordinate fields are mean actual fiber positions.','input_sha256':{str(source.relative_to(ROOT)):sha(source)},'code_sha256':{str(p.relative_to(ROOT)):sha(p) for p in [CODE/'design.json',Path(__file__)]},'query_files_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(queries.glob('*'))},'output_sha256':{str((WORK/p).relative_to(ROOT)):sha(WORK/p) for p in outputs}}
    (RESULT/'acquisition.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({k:record[k] for k in ['unique_hosts','clean_hosts','status_counts','clean_z_range']},indent=2))


if __name__=='__main__':
    main()
