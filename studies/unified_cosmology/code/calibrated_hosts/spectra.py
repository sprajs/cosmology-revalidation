"""Read only matched native DESI HEALPix FITS rows using HTTP byte ranges."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import datetime, json, time
import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from acquire import ROOT, WORK, RESULT, CODE, sha


def main():
    clean=pd.read_csv(WORK/'clean-matches.csv',dtype={'targetid':str})
    hosts=pd.read_csv(WORK/'input-hosts.csv').set_index('oz_OzDES_ID')
    out=WORK/'spectra';out.mkdir(exist_ok=True)
    clean['source_url']=[f'https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/{r.survey}/{r.program}/{r.healpix//100}/{r.healpix}/coadd-{r.survey}-{r.program}-{r.healpix}.fits' for r in clean.itertuples()]

    def group(item):
        url,rows=item;pending=[r for r in rows.to_dict('records') if not (out/f"{r['targetid']}.json").exists()]
        if not pending:return [json.loads((out/f'{x}.json').read_text())for x in rows.targetid]
        head=requests.head(url,timeout=45);head.raise_for_status();assert head.headers.get('Accept-Ranges')=='bytes'
        started=time.monotonic()
        with fits.open(url,use_fsspec=True,fsspec_kwargs={'block_size':65536,'cache_type':'readahead'},memmap=False) as hdus:
            targetids=hdus['FIBERMAP'].data['TARGETID']
            for r in pending:
                target=r['targetid'];ix=np.flatnonzero(targetids==int(target));assert len(ix)==1
                index=int(ix[0]);native=hdus['FIBERMAP'].data[index]
                meta={k:(native[k].item() if isinstance(native[k],np.generic) else native[k]) for k in hdus['FIBERMAP'].data.names}
                meta={k:(v.decode().strip() if isinstance(v,bytes) else v) for k,v in meta.items()}
                h=hosts.loc[r['oz_OzDES_ID']]
                separation=SkyCoord(float(meta['TARGET_RA'])*u.deg,float(meta['TARGET_DEC'])*u.deg).separation(SkyCoord(h.oz_RA*u.deg,h.oz_DEC*u.deg)).arcsec
                arrays={};units={};headers={}
                for b in 'BRZ':
                    arrays[b+'_WAVELENGTH']=np.array(hdus[b+'_WAVELENGTH'].data)
                    for field in ['FLUX','IVAR','MASK']:
                        arrays[b+'_'+field]=np.array(hdus[b+'_'+field].section[index,:])
                    arrays[b+'_RESOLUTION']=np.array(hdus[b+'_RESOLUTION'].section[index,:,:])
                    n=len(arrays[b+'_WAVELENGTH'])
                    assert all(arrays[b+'_'+x].shape==(n,)for x in ['FLUX','IVAR','MASK'])
                    assert arrays[b+'_RESOLUTION'].shape[-1]==n and np.all(np.diff(arrays[b+'_WAVELENGTH'])>0)
                    units[b]=hdus[b+'_FLUX'].header.get('BUNIT')
                    headers[b]={x:hdus[b+'_'+x].header.tostring()for x in ['FLUX','IVAR','MASK','RESOLUTION','WAVELENGTH']}
                npz=out/f'{target}.npz';np.savez_compressed(npz,**arrays)
                record={'targetid':target,'oz_OzDES_ID':r['oz_OzDES_ID'],'desi_z':r['z'],'ozdes_z':r['oz_z'],'source_url':url,'source_headers':{k:head.headers.get(k)for k in ['Content-Length','ETag','Last-Modified','Accept-Ranges']},'native_row_index':index,'native_target_coordinate_separation_arcsec':float(separation),'native_coordinate_gate':bool(separation<=1),'native_objtype':str(meta['OBJTYPE']),'native_coadd_fiberstatus':int(meta['COADD_FIBERSTATUS']),'fibermap':meta,'flux_units':units,'extension_headers':headers,'source_hash_scope':'Only retrieved native row arrays and headers are hashed. Remote full FITS not downloaded or claimed fully hashed.','array_file':str(npz.relative_to(ROOT)),'array_sha256':sha(npz),'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'group_elapsed_seconds_at_write':time.monotonic()-started}
                (out/f'{target}.json').write_text(json.dumps(record,indent=2)+'\n')
        return [json.loads((out/f'{x}.json').read_text())for x in rows.targetid]

    records=sum(list(ThreadPoolExecutor(4).map(group,clean.groupby('source_url'))),[])
    assert len(records)==len(clean)
    failed=[r['targetid']for r in records if not r['native_coordinate_gate'] or r['native_objtype']!='TGT' or r['native_coadd_fiberstatus']!=0]
    summary={'schema':'calibrated-host-native-spectra-v1','status':'complete','spectra':len(records),'unique_remote_healpix_files':clean.source_url.nunique(),'native_identity_or_quality_failures':failed,'maximum_native_target_offset_arcsec':max(r['native_target_coordinate_separation_arcsec']for r in records),'total_saved_array_bytes':sum((ROOT/r['array_file']).stat().st_size for r in records),'input_sha256':{str((WORK/'clean-matches.csv').relative_to(ROOT)):sha(WORK/'clean-matches.csv')},'code_sha256':{str(Path(__file__).relative_to(ROOT)):sha(__file__)},'retrieved_files_sha256':{str(p.relative_to(ROOT)):sha(p)for p in sorted(out.glob('*'))},'quality_note':'FIBERMAP coadded metadata have known DR1 caveats; source TAP v1-style summary quantities remain separate. Full native resolution retained despite documented DR1 coadd mis-weighting issue.','remote_files_fully_downloaded':False}
    (RESULT/'spectra.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps({k:v for k,v in summary.items()if not k.endswith('sha256')},indent=2))


if __name__=='__main__':main()
