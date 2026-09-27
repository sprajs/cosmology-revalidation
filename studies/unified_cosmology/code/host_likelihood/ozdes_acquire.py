"""Acquire every distinct publicly linked high-redshift host spectrum."""
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import hashlib,json,urllib.request
from pathlib import Path
import pandas as pd
from astropy.io import fits
from acquire import ROOT,WORK,OUT,sha

def main():
    source=ROOT/'.work/unified-cosmology/survey-selection/followup/dovekie-ozdes-clean.csv'
    parent=ROOT/'studies/unified_cosmology/results/survey_selection/ozdes-crossmatch.json'
    rows=pd.read_csv(source,dtype={'SNID':str,'oz_OzDES_ID':str})
    columns=['oz_OzDES_ID','oz_spectrum_url','oz_RA','oz_DEC','oz_z']
    assert rows.groupby('oz_OzDES_ID')[columns[1:]].nunique().to_numpy().max()==1
    targets=rows[columns].drop_duplicates('oz_OzDES_ID')
    folder=WORK/'ozdes';folder.mkdir(exist_ok=True)
    def get(row):
        path=folder/row.oz_spectrum_url.rsplit('/',1)[-1]
        if not path.exists():
            failures=[]
            for trial in range(3):
                try:
                    data=urllib.request.urlopen(row.oz_spectrum_url,timeout=45).read()
                    temp=path.with_suffix('.fits.partial');temp.write_bytes(data)
                    with fits.open(temp) as h:
                        assert h[0].data.ndim==1 and h['VARIANCE'].data.shape==h[0].data.shape
                    temp.replace(path);break
                except Exception as exc:failures.append(type(exc).__name__+': '+str(exc))
            else:return dict(target=row.oz_OzDES_ID,url=row.oz_spectrum_url,status='download_failed',attempts=failures)
        with fits.open(path) as h:
            header=h[0].header
            assert header['SOURCE']==row.oz_OzDES_ID
            assert abs(header['Z']-row.oz_z)<1e-5
            units=header.get('BUNIT');ext=len(h)
        return dict(target=row.oz_OzDES_ID,url=row.oz_spectrum_url,path=str(path.relative_to(ROOT)),
                    sha256=sha(path),bytes=path.stat().st_size,status='downloaded',BUNIT=units,HDUs=ext)
    results=[]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for future in as_completed([pool.submit(get,row) for row in targets.itertuples(index=False)]):
            results.append(future.result())
            if len(results)%100==0:print('spectra',len(results),'/',len(targets),flush=True)
    records=dict(completed_utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),
       design_sha256=sha(Path(__file__).with_name('ozdes-design.json')),
       inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in [source,parent]},
       SN_rows=len(rows),distinct_targets=len(targets),shared_targets=rows.groupby('oz_OzDES_ID').size().loc[lambda n:n>1].to_dict(),
       records=sorted(results,key=lambda x:x['target']))
    (OUT/'ozdes-acquisition.json').write_text(json.dumps(records,indent=2)+'\n')
    print(json.dumps({'complete':sum(x['status']=='downloaded' for x in results),'requested':len(targets)}))

if __name__=='__main__':main()
