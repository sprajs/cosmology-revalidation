"""Cache official documentation and service schema used in DESI host recovery."""
import datetime,json
from pathlib import Path
import requests
from acquire import ROOT,WORK,RESULT,sha

URLS={
    'noirlab-desi.html':'https://datalab.noirlab.edu/data/desi',
    'known-issues.html':'https://data.desi.lbl.gov/doc/releases/dr1/known-issues/',
    'desi-access.html':'https://data.desi.lbl.gov/doc/access/',
    'desi-acknowledgments.html':'https://data.desi.lbl.gov/doc/acknowledgments/',
    'fastspecfit-Dn.html':'https://fastspecfit.readthedocs.io/en/3.0.0/_modules/fastspecfit/photometry.html'
}


def main():
    WORK.mkdir(parents=True,exist_ok=True);RESULT.mkdir(parents=True,exist_ok=True)
    records=[]
    for name,url in URLS.items():
        p=WORK/name
        if not p.exists():
            response=requests.get(url,timeout=60);response.raise_for_status();p.write_bytes(response.content)
        records.append({'url':url,'path':str(p.relative_to(ROOT)),'sha256':sha(p)})
    for name,query in {
        'zpix-schema.csv':"SELECT column_name,datatype,description FROM TAP_SCHEMA.columns WHERE table_name='desi_dr1.zpix'",
        'table-schema.csv':"SELECT table_name,description FROM TAP_SCHEMA.tables WHERE schema_name='desi_dr1'"
    }.items():
        p=WORK/name
        if not p.exists():
            response=requests.post('https://datalab.noirlab.edu/tap/sync',data={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':query},timeout=60)
            response.raise_for_status();p.write_bytes(response.content)
        assert not p.read_text().startswith('<?xml'), 'Service error is not a valid catalogue schema.'
    # Keep the actual installed-client discovery failure distinct from native recovery.
    from sparcl.client import SparclClient
    try:
        c=SparclClient(connect_timeout=20,read_timeout=45)
        sparcl={'status':'initialized','client':str(c)}
    except Exception as e:
        sparcl={'status':'failed','exception_type':type(e).__name__,'message':str(e),'action':'Native NERSC byte-range retrieval is independently executable; no claim that public spectra are unavailable.'}
    doc={'schema':'calibrated-host-source-registry-v1','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'recorded','sources':records,'sparcl_client_probe':sparcl,'source_sha256':{str(Path(__file__).relative_to(ROOT)):sha(__file__)},'discovery_files_sha256':{str(p.relative_to(ROOT)):sha(p)for p in WORK.glob('test-*')},'catalogue_schema_sha256':{str(p.relative_to(ROOT)):sha(p)for p in WORK.glob('*schema.csv')},'known_issue_interpretation':'DR1 inverse variances understate uncertainty at high per-pixel SNR; no invented universal multiplicative repair. Multi-exposure native resolution matrices have documented mis-weighting and are preserved as released.'}
    (RESULT/'sources.json').write_text(json.dumps(doc,indent=2)+'\n')
    print(json.dumps({'documents':len(records),'sparcl':sparcl['status']},indent=2))


if __name__=='__main__':main()
