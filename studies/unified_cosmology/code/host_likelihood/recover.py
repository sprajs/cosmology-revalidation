"""Recover primary documentation, pinned author code and unblinded YSE fluxes."""
from __future__ import annotations
from datetime import datetime, timezone
import hashlib,json,subprocess,tarfile,urllib.request
from pathlib import Path
from acquire import ROOT,WORK,OUT,sha

def fetch(url,path):
    if not path.exists():path.write_bytes(urllib.request.urlopen(url,timeout=120).read())
    return dict(url=url,path=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=sha(path))

def main():
    records=[]
    meta=WORK/'yse-record.json'
    records.append(fetch('https://zenodo.org/api/records/7317476',meta))
    for item in json.loads(meta.read_text())['files']:
        if item['key'] not in ['yse_dr1_zenodo.tar.gz','README.md']:continue
        path=WORK/('yse-README.md' if item['key']=='README.md' else item['key'])
        record=fetch(item['links']['self'],path)
        algorithm,digest=item['checksum'].split(':')
        with path.open('rb') as stream:assert hashlib.file_digest(stream,algorithm).hexdigest()==digest
        record['published_checksum']=item['checksum'];records.append(record)
    with tarfile.open(WORK/'yse_dr1_zenodo.tar.gz') as archive:
        archive.extractall(WORK/'yse',filter='data')
    for name,url in [('frankenblast','https://arxiv.org/pdf/2509.08874'),('cigars','https://arxiv.org/pdf/2508.15899'),
                     ('yse','https://inspirehep.net/files/32e400834761c5a3597e3df73856cf1a')]:
        records.append(fetch(url,WORK/(name+'.pdf')))
        subprocess.run(['pdftotext',str(WORK/(name+'.pdf')),str(WORK/(name+'.txt'))],check=True)
    repositories=[]
    for repo,commit in [('anugent96/frankenblast-host','94f81ed6d27a190832c891bb8452e54b0ddbff2e'),
                        ('snai-analysis/cigars','160f823ea6e0bd0d35b0ee5428e66927db537d52')]:
        folder=WORK/repo.split('/')[-1]
        if not folder.exists():subprocess.run(['git','clone','https://github.com/'+repo+'.git',str(folder)],check=True)
        head=subprocess.check_output(['git','-C',str(folder),'rev-parse','HEAD'],text=True).strip()
        if head!=commit:subprocess.run(['git','-C',str(folder),'checkout',commit],check=True)
        tracked=subprocess.check_output(['git','-C',str(folder),'ls-files'],text=True).splitlines()
        sources={str((folder/name).relative_to(ROOT)):sha(folder/name) for name in tracked
                 if (folder/name).is_file() and (name.endswith(('.py','.md','.yaml','.json','.dvc')))}
        repositories.append(dict(url='https://github.com/'+repo,commit=commit,source_sha256=sources))
    record=dict(completed_utc=datetime.now(timezone.utc).isoformat(),code_sha256=sha(__file__),files=records,repositories=repositories,
      cigars_scope='Released target is simulated: paper mock D0 Nsel=1578, LSST-inspired selection; not an observed DES host likelihood. Bank-generation code omits later reweighting; provided bank required for exact author demonstration.',
      yse_scope='Calibrated AB difference flux at zeropoint27.5. Even the noSNRcut file drops negative flux and MAGERR>1.5; per-object redshift_err=.005 is a placeholder. These are not selection-complete untruncated light curves.')
    (OUT/'additional-acquisition.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps({'downloaded_sources':len(records),'repositories':len(repositories)}))

if __name__=='__main__':main()
