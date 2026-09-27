#!/usr/bin/env python3
"""Pin the documented Pippin preprocessing source used for the closure audit."""
import json,hashlib
from common import ROOT,WORK,RESULTS,sha
from acquire import get
COMMIT='e15ee2161292821d4f7ecf56fc010cde4080bff6'
FILES=['pippin/dataprep.py','pippin/classifiers/supernnova.py','pippin/tasks/dataprep','pippin/tasks/supernnova']

def main():
    folder=WORK/'pippin-source';folder.mkdir(exist_ok=True)
    tree=get(f'https://api.github.com/repos/dessn/Pippin/git/trees/{COMMIT}?recursive=1',folder/'tree.json');t=json.loads((folder/'tree.json').read_text());assert t['sha']==COMMIT and not t.get('truncated')
    byname={v['path']:v for v in t['tree']};records=[]
    for n in FILES:
        r=get(f'https://raw.githubusercontent.com/dessn/Pippin/{COMMIT}/{n}',folder/n.replace('/','_'));data=open(r['path'],'rb').read();digest=hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest();assert digest==byname[n]['sha'];r['git_blob_sha1']=digest;records.append(r)
    result={'code_sha256':sha(__file__),'repository':'dessn/Pippin','commit':COMMIT,'assets':records,'semantics':{'clump_estimator':'OPT_SETPKMJD=16; SNANA(text:key) PKMJDINI','phot_rejection':'PHOTFLAG bits8,16,32,64,128,256,512','redshift_override':'REDSHIFT_FINAL and its error','window_source':'DATAPREP clump_file passed with --photo_window_files'},'limitation':'Pinned publicly available pipeline documents the reconstruction. Exact historical production executable and complete preprocessing output are not presumed identical.'}
    (RESULTS/'processing-source.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
