from pathlib import Path
import urllib.request,json,hashlib,datetime
P=Path(__file__).parent;O=P/'author_20211111';O.mkdir(exist_ok=True)
commit='aaa709ead7a7339d56a2a4604d78991329d9368f';tree=json.loads((P/'public_metadata/author-20211111-tree.json').read_text());by={r['path']:r for r in tree['tree']}
paths=[('nir.FITRES.gz','output/fit_nir/DES_RAISIN_NIR_SIM/DES_RAISIN_SIM/FITOPT000.FITRES.gz'),('optnir.FITRES.gz','output/fit_all/DES_RAISIN_OPTNIR_SIM/DES_RAISIN_SIM/FITOPT000.FITRES.gz')]
rows=[]
for name,path in paths:
 meta=by[path];assert meta['size']<=30_000_000;url=f'https://raw.githubusercontent.com/djones1040/RAISIN_cosmo/{commit}/{path}'
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'supernova-provenance-audit'}),timeout=45) as f:b=f.read(30_000_001);headers=dict(f.headers)
 assert len(b)==meta['size'];assert hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()==meta['sha'];(O/name).write_bytes(b)
 rows.append({'commit':commit,'source_path':path,'path':str(O/name),'url':url,'size':len(b),'git_blob':meta['sha'],'sha256':hashlib.sha256(b).hexdigest(),'headers':headers,'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
(P/'historical-tables-acquisition.json').write_text(json.dumps(rows,indent=2)+'\n');print([{k:v for k,v in r.items() if k in ['path','size','git_blob','sha256']} for r in rows])
