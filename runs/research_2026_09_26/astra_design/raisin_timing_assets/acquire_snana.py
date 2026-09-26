from pathlib import Path
import urllib.request,json,hashlib,tarfile,io,shutil,datetime,subprocess
P=Path(__file__).parent;tag='f8c25e0de3c68f818406974467c9736d44d65784';commit='10ec91297e4482d593cb5d3d055b10d4aa915071';meta=P/'snana_v11_04d';meta.mkdir(exist_ok=True)
def get(url,cap=30_000_000):
 with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'supernova-provenance-audit'}),timeout=40) as f:b=f.read(cap+1);assert len(b)<=cap;return b,dict(f.headers)
link=subprocess.check_output(['git','ls-remote','https://github.com/RickKessler/SNANA.git','refs/tags/v11_04d*'],text=True);(meta/'ls-remote.txt').write_text(link);assert tag in link and commit in link
b,_=get('https://api.github.com/repos/RickKessler/SNANA/git/tags/'+tag);(meta/'tag.json').write_bytes(b);assert json.loads(b)['object']['sha']==commit
b,_=get('https://api.github.com/repos/RickKessler/SNANA/git/trees/'+commit+'?recursive=1');(meta/'tree.json').write_bytes(b);tree=json.loads(b);assert not tree['truncated']
url='https://codeload.github.com/RickKessler/SNANA/tar.gz/'+commit;b,headers=get(url);(meta/'source.tar.gz').write_bytes(b)
source=meta/'source';source.mkdir(exist_ok=False)
with tarfile.open(fileobj=io.BytesIO(b),mode='r:gz') as t:
 members=t.getmembers();prefix=members[0].name.split('/')[0]+'/'
 for m in members:
  if not m.isfile():continue
  assert m.name.startswith(prefix);relative=m.name[len(prefix):];assert '..' not in Path(relative).parts and not Path(relative).is_absolute()
  dst=source/relative;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_bytes(t.extractfile(m).read());dst.chmod(m.mode)
checks=[]
for r in tree['tree']:
 if r['type']!='blob':continue
 p=source/r['path'];d=p.read_bytes();assert hashlib.sha1(f'blob {len(d)}\0'.encode()+d).hexdigest()==r['sha'];checks.append({'path':r['path'],'git_blob':r['sha'],'sha256':hashlib.sha256(d).hexdigest(),'bytes':len(d)})
assert len(checks)==sum(p.is_file() for p in source.rglob('*'))
shutil.copytree(source,meta/'build')
out={'tag':'v11_04d','tag_object':tag,'commit':commit,'url':url,'source_archive_bytes':len(b),'source_archive_sha256':hashlib.sha256(b).hexdigest(),'files':checks,'verified_files':len(checks),'source_selection':'coherent November2021 archived FITRES headers; no fitting or outcome-based source selection','utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};(meta/'acquisition.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['files']},indent=2))
