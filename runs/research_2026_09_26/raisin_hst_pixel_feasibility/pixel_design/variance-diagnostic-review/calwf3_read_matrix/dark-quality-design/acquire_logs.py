from pathlib import Path
import csv,json,hashlib,sys,time,urllib.parse
import urllib.request
R=Path.cwd();O=Path(__file__).resolve().parent
H=R/'runs/research_2026_09_26/raisin_hst_pixel_feasibility'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
if sys.argv[1]=='prepare':
 roots=list(csv.DictReader((H/'quality_reference_audit/candidate-quality-ledger.csv').open()))
 remaining={r['root']:r for r in roots if r['root']!='idp247tnq'}
 assert len(remaining)==13
 records={};used={}
 for q in sorted((H/'dark_ramp_metadata/queries').glob('*products.json')):
  d=json.loads(q.read_text())
  for x in d.get('data',[]):
   name=x.get('productFilename','');root=name.removesuffix('_log.txt')
   if name.endswith('_log.txt') and root in remaining:
    assert x['description']=='DADS LOG file' and x['dataRights']=='PUBLIC' and x['productSubGroupDescription']=='LOG'
    v={'root':root,'filename':name,'uri':x['dataURI'],'size':int(x['size']),'visit':remaining[root]['visit'],'EXPFLAG':remaining[root]['EXPFLAG']}
    if root in records:assert records[root]==v
    records[root]=v;used[str(q.relative_to(R))]=sha(q)
 assert set(records)==set(remaining)
 rec=sorted(records.values(),key=lambda x:(x['visit'],x['root']))
 p={'status':'frozen before network','scope':'metadata-only listed DADS text logs; no arrays or image products','authority':'root explicit at most13 logs,10MB,90s','max_files':13,'max_bytes':10000000,'max_seconds':90,'per_request_max_seconds':12,'selection':'All14 frozen candidates except already acquired idp247tnq, sorted visit/root; no replacements on failure','records':rec,'input_hashes':used|{'runs/research_2026_09_26/raisin_hst_pixel_feasibility/quality_reference_audit/candidate-quality-ledger.csv':sha(H/'quality_reference_audit/candidate-quality-ledger.csv')},'script_sha256':sha(Path(__file__))}
 (O/'protocol.json').write_text(json.dumps(p,indent=2)+'\n')
 print(json.dumps({'protocol_sha256':sha(O/'protocol.json'),'files':len(rec),'listed_bytes':sum(r['size'] for r in rec)}))
elif sys.argv[1]=='acquire':
 p=json.loads((O/'protocol.json').read_text());assert sha(Path(__file__))==p['script_sha256']
 for name,s in p['input_hashes'].items():assert sha(R/name)==s
 assert not (O/'acquisition.json').exists()
 (O/'logs').mkdir(exist_ok=True);start=time.monotonic();total=0;out=[]
 for rec in p['records']:
  left=p['max_seconds']-(time.monotonic()-start)
  if left<1:break
  if total+rec['size']>p['max_bytes']:break
  entry=dict(rec);url='https://mast.stsci.edu/api/v0.1/Download/file?'+urllib.parse.urlencode({'uri':rec['uri']});entry['url']=url
  try:
   resp=urllib.request.urlopen(url,timeout=min(12,left))
   entry['status_code']=resp.status;raw=bytearray()
   while True:
    block=resp.read(16384)
    if not block:break
    if time.monotonic()-start>p['max_seconds']:raise TimeoutError('frozen wall cap')
    if total+len(raw)+len(block)>p['max_bytes']:raise RuntimeError('frozen byte cap')
    raw.extend(block)
   total+=len(raw);entry['bytes']=len(raw)
   entry['expected_size_match']=len(raw)==rec['size']
   dest=O/'logs'/rec['filename'];dest.write_bytes(raw);entry['sha256']=sha(dest)
   entry['utf8_valid']=True;raw.decode('utf8');entry['status']='saved' if entry['expected_size_match'] else 'size_mismatch'
  except Exception as e:entry['status']='failure';entry['error']=str(e)
  finally:
   if 'resp' in locals():resp.close()
  out.append(entry)
  (O/'acquisition.json').write_text(json.dumps({'records':out,'elapsed_seconds':time.monotonic()-start,'total_bytes':total,'protocol_sha256':sha(O/'protocol.json'),'unattempted':[x['root'] for x in p['records'][len(out):]],'no_image_arrays':True},indent=2)+'\n')
 print(json.dumps({'saved':sum(r['status']=='saved' for r in out),'attempted':len(out),'bytes':total,'seconds':time.monotonic()-start}))
else:raise RuntimeError('unknown command')
