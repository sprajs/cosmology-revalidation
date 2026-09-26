from pathlib import Path
import urllib.request,json,hashlib,zipfile,struct,zlib,time
P=Path(__file__).parent
amend={'scope':'Resource-only amendment before full spectra acquisition; no scientific scores inspected','original_download_cap_bytes':30000000,'new_spectra_download_cap_bytes':40000000,'authorization':'Parent explicitly authorized the 31,173,962-byte primary spectra archive, preserving HEAD and range metadata. The 153,724,421-byte template stays metadata-only.','spectra_url':'https://csp.obs.carnegiescience.edu/data/CSPII_NIR_Ia_spectra.zip'}
(P/'resource-amendment.json').write_text(json.dumps(amend,indent=2)+'\n')
url=amend['spectra_url']; out=P/'CSPII_NIR_Ia_spectra.zip'
if not out.exists():
    req=urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'})
    with urllib.request.urlopen(req,timeout=60) as r, out.with_suffix('.partial').open('wb') as f:
        n=0
        while chunk:=r.read(1024*1024):
            n+=len(chunk)
            if n>40000000: raise RuntimeError('Resource cap exceeded')
            f.write(chunk)
        headers=dict(r.headers); status=r.status
    out.with_suffix('.partial').rename(out)
else: headers={};status='already present'
rows=[]
with zipfile.ZipFile(out) as z:
    for i in z.infolist():
        b=z.read(i)
        rows.append({'name':i.filename,'size':i.file_size,'compressed':i.compress_size,'CRC':i.CRC,'sha256':hashlib.sha256(b).hexdigest()})
        # Extract only within owned path, no symlink or parent path escape.
        rel=Path(i.filename)
        if rel.is_absolute() or '..' in rel.parts: raise RuntimeError(i.filename)
        if i.is_dir():continue
        p=P/'spectra'/rel;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
(P/'spectra-acquisition.json').write_text(json.dumps({'url':url,'status':status,'headers':headers,'bytes':out.stat().st_size,'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'members':rows},indent=2)+'\n')
# Only the small first local member of the deferred template archive.
u='https://csp.obs.carnegiescience.edu/data/NIR_Ia_template.zip'
req=urllib.request.Request(u,headers={'User-Agent':'Mozilla/5.0','Range':'bytes=0-767'})
with urllib.request.urlopen(req,timeout=60) as r:
    b=r.read(4096); h=dict(r.headers);status=r.status
if status!=206 or len(b)!=768:raise RuntimeError((status,len(b)))
(P/'template-readme-local-range').write_bytes(b)
a=struct.unpack('<4s5H3I2H',b[:30]); sig,version,flags,method,mt,md,crc,cs,us,nl,xl=a
assert sig==b'PK\x03\x04' and method==8
name=b[30:30+nl].decode();dat=zlib.decompress(b[30+nl+xl:30+nl+xl+cs],-15)
assert zlib.crc32(dat)==crc and len(dat)==us
(P/name).write_bytes(dat)
(P/'template-readme-acquisition.json').write_text(json.dumps({'url':u,'status':status,'headers':h,'range_bytes':len(b),'member':name,'bytes':len(dat),'sha256':hashlib.sha256(dat).hexdigest(),'archive_downloaded':False},indent=2)+'\n')
print(json.dumps({'spectra_bytes':out.stat().st_size,'members':len(rows),'spectra_sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'template_readme':dat.decode()},indent=2))
