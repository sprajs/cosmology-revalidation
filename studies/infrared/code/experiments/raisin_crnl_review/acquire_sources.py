from pathlib import Path
import urllib.request, json, hashlib, datetime, subprocess
P=Path(__file__).resolve().parent
sources={
 'WFC3-2019-01.pdf':'https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs/_documents/2019/WFC3-2019-01.pdf',
 'WFC3-2010-07.pdf':'https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs/_documents/2010/WFC3-2010-07.pdf',
 'WFC3-ISR-2020-10.pdf':'https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs/_documents/2020/WFC3-ISR-2020-10.pdf',
 'handbook-7-7.html':'https://hst-docs.stsci.edu/wfc3dhb/chapter-7-wfc3-ir-sources-of-error/7-7-count-rate-non-linearity',
 'handbook-9-1.html':'https://hst-docs.stsci.edu/wfc3dhb/chapter-9-wfc3-data-analysis/9-1-photometry',
 'handbook-3-3.html':'https://hst-docs.stsci.edu/wfc3dhb/chapter-3-wfc3-data-calibration/3-3-ir-data-calibration-steps',
 'WFC3-2025-09-arxiv2602.12110v1.html':'https://arxiv.org/html/2602.12110v1',
}
records=[]
for name,url in sources.items():
 out=P/'sources'/name
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=30) as r:
   data=r.read(5_000_001)
   assert len(data)<=5_000_000
   final=r.geturl(); ctype=r.headers.get('Content-Type')
  assert not out.exists()
  out.write_bytes(data)
  row=dict(url=url,final_url=final,path=str(out.relative_to(P)),bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),content_type=ctype)
  if name.endswith('.pdf'):
   txt=out.with_suffix('.txt'); subprocess.run(['pdftotext','-layout',str(out),str(txt)],check=True)
   row['text_sha256']=hashlib.sha256(txt.read_bytes()).hexdigest()
 except Exception as e: row=dict(url=url,error=repr(e))
 records.append(row)
(P/'acquisition.json').write_text(json.dumps(dict(acquired_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),scope='Primary text documentation only; no image pixels or photometry acquired',per_source_cap_bytes=5_000_000,sources=records),indent=2)+'\n')
print(json.dumps(records,indent=2))
