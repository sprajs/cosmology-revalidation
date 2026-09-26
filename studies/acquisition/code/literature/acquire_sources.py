#!/usr/bin/env python3
"""Pin supplementary primary papers without updating phase1 or overwriting sources."""
from pathlib import Path
import concurrent.futures, subprocess, json, hashlib, re, datetime
ROOT=Path(__file__).resolve().parents[3]
OUT=ROOT/'phase2/literature/sources'; OUT.mkdir(parents=True,exist_ok=True)
IDS=['1811.02379','2012.07180','2307.13696','2111.10382','2201.11142','2402.18690','2104.07795','2112.03864','2301.10644','1710.00845','1710.00846','1506.01354','1610.08972','1912.02191','1912.04257','2408.07175','2501.06664','1807.06209','1610.04677','1507.01602']
def acquire(arxiv):
 result={'arxiv_id':arxiv,'retrieved_at':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':[]}
 def get(url,path):
  if not path.exists():
   tmp=path.with_suffix(path.suffix+'.part')
   subprocess.run(['curl','--fail','-sSL','--max-time','45','--retry','1',url,'-o',str(tmp)],check=True)
   tmp.rename(path)
  b=path.read_bytes();result['files'].append({'url':url,'path':str(path.relative_to(ROOT)),'sha256':hashlib.sha256(b).hexdigest(),'bytes':len(b)})
 try:
  landing=OUT/(arxiv+'.html'); get('https://arxiv.org/abs/'+arxiv,landing)
  t=landing.read_text(); versions=re.findall(r'\[v(\d+)\]',t)
  v=max(map(int,versions)) if versions else int(re.findall(arxiv.replace('.','\\.')+r'v(\d+)',t)[-1])
  result['version']=f'arXiv:{arxiv}v{v}'
  title=re.search(r'<meta name="citation_title" content="([^"]+)"',t); result['title']=title.group(1) if title else ''
  pdf=OUT/f'{arxiv}v{v}.pdf';get(f'https://arxiv.org/pdf/{arxiv}v{v}',pdf)
  txt=pdf.with_suffix('.txt')
  if not txt.exists():subprocess.run(['pdftotext','-layout',str(pdf),str(txt)],check=True)
  result['files'].append({'path':str(txt.relative_to(ROOT)),'sha256':hashlib.sha256(txt.read_bytes()).hexdigest(),'bytes':txt.stat().st_size,'derived_from':str(pdf.relative_to(ROOT))})
 except Exception as e:result['error']=str(e)
 return result
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:results=list(pool.map(acquire,IDS))
 out=OUT/'acquisition.json';out.write_text(json.dumps(results,indent=2)+'\n')
 for r in results:print(r['arxiv_id'],r.get('version'),r.get('error','OK'))
