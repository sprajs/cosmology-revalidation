from pathlib import Path
import urllib.request,hashlib,json
from html.parser import HTMLParser
P=Path(__file__).parent
class Text(HTMLParser):
 def __init__(self):super().__init__();self.parts=[]
 def handle_data(self,x):
  if x.strip():self.parts.append(x.strip())
rows=[]
for name,url in [('paper','https://arxiv.org/html/1709.05146v2'),('filters','https://csp.obs.carnegiescience.edu/data/filters')]:
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=45) as r:b=r.read(12000000);h=dict(r.headers)
  assert len(b)<12000000
  (P/(name+'.html')).write_bytes(b);p=Text();p.feed(b.decode());(P/(name+'.txt')).write_text('\n'.join(p.parts)+'\n');rows.append({'url':url,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest(),'headers':h})
 except Exception as e:rows.append({'url':url,'error':str(e)})
(P/'acquisition.json').write_text(json.dumps(rows,indent=2)+'\n');print(json.dumps(rows,indent=2))
