from pathlib import Path
import urllib.request,json,hashlib,datetime
from types import SimpleNamespace
from astropy.io.votable import parse_single_table
R=Path(__file__).resolve().parents[3];O=R/'phase2/assumptions/image_access';O.mkdir(parents=True,exist_ok=True)
targets=[('1246275',54.647026,-26.401205),('1246281',53.725414,-27.622061),('1246314',54.836567,-26.640186)]
records=[]
urls={'sia_docs':'https://datalab.noirlab.edu/docs/manual/UsingAstroDataLab/DataAccessInterfaces/SimpleImageAccessSIA/SimpleImageAccessSIA.html','des_access':'https://www.darkenergysurvey.org/the-des-project/data-access/'}
for cid,ra,dec in targets:
 urls['sia_des_dr2_se_'+cid]=f'https://datalab.noirlab.edu/sia/des_dr2_se?POS={ra},{dec}&SIZE=0.003&VERB=3'
for name,url in urls.items():
 p=O/(name+('.xml' if name.startswith('sia_des') else '.html'))
 try:
  r=urllib.request.urlopen(url,timeout=50);response=SimpleNamespace(content=r.read(),url=r.url,status_code=r.status,ok=r.status==200);p.write_bytes(response.content);rec={'name':name,'url':url,'final_url':response.url,'status':response.status_code,'bytes':len(response.content),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.file_digest(p.open('rb'),'sha256').hexdigest(),'path':str(p.relative_to(R))}
  if name.startswith('sia_des') and response.ok:
   try:
    t=parse_single_table(p).to_table();t.write(O/(name+'.csv'),format='ascii.csv',overwrite=True);rec['rows']=len(t);rec['columns']=t.colnames;print(name,len(t),t.colnames,flush=True)
   except Exception as e:rec['parse_error']=str(e)
  records.append(rec)
 except Exception as e:records.append({'name':name,'url':url,'error':str(e)})
 (O/'access_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps(records,indent=2))
