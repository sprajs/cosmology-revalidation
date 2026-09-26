from pathlib import Path
import json,re,csv
R=Path.cwd();P=Path(__file__).parent
sources={'author':'runs/research_2026_09_26/raisin_differential/mass-threshold/code-tree.json','release':'sources/updates/2026-09-26-raisin/tree.json'}
summary={};hits=[]
for repo,path in sources.items():
 d=json.loads((R/path).read_text());assert d['truncated'] is False
 tt=d['tree'];files=[r for r in tt if r['type']=='blob']
 flagged=[]
 for r in files:
  n=r['path'];low=n.lower();why=[]
  if re.search(r'(^|[/_.])(head|phot)([/_.]|$)',low):why.append('head_phot_token')
  if ('sim' in low and low.endswith(('.fits','.fits.gz','.fit','.fit.gz','.tar','.tgz','.zip','.tar.gz'))):why.append('simulation_binary_or_archive')
  if ('des_raisin_sim' in low and not low.endswith(('.fitres','.fitres.gz','.lcplot','.lcplot.gz','.pdf','.png','.jpg','.svg','.root','.hbook'))):why.append('target_version_noncompressedtable')
  if re.search(r'(^|/)(readme[^/]*|\.gitignore|\.gitattributes)$',low):why.append('repository_documentation')
  if why:flagged.append({'repo':repo,'path':n,'size':r.get('size'),'blob':r['sha'],'reason':'|'.join(why)})
 hits+=flagged;summary[repo]={'commit':d['sha'],'tree_entries':len(tt),'blobs':len(files),'truncated':d['truncated'],'flagged':len(flagged),'head_phot_token':sum('head_phot_token' in r['reason'] for r in flagged),'simulation_binary_or_archive':sum('simulation_binary_or_archive' in r['reason'] for r in flagged)}
with (P/'tree-candidates.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(hits[0]));w.writeheader();w.writerows(hits)
(P/'tree-summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
for r in hits:
 if 'head_phot_token' in r['reason'] or 'simulation_binary_or_archive' in r['reason']:print(r)
