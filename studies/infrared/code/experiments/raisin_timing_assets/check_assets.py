from pathlib import Path
import json,hashlib,gzip,re,collections
R=Path.cwd();P=Path(__file__).parent
old=json.loads((P/'public_metadata/author-20211111-tree.json').read_text());new=json.loads((R/'runs/research_2026_09_26/raisin_differential/mass-threshold/code-tree.json').read_text());oldby={r['path']:r for r in old['tree']};newby={r['path']:r for r in new['tree']}
paths=['fit/sim/REFAC_DES_RAISIN_NIR.nml','fit/sim/REFAC_DES_RAISIN_optnir.nml','snoopy.B18/SNooPy_test.fits','snoopy.B18/snoopy.info','kcor/kcor_DES_NIR.fits','sim/inputs/DES/sim_DES_SNOOPY.input','sim/simlibs/DES_RAISIN.simlib']
compare=[]
for n in paths:
 a,b=oldby[n],newby[n];compare.append({'source_path':n,'historical_blob':a['sha'],'latest_blob':b['sha'],'equal':a['sha']==b['sha'],'bytes':a.get('size')})
headers={}
for n in ['nir','optnir']:
 lines=[]
 with gzip.open(P/f'author_20211111/{n}.FITRES.gz','rt') as f:
  for line in f:
   if line.startswith('VARNAMES:'):break
   lines.append(line.rstrip())
 headers[n]=lines
summary={}
for label,d in [('author2021',old),('authorlatest',new),('releasev1p0',json.loads((P/'public_metadata/release-v1p0-tree.json').read_text())),('releasev1p1',json.loads((R/'sources/updates/2026-09-26-raisin/tree.json').read_text()))]:
 assert d['truncated'] is False
 raw=[x for x in d['tree'] if x['type']=='blob' and re.search(r'(^|[/_.])(HEAD|PHOT)(\.|$)',x['path'])]
 summary[label]={'commit':d['sha'],'entries':len(d['tree']),'blobs':sum(x['type']=='blob' for x in d['tree']),'HEAD_PHOT_candidates':raw}
res={'source_asset_comparison':compare,'table_headers':headers,'trees':summary,'version_amendment_protocol_sha256':hashlib.sha256((P/'version-amendment-protocol.json').read_bytes()).hexdigest()}
(P/'asset-check-result.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps({'allseven_source_assets_unchanged':all(r['equal'] for r in compare),'headers':headers,'trees':summary},indent=2))
