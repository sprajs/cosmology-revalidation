"""CID-only metadata gate; does not parse distance, peak or chi-square outcomes."""
from pathlib import Path
import gzip,json,hashlib
R=Path.cwd();O=Path(__file__).resolve().parent;A=R/'runs/research_2026_09_26/astra_design/raisin_timing_assets/author_20211111'
def read_cids(p):
 ids=[]
 with gzip.open(p,'rt') as f:
  for line in f:
   if line.startswith('VARNAMES:'):assert line.split()[1]=='CID'
   elif line.startswith('SN:'):ids.append(line.split(maxsplit=2)[1])
 assert len(ids)==len(set(ids));return ids
n=read_cids(A/'nir.FITRES.gz');j=read_cids(A/'optnir.FITRES.gz');first=[str(i) for i in range(1,501)];missing=sorted(set(n)-set(j),key=int)
d={'status':'Metadata only, fixed future domains; no paired outcomes calculated.','NIR_count':len(n),'joint_count':len(j),'first500_NIR_present':all(c in n for c in first),'first500_joint_missing':[c for c in first if c not in j],'global_joint_missing':missing,'first8_joint_present':all(str(c) in j for c in range(1,9)),'future_expansion_ledger':first,'inputs_sha256':{str((A/name).relative_to(R)):hashlib.sha256((A/name).read_bytes()).hexdigest() for name in ['nir.FITRES.gz','optnir.FITRES.gz']}}
(O/'membership-only.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({k:v for k,v in d.items() if k!='future_expansion_ledger'},indent=2))
