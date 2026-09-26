from pathlib import Path
import re,csv,tarfile,collections,json
from decimal import Decimal
P=Path(__file__).parent;R=Path.cwd();AR=R/'runs/research_2026_09_26/csp_dr3_provenance/CSP_Photometry_DR3.tgz'
text=(P/'Krisciunas2017.txt').read_text();rows=[]
for line in text.splitlines():
 m=re.match(r'^2006 Jan\s+(\d+)\s+([\d,]+\.\d+)\s+SN\s+(\w+)',line)
 if m:
  day,jd,name=m.groups();rows.append({'CID':name.lower(),'day_UT':int(day),'JD':jd.replace(',',''),'MJD':str(Decimal(jd.replace(',',''))-Decimal('2400000.5'))})
assert len(rows)==21
with tarfile.open(AR) as t:
 raw=[line.split() for line in t.extractfile('DR3/SN_photo.dat').read().decode().splitlines()]
 actual=[]
 for m in t:
  if not(m.isfile() and m.name.endswith('_snpy.txt')):continue
  lines=t.extractfile(m).read().decode().splitlines();cid=lines[0].split()[0].lower().removeprefix('sn');band=None
  for line in lines[1:]:
   a=line.split()
   if not a or a[0].startswith('#'):continue
   if a[0]=='filter':band=a[1];continue
   actual.append((cid,band,Decimal(a[0])+53000))
 rawmeta=[(a[0].lower().removeprefix('sn'),a[1],Decimal(a[2])+53000) for a in raw]
 for r in rows:
  dt=Decimal(r['MJD']);r['raw_V_candidates_0p01d']=sum(n==r['CID'] and b in ['V','V0','V1'] and abs(v-dt)<=Decimal('.01') for n,b,v in rawmeta)
  r['SNpy_V_candidates_0p01d']=sum(n==r['CID'] and b in ['V','V0','V1'] and abs(v-dt)<=Decimal('.01') for n,b,v in actual)
with (P/'table9-missing-metadata.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
res={'paper_table9_rows':len(rows),'objects':len({r['CID'] for r in rows}),'raw_V_matched':sum(r['raw_V_candidates_0p01d']>0 for r in rows),'SNpy_V_matched':sum(r['SNpy_V_candidates_0p01d']>0 for r in rows),'scope':'Source completeness only; dates/names parsed, not paper magnitudes/errors. No restoration, native fit or likelihood claim.'}
(P/'table9-result.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
