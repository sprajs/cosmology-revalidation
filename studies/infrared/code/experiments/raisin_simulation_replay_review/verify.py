from pathlib import Path
from decimal import Decimal as D
import gzip,csv,json,hashlib,collections
R=Path.cwd();P=Path(__file__).parent;S=R/'runs/research_2026_09_26/raisin_simulation_execution_review/nominal_fitres'
def read(p,tag):
 with gzip.open(p,'rt') as f:
  out=[]
  for line in f:
   a=line.split()
   if not a:continue
   if a[0]=='VARNAMES:':keys=a[1:]
   if a[0]==tag:assert len(a)==len(keys)+1;out.append(dict(zip(keys,a[1:])))
 return out
f=read(S/'nir.FITRES.gz','SN:');fd={int(r['CID']):r for r in f};p=read(S/'lcplot-feasibility/FITOPT000.LCPLOT.gz','OBS:');groups=collections.defaultdict(list)
for r in p:groups[int(r['CID'])].append(r)
source={};cid=None
for l in (R/'runs/research_2026_09_26/raisin_sign_source/sim/simlibs/DES_RAISIN.simlib').read_text().splitlines():
 a=l.split()
 if not a:continue
 if a[0]=='LIBID:':cid=int(a[1]);source[cid]=collections.defaultdict(list)
 if a[0]=='S:' and a[3] in ['J','H']:source[cid][a[3]].append(D(a[1]))
saved={int(r['plot_CID']):r for r in csv.DictReader((P/'metadata-candidates-500.csv').open())};stats=collections.Counter()
for cid,rr in groups.items():
 clocks={D(r['MJD'])-D(r['Tobs']) for r in rr if r['DATAFLAG']=='0'};assert len(clocks)==1;clock=clocks.pop()
 d=[r for r in rr if r['DATAFLAG']=='1'];times=collections.defaultdict(list)
 for r in d:times[r['BAND']].append(D(r['MJD']))
 # Different, exhaustive one-to-one subset check (permutations; at most six source epochs).
 from itertools import permutations
 found=[]
 for lid,bands in source.items():
  ok=True
  for b,ts in times.items():
   if not any(all(abs(t-v)<=D('.0056') for t,v in zip(ts,seq)) for seq in permutations(bands[b],len(ts))):ok=False;break
  if ok:found.append(lid)
 assert len(found)==1
 candidates=[int(r['CID']) for r in f if int(r['SIM_LIBID'])==found[0] and abs(D(r['PKMJD'])-clock)<=D('.0056')]
 savedids=[int(x) for x in saved[cid]['candidate_CIDs'].split('|')];assert candidates==savedids
 assert int(saved[cid]['identified_SIMLIB'])==found[0]
 stats['objects']+=1;stats['sameCID_clock']+=abs(D(fd[cid]['PKMJD'])-clock)<=D('.0056');stats['sameCID_cadence']+=int(fd[cid]['SIM_LIBID'])==found[0];stats['sameCID_count']+=int(fd[cid]['NDOF'])+1==len(d)
 stats['sameCID_combined']+=cid in candidates;stats['multiple_candidates']+=len(candidates)>1
for name,blob,sha in [('nir.FITRES.gz','7dd16e0978e2be4f09ebdc14698aaeca902df26e','8b801d43ad37200ee4ccaff84241fe4c4a23539aaa20b531cf8570e069056a59'),('lcplot-feasibility/FITOPT000.LCPLOT.gz','d81d64d3f1f395ad9dbbab71837d89f4811320af','6aeb02bab3648daad6435405bbd1c9ab7f8f4db6949bf24ab00cba7296de6f97')]:
 b=(S/name).read_bytes();assert hashlib.sha256(b).hexdigest()==sha;assert hashlib.sha1(f'blob {len(b)}\0'.encode()+b).hexdigest()==blob
res={'pass':True,'method':'Separate stdlib Decimal clock parsing and exhaustive same-band subset permutations; no imports from primary scripts or collaborator code. Source compressed Git blobs independently verified.', 'stats':dict(stats)}
(P/'verification.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
