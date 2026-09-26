"""Resolve NIR physical-filter metadata after discovering archive label collapse."""
from pathlib import Path
import json,csv,tarfile,collections,sys,hashlib
from decimal import Decimal
O=Path(__file__).resolve().parent;BASE=O.parent;ROOT=O.parents[4];sys.path.insert(0,str(BASE));from inventory import table
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
def main():
 p=json.loads((O/'refinement-protocol.json').read_text());assert sha(__file__)==p['source_sha256'];a=ROOT/p['archive'];assert sha(a)==p['archive_sha256'];raw=[];snpy=[]
 with tarfile.open(a) as t:
  for line in t.extractfile('DR3/SN_photo.dat').read().decode().splitlines():
   x=line.split();raw.append((x[0].lower(),x[1],str(Decimal(x[2]).normalize())))
  for m in t:
   if not(m.isfile() and m.name.endswith('_snpy.txt')):continue
   lines=t.extractfile(m).read().decode().splitlines();cid=lines[0].split()[0].lower();b=None
   for l in lines[1:]:
    x=l.split()
    if not x or x[0].startswith('#'):continue
    if x[0]=='filter':b=x[1];continue
    snpy.append((cid,b,str(Decimal(x[0]).normalize())))
 mapping={'Jdw':'J','Hdw':'H'};pred=collections.Counter((cid,mapping.get(b,b),t) for cid,b,t in raw);actual=collections.Counter(snpy)
 missing=pred-actual;extra=actual-pred
 nir={'Y','Ydw','J','Jrc2','Jdw','H','Hdw'};prednir=collections.Counter((cid,mapping.get(b,b),t) for cid,b,t in raw if b in nir);actnir=collections.Counter(x for x in snpy if x[1] in nir)
 assert prednir==actnir
 DATA=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b/photometry/RAISIN/CSPDR3_RAISIN';refs={}
 for name in (DATA/'CSPDR3_RAISIN.LIST').read_text().splitlines():
  h,r=table(DATA/name);refs['sn'+h['SNID'][0].lower()]=h
 good={r['CID'].lower() for r in csv.DictReader((BASE/'dr3_extension/objects.csv').open()) if r.get('identity_gate')=='True'}
 rr=[]
 for cid in sorted(set(x[0] for x in raw)):
  if cid[2:] not in good:continue
  h=refs[cid];peak=float(h['PEAKMJD'][0]);z=float(h['REDSHIFT_HELIO'][0])
  for b in sorted(nir):
   times=[float(Decimal(t))+53000 for c,bb,t in raw if c==cid and bb==b];row=dict(CID=h['SNID'][0],physical_raw_filter=b,epochs=len(times))
   for label,lo,hi in [('early',-7,7),('late',10,20)]:
    ts=[t for t in times if lo<=(t-peak)/(1+z)<=hi];ps=[(t-peak)/(1+z) for t in ts]
    row[label+'_nights']=len(set(int(t) for t in ts));row[label+'_span']=max(ps)-min(ps) if ps else 0.
   row['eligible']=all(row[l+'_nights']>=3 and row[l+'_span']>=2 for l in ['early','late']);rr.append(row)
 with(O/'physical-filter-cadence.csv').open('w') as f:w=csv.DictWriter(f,rr[0].keys());w.writeheader();w.writerows(rr)
 yy=[r['CID'] for r in rr if r['physical_raw_filter']=='Y' and r['eligible']];paired=[c for c in yy if any(r['CID']==c and r['physical_raw_filter'] in ['J','Jrc2','Jdw','H','Hdw'] and r['eligible'] for r in rr)]
 result=dict(exact_NIR_metadata_multiset_mapping=True,NIR_rows=sum(prednir.values()),Jdw_to_J_rows=sum(b=='Jdw' for _,b,_ in raw),Hdw_to_H_rows=sum(b=='Hdw' for _,b,_ in raw),remaining_global_mismatch_source_labels=dict(collections.Counter({b:sum(n for(c,bb,t),n in missing.items() if bb==b) for b in set(k[1] for k in missing)})),remaining_global_mismatch_target_labels=dict(collections.Counter({b:sum(n for(c,bb,t),n in extra.items() if bb==b) for b in set(k[1] for k in extra)})),
   physical_Y_eligible=yy,physical_Y_second_NIR_eligible=paired,physical_exact_filter_eligible_counts={b:sum(r['physical_raw_filter']==b and r['eligible'] for r in rr) for b in sorted(nir)},
   interpretation='NIR metadata label collapse is exact; no equality of transmission curves, SN synthetic magnitudes, or executed correction asserted. Rawphysical filter labels are retained for design.')
 (O/'refinement-result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
