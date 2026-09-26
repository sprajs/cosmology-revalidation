"""Metadata multiset test, no magnitude/error parsing."""
from pathlib import Path
import tarfile,collections,json,csv,hashlib
from decimal import Decimal
O=Path(__file__).resolve().parent;ROOT=O.parents[4];sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
def main():
 p=json.loads((O/'protocol.json').read_text());archive=ROOT/p['archive'];assert sha(archive)==p['archive_sha256'] and sha(__file__)==p['source_sha256']
 raw=[];snpy=[]
 with tarfile.open(archive) as t:
  for line in t.extractfile('DR3/SN_photo.dat').read().decode().splitlines():
   a=line.split();raw.append((a[0].lower(),a[1],str(Decimal(a[2]).normalize())))
  for member in t:
   if not(member.isfile() and member.name.endswith('_snpy.txt')):continue
   lines=t.extractfile(member).read().decode().splitlines();cid=lines[0].split()[0].lower();band=None
   for line in lines[1:]:
    a=line.split()
    if not a or a[0].startswith('#'):continue
    if a[0]=='filter':band=a[1];continue
    snpy.append((cid,band,str(Decimal(a[0]).normalize())))
 actual=collections.Counter(snpy);identity=collections.Counter(raw);mapping={'Jdw':'J','Hdw':'H'};merged=collections.Counter((cid,mapping.get(b,b),t) for cid,b,t in raw)
 comparisons={}
 for label,predicted in [('identity',identity),('Jdw_to_J_Hdw_to_H',merged)]:
  miss=predicted-actual;extra=actual-predicted;comparisons[label]=dict(exact_multiset_equality=predicted==actual,missing_count=sum(miss.values()),extra_count=sum(extra.values()))
 ledger=[dict(CID=cid[2:],raw_filter=b,SNpy_candidate=mapping.get(b,b),time_token=t,multiplicity=n,merged_target_multiplicity=actual[(cid,mapping.get(b,b),t)]) for (cid,b,t),n in sorted(identity.items()) if b in mapping]
 with(O/'changed-filter-ledger.csv').open('w') as f:w=csv.DictWriter(f,ledger[0].keys());w.writeheader();w.writerows(ledger)
 result=dict(raw_rows=len(raw),snpy_rows=len(snpy),raw_filter_counts=dict(collections.Counter(b for _,b,_ in raw)),SNpy_filter_counts=dict(collections.Counter(b for _,b,_ in snpy)),comparisons=comparisons,
  Jdw_objects=sorted(set(cid[2:] for cid,b,t in raw if b=='Jdw')),Hdw_objects=sorted(set(cid[2:] for cid,b,t in raw if b=='Hdw')),
  ignored_fields='Every magnitude and error token; no flux values converted, compared or scored',conclusion='A metadata mapping identifies labels, not physical passband equivalence or the executed photometric transformation.')
 (O/'result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
