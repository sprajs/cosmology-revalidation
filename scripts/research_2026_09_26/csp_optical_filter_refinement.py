"""Post-universal-map metadata refinement, keeping magnitude fields out of joins."""
from pathlib import Path
from decimal import Decimal as D
from collections import Counter,defaultdict
import csv,json,tarfile,hashlib
ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'runs/research_2026_09_26/csp_optical_magnitude_lineage'
OUT=BASE/'metadata-refinement'
ARC=ROOT/'runs/research_2026_09_26/csp_dr3_provenance/CSP_Photometry_DR3.tgz'
NORM=lambda s:str(D(s).normalize())
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    OUT.mkdir(exist_ok=False)
    (OUT/'design.json').write_text(json.dumps({'status':'Post-outcome refinement after frozen universalV0→V map failed760rows. Resolve per-row labels only from name/time/allowedV0→V orV0 metadata; preserve ambiguous keys, then test magnitude/error.','parent_protocol_sha256':sha(BASE/'protocol.json'),'archive_sha256':sha(ARC),'source_sha256':sha(Path(__file__))},indent=2)+'\n')
    raw=[];actual=[];opt={'u','g','r','i','B','V','V0','V1'}
    with tarfile.open(ARC) as t:
        for line in t.extractfile('DR3/SN_photo.dat').read().decode().splitlines():
            c,b,tt,m,e=line.split()
            if b in opt:raw.append((c.lower(),b,NORM(tt),NORM(m),NORM(e)))
        for mem in t:
            if not(mem.isfile() and mem.name.endswith('_snpy.txt')):continue
            ls=t.extractfile(mem).read().decode().splitlines();c=ls[0].split()[0].lower();b=None
            for line in ls[1:]:
                v=line.split()
                if not v or v[0].startswith('#'):continue
                if v[0]=='filter':b=v[1];continue
                if b in opt:actual.append((c,b,*map(NORM,v)))
    metadata=Counter(r[:3] for r in actual);predicted=[];ledger=[];groups=defaultdict(list)
    for c,b,t,m,e in raw:
        cand=[bb for bb in (['V0','V'] if b=='V0' else [b]) if metadata[(c,bb,t)]>0]
        # Do not decrement to choose among branches based on iteration order.
        if len(cand)==1:predicted.append((c,cand[0],t,m,e))
        if b=='V0':
            row=dict(name=c,raw_band=b,time=t,target_candidates='|'.join(cand),metadata_unique=len(cand)==1)
            ledger.append(row)
            if len(cand)==1:groups[(c,cand[0])].append((c,cand[0],t,m,e))
    pc=Counter(predicted);ac=Counter(actual)
    summary=[]
    for (c,b),rr in sorted(groups.items()):
        counter=Counter(rr);tt=[float(x[2])+53000 for x in rr]
        summary.append(dict(name=c,SNpy_band=b,raw_V0_rows=len(rr),unmatched_magnitude_error=sum((counter-ac).values()),MJD_min=min(tt),MJD_max=max(tt)))
    for n,rows in [('V0-row-metadata.csv',ledger),('V0-object-groups.csv',summary)]:
        with (OUT/n).open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    result=dict(raw_rows=len(raw),actual_rows=len(actual),metadata_unresolved_rows=len(raw)-len(predicted),complete_exact_mag_error_multiset=pc==ac,expected_unmatched=sum((pc-ac).values()),actual_unmatched=sum((ac-pc).values()),V0_targets={b:sum(r['raw_V0_rows'] for r in summary if r['SNpy_band']==b) for b in ['V0','V']},V0_objects_by_target={b:sum(r['SNpy_band']==b for r in summary) for b in ['V0','V']},interpretation='Metadata-specificV0mapping observed exactly; no physical equivalence or earlier correction inferred. Original universal-map failure retained.')
    (OUT/'result.json').write_text(json.dumps(result,indent=2)+'\n');(OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes());(OUT/'manifest.json').write_text(json.dumps({p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'},indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
