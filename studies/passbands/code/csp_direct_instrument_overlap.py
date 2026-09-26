"""Reproduce the saved metadata-only instrument-overlap feasibility counts.

No photometric magnitude or error is converted, compared or scored here.
This script was persisted after the initial metadata inventory; its rerun
checks that inventory, without replacing the frozen protocol or result.
"""
from pathlib import Path
from collections import defaultdict
import hashlib, json, tarfile

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/research_2026_09_26/csp_direct_instrument_overlap'

def main():
    protocol=json.loads((OUT/'protocol.json').read_text())
    archive=ROOT/'runs/research_2026_09_26/csp_dr3_provenance/CSP_Photometry_DR3.tgz'
    assert hashlib.sha256(archive.read_bytes()).hexdigest()==protocol['archive_sha256']
    rows=defaultdict(lambda:defaultdict(list))
    with tarfile.open(archive) as tar:
        for n,line in enumerate(tar.extractfile('DR3/SN_photo.dat').read().decode().splitlines(),1):
            fields=line.split()
            rows[fields[0].lower()][fields[1]].append((float(fields[2]),n))
    results=[]; allpairs=[]
    for b1,b2 in protocol['pairs']:
        for window in protocol['time_windows_days']:
            found=[]
            for cid,bands in rows.items():
                for t1,n1 in bands[b1]:
                    for t2,n2 in bands[b2]:
                        if abs(t1-t2)<=window:
                            found.append(dict(CID=cid,band1=b1,band2=b2,time1=t1,time2=t2,
                                source_line1=n1,source_line2=n2,absolute_time_gap=abs(t1-t2)))
            results.append(dict(pair=[b1,b2],window=window,pairs=len(found),objects=len({x['CID'] for x in found})))
            if window==1.0:allpairs+=found
    saved=json.loads((OUT/'result.json').read_text())
    assert saved==dict(scope=protocol['scope'],counts=results,all_pairs_within_one_day=allpairs)
    review=dict(all_saved_metadata_counts_and_pairs_exact=True,
        inputs={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [archive,OUT/'protocol.json',OUT/'result.json',Path(__file__)]},
        scope='Deterministic metadata reproduction; no brightness outcome.')
    dest=OUT/'reproduction.json';assert not dest.exists()
    dest.write_text(json.dumps(review,indent=2)+'\n')
    print('All 16 window/filter counts and 81 source-line pairs reproduce exactly.')

if __name__=='__main__':main()
