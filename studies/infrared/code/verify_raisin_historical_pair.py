"""Root metadata-only independent check of source-selected historical pairing.

No fitted distance, chi-square, or luminosity residual is used to pair events.
"""
from pathlib import Path
from collections import defaultdict, Counter
from decimal import Decimal
import gzip, hashlib, json

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'runs/research_2026_09_26'
SRC = BASE/'astra_design/raisin_timing_assets'
OUT = BASE/'raisin_historical_pair_root_review'
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()

def table(path, marker):
    with gzip.open(path, 'rt') as f:
        rows=[]
        for line in f:
            t=line.split()
            if not t: continue
            if t[0]=='VARNAMES:': keys=t[1:]
            elif t[0]==marker:
                assert len(t)==len(keys)+1
                rows.append(dict(zip(keys,t[1:])))
    return rows

def main():
    OUT.mkdir(exist_ok=False)
    acquisition=json.loads((SRC/'historical-tables-acquisition.json').read_text())
    for item in acquisition:
        p=Path(item['path']); raw=p.read_bytes()
        assert sha(p)==item['sha256']
        assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()==item['git_blob']
    nir=table(SRC/'author_20211111/nir.FITRES.gz','SN:')
    joint=table(SRC/'author_20211111/optnir.FITRES.gz','SN:')
    nd={r['CID']:r for r in nir}; jd={r['CID']:r for r in joint}
    assert len(nd)==len(nir)==30000 and len(jd)==len(joint)==29995
    keys=[k for k in nir[0] if k.startswith('SIM_') and k in joint[0]]+['PKMJDINI','zHEL','zCMB','zHD','MWEBV']
    mismatch={k:sum(nd[c][k]!=jd[c][k] for c in nd.keys()&jd.keys()) for k in keys}
    assert not any(mismatch.values())
    plotpath=BASE/'raisin_simulation_execution_review/nominal_fitres/lcplot-feasibility/FITOPT000.LCPLOT.gz'
    raw=plotpath.read_bytes();assert hashlib.sha1(f'blob {len(raw)}\0'.encode()+raw).hexdigest()=='d81d64d3f1f395ad9dbbab71837d89f4811320af'
    groups=defaultdict(list)
    for r in table(plotpath,'OBS:'):groups[r['CID']].append(r)
    cadencepath=BASE/'raisin_sign_source/sim/simlibs/DES_RAISIN.simlib'
    cadence={}
    for line in cadencepath.read_text().splitlines():
        t=line.split()
        if not t:continue
        if t[0]=='LIBID:':lid=t[1];cadence[lid]=defaultdict(list)
        if t[0]=='S:' and t[3] in ['J','H']:cadence[lid][t[3]].append(Decimal(t[1]))
    # Same-CID cadence is checked directly; no nearest outcome or CID remapping.
    results=[]
    for cid, rs in groups.items():
        times={Decimal(r['MJD'])-Decimal(r['Tobs']) for r in rs if r['DATAFLAG']=='0'}
        assert len(times)==1
        clock=times.pop();delta=abs(clock-Decimal(nd[cid]['PKMJD']))
        observations=[r for r in rs if r['DATAFLAG']=='1']
        expected=cadence[nd[cid]['SIM_LIBID']]
        used=defaultdict(set)
        for r in observations:
            # Demand unique per-band matches and retain duplicate-epoch multiplicity.
            candidates=[i for i,t in enumerate(expected[r['BAND']]) if i not in used[r['BAND']] and abs(t-Decimal(r['MJD']))<=Decimal('.0056')]
            assert len(candidates)==1,(cid,r['BAND'],candidates)
            used[r['BAND']].add(candidates[0])
        assert delta<=Decimal('.0056') and len(observations)==int(nd[cid]['NDOF'])+1
        results.append(dict(CID=cid,clock_gap_day=str(delta),observed_rows=len(observations),SIM_LIBID=nd[cid]['SIM_LIBID']))
    assert len(results)==500
    answer=dict(pass_metadata_pairing=True,scope='Original HEAD/PHOT, rounding and exact native execution remain separate gates.',
                method='Independent stdlib token/Decimal parse and direct same-CID one-to-one cadence match; no collaborator imports.',
                historical_commit='aaa709ead7a7339d56a2a4604d78991329d9368f',
                common_rows=len(nd.keys()&jd.keys()),metadata_mismatches=mismatch,
                coherent_plot_objects=len(results),max_clock_gap_day=str(max(Decimal(r['clock_gap_day']) for r in results)),objects=results)
    (OUT/'result.json').write_text(json.dumps(answer,indent=2)+'\n')
    (OUT/'executed-source.py').write_bytes(Path(__file__).read_bytes())
    paths=[plotpath,cadencepath,SRC/'historical-tables-acquisition.json']+[Path(i['path']) for i in acquisition]
    (OUT/'manifest.json').write_text(json.dumps(dict(inputs={str(p.relative_to(ROOT)):sha(p) for p in paths},outputs={p.name:sha(p) for p in OUT.iterdir() if p.is_file() and p.name!='manifest.json'}),indent=2)+'\n')
    print(json.dumps({k:v for k,v in answer.items() if k!='objects'},indent=2))

if __name__=='__main__':main()
