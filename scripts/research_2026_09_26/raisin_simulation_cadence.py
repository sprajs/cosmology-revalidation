#!/usr/bin/env python3
"""Does the archived RAISIN simulation cadence include omitted signed epochs?

Source/cadence accounting only; this does not execute a bias correction.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'runs/research_2026_09_26/raisin_sign_source'
LEDGER=ROOT/'runs/research_2026_09_26/raisin_flux_sign_audit/original-retention.csv'
SIM=OUT/'sim/simlibs/DES_RAISIN.simlib'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def freeze():
    p=OUT/'cadence-protocol.json'
    if p.exists():raise FileExistsError(p)
    b=LEDGER.read_bytes();(OUT/'frozen-original-retention.csv').write_bytes(b)
    write(p,dict(created_utc=datetime.now(timezone.utc).isoformat(),
        stage='After data sign-retention finding, before simulation cadence matching.',
        source_commit='b888214a5cbae38ac0bf886488ce7734f5b77e87',
        rule='For every original DIFFIMG row in the frozen17-object ledger, match same RAISIN FIELD and band in archived DES_RAISIN.simlib with abs(MJD difference)<=.00055day. Retain0,1,andmultiplematches. No flux rescaling, weights or new selection.',
        summaries='Counts by original flux sign, RAISIN retention and SIMLIB presence, allrows and previously frozen restphase[-7,45],M20[-10,40], andoffseason absobserverday>180. These are provenance counts, not physical noise outcomes.',
        interpretation='A cadence derived from retained observations does not itself impose a positive-flux cut on newly generated noise. Exact generator/fitter execution and any preprocessing must be reproduced before deciding whether published bias corrections account for this effect.',
        inputs_sha256={str(p.relative_to(ROOT)):sha(p) for p in [LEDGER,SIM,Path(__file__).resolve()]},
        frozen_ledger_sha256=hashlib.sha256(b).hexdigest()))
    print(json.dumps({'protocol_sha256':sha(p)}))

def score():
    p=OUT/'cadence-protocol.json';q=json.loads(p.read_text())
    assert sha(OUT/'frozen-original-retention.csv')==q['frozen_ledger_sha256']
    for k,h in q['inputs_sha256'].items():
        if ROOT/k!=LEDGER:assert sha(ROOT/k)==h,k
    if (OUT/'cadence-result.json').exists():raise FileExistsError('Keep original result')
    rows=[];lib=None;field=None
    for lineno,line in enumerate(SIM.read_text().splitlines(),1):
        s=line.split()
        if not s:continue
        if s[0]=='LIBID:':lib=int(s[1]);field=None
        elif s[0]=='FIELD:':field=s[1]
        elif s[0]=='S:':
            assert field is not None
            rows.append(dict(LIBID=lib,field=field,MJD=float(s[1]),band=s[3],source_line=lineno))
    sim=pd.DataFrame(rows)
    data=pd.read_csv(OUT/'frozen-original-retention.csv',dtype={'raisin_CID':str,'original_CID':str})
    assert len(data)==9843 and data.raisin_CID.nunique()==17
    assert set(data.raisin_CID)<=set(sim.field)
    results=[]
    for (cid,band),group in data.groupby(['raisin_CID','band'],sort=False):
        ss=sim[(sim.field==cid)&(sim.band==band)]
        for idx,row in group.iterrows():
            m=ss[abs(ss.MJD-row.MJD)<=.00055]
            results.append(dict(source_ledger_row=idx,simlib_matches=len(m),
                simlib_line=';'.join(str(i)for i in m.source_line),simlib_LIBID=';'.join(str(i)for i in m.LIBID)))
    match=pd.DataFrame(results).set_index('source_ledger_row').sort_index()
    assert np.array_equal(match.index,np.arange(len(data)))
    data=pd.concat([data,match],axis=1);data['simlib_present']=data.simlib_matches>0
    counts=[]
    for label,mask in [('all',np.ones(len(data),dtype=bool)),('optical_window',data.original_optical_fit_window),
                       ('m20_window',data.m20_window),('offseason',data.offseason)]:
        sub=data.loc[mask]
        for sign in sorted(data.original_sign.unique()):
            ss=sub[sub.original_sign==sign]
            counts.append(dict(window=label,sign=sign,original_rows=len(ss),raisin_retained=int(ss.retained.sum()),
                simlib_present=int(ss.simlib_present.sum()),simlib_unique=int((ss.simlib_matches==1).sum()),
                simlib_ambiguous=int((ss.simlib_matches>1).sum()),
                missing_raisin_but_simlib=int((~ss.retained&ss.simlib_present).sum()),
                in_raisin_missing_simlib=int((ss.retained&~ss.simlib_present).sum())))
    data.to_csv(OUT/'cadence-matches.csv',index=False)
    pd.DataFrame(counts).to_csv(OUT/'cadence-counts.csv',index=False)
    write(OUT/'cadence-result.json',dict(protocol_sha256=sha(p),source_ledger_unchanged_now=sha(LEDGER)==q['frozen_ledger_sha256'],
        simlib_blocks=sim.LIBID.nunique(),simlib_rows=len(sim),matched_objects=data.raisin_CID.nunique(),counts=counts,
        no_bias_correction_or_fit_executed=True,
        outputs_sha256={p.name:sha(p)for p in [OUT/'cadence-matches.csv',OUT/'cadence-counts.csv']}))
    (OUT/'cadence-executed-source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(counts,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','score']);a=p.parse_args();globals()[a.mode]()
