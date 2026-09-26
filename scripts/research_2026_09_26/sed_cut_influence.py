#!/usr/bin/env python3
"""Post-cut-audit fixed-ID influence check; primary results are untouched."""
from pathlib import Path
from datetime import datetime, timezone
import argparse, hashlib, json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
RUN=ROOT/'runs/research_2026_09_26'
OUT=RUN/'sed_cut_influence'
P=[RUN/'sed_quality_cut_probe/per-object.csv',
   RUN/'sed_nonlinear_validation/resolved/paired-responses.csv',
   RUN/'sed_nonlinear_validation/contrast-membership.csv',
   RUN/'astra_design/validation1020/analysis/object-scores.csv']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def freeze():
    OUT.mkdir(exist_ok=True)
    q=OUT/'protocol.json'
    if q.exists():raise FileExistsError(q)
    cuts=pd.read_csv(P[0],dtype={'CID':str})
    keep=cuts[['native_joint_pass','nominal_observed_joint_pass','observer_joint_pass','sed_joint_pass']].all(axis=1)
    ids=cuts.loc[~keep,'CID'].tolist()
    write(q,dict(created_utc=datetime.now(timezone.utc).isoformat(),stage='Follow-up after coordinate-cut outcomes; before this fixed-ID influence score.',
        removed_ids=ids,rule='Remove union of failures of native, unperturbed observed, observer and SED colour/shape eligibility from BOTH arms; no new threshold search.',
        estimands='Sum original published-mask M,I,G on retained objects without refitting or changing modes. Recompute full-nonlinear high-minus-low responses within ORIGINAL frozen quartile groups, renormalizing means to retained members. Report untouched full-cohort results alongside.',
        limits='Post hoc descriptive influence only; not actual re-selection, significance recalibration or corrected cosmology. Classifier/errors/detection/BBC unchanged.',
        inputs_sha256={str(p.relative_to(ROOT)):sha(p)for p in P+[Path(__file__).resolve()]}))
    print(json.dumps({'removed':ids,'protocol_sha256':sha(q)}))

def score():
    q=json.loads((OUT/'protocol.json').read_text())
    for p,h in q['inputs_sha256'].items():assert sha(ROOT/p)==h
    if (OUT/'result.json').exists():raise FileExistsError('Keep original result')
    rm=set(q['removed_ids'])
    rows=pd.read_csv(P[1],dtype={'CID':str})
    mem=pd.read_csv(P[2],dtype={'CID':str})
    scores=pd.read_csv(P[3],dtype={'CID':str})
    scores=scores[scores.arm=='published_mask'].copy()
    assert len(scores)==1020 and scores.CID.nunique()==1020
    assert set(mem.CID)==set(scores.CID)==set(rows.CID)
    assert rm<=set(mem.CID)
    low=set(mem.loc[mem.low_quartile,'CID']);high=set(mem.loc[mem.high_quartile,'CID'])
    out={}
    for label,exclude in [('original',set()),('common_cut',rm)]:
        s=scores[~scores.CID.isin(exclude)]
        M=float(s.matched_filter.sum());I=float(s.information.sum());G=float(s.fixed_prediction_gain.sum())
        assert abs(G-(M-.5*I))<1e-10
        details={}
        for target in ['native_noiseless','observed_flux']:
            d=rows[(rows.target==target)&~rows.CID.isin(exclude)]
            arm={}
            for mode in ['observer','sed']:
                m=d[d['mode']==mode].set_index('CID').nonlinear_standardized_mag
                lo=sorted(low-exclude);hi=sorted(high-exclude)
                arm[mode]=float(m.loc[hi].mean()-m.loc[lo].mean())
            details[target]={**arm,'sed_minus_observer_mag':arm['sed']-arm['observer']}
        out[label]=dict(n=len(s),n_low=len(low-exclude),n_high=len(high-exclude),M=M,I=I,G=G,
            descriptive_amplitude=M/I,contrasts=details)
    scores[scores.CID.isin(rm)].to_csv(OUT/'removed-object-scores.csv',index=False)
    write(OUT/'result.json',dict(protocol_sha256=sha(OUT/'protocol.json'),results=out))
    (OUT/'executed_source.py').write_bytes(Path(__file__).read_bytes())
    print(json.dumps(out,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['freeze','score']);a=p.parse_args();globals()[a.mode]()
