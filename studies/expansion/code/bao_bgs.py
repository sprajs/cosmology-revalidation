"""Post-pilot exploratory BAO inequalities with BGS; preserve chronology."""
from pathlib import Path
import hashlib
import json
from datetime import datetime, timezone

import numpy as np
import pandas as pd
from scipy.stats import norm
from bao_shape import ROOT, OUT, load_data, cone_matrices, tail

SOURCE=Path(__file__).read_bytes()
PROTOCOL=ROOT/'docs/research-2026-09-26/bao-bgs-amendment.md'
PROTOCOL_BYTES=PROTOCOL.read_bytes()


def main():
    out=OUT/'bgs_exploratory'
    if out.exists():raise FileExistsError(out)
    out.mkdir()
    (out/'executed_source.py').write_bytes(SOURCE)
    (out/'protocol.md').write_bytes(PROTOCOL_BYTES)
    z,y,C,hashes,order,raw,rcov=load_data()
    n=len(z);t=np.log1p(z)
    idx=int(np.flatnonzero(raw.kind=='DV_over_rs')[0]);b=float(raw.z.iloc[idx]);tb=np.log1p(b)
    factor=(b/(1+b)*tb**2)**(1/3)
    scale=np.r_[np.ones(n),1+z,1/factor]
    indices=order+[idx]
    yy=raw.value.to_numpy()[indices]*scale
    CC=rcov[np.ix_(indices,indices)]*np.outer(scale,scale)
    rows=[];names=[];families=[]
    for j in range(n):
        row=np.zeros(2*n+1);row[-1]=1;row[n+j]=-1
        rows.append(row);names.append(f'BGS_to_radial_z{z[j]}');families.append('BGS_radial')
    for j in range(n):
        row=np.zeros(2*n+1);row[j]=1;row[n+j]=-t[j]
        rows.append(row);names.append(f'AP_z{z[j]}');families.append('AP')
    A,_,labels=cone_matrices(t)
    for row,name in zip(A,labels):
        if name.startswith('positive'):continue
        rows.append(np.r_[row,0]);names.append(name);families.append('interval_or_radial')
    B=np.array(rows)
    assert len(rows)==28
    coast=np.r_[t,np.ones(n),1.]
    assert np.max(np.abs(B@coast))<1e-14
    # Independent direct checks on generated finite-q histories, including q<0.
    generated={}
    for q in [-.5,0,.5,1.]:
        g=30*np.exp(-q*t);gb=30*np.exp(-q*tb)
        d=30*t if q==0 else -30*np.expm1(-q*t)/q
        db=30*tb if q==0 else -30*np.expm1(-q*tb)/q
        sb=((db/tb)**2*gb)**(1/3)
        vals=B@np.r_[d,g,sb]
        generated[str(q)]=float(vals.min())
        if q>=0:assert vals.min()>-1e-10
    vc=B@CC@B.T;sd=np.sqrt(np.diag(vc));contrasts=B@yy;standardized=contrasts/sd
    obs=float(np.maximum(0,-standardized.min()))
    bgsobs=float(np.maximum(0,-standardized[:n].min()))
    L=np.linalg.cholesky(CC);response=B@L/sd[:,None]
    rng=np.random.default_rng(260928);N=1000000
    counts=0;bgscounts=0
    for start in range(0,N,20000):
        values=rng.normal(size=(min(20000,N-start),len(yy)))@response.T
        counts+=int(np.sum(-values.min(axis=1)>=obs))
        bgscounts+=int(np.sum(-values[:,:n].min(axis=1)>=bgsobs))
    # Reuse binomial summary without retaining the full 28-million contrast array.
    calibrated=lambda count:tail(np.r_[np.ones(count),np.zeros(N-count)],.5)
    table=pd.DataFrame({'name':names,'family':families,'contrast':contrasts,'sd':sd,
                        'standardized':standardized,'one_sided_marginal_p':norm.cdf(standardized)})
    table.to_csv(out/'all_contrasts.csv',index=False)
    results={'status':'exploratory after cone and BGS plug-in inspection; no replacement of primary',
             'BGS_z':b,'BGS_conversion_factor':factor,'BGS_s':float(yy[-1]),
             'BGS_s_sd':float(np.sqrt(CC[-1,-1])),
             'joint_max_deficit_sigma':obs,'joint_tail':calibrated(counts),
             'BGS_only_max_deficit_sigma':bgsobs,'BGS_family_tail':calibrated(bgscounts),
             'worst_contrast':names[int(np.argmin(standardized))],
             'constant_q_validation':generated,'seed':260928,'null_draws':N,
             'interpretation':'Conditional evidence about some past acceleration under flat constant-ruler compressed-BAO assumptions; not present q0 or new physics.',
             'multiplicity':'Joint Gaussian maximum across the 28 listed contrasts; does not correct all possible adaptive research paths.'}
    (out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
    np.savez_compressed(out/'matrices.npz',B=B,data=yy,covariance=CC,contrast_covariance=vc)
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    outputs=[out/'results.json',out/'matrices.npz',out/'all_contrasts.csv']
    (out/'manifest.json').write_text(json.dumps({'utc':datetime.now(timezone.utc).isoformat(),
        'input_hashes':hashes,'source_sha256':hashlib.sha256(SOURCE).hexdigest(),
        'protocol_sha256':hashlib.sha256(PROTOCOL_BYTES).hexdigest(),
        'bao_shape_source_sha256':sha(Path(__file__).with_name('bao_shape.py')),
        'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in outputs}},indent=2)+'\n')
    print(json.dumps(results,indent=2))


if __name__=='__main__':main()
