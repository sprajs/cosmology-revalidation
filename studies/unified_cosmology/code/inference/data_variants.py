"""Acquire the pinned primary Pantheon+ release as an alternative SN likelihood.

Never multiply this sample with Dovekie: events and calibration overlap.
The magnitude-column name MU below is a common adapter interface only;
m_b_corr is used with a free zero point, without SH0ES or Cepheid anchoring.
"""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import quote
import numpy as np
import pandas as pd
from scipy.linalg import cho_factor,cho_solve

ROOT = Path(__file__).resolve().parents[4]
WORK = ROOT/'.work/unified-cosmology/inference/pantheon'
COMMIT = 'c447f0fea703fcd0fff57de5000947b5ca81286b'


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    WORK.mkdir(parents=True,exist_ok=True)
    records = []
    for name in ['Pantheon+SH0ES.dat','Pantheon+SH0ES_STAT+SYS.cov','README']:
        url = f'https://raw.githubusercontent.com/PantheonPlusSH0ES/DataRelease/{COMMIT}/'+quote('Pantheon+_Data/4_DISTANCES_AND_COVAR/'+name)
        path = WORK/name
        if not path.exists():
            with urlopen(url,timeout=120) as response:
                body = response.read()
            temporary = path.with_suffix(path.suffix+'.part')
            temporary.write_bytes(body);temporary.replace(path)
        old = ROOT/'data/distances'/name
        records.append({'url':url,'path':str(path.relative_to(ROOT)),'bytes':path.stat().st_size,
                        'sha256':digest(path),'identical_to_preexisting_input':digest(old)==digest(path) if old.exists() else None})
    table = pd.read_csv(WORK/'Pantheon+SH0ES.dat',sep=r'\s+')
    raw = np.loadtxt(WORK/'Pantheon+SH0ES_STAT+SYS.cov')
    n = int(raw[0]);assert raw.size==n*n+1 and n==len(table)
    covariance = raw[1:].reshape(n,n)
    asymmetry = float(np.max(abs(covariance-covariance.T)))
    assert asymmetry<1e-7
    covariance = (covariance+covariance.T)/2
    mask = table.zHD.gt(.01).to_numpy()
    retained = table.loc[mask]
    covariance = covariance[np.ix_(mask,mask)]
    p = cho_solve(cho_factor(covariance,lower=True),np.eye(mask.sum()))
    identity = float(np.max(abs(p@covariance-np.eye(len(p)))))
    assert identity<1e-10
    out = WORK/'pantheon-total.npz'
    np.savez_compressed(out,CID=retained.CID.to_numpy(dtype=str),IDSURVEY=retained.IDSURVEY.to_numpy(),
                        zHD=retained.zHD.to_numpy(),zHEL=retained.zHEL.to_numpy(),
                        MU=retained.m_b_corr.to_numpy(),covariance=covariance,precision=p)
    result = {'status':'verified','primary_release_commit':COMMIT,'files':records,
              'raw_rows':n,'retained_rows':len(retained),'distinct_supernova_names':int(retained.CID.nunique()),
              'selection':'zHD>0.01, noCepheidabsoluteanchor; duplicate measurements retained with released covariance',
              'original_covariance_max_asymmetry':asymmetry,'inverse_identity_max_error':identity,
              'output':{'path':str(out.relative_to(ROOT)),'sha256':digest(out)},
              'code_sha256':digest(Path(__file__)),
              'scope':'Alternative likelihood only; m_b_corr is integrated over a global absolute-magnitude offset; never union with Dovekie without event/calibration crosscovariance.'}
    dest = ROOT/'studies/unified_cosmology/results/inference/pantheon-interface.json'
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
