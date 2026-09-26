"""Independent explicit predictive-covariance check of 36 implementation scores."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy.linalg import solve_triangular
P=Path(__file__).resolve().parent/'exact43'/'comparison';a=np.load(P/'matched-matrices.npz');r=json.loads((P/'result.json').read_text());m=pd.read_csv(P/'object-checks.csv',dtype={'CID':str});errors={}
for arm in r['arms']:
    rs=[a[c+'__'+arm+'_r'] for c in m.CID];ts=[a[c+'__'+arm+'_T'] for c in m.CID]
    for family,ii in [('observer',slice(0,3)),('rest_phase',slice(3,6))]:
        F=sum(t[:,ii].T@t[:,ii] for t in ts);u=sum(t[:,ii].T@y for t,y in zip(ts,rs))
        for sigma in [.01,.02,.05]:
            for cv,groups in [('object',[[i] for i in range(len(m))]),('field',[list(np.flatnonzero(m.field==f)) for f in sorted(m.field.unique())])]:
                total=0.
                for ids in groups:
                    T=np.vstack([ts[i][:,ii] for i in ids]);y=np.concatenate([rs[i] for i in ids]);V=np.linalg.inv(np.eye(3)/sigma**2+F-T.T@T);mean=T@V@(u-T.T@y)
                    predictive=np.eye(len(y))+T@V@T.T;L=np.linalg.cholesky(predictive);z=solve_triangular(L,y-mean,lower=True)
                    total+=-.5*z@z-np.log(np.diag(L)).sum()+.5*y@y
                name=f'{cv}_{family}_{sigma:g}';errors[arm+'_'+name]=float(total-r['arms'][arm]['candidate_scores'][name])
assert max(abs(e) for e in errors.values())<1e-9
out={'method':'Reconstruct full heldout predictive covariance I+T Vtrain T^T and use Cholesky Gaussian log density; independent of saved sufficient-stat score implementation.','max_absolute_error':max(abs(e) for e in errors.values()),'errors':errors,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(P/'independent-score-check.json').write_text(json.dumps(out,indent=2)+'\n');print('maximum direct predictive score error',out['max_absolute_error'])
