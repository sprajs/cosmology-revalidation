"""Independent Gaussian marginal-ratio check of leave-object-out scores."""
from pathlib import Path
import numpy as np,pandas as pd,json,hashlib
from scipy.linalg import cho_factor,cho_solve
p=Path(__file__).resolve().parent
a=np.load(p/'residual-transfer-arrays.npz');r=a['residual']
info=pd.read_csv(p/'residual-transfer-objects.csv')
reported=json.loads((p/'residual-transfer.json').read_text())
def lp(y,C):
    f=cho_factor(C,lower=True)
    return -.5*(y@cho_solve(f,y)+2*np.log(np.diag(f[0])).sum()+len(y)*np.log(2*np.pi))
checks={}
for family in ['observer','rest_phase']:
    T=a[family]
    for sigma in [.01,.02,.05]:
        C=np.eye(len(r))+sigma**2*T@T.T
        full=lp(r,C);total=0.
        for row in info.itertuples():
            test=np.arange(row.row_start,row.row_end)
            train=np.setdiff1d(np.arange(len(r)),test)
            score=full-lp(r[train],C[np.ix_(train,train)])-lp(r[test],np.eye(len(test)))
            total+=score
        name=f'{family}_{sigma:g}'
        delta=total-reported['results'][name]['sum_delta_log_score']
        checks[name]={'direct_joint_minus_training_score':float(total),'difference_from_operator_score':float(delta)}
assert max(abs(v['difference_from_operator_score']) for v in checks.values())<1e-9
out={'checks':checks,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(p/'residual-score-check.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
