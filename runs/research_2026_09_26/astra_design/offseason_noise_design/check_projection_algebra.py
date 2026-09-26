"""Synthetic-only validation of heterogeneous-error baseline projection."""
from pathlib import Path
import json
import numpy as np
P=Path(__file__).resolve().parent
s=np.array([.5,1,2,4.]);a=1/s
U,_,_=np.linalg.svd(a[:,None],full_matrices=True);Q=U[:,1:].T
rng=np.random.default_rng(26092698);n=100000
F=17+rng.normal(size=(n,4))*s;y=(F/s)@Q.T
H=np.zeros((4,4));H[0,1]=H[1,0]=.5;H[1,2]=H[2,1]=.5
A=Q@H@Q.T;q=np.einsum('ij,jk,ik->i',y,A,y)-np.trace(A)
assert np.max(abs(Q@a))<1e-12
assert np.max(abs(y-(F/s-17/s)@Q.T))<1e-12
r=dict(synthetic_only=True,baseline_projection_max=float(np.max(abs(Q@a))),projection_noise_cov_error=float(np.max(abs(np.cov(y,rowvar=False)-np.eye(3)))),lag_quadratic_mean=float(q.mean()),lag_quadratic_mcse=float(q.std()/np.sqrt(n)),mean_chi2_per_dof=float(np.sum(y*y)/n/3))
assert abs(r['lag_quadratic_mean'])<5*r['lag_quadratic_mcse']
(P/'synthetic-algebra-check.json').write_text(json.dumps(r,indent=2)+'\n')
print(json.dumps(r,indent=2))
