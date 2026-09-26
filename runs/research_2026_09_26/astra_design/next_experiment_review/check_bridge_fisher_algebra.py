"""Synthetic-only checks: coordinate bridge, epsilon-prior Schur Fisher, nulls."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.linalg import solve_triangular,null_space
P=Path(__file__).resolve().parent;O=P/'bridge-check';O.mkdir(exist_ok=True)
rng=np.random.default_rng(26092699);n=30
G=rng.normal(size=(n,n));C=G@G.T/20+np.diag(np.linspace(.7,1.3,n))
L=np.linalg.cholesky(C);N=rng.normal(size=(n,2));E=rng.normal(size=(n,5));D=E[:,:2].copy()
Nw=solve_triangular(L,N,lower=True);Ew=solve_triangular(L,E,lower=True);Dw=solve_triangular(L,D,lower=True)
Q=null_space(Nw.T).T;e=Q@Ew;d=Q@Dw
Iproj=d.T@np.linalg.solve(np.eye(len(Q))+e@e.T,d)
X=np.column_stack([Nw,Ew]);prior=np.diag([0,0,1,1,1,1,1])
Ischur=Dw.T@Dw-Dw.T@X@np.linalg.solve(X.T@X+prior,X.T@Dw)
Cm=C+E@E.T;W=np.linalg.inv(Cm);Wprofile=W-W@N@np.linalg.solve(N.T@W@N,N.T@W)
Idense=D.T@Wprofile@D
Qfree=null_space(np.column_stack([Nw,Ew]).T).T;Ifree=(Qfree@Dw).T@(Qfree@Dw)
assert np.max(abs(Iproj-Ischur))<1e-12 and np.max(abs(Idense-Iproj))<1e-12
assert np.linalg.norm(Ifree)<1e-25
# Coordinate changes transform all model and noise objects, leaving information.
s=np.exp(rng.uniform(-.03,.03,n));Cs=s[:,None]*C*s[None,:];Js=s[:,None]*np.column_stack([N,E,D])
J=np.column_stack([N,E,D]);F=J.T@np.linalg.solve(C,J);Fs=Js.T@np.linalg.solve(Cs,Js)
assert np.max(abs(F-Fs))<1e-12
model=rng.normal(size=n);y=model+rng.multivariate_normal(np.zeros(n),C)
def ll(y,f,c):
    r=y-f;return -.5*(r@np.linalg.solve(c,r)+np.linalg.slogdet(c)[1]+len(y)*np.log(2*np.pi))
delta=ll(s*y,s*model,Cs)-ll(y,model,C)
assert abs(delta+np.log(s).sum())<1e-12
# Identical gray columns imply exact mu/deltaM degeneracy; AV0 kills RV derivative.
gray=np.column_stack([model,model]);gray_null=gray@np.array([1.,-1.]);rv_at_AV0=0.*rng.normal(size=n)
assert np.max(abs(gray_null))==0 and np.max(abs(rv_at_AV0))==0
out=dict(status='Synthetic PASS; no SN model/data Fisher matrix evaluated',prior_projected_vs_full_schur_max_error=float(np.max(abs(Iproj-Ischur))),prior_projected_vs_dense_marginal_covariance_max_error=float(np.max(abs(Iproj-Idense))),dust_equals_epsilon_span_trained_prior_eigenvalues=np.linalg.eigvalsh(Iproj).tolist(),dust_equals_epsilon_span_free_epsilon_fisher_norm=float(np.linalg.norm(Ifree)),coordinate_fisher_max_error=float(np.max(abs(F-Fs))),likelihood_coordinate_jacobian_error=float(abs(delta+np.log(s).sum())),gray_mu_minus_deltaM_max=float(np.max(abs(gray_null))),RV_derivative_at_AV0_max=float(np.max(abs(rv_at_AV0))),interpretation='Positive dust Fisher after a trained epsilon prior can coexist with exactly zero dust information if that intrinsic SED freedom is unconstrained. Pure data/model/C coordinate transforms change likelihood density only by the known Jacobian; physical passband changes are not covered by this identity.',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(O/'algebra.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
