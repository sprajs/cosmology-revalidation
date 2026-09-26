"""Outcome-free algebra checks. Synthetic vectors only; no observed files loaded."""
from pathlib import Path
import numpy as np,json,hashlib,itertools
P=Path(__file__).resolve().parent;K=2.5/np.log(10);rng=np.random.default_rng(260927)

def ratio_box_bounds(v0,v1,q0,q1,low,high):
 # a1/a0=(q0/q1)*(v1.y)/(v0.y); positivity over entire box required.
 assert np.sum(np.where(v0>=0,v0*low,v0*high))>0
 assert np.sum(np.where(v1>=0,v1*low,v1*high))>0
 assert q0>0 and q1>0
 vals=[]
 for bits in itertools.product([False,True],repeat=len(low)):
  y=np.where(bits,high,low);vals.append((q0/q1)*(v1@y)/(v0@y))
 return min(vals),max(vals)

h=np.array([.8,1.,.7,1.2]);h1=np.array([.025,-.012,.018,-.032]);h2=np.array([-.002,.003,-.001,.004]);C=np.diag([.05,.08,.07,.04])**2+.0001*np.ones((4,4));W=np.linalg.inv(C);A=1.3;y=A*h+np.array([.01,-.005,.002,-.007]);S=h@W@h;T=h1@W@h;U=h1@W@h1;V=h2@W@h;B=h@W@y
f=lambda t:h+t*h1+t*t*h2/2
D=lambda t:-K*np.log((f(t)@W@y)/(f(t)@W@f(t)))
d1=K*(2*T/S-h1@W@y/B);d2=K*(2*(U+V)/S-4*T*T/S**2-h2@W@y/B+(h1@W@y/B)**2)
step=.001;num1=(D(step)-D(-step))/(2*step);num2=(D(step)-2*D(0)+D(-step))/step**2;assert abs(num1-d1)<1e-8 and abs(num2-d2)<1e-8
noiseless=[]
a=T/S;b=U/S;c=V/S
for eta in [.002,.001,.0005]:
 f1=f(eta);actual=-K*np.log((f1@W@(A*h))/(f1@W@f1)/A);approx=K*(a*eta+(b+c/2-1.5*a*a)*eta**2);noiseless.append({'eta':eta,'exact':actual,'second_order':approx,'remainder':actual-approx})
# Fixed matrices, positive denominator: extrema of linear-fractional ratio lie at boxvertices.
v0=W@h;hh=f(.2);v1=W@hh;q0=h@W@h;q1=hh@W@hh;low=y-.0005;high=y+.0005;lo,hi=ratio_box_bounds(v0,v1,q0,q1,low,high)
pts=rng.uniform(low,high,(10000,4));ratios=(q0/q1)*(pts@v1)/(pts@v0);assert ratios.min()>=lo and ratios.max()<=hi
# Explicit covariance term for a linear noisy timing estimator eta=l.n.
ell=np.array([.25,-.1,.15,-.2]);var_eta=float(ell@C@ell);cov_eta_e0=float(ell@C@W@h/(A*S));cov_eta_e1=float(ell@C@W@h1/(A*S));ind=K*(b+c/2-1.5*a*a)*var_eta;correction=K*(-cov_eta_e1+a*cov_eta_e0)
res={'scope':'Synthetic algebra/unit checks only; no real photons, archived timing outcomes or native fits.','exact_derivative_errors':{'first':float(num1-d1),'second':float(num2-d2)},'noiseless_taylor':noiseless,'box_ratio_bounds':[float(lo),float(hi)],'random_interior_containment':True,'illustrative_same_noise_terms':{'independent_Gaussian_second_order_mean':float(ind),'same_noise_covariance_extra':float(correction),'both':float(ind+correction)},'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()};(P/'algebra-check.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res,indent=2))
