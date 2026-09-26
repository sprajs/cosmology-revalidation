from pathlib import Path
import json,hashlib
import numpy as np
O=Path(__file__).resolve().parent
op=np.load(O/'fixed-operators.npz');t=op['times_seconds'];h=op['h'];checks={}
# Unsigned subtraction must not wrap; test FITS-BZERO-equivalent DN endpoints.
a=np.array([65535,0,32768],dtype=np.uint16);b=np.array([0,65535,32767],dtype=np.uint16)
d=b.astype(np.float64)-a.astype(np.float64)
assert np.array_equal(d,[-65535,65535,-1]);checks['unsigned_promotion']=True
# A deterministic algebra check of full read covariance and direct contractions.
rng=np.random.default_rng(260926);dy=rng.normal(size=(23,7))*np.arange(1,8)[None,:]+rng.normal(size=(23,1))*100
G=dy.T@dy/(2*len(dy));direct=np.mean((dy@h.T)**2,axis=0)/2;matrix=np.einsum('ki,ij,kj->k',h,G,h)
err=float(np.max(abs(direct-matrix)));assert err<1e-13;checks['direct_vs_matrix_max_error']=err
assert np.max(abs(h.sum(axis=1)))<1e-15 and np.max(abs(h@t-1))<1e-13
checks['slope_intercept_and_scale']=True
# A common-zero vector cancels from fixed slopes, without subtracting fitted slope.
z=rng.normal(size=(23,1))*1000;assert np.max(abs((dy-z)@h.T-dy@h.T))<1e-13
checks['zero_read_common_term_cancels']=True
# Gain is a deterministic coordinate transformation, not an inferred variance.
for g in [2.34,2.37,2.31,2.38]:assert np.allclose(np.mean(((g*dy)@h.T)**2,axis=0)/2,g*g*direct,rtol=1e-13,atol=1e-14)
checks['fixed_gain_squared_scaling']=True
# No event rejection: a step event remains in the total second moment.
step=np.array([0,0,0,100,100,100,100],float);s=(h@step)**2/2
assert np.all(s>0);checks['cosmic_step_retained_slope_second_moments']=s.tolist()
# Detrending loses the desired slope direction; it cannot be the primary statistic.
A=np.column_stack([np.ones(7),t]);P=np.eye(7)-A@np.linalg.inv(A.T@A)@A.T
assert np.linalg.norm(P@t)<1e-10 and abs(h[0]@t-1)<1e-13
checks['detrending_erases_slope_mode']=True
# Masked geometric aperture: preserve zero-sum weights but do not normalize missing flux.
a=np.array([.2,1,.6,0,0,0]);b=np.array([0,0,0,.3,1,.7]);m=np.array([1,1,0,1,1,1])
w=m*a-(m*a).sum()/(m*b).sum()*m*b
assert abs(w.sum())<1e-15;assert abs(w@np.ones(6))<1e-15
checks['fixed_mask_background_cancellation']=True
# Nominal clock scaling: a uniform rate-scale error changes h and second moment exactly.
for e in [-.01,.01]:
 tt=(1+e)*t;u=np.ones(7);hh=(tt-tt.mean())/np.sum((tt-tt.mean())**2)
 assert np.allclose(hh,h[0]/(1+e),atol=1e-17)
checks['clock_scale_relation']=True
res={'status':'PASS_SYNTHETIC_ONLY','checks':checks,'no_observed_arrays':True,'seed':260926,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'operator_sha256':hashlib.sha256((O/'fixed-operators.npz').read_bytes()).hexdigest()}
(O/'synthetic-result.json').write_text(json.dumps(res,indent=2)+'\n')
print(json.dumps({'status':res['status'],'direct_error':err}))
