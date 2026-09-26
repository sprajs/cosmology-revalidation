"""Independent covariance algebra, no observed pixels or peer numerical imports."""
from pathlib import Path
import numpy as np,json,hashlib
P=Path(__file__).resolve().parent
t=np.arange(7,dtype=float)*50; tau=t-t.mean(); one=np.ones(7); sigma2=21**2/2
D=np.eye(7)-np.outer(one,np.eye(7)[0])
X=np.column_stack([one,t]); R=np.eye(7)-X@np.linalg.solve(X.T@X,X.T)
rows=[]
for p in [0,.4,1,3,6,10]:
 u=np.abs((np.arange(7)-3)/3)**p; tt=t-u@t/u.sum(); Q=u@(tt*tt);h=u*tt/Q
 reported=441/Q;actual=sigma2*(h@h)
 poisson=np.minimum(t[:,None],t[None,:]);increment=sum(50*h[k:].sum()**2 for k in range(1,7))
 K=sigma2*np.eye(7)+93*np.outer(one,one)+.007*np.outer(t,t)
 # A common zero-read subtraction leaves any fixed-intercept slope unchanged,
 # even for correlated reads. The covariance need not be white.
 common_error=abs(h@D@K@D.T@h-h@K@h)
 # A random linear drift is invisible after fitting that ramp's own line.
 drift=.002*np.outer(t,t)
 projected=abs(R@drift@R.T).max()
 drift_slope=float(h@drift@h)
 alternatives=[]
 for rho in [0,.25,.5,.8]:
  # Illustrative AR(1) family with same adjacent-read difference variance441.
  Kwhite=441/(2*(1-rho))*rho**np.abs(np.arange(7)[:,None]-np.arange(7)[None,:])
  alternatives.append({'rho':rho,'adjacent_difference_variance':float(Kwhite[0,0]+Kwhite[1,1]-2*Kwhite[0,1]),'slope_read_variance':float(h@Kwhite@h),'reported_over_this_model':float(reported/(h@Kwhite@h))})
 assert common_error<1e-12 and projected<1e-12 and abs(drift_slope-.002)<1e-15
 assert abs(h.sum())<1e-17 and abs(h@t-1)<1e-14
 assert abs(h@poisson@h-increment)<1e-16
 rows.append({'power':p,'sum_h':float(h.sum()),'h_dot_t':float(h@t),'reported_read_var':float(reported),'independent_single_read_var':float(actual),'reported_over_independent':float(reported/actual),'reported_Poisson_coefficient':1/300,'cumulative_Poisson_coefficient':float(h@poisson@h),'common_reference_covariance_identity_error':float(common_error),'intraramp_detrended_drift_cov_max':float(projected),'unseen_drift_slope_variance':drift_slope,'same_adjacent_CDS_covariance_examples':alternatives})
(P/'result.json').write_text(json.dumps({'scope':'Synthetic algebra only. AR1 models are counterexamples, not fitted detector noise.','n':7,'times_seconds':t.tolist(),'rows':rows},indent=2)+'\n')
print(json.dumps({'independent_read_ratio_range':[min(r['reported_over_independent'] for r in rows),max(r['reported_over_independent'] for r in rows)],'equal_weight_same_CDS_examples':rows[0]['same_adjacent_CDS_covariance_examples'],'detrending_counterexample':{'residual_cov_max':rows[0]['intraramp_detrended_drift_cov_max'],'slope_variance':rows[0]['unseen_drift_slope_variance']}},indent=2))
