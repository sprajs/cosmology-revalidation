"""Independent central-value stress budget for the exploratory BGS inequality.

This is a necessary-contrast repair scale, not a fitted systematic prior.
"""
from pathlib import Path
import numpy as np, pandas as pd, json, hashlib
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[3]
p=ROOT/'sources/repos/CobayaSampler__bao_data/desi_bao_dr2'
mean=p/'desi_gaussian_bao_ALL_GCcomb_mean.txt'
covpath=p/'desi_gaussian_bao_ALL_GCcomb_cov.txt'
d=pd.read_csv(mean,sep=r'\s+',comment='#',names=['z','v','kind'])
C=np.loadtxt(covpath);y=d.v.to_numpy()
ib=int(np.flatnonzero(d.kind=='DV_over_rs')[0]);b=float(d.z.iloc[ib]);tb=np.log1p(b)
factor=(b/(1+b)*tb**2)**(1/3)
j=int(np.flatnonzero((d.z==.934)&(d.kind=='DH_over_rs'))[0])
v=np.zeros(len(y));v[ib]=1/factor;v[j]=-(1+d.z.iloc[j])
contrast=float(v@y);variance=float(v@C@v)
correction=-contrast*(C@v)/variance
sb=float(y[ib]/factor);gj=float((1+d.z.iloc[j])*y[j]);ratio=sb/gj
u=brentq(lambda x:np.sin(x)/x-ratio**1.5,1e-9,np.pi/2)
radius_max=tb*gj/u
dmmax=float(d.loc[d.kind=='DM_over_rs','v'].max())
maxbending=(np.sin(tb*gj/dmmax)/(tb*gj/dmmax))**(2/3)
out={
 'scope':'Central-value necessary-inequality repair, not a global fit or empirical systematic prior.',
 'BGS_z':b,'radial_z':float(d.z.iloc[j]),'s_BGS':sb,'g_radial':gj,
 'contrast':contrast,'contrast_sd':float(np.sqrt(variance)),'deficit_sigma':float(-contrast/np.sqrt(variance)),
 'BGS_DV_relative_increase_if_alone':1/ratio-1,
 'radial_DH_relative_decrease_if_alone':1-ratio,
 'ruler_ratio_highz_to_BGS_at_coasting_equality':ratio,
 'minimum_one_contrast_covariance_metric_shift':{
   'chi2':float(correction@np.linalg.solve(C,correction)),
   'BGS_DV_relative_change':float(correction[ib]/y[ib]),
   'radial_DH_relative_change':float(correction[j]/y[j]),
   'all_observable_shifts':correction.tolist(),
   'repaired_contrast':float(v@(y+correction))},
 'closed_curvature_on_increasing_distance_branch':{
   'required_sinc_argument':float(u),
   'required_curvature_radius_over_rd_at_most':float(radius_max),
   'required_abs_Omega_k_at_least_if_g0_ge_gj':float((u/tb)**2),
   'largest_observed_DM_over_rd':dmmax,
   'radius_over_rd_at_least_to_fit_that_central_DM':dmmax,
   'maximum_fractional_BGS_s_reduction_at_that_radius_and_constant_gj':float(1-maxbending),
   'qualification':'Flat coasting minimizes BGS chi at fixed later g; for positive monotone g and transverse distance on the increasing sin branch, the closed-curvature bound follows. Does not admit antipodal/conjugate-point histories or propagate measurement uncertainty. Such extreme histories also challenge BAO compression.'},
 'inputs_sha256':{str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in (mean,covpath)},
 'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
(Path(__file__).parent/'bao-scale-budget.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
