"""Reproduction of Ray v2's simplified estimator, not Sah's C2 model."""
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.interpolate import CubicSpline
from astropy.coordinates import SkyCoord
import astropy.units as u
from audit import SRC,OUT,cosine,cubic_mu,write

t=pd.read_csv(SRC/'Analysis C1/Pantheon+SH0ES.dat',sep=r'\s+')
t=t[(t.zHEL>.00937)&(t.zHEL<=.8)]
a=pd.read_csv(SRC/'median_deltaage.csv',header=None).sort_values(0);sp=CubicSpline(a[0],a[1])
rows=[]
for frame in ['zHEL','zCMB','zHD']:
    for age in [False,True]:
        for subset in ['all','dipole','antidipole']:
            co=cosine(t.RA.to_numpy(),t.DEC.to_numpy(),167.8,-7.1)
            tt=t if subset=='all' else t[co>=0] if subset=='dipole' else t[co<0]
            y=tt.mB.to_numpy()-(.03*sp(tt.zHEL.to_numpy()) if age else 0)
            A=np.column_stack((np.ones(len(tt)),-tt.x1.to_numpy(),tt.c.to_numpy()))
            def fun(q,detail=False):
                mu,_=cubic_mu(tt[frame].to_numpy(),tt.zHEL.to_numpy(),q,1)
                if mu is None:return 1e40
                mean=np.linalg.lstsq(A,y-mu,rcond=None)[0]
                chi=float(np.sum(((y-mu-A@mean)/.15)**2))
                return (chi,mean) if detail else chi
            r=minimize_scalar(fun,bounds=(-2,1.5),method='bounded',options={'xatol':1e-12})
            chi,mean=fun(r.x,True)
            rows.append(dict(frame=frame,age_correction=age,subset=subset,n=len(tt),q=float(r.x),j=1,M=float(mean[0]),alpha=float(mean[1]),beta=float(mean[2]),chi2=chi))
write('ray-v2-simple-fit.json',rows)

c=SkyCoord(l=264*u.deg,b=48*u.deg,frame='galactic').icrs
out=dict(n=len(t),actual_transformed_icrs=[c.ra.deg,c.dec.deg],counts={})
for tag,ra,dec in [('wrong_equatorial',264,48),('response_rounded',167.8,-7.1),('astropy_transformed',c.ra.deg,c.dec.deg),('sah_code_direction',168,-7),('planck_like_rounded',167.94,-6.94)]:
    ca=cosine(t.RA.to_numpy(),t.DEC.to_numpy(),ra,dec);out['counts'][tag]=dict(positive=int(sum(ca>=0)),negative=int(sum(ca<0)))
write('response-coordinate-check.json',out)
print(rows)
