"""Prespecified boundary-response rerun; original low-z results preserved."""
import numpy as np
import pandas as pd
from audit import load_c1,Fit,CONST,exact,OUT,write

rng=np.random.default_rng(260609650)
rows=[];summ=[]
for zmax in [.1,.2]:
    tab,z,zh,cd,c=load_c1(zmin=.01,zmax=zmax,frame='zHD');f=Fit(z,zh,cd,c)
    noise=(f.u*np.sqrt(f.e))@rng.normal(size=(len(z),100))
    for kind,(q,j) in {'lcdm':(-.55,1),'eds':(.5,1),'coasting':(0,0)}.items():
        y=CONST+5*np.log10(exact(z,kind)*(1+zh)/(1+z))-19.3
        fit=f.fit(y,start=[q,j],wide=True);temp=[]
        for i in range(100):
            a=f.fit(y+noise[:,i],start=[fit['q'],fit['j']],wide=True)
            temp.append(a['q']);rows.append(dict(model=kind,zmax=zmax,realization=i,**a))
        summ.append(dict(model=kind,zmax=zmax,noiseless=fit,mean_q=float(np.mean(temp)),sd_q=float(np.std(temp,ddof=1)),quantiles_q=np.quantile(temp,[.025,.5,.975]).tolist(),jerk_bound_hits=sum(abs(r['j'])>1999.99 for r in rows[-100:]),q_bound_hits=sum(abs(r['q'])>14.999 for r in rows[-100:]),successes=sum(r['success'] for r in rows[-100:])))
        write('lowz-wide-summary.json',summ);write('lowz-wide-realizations.json',rows)
        print(summ[-1],flush=True)
