"""Independent optimizer reruns only for non-success synthetic terminations."""
import numpy as np
import pandas as pd
import json
from scipy.optimize import minimize
from audit import load_c1,Fit,CONST,exact,OUT,write

out=[]
for wide in [False,True]:
    rng=np.random.default_rng(260609650)
    if wide:table=pd.DataFrame(json.loads((OUT/'lowz-wide-realizations.json').read_text()))
    else:table=pd.read_csv(OUT/'synthetic-realizations.csv')
    for upper in ([.1,.2] if wide else [.1,.2,.4,.8]):
        tab,z,zh,cd,c=load_c1(zmin=.01,zmax=upper,frame='zHD');f=Fit(z,zh,cd,c)
        noise=(f.u*np.sqrt(f.e))@rng.normal(size=(len(z),100))
        failed=table[(table.zmax==upper)&(~table.success)]
        for r in failed.itertuples():
            y=CONST+5*np.log10(exact(z,r.model)*(1+zh)/(1+z))-19.3+noise[:,r.realization]
            old=f.objective(np.array([r.q,r.j]),y,False,False)[0]
            a=minimize(lambda p:f.objective(p,y,False,False)[0],[r.q,r.j],method='Nelder-Mead',bounds=[(-15,15),(-2000,2000)] if wide else [(-3,3),(-20,20)],options={'xatol':1e-8,'fatol':1e-9,'maxiter':3000})
            out.append(dict(wide=wide,model=r.model,zmax=upper,realization=int(r.realization),original_q=float(r.q),original_j=float(r.j),rerun_q=float(a.x[0]),rerun_j=float(a.x[1]),nll_improvement=float(old-a.fun),success=bool(a.success)))
write('synthetic-optimizer-flag-check.json',out)
print(out)
