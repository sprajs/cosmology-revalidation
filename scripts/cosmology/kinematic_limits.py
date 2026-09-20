"""Profile constrained histories; no Wilks significance assigned to boundary null."""
import json
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.interpolate import PchipInterpolator
from core import ROOT,Pantheon,mu,manifest

out=ROOT/'runs/cosmology/kinematic-limits';out.mkdir(parents=True,exist_ok=True)
sn=Pantheon();p=ROOT/'runs/cosmology/templates/c14-cpl63-median.csv'
d=pd.read_csv(p);sn.reset_template(PchipInterpolator(d.z,d.delta_mu)(sn.z))
rng=np.random.default_rng(41012);results={}
for amp in [0.,1.]:
    case={}
    for label,bounds in [('unconstrained',[(-3.,2.)]*5),('q_nonnegative',[(0.,2.)]*5)]:
        opts=[]
        for k in range(8):
            x=np.array([rng.uniform(a,b) for a,b in bounds]) if k else np.full(5,.1)
            r=minimize(lambda q:sn.chisq(q,'qbins',amp)[0],x,bounds=bounds,method='L-BFGS-B',
                       options={'maxiter':5000,'ftol':1e-12,'gtol':1e-7})
            opts.append({'q':r.x.tolist(),'chisq':float(r.fun),'success':bool(r.success),'message':r.message})
        best=min(opts,key=lambda r:r['chisq'])
        cross=minimize(lambda q:sn.chisq(q,'qbins',amp)[0],best['q'],bounds=bounds,method='Powell',
                       options={'maxiter':5000,'xtol':1e-9,'ftol':1e-12})
        opts.append({'q':cross.x.tolist(),'chisq':float(cross.fun),'success':bool(cross.success),'message':cross.message,'method':'independent Powell'})
        # Prefer successful equivalent solutions over an L-BFGS line-search warning.
        successful=[r for r in opts if r['success'] and r['chisq']<best['chisq']+1e-5]
        if successful: best=min(successful,key=lambda r:r['chisq'])
        case[label]={'best':best,'all_starts':opts,'full_chisq':float(sn.chisq_full(best['q'],'qbins',amp)[0])}
    case['delta_chisq_nonnegative']=case['q_nonnegative']['best']['chisq']-case['unconstrained']['best']['chisq']
    case['coasting_chisq']=float(sn.chisq([0.]*5,'qbins',amp)[0])
    results[f'fixed_age_amplitude_{amp}']=case
# Local information about separation of luminosity evolution from expansion parameters.
theta=np.array([.31,-.84,-.6]); step=np.array([1e-4,1e-4,1e-4])
cols=[]
for k in range(3):
    a=theta.copy();b=theta.copy();a[k]+=step[k];b[k]-=step[k]
    cols.append(((mu(sn.z,a,'cpl',sn.zhel)-mu(sn.z,b,'cpl',sn.zhel))/(2*step[k]))[0])
J=np.column_stack(cols);f=sn.template
F=J.T@sn.A@J;v=J.T@sn.A@f
ftotal=float(f@sn.A@f);remaining=float(ftotal-v@np.linalg.solve(F,v))
results['local_evolution_cosmology_degeneracy']={'fiducial_CPL':theta.tolist(),'evolution_total_information':ftotal,
    'information_after_profiling_CPL':remaining,'fraction_absorbed_by_linearized_CPL':1-remaining/ftotal,
    'caveat':'local tangent geometry, not global posterior or proof of physical evolution'}
(out/'results.json').write_text(json.dumps(results,indent=2)+'\n')
manifest(out,'Finite-bin nonacceleration profile test and local nuisance degeneracy',sn.inputs+[p],
    {'seed':41012,'starts':8,'smoothing_prior':None,'q_bounds':[-3,2],'nonaccelerating_bounds':[0,2],
     'age_amplitudes':[0,1],'boundary_likelihood_ratio_calibration':'none; do not call sigma'},[out/'results.json'])
print(json.dumps({k:v.get('delta_chisq_nonnegative') for k,v in results.items()},indent=2))
