"""Independent finite-distribution and nuisance-fit falsifiers."""
import sys,json
from pathlib import Path
import numpy as np,pandas as pd
sys.path.insert(0,str(Path(__file__).resolve().parent))
from compare import Forecast,raw_groups,OUT
from selection import contrast

def fixture():
    rng=np.random.default_rng(314159);n=240;ns=350
    r=pd.DataFrame({'z':np.linspace(.03,1.1,n),'field':['C1']*n,'host':np.zeros(n,dtype=int),'fold':np.arange(n)%5,'u':rng.normal(-19.2,.2,n),'c':rng.normal(0,.05,n),'x':rng.normal(0,1,n)})
    s=pd.DataFrame({'z':np.linspace(.025,1.2,ns),'field':['C1']*ns,'host':np.zeros(ns,dtype=int),'u':rng.normal(-19.3,.2,ns),'c':rng.normal(0,.05,ns),'x':rng.normal(0,1,ns),'generated_attempt_index':np.arange(1,ns+1)})
    for f in [r,s]:
        for col in ['sigma_m','sigma_x','sigma_c','fitprob']:f[col]=.1
    return r,s

def main():
    r,s=fixture();f=Forecast(r,s,.1);a,al=f.run(details=True)
    # Scalar CRPS by direct double summation, separate from vectorized scorer.
    target=f.test[11];group=next(g for g in f.groups if target in g[0]);ri,si,w=group;j=np.flatnonzero(ri==target)[0];p=w[j]/w[j].sum();obs=r.iloc[target].c;v=s.c.to_numpy()[si]
    exact=sum(pi*abs(vi-obs) for pi,vi in zip(p,v))-.5*sum(pi*pj*abs(vi-vj) for pi,vi in zip(p,v) for pj,vj in zip(p,v))
    np.testing.assert_allclose(a['crps_c'][11],exact,rtol=1e-13,atol=1e-14)
    shifted=r.copy();shifted.u+=42;b,bl=Forecast(shifted,s,.1).run()
    for k in b:np.testing.assert_allclose(a[k],b[k],atol=4e-13)
    # A heldout-only magnitude change cannot be absorbed by train alignment.
    poisoned=r.copy();poisoned.loc[poisoned.fold==0,'u']+=1;c,cl=Forecast(poisoned,s,.1).run()
    np.testing.assert_allclose(al['offsets'],cl['offsets'],atol=1e-14);np.testing.assert_allclose(a['energy_cx'],c['energy_cx'],atol=1e-14);assert c['crps_m'].mean()>a['crps_m'].mean()+.5
    # Sparse exact-field support is rejected, not borrowed from another field.
    absent=r.copy();absent.loc[0,'field']='absent';_,sup=raw_groups(absent,s,.1);assert np.array_equal(sup[0],[0,0])
    # Same marginal counts can hide individual flips; paired variance detects them.
    sel=contrast([1,0,1,0],[0,1,1,0],[1,2,3,4]);assert sel['delta_vs_P21']==0 and sel['newly_selected_vs_P21']==sel['lost_vs_P21']==1 and sel['paired_standard_error']>0
    assert al['rank']==7 and all((a[k]>=-1e-12).all() for k in ['crps_m','crps_c','crps_x','energy_cx','energy_mcx'])
    out={'all_pass':True,'tests':['direct double-sum CRPS agreement','global magnitude shift invariance','heldout magnitude cannot contaminate train offset','exact-field unsupported rejection','paired selection flips and variance','full nuisance rank and nonnegative proper scores'],'primary_crps_check_absolute_error':float(abs(a['crps_c'][11]-exact))}
    (OUT/'numerical-validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out))

if __name__=='__main__':main()
