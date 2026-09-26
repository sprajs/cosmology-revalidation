#!/usr/bin/env python3
"""Independent block-normal check of the conditional integrated likelihood.
Verifies the algebra used in fullcov_eiv.py against one 2N block Gaussian,
and the zero-age-error GLS limit; neither uses LINMIX.
"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.linalg import cholesky,solve_triangular
sys.path.insert(0,str(Path(__file__).resolve().parent))
from analyse import ROOT,OUT,save_manifest,wls

def gaussian_nll(res,cov):
    L=cholesky(cov,lower=True);u=solve_triangular(L,res,lower=True);return np.log(np.diag(L)).sum()+.5*u@u

def main():
    rng=np.random.default_rng(20260921);n=13; q=rng.normal(size=(n,n));C=q@q.T*.0002+np.eye(n)*.01;sa=rng.uniform(.3,2,n);m=rng.uniform(3,7,n);tau=1.3;alpha=.1;b=-.031;s=.025;age=m+rng.normal(0,2,n);y=rng.normal(0,.1,n)
    av=tau**2+sa**2;am=m+tau**2/av*(age-m);apv=tau**2*sa**2/av
    conditional=.5*np.sum(np.log(av)+(age-m)**2/av)+gaussian_nll(y-alpha-b*am,C+np.diag(s*s+b*b*apv))
    block=np.block([[np.diag(av),np.eye(n)*b*tau**2],[np.eye(n)*b*tau**2,C+np.eye(n)*(b*b*tau**2+s*s)]])
    joint=gaussian_nll(np.r_[age-m,y-alpha-b*m],block)
    assert abs(conditional-joint)<1e-10
    # For exactly measured ages, dependence on m/tau separates and GLS is exact.
    vv=C+np.eye(n)*s*s; coef,cov,chi=wls(age,y,np.ones(n),cov=vv);X=np.c_[np.ones(n),age]
    from scipy.optimize import minimize
    rr=minimize(lambda par:gaussian_nll(y-X@par,vv),np.array([0.,0.]),jac=lambda par: -X.T@np.linalg.solve(vv,y-X@par),method='BFGS',tol=1e-10)
    discrepancy=float(np.max(abs(rr.x-coef))); assert discrepancy<1e-7
    result={'seed':20260921,'block_vs_conditional_nll_abs_difference':abs(float(conditional-joint)),'zero_age_error_gls_max_parameter_difference':discrepancy,'status':'pass','limits':'Validates Gaussian formula and zero-age-error limit; does not establish Gaussian age summaries or empirical selection model are correct.'}
    dest=OUT/'algebra-validation.json';dest.write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));save_manifest('A2-algebra-validation',[Path(__file__),ROOT/'scripts/age_signal/fullcov_eiv.py'],{'seed':20260921},[dest])
if __name__=='__main__':main()
