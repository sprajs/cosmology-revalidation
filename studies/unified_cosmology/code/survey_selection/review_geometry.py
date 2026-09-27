#!/usr/bin/env python3
"""Independent bounded review of root-owned shared geometry and SN nuisance algebra."""
import json,sys
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from common import ROOT,WORK,RESULTS,sha
sys.path.insert(0,str(ROOT/'studies/unified_cosmology/code/inference'))
from likelihood import ReleasedDistances,expansion_diagnostics
from late_geometry import Geometry

class Provider:
    def __init__(self,om,w,wa,H0=70):self.om,self.w,self.wa,self.H0=om,w,wa,H0
    def E(self,z):return np.sqrt(self.om*(1+z)**3+(1-self.om)*(1+z)**(3*(1+self.w+self.wa))*np.exp(-3*self.wa*z/(1+z)))
    def get_angular_diameter_distance(self,z):return np.array([299792.458/self.H0*quad(lambda x:1/self.E(x),0,float(v),epsabs=1e-12,epsrel=1e-12)[0]/(1+v) for v in z])

def main():
    paths=[ROOT/'studies/unified_cosmology/code/inference'/n for n in ['likelihood.py','run.py','late_geometry.py']]
    reviewed_hashes={str(p.relative_to(ROOT)):sha(p) for p in paths}
    snapshot=WORK/'geometry-reviewed-source';snapshot.mkdir(parents=True,exist_ok=True)
    snapshots={}
    for p in paths:
        dest=snapshot/sha(p)/p.name;dest.parent.mkdir(parents=True,exist_ok=True)
        if dest.exists():assert sha(dest)==sha(p)
        else:dest.write_bytes(p.read_bytes())
        snapshots[str(dest.relative_to(ROOT))]=sha(dest)
    rows=[]
    for evo,sigma in [('none',0.),('linear',0.),('smooth01',.1),('smooth03',.3)]:
        g=Geometry(model='cpl',evolution=evo)
        l=ReleasedDistances.__new__(ReleasedDistances);l.data_file=str(WORK/'normalized/dovekie-total.npz');l.smooth_sigma=sigma;l.initialize()
        for om,w,wa,eps in [(.31,-1,0,.07),(.27,-.7,-1.1,-.09),(.55,-.3,.8,.15)]:
            if evo!='linear':eps=0.
            l.provider=Provider(om,w,wa);derived={};ll=l.logp(epsilon=eps,_derived=derived)
            theta=[om,9900,w,wa]+([eps] if evo=='linear' else [])
            cc=g.components(theta)[0][0];direct=g.direct_chi2(theta)
            l.provider=Provider(om,w,wa,93.);dd={};ll2=l.logp(epsilon=eps,_derived=dd)
            rows.append({'evolution':evo,'Omega_m':om,'w0':w,'wa':wa,'epsilon':eps,'shared_vs_independent_direct_chi2':abs(derived['sn_chi2']-direct),'compressed_vs_independent_chi2':abs(cc-direct),'H0_global_offset_logp_invariance':abs(ll-ll2)})
    assert max(r['shared_vs_independent_direct_chi2'] for r in rows)<1e-7
    assert max(r['compressed_vs_independent_chi2'] for r in rows)<1e-4
    assert max(r['H0_global_offset_logp_invariance'] for r in rows)<1e-7
    diagnostics=[]
    zz=np.concatenate([np.arange(5)*.001+c for c in [0.,.5,1.]])
    for om,w0,wa in [(.31,-1,0),(.27,-.7,-1.1),(.55,-.3,.8)]:
        p=Provider(om,w0,wa);calc=expansion_diagnostics(zz,p.E(zz)*70)
        for z,label in [(0.,'0'),(.5,'05'),(1.,'1')]:
            w=w0+wa*z/(1+z);f=(1-om)*(1+z)**(3*(1+w0+wa))*np.exp(-3*wa*z/(1+z))/p.E(z)**2
            q=.5+1.5*w*f;qp=1.5*(wa/(1+z)**2*f+w*3*w*f*(1-f)/(1+z));j=q*(2*q+1)+(1+z)*qp
            diagnostics.append({'z':z,'Omega_m':om,'w0':w0,'wa':wa,'q_error':abs(calc['q'+label]-q),'j_error':abs(calc['j'+label]-j)})
    assert max(x['q_error'] for x in diagnostics)<1e-7 and max(x['j_error'] for x in diagnostics)<1e-5
    assert reviewed_hashes=={str(p.relative_to(ROOT)):sha(p) for p in paths},'Reviewed sources changed during calculation'
    result={'status':'passed_bounded_math_review','code_sha256':sha(__file__),'reviewed_source_sha256':snapshots,'reviewed_checkout_source_sha256':reviewed_hashes,'snapshot_scope':'This bounded review applies only to the recorded checkout hashes and preserved source snapshots. Later checkout changes are not automatically reviewed.','SN_geometry':rows,'acceleration_diagnostics':diagnostics,'scope':'Independent FLRW distance integration, samefull1820covariance, properglobaloffset projection, positiveepsilonmeansfainterSNe, Gaussianluminositymodes integratedonce; no CMBspectrallikelihoodreview or convergencecertification.','findings':[],'required_labels':['CPL+flatFLRW and GR/earlyphysics are modelconditions, not nonparametric evidence.','Lategeometry freeH0rdrag removesCMBsoundhorizon information but omitsradiation and retainsmatter+CPLparametrization.','Smoothluminosity modes have fixeddeclaredGaussian priors; they are sensitivityassumptions, not measuredagecorrections.','Likelihood offsetpriorisimproper; no absoluteevidence claim.']}
    (RESULTS/'geometry-review.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'status':result['status'],'max_direct_chi2_error':max(r['shared_vs_independent_direct_chi2'] for r in rows),'max_j_error':max(x['j_error'] for x in diagnostics)}))
if __name__=='__main__':main()
