"""Separately validate a lower spectral lmax; never change the declared target.

Pad only multipoles beyond every nonzero likelihood window/response. Preserve
max_eta_k72000 and lens-potential accuracy4 so k coverage is not reduced.
"""
import os
os.environ.setdefault('CLIPY_NOJAX','1')
import json,time
import numpy as np
import camb
from cobaya.model import get_model
from modern_adapter import modern_info
from adapter import reference_point,WORK
from acquire import RESULTS,HERE,sha

class Provider:
    def __init__(self,r,cls):self.r=r;self.cls=cls
    def get_Cl(self,ell_factor=False,units='FIRASmuK2'):
        d={k:v.copy() for k,v in self.cls.items()};ell=d['ell'];ll=ell*(ell+1)
        assert units in ['FIRASmuK2','muK2']
        if not ell_factor:
            for k in ['tt','ee','bb','te']:d[k][1:]*=2*np.pi/ll[1:]
            d['pp'][1:]*=2*np.pi/ll[1:]**2
        return d
    def get_Hubble(self,z,units='km/s/Mpc'):
        v=np.atleast_1d(self.r.hubble_parameter(z));return v/299792.458 if units=='1/Mpc' else v
    def get_angular_diameter_distance(self,z):return np.atleast_1d(self.r.angular_diameter_distance(z))
    def get_comoving_radial_distance(self,z):return np.atleast_1d(self.r.comoving_radial_distance(z))
    def get_param(self,k):return self.r.get_derived_params()[k]


def main():
    m=get_model(modern_info());m.add_requirements({'CAMBdata':None});ref=reference_point(m)
    points=[('reference',ref)]
    for i in [0,32,64,96]:
        rec=json.loads((WORK/f'spectral-training/train/{i:04d}.json').read_text())
        assert rec['status']=='finite_exact'
        points.append((f'train{i:04d}',ref|rec['physical_point']))
    out={'status':'running','target_unchanged':True,'code_sha256':sha(__file__),'adapter_sha256':sha(HERE/'modern_adapter.py'),
         'scope':'Deterministic reference and four already frozen exact training points; precision/runtime test, not posterior result. Lower-lmax spectra padded only beyond all data window support.', 'rows':[]}
    # Native ACT asks9001 and indexes window columns to8501; their weights above6325 are exactlyzero.
    act=m.likelihood['act_dr6_cmbonly.ACTDR6CMBonly']
    maxused=max(int(meta['window'].values[np.where(np.any(meta['window'].weight!=0,axis=1))[0]].max()) for meta in act.spec_meta)
    assert maxused==6325
    out['ACT_last_nonzero_window_ell']=maxused
    for label,point in points:
        t=time.monotonic();baseline=m.logposterior(point);reference_seconds=time.monotonic()-t
        original={k:v.provider for k,v in m.likelihood.items()}
        extra=dict(modern_info()['theory']['camb']['extra_args']);extra.update(m.theory['camb'].extra_args)
        physical={k:point[k] for k in ['H0','ombh2','omch2','ns','tau','w','wa']};physical['As']=1e-10*np.exp(point['logA'])
        params=m.parameterization.to_input(point)
        n=len(m.provider.get_Cl(ell_factor=True)['ell'])
        for lmax in [9001,6500,6325]:
            t=time.monotonic();pars=camb.set_params(**physical,**(extra|{'lmax':lmax,'max_eta_k':72000.}))
            r=camb.get_results(pars);native=r.get_cmb_power_spectra(lmax=min(n-1,lmax+100),CMB_unit='muK')
            cls={k:np.zeros(n) for k in ['tt','ee','bb','te','pp']};cls['ell']=np.arange(n)
            for k,j in [('tt',0),('ee',1),('bb',2),('te',3)]:cls[k][:len(native['total'])]=native['total'][:,j]
            cls['pp'][:len(native['lens_potential'])]=native['lens_potential'][:,0]
            provider=Provider(r,cls);values={}
            try:
                for name,like in m.likelihood.items():
                    like.provider=provider;values[name]=float(like.logp(**{k:params[k] for k in like.input_params}))
            finally:
                for name,like in m.likelihood.items():like.provider=original[name]
            delta={k:values[k]-float(baseline.loglikes[j]) for j,k in enumerate(m.likelihood)}
            row={'point':label,'requested_lmax':lmax,'CAMB_max_l':pars.max_l,'max_eta_k':pars.max_eta_k,
                 'seconds':time.monotonic()-t,'native_baseline_seconds':reference_seconds,'component_differences':delta,
                 'total_loglike_difference':sum(delta.values())}
            if lmax==9001:assert max(abs(v) for v in delta.values())<1e-6,delta
            out['rows'].append(row);(RESULTS/'theory-trim-check.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(row),flush=True)
    out['status']='completed_not_adopted';(RESULTS/'theory-trim-check.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
