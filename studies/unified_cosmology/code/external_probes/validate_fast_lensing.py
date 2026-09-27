"""Native-vs-prebinned algebra and stored exact-spectrum closure; no CAMB rerun."""
import os
os.environ.setdefault('CLIPY_NOJAX','1')
import json,time
import numpy as np
from act_dr6_lenslike.act_dr6_lenslike import ACTDR6LensLike,generic_lnlike,standardize,pp_to_kk
from modern_adapter import modern_info
from fast_lensing import FastACTDR6LensLike
from acquire import RESULTS,WORK,PACKAGES,HERE,sha


def raw_cls(path):
    with np.load(path) as f:ell=f['ell'];dls=f['spectra'];expected=f['loglikes']
    ll=ell*(ell+1);f=np.zeros(len(ell));f[1:]=2*np.pi/ll[1:]
    cls={s:v*f for s,v in zip(['tt','ee','bb','te'],dls[:4])};cls['pp']=np.zeros(len(ell));cls['pp'][1:]=dls[4,1:]*2*np.pi/ll[1:]**2;cls['ell']=ell
    return cls,expected


def main():
    key='act_dr6_lenslike.ACTDR6LensLike';opts=modern_info()['likelihood'][key]
    start=time.monotonic();native=ACTDR6LensLike(opts,packages_path=str(PACKAGES));fast=FastACTDR6LensLike(opts,packages_path=str(PACKAGES))
    construction=time.monotonic()-start
    assert np.array_equal(native.data['cov'],fast.data['cov']) and np.array_equal(native.data['cinv'],fast.data['cinv'])
    index=[0,1,2,3,4,32,64,96,128,160,192,224,256,288,320,352]
    rows=[];last=None
    for i in index:
        file=WORK/f'spectral-training/train/{i:04d}.npz';meta=json.loads(file.with_suffix('.json').read_text());cls,logs=raw_cls(file)
        ln,b=generic_lnlike(native.data,cls['ell'],pp_to_kk(cls['pp'],cls['ell']),cls['ell'],cls['tt'],cls['ee'],cls['te'],cls['bb'],native.trim_lmax,return_theory=True)
        kk=standardize(cls['ell'],pp_to_kk(cls['pp'],cls['ell']),native.trim_lmax);cc={s:standardize(cls['ell'],cls[s],native.trim_lmax) for s in ['tt','ee','bb','te']}
        ll,bb=fast.binned_response.loglike(kk,cc,True);wrapper=fast.loglike(cls)
        stored=float(logs[meta['likelihood_names'].index(key)])
        rec={'index':i,'source_sha256':sha(file),'loglike_difference':float(ll-ln),'wrapper_difference':float(wrapper-ln),
            'stored_exact_loglike_difference':float(ll-stored),'max_binned_abs_difference':float(np.max(abs(bb-b))),
            'max_binned_relative_difference':float(np.max(abs((bb-b)/b)))}
        assert abs(ll-ln)<1e-8 and abs(ll-stored)<1e-8 and abs(wrapper-ln)<1e-8,rec
        rows.append(rec);last=cls
    rng=np.random.default_rng(2727751);stress=[]
    for j in range(16):
        cl={k:v.copy() for k,v in last.items()};ell=cl['ell']
        for k in ['tt','ee','bb','te','pp']:
            cl[k]*=1+rng.normal(0,.03)+rng.normal(0,.03)*np.sin(ell/rng.uniform(50,2000))
        a=float(native.loglike(cl));b=float(fast.loglike(cl));assert abs(a-b)<1e-8
        stress.append(b-a)
    timing={}
    for name,like in [('native',native),('prebinned',fast)]:
        like.loglike(last);times=[]
        for _ in range(20):
            t=time.perf_counter();like.loglike(last);times.append(time.perf_counter()-t)
        timing[name]={'median_seconds':float(np.median(times)),'seconds':times}
    out={'status':'passed','active_adapter_changed':False,'scientific_inputs_changed':False,
        'code_sha256':{f:sha(HERE/f) for f in ['fast_lensing.py','validate_fast_lensing.py']},
        'acquisition_sha256':sha(RESULTS/'modern-acquisition.json'),'construction_seconds':construction,
        'actual_frozen_training_spectra':rows,'algebraic_stress_seed':2727751,'algebraic_stress_loglike_differences':stress,
        'native_joint_covariance_and_precision_identical':True,'likelihood_tolerance':1e-8,'timing':timing,
        'scope':'Exact matrix reassociation for fixed native linear responses, checked against unmodified native code and stored exact likelihoods. Stress curves are algebra tests, not physical cosmologies. Unsupported nonlinear/configuration variants are rejected. No approximation or posterior result.'}
    (RESULTS/'fast-lensing-validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
