"""Independent observable arithmetic and cross-implementation CMB checks."""
import os
os.environ.setdefault('CLIPY_NOJAX', '1')
import json
import time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.linalg import solve_triangular
from cobaya.model import get_model
import clipy
from adapter import external_info, reference_point, WORK, PACKAGES, ROOT
from acquire import sha, HERE, RESULTS
from lensing_precision import audit as audit_lensing_precision


def vector(likelihood, cls, params):
    names = ['tt','ee','bb','te','tb','eb']
    if len(likelihood.lmax) == 7:
        names = ['pp'] + names
    return np.concatenate([cls[name][:int(lmax)+1]
                           for name,lmax in zip(names,likelihood.lmax) if lmax>=0]
                          + [np.array([params[name] for name in likelihood.extra_parameter_names])])


def main():
    started = time.monotonic()
    out = {'status':'passed', 'comparison_scope':'Actual CMB spectra/lensing, not a compressed distance prior.',
           'reference_is_not_MAP':True, 'official_check_vectors':{}, 'cross_implementation':[],
           'bao':{}, 'theory_probes':[], 'benchmarks':{}}
    baseline = PACKAGES/'data/planck_2018/baseline/plc_3.0'
    paths = {
        'full':'hi_l/plik/plik_rd12_HM_v22b_TTTEEE.clik',
        'lite':'hi_l/plik_lite/plik_lite_v22_TTTEEE.clik',
        'TT':'low_l/commander/commander_dx12_v3_2_29.clik',
        'EE':'low_l/simall/simall_100x143_offlike5_EE_Aplanck_B.clik',
        'lensing':'lensing/smicadx12_Dec5_ftl_mv2_ndclpp_p_teb_consext8.clik_lensing',
    }
    clik = {}
    for name,rel in paths.items():
        p = baseline/rel
        like = clipy.clik(str(p))
        clik[name] = like
        raw = clipy.cldf.open(str(p))
        if name == 'lensing':
            # Lensing uses its own initialization checksum rather than CMB check_param.
            out['official_check_vectors'][name] = {'scope':'Compared to independent native likelihood below.'}
            continue
        point = np.array(raw['clik']['check_param'])
        expected = float(raw['clik']['check_value'])
        actual = float(like(point))
        difference = actual-expected
        assert abs(difference)<1e-4, (name,difference)
        out['official_check_vectors'][name] = {'expected_loglike':expected,
            'actual_loglike':actual, 'difference':difference}
    info = external_info('lite','cpl')
    model = get_model(info)
    model.add_requirements({'CAMBdata':None,'omegam':None})
    point = reference_point(model)
    point['A_planck'] = 1.
    result = model.logposterior(point)
    assert np.isfinite(result.logpost)
    cls = model.provider.get_Cl(units='FIRASmuK2')
    dls = model.provider.get_Cl(ell_factor=True)
    for name,key in [('lite','planck_2018_highl_plik.TTTEEE_lite_native'),
                     ('TT','planck_2018_lowl.TT'), ('EE','planck_2018_lowl.EE'),
                     ('lensing','planck_2018_lensing.native')]:
        for calibration in [.997,1.,1.003]:
            parameters = {'A_planck':calibration}
            direct = float(clik[name](vector(clik[name],cls,parameters)))
            like = model.likelihood[key]
            if name == 'lite':
                native = -.5*like.get_chi_squared(0,dls['tt'],dls['te'],dls['ee'],calibration)
            elif name in ('TT','EE'):
                native = like.log_likelihood(dls[name.lower()],calibration)
            else:
                native = like.log_likelihood(dls,**parameters)
            difference = float(direct-native)
            precision = None
            if name == 'lensing':
                precision = audit_lensing_precision(like,clik[name],dls,calibration,direct)
            else:
                assert abs(difference)<1e-4, (name,calibration,difference)
            out['cross_implementation'].append({'likelihood':name,'A_planck':calibration,
                'clipy_loglike':direct,'native_loglike':float(native),'difference':difference,
                'finite_export_precision_audit':precision})
    # Rebuild DESI observables using an independent dense Cholesky calculation.
    base = PACKAGES/'data/bao_data/desi_bao_dr2'
    mean = base/'desi_gaussian_bao_ALL_GCcomb_mean.txt'
    covariance = base/'desi_gaussian_bao_ALL_GCcomb_cov.txt'
    data = pd.read_csv(mean,sep=r'\s+',comment='#',names=['z','value','kind'])
    cov = np.loadtxt(covariance)
    assert cov.shape==(13,13) and np.array_equal(cov,cov.T)
    eig = np.linalg.eigvalsh(cov)
    assert eig.min()>0
    cmb = model.provider.get_CAMBdata()
    rd = float(model.provider.get_param('rdrag'))
    prediction=[]
    errors=[]
    for r in data.itertuples():
        dh = 299792.458/cmb.hubble_parameter(r.z)
        dm = (1+r.z)*cmb.angular_diameter_distance(r.z)
        independent_dm = quad(lambda z:299792.458/cmb.hubble_parameter(z),0,r.z,
                              epsabs=1e-9,epsrel=1e-11)[0]
        errors.append(abs(dm-independent_dm))
        prediction.append({'DM_over_rs':dm/rd,'DH_over_rs':dh/rd,
                           'DV_over_rs':(r.z*dm**2*dh)**(1/3)/rd}[r.kind])
    residual = data.value.to_numpy()-prediction
    whitened = solve_triangular(np.linalg.cholesky(cov),residual,lower=True)
    chi2 = float(whitened@whitened)
    released = float(-2*result.loglikes[list(model.likelihood).index('bao.desi_dr2')])
    # Scalar CAMB background calls and Cobaya's vector distances round slightly
    # differently; require 1e-8 relative precision in this independent chi2.
    assert abs(chi2-released)<1e-8*max(1,abs(released)), (chi2,released)
    assert max(errors)<2e-5
    block_errors=[]
    for p in sorted(base.glob('*_mean.txt')):
        if p==mean: continue
        d=pd.read_csv(p,sep=r'\s+',comment='#',names=['z','value','kind'])
        inds=[int(np.flatnonzero((data.z==r.z)&(data.kind==r.kind))[0]) for r in d.itertuples()]
        c=np.atleast_2d(np.loadtxt(p.with_name(p.name.replace('_mean','_cov'))))
        block_errors.append(float(np.max(abs(c-cov[np.ix_(inds,inds)]))))
        assert np.array_equal(d.value.to_numpy(),data.value.to_numpy()[inds])
    assert max(block_errors)<1e-8
    out['bao']={'dimension':13,'mean_sha256':sha(mean),'covariance_sha256':sha(covariance),
        'covariance_min_eigenvalue':float(eig.min()),'rdrag_Mpc':rd,
        'H0_km_s_Mpc':float(point['H0']),'omegam':float(model.provider.get_param('omegam')),
        'prediction':prediction,'chi2_independent':chi2,'chi2_cobaya':released,
        'chi2_difference':chi2-released,
        'max_distance_quadrature_difference_Mpc':max(errors),
        'max_individual_block_covariance_difference':max(block_errors),
        'cross_redshift_covariance_max_abs':float(np.max(abs(cov[data.z.to_numpy()[:,None]!=data.z.to_numpy()[None,:]])))}
    # Explicitly probe crossing and early-DE domains without adding a hidden prior.
    for w,wa in [(-1.,0.),(-.8,-.6),(-1.2,.4),(-.5,.7)]:
        p=dict(point,w=w,wa=wa)
        t=time.monotonic()
        try:
            z=model.logposterior(p)
            valid=bool(np.isfinite(z.logpost))
            out['theory_probes'].append({'w':w,'wa':wa,'status':'finite' if valid else 'nonfinite_theory_or_likelihood',
                'loglikes':[float(v) if np.isfinite(v) else None for v in z.loglikes],
                'seconds':time.monotonic()-t})
        except Exception as e:
            out['theory_probes'].append({'w':w,'wa':wa,'status':'theory_error',
                                        'exception':type(e).__name__,'message':str(e)})
    for mode in ['full','lite']:
        m=get_model(external_info(mode,'cpl'))
        p=reference_point(m)
        r=m.logposterior(p)
        assert np.isfinite(r.logpost)
        pieces={}
        for key in ['H0','A_planck']+(['A_cib_217'] if mode=='full' else []):
            m.logposterior(p) # exclude changing back to a different cosmology
            times=[]
            for i in range(5):
                q=dict(p)
                q[key] += (i+1)*(.0001 if key=='A_planck' else .1)
                t=time.perf_counter(); v=m.logposterior(q); times.append(time.perf_counter()-t)
                assert np.isfinite(v.logpost)
            pieces[key]={'median_seconds':float(np.median(times)),'seconds':times}
        again=m.logposterior(p)
        repeat=m.logposterior(p)
        assert np.array_equal(r.loglikes,again.loglikes) and np.array_equal(r.loglikes,repeat.loglikes)
        out['benchmarks'][mode]={'timing':pieces,'reference_point':{k:float(v) for k,v in p.items()},
            'reference_loglikes':dict(zip(m.likelihood,[float(v) for v in r.loglikes])),
            'exact_repeat':True,'sampled_parameter_count':len(p)}
    out['seconds']=time.monotonic()-started
    out['code_sha256']={p.name:sha(p) for p in HERE.glob('*.py')}
    out['design_sha256']=sha(HERE/'design.json')
    out['acquisition_sha256']=sha(RESULTS/'acquisition.json')
    (RESULTS/'validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
