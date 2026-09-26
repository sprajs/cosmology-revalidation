#!/usr/bin/env python3
"""Fresh checks of physical extinction identities and selected local implementations.

Run: phase2/env-official/bin/python scripts/physics_audit/extinction_check.py
Independent equations below are compared with freshly compiled, unmodified
SNANA function bodies and extinction.py. This is not a simulation/refit.
"""
from pathlib import Path
import hashlib
import json
import sys
import tempfile
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.integrate import quad
from scipy.optimize import brentq
from astropy.io import fits
import extinction

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts/salt_dust_audit'))
from snana_extinction import SnanaExtinction, build_library, CURRENT, HISTORICAL
sys.path.insert(0, str(ROOT / 'scripts/standardization'))
from dust_identifiability import slab


def f99(wave, rv, ebv=1.):
    """F99/FM_UNRED convention; coefficient normalization is E times (R+k)."""
    wave = np.asarray(wave, dtype=float)
    x = 10000 / wave
    c2 = -.824 + 4.717 / rv
    c1 = 2.03 - 3.007 * c2

    def uv(q):
        y = np.maximum(q - 5.9, 0.)
        return c1 + c2*q + 3.23*q*q/((q*q-4.596**2)**2 + q*q*.99**2) + .41*(.5392*y*y+.05644*y**3)

    knots = np.array([0, 1/2.65, 1/1.22, 1/.6, 1/.547, 1/.467, 1/.411, 1/.270, 1/.260])
    # Use full optical coefficients in the locally implemented FM_UNRED family.
    values = np.array([-rv, -.914616129*rv, -.7325*rv,
        -.422809+.00270*rv+2.13572e-4*rv**2,
        -.051354+.00216*rv-7.35778e-5*rv**2,
        .700127+.00184*rv-3.32598e-5*rv**2,
        1.19456+.01707*rv-5.46959e-3*rv**2+7.97809e-4*rv**3-4.45636e-5*rv**4,
        uv(knots[-2]), uv(knots[-1])])
    k = np.where(wave <= 2700, uv(x), CubicSpline(knots, values, bc_type='natural')(x))
    return ebv*(rv+k)


def ccm(wave, rv, ebv=1., odonnell=False):
    x = 10000/np.asarray(wave, float)
    a, b = np.zeros_like(x), np.zeros_like(x)
    ir = (x >= .3) & (x < 1.1)
    a[ir], b[ir] = .574*x[ir]**1.61, -.527*x[ir]**1.61
    op = (x >= 1.1) & (x < 3.3)
    y = x[op]-1.82
    ca = [1,.104,-.609,.701,1.137,-1.718,-.827,1.647,-.505] if odonnell else [1,.17699,-.50447,-.02427,.72085,.01979,-.77530,.32999]
    cb = [0,1.952,2.908,-3.989,-7.985,11.102,5.491,-10.805,3.347] if odonnell else [0,1.41338,2.28305,1.07233,-5.38434,-.62251,5.30260,-2.09002]
    a[op], b[op] = np.polynomial.polynomial.polyval(y,ca), np.polynomial.polynomial.polyval(y,cb)
    uv = (x >= 3.3) & (x < 8)
    q = x[uv]; y = np.maximum(q-5.9, 0)
    a[uv] = 1.752-.316*q-.104/((q-4.67)**2+.341)-.04473*y*y-.009779*y**3
    b[uv] = -3.09+1.825*q+1.206/((q-4.62)**2+.263)+.213*y*y+.1207*y**3
    far = (x >= 8) & (x <= 10); y=x[far]-8
    a[far] = np.polynomial.polynomial.polyval(y,[-1.073,-.628,.137,-.070])
    b[far] = np.polynomial.polynomial.polyval(y,[13.670,4.257,-.420,.374])
    return ebv*(rv*a+b)


def old99(wave, rv, ebv=1.):
    coeff=[.0855929205,1.91547833,-1.65101945,.750611119,-.200041118,
           .0330155576,-.00346344458,.000230741420,-.00000943018242,
           .000000214917977,-.00000000208276810]
    return ccm(wave,rv,ebv,odonnell=True)*np.polynomial.polynomial.polyval(np.asarray(wave)/1000,coeff)


def main():
    errors={'fresh_C_F99_vs_equation':0., 'extinction_py_F99_vs_equation':0.,
            'fresh_C_CCM_vs_equation':0., 'fresh_C_O94_vs_equation':0.,
            'fresh_C_historical99_vs_equation':0., 'current_minus99_vs_historical99':0.}
    waves=np.unique(np.r_[np.linspace(1000,10000,1801),2700.,10000/3.3,10000/1.1])
    rvs=[.4,.7,1.,1.5,2.,3.1,4.,6.]
    with tempfile.TemporaryDirectory(prefix='supernova-extinction-') as temp:
        lib=SnanaExtinction(build_library(Path(temp)))
        for rv in rvs:
            own=f99(waves,rv)
            comparisons={
              'fresh_C_F99_vs_equation':lib(waves,rv)-own,
              'extinction_py_F99_vs_equation':extinction.fitzpatrick99(waves,rv,rv)-own,
              'fresh_C_CCM_vs_equation':lib(waves,rv,option=89)-ccm(waves,rv),
              'fresh_C_O94_vs_equation':lib(waves,rv,option=94)-ccm(waves,rv,odonnell=True),
              'fresh_C_historical99_vs_equation':lib(waves,rv,historical=True)-old99(waves,rv),
              'current_minus99_vs_historical99':lib(waves,rv,option=-99)-lib(waves,rv,historical=True)}
            for key,value in comparisons.items(): errors[key]=max(errors[key],float(np.max(abs(value))))
        assert max(errors.values())<1e-7, errors
    exactroot=brentq(lambda r:float(f99(8000,r)),.1,2)
    oldroot=brentq(lambda r:float(old99(np.array([8000.]),r)[0]),.1,2)
    # Derive a local broadband derivative from a positive test SED and throughput.
    w=np.linspace(4000,8000,20001); weight=w*(w/5500)**-2
    k=f99(w,3.1); normalize=np.trapezoid(weight,w)
    def amag(e):return -2.5*np.log10(np.trapezoid(weight*10**(-.4*e*k),w)/normalize)
    epsilon=1e-4
    slope=(amag(epsilon)-amag(-epsilon))/(2*epsilon)
    mean=float(np.trapezoid(weight*k,w)/normalize)
    second=(amag(epsilon)+amag(-epsilon)-2*amag(0))/epsilon**2
    var=float(np.trapezoid(weight*(k-mean)**2,w)/normalize)
    assert abs(slope-mean)<1e-8
    assert abs(second+.4*np.log(10)*var)<1e-6
    tau=np.array([0.,1e-8,.1,1.,10.])
    transmission=10**(-.4*slab(tau))
    errslab=max(abs(quad(lambda u:np.exp(-t*u),0,1)[0]-s) for t,s in zip(tau,transmission))
    assert errslab<1e-12
    assert slab(0.)==0.
    tiny=np.array([1e-16,1e-12,1e-8])
    thin_error=float(np.max(abs(slab(tiny)/(2.5/np.log(10)*tiny/2)-1)))
    assert thin_error<1e-8
    positive=np.array([.01,.1,.3,1.,3.,10.])
    old_positive=-2.5*np.log10(-np.expm1(-positive)/positive)
    assert np.array_equal(slab(positive),old_positive)
    rejected=[]
    for bad in [-1.,np.inf,-np.inf,np.nan]:
        try:slab(bad)
        except ValueError:rejected.append(str(bad))
        else:raise AssertionError(f'Invalid optical depth accepted: {bad}')
    passive_wave=np.linspace(1000,35000,3401)
    minphysical=min(float(f99(passive_wave,r,.1).min()) for r in np.linspace(2,6,17))
    assert minphysical>0
    # Re-read stored simulation truth; these are simulated objects, not dust data.
    base=ROOT/'sources/repos/des-science__DES-SN5YR@1.3/1_SIMULATIONS/SNIa_SIMULATIONS'
    files=sorted(base.glob('*/*_HEAD.FITS.gz'))
    totals={'files':len(files),'rows':0,'RV_lt2':0,'old_A8000_negative':0,'exact_A8000_negative':0}
    example=None
    for p in files:
        with fits.open(p,memmap=False) as hd:
            d=hd[1].data;rv=np.asarray(d['SIM_RV'],float);av=np.asarray(d['SIM_AV'],float)
            assert np.all(rv>0) and np.all(av>=0)
            totals['rows']+=len(rv);totals['RV_lt2']+=int(np.sum(rv<2))
            old_values=old99(8000.,rv,av/rv)
            exact_values=f99(8000.,rv,av/rv)
            assert np.array_equal(old_values<0,(rv<oldroot)&(av>0))
            assert np.array_equal(exact_values<0,(rv<exactroot)&(av>0))
            # Direct equations evaluated on every stored truth pair; no full
            # SNANA light curves or stored broadband fluxes are asserted.
            totals['old_A8000_negative']+=int(np.sum(old_values<0))
            totals['exact_A8000_negative']+=int(np.sum(exact_values<0))
            if example is None and np.any((rv<oldroot)&(av>0)):
                j=np.flatnonzero((rv<oldroot)&(av>0))[0]
                example={'path':str(p.relative_to(ROOT)),'SNID':str(d['SNID'][j]).strip(),
                  'RV':float(rv[j]),'AV':float(av[j]),
                  'old_A8000_mag':float(old99(np.array([8000.]),rv[j],av[j]/rv[j])[0]),
                  'exact_A8000_mag':float(f99(8000.,rv[j],av[j]/rv[j]))}
    result={'scope':'Equation checks and sampled extinction routines only; no full vendor, data-fit, selection, BBC or cosmology verification.',
      'equation_max_absolute_errors_mag_per_E':errors,'comparison_grid':{'wavelength_min_A':1000,'wavelength_max_A':10000,'n_wavelengths':len(waves),'RV':rvs},
      'A8000_zero_RV':{'historical':oldroot,'F99_spline':exactroot},
      'example_E0p1_RV0p4':{'A6500_mag':float(f99(6500,.4,.1)),'A8000_mag':float(f99(8000,.4,.1)),
          'flux_multiplier_8000':float(10**(-.4*f99(8000,.4,.1)))},
      'F99_RV3p1_A5495_per_input_AV':float(f99(5495.,3.1)/3.1),
      'broadband_derivative':{'finite_difference':slope,'weighted_k':mean,'second_derivative':second,
          'minus_K_variance_k':-.4*np.log(10)*var,'A_E0p1':amag(.1),'linear_A_E0p1':.1*mean},
      'host_vs_MW_at_obs8000_z0p5_E0p1_RV3p1':{'host_mag':float(f99(8000/1.5,3.1,.1)),'MW_mag':float(f99(8000,3.1,.1))},
      'slab_quadrature_max_error':float(errslab),'slab_zero_transmission_limit':float(transmission[0]),
      'slab_regression':{'zero_magnitude':float(slab(0.)),'thin_limit_max_relative_error':thin_error,
                         'old_positive_grid_bitwise_unchanged':True,'rejected_nonphysical_inputs':rejected},
      'positive_grid_min_A_mag':minphysical,'simulation_truth_counts':totals,'simulation_truth_example':example,
      'versions':{'numpy':np.__version__,'extinction':extinction.__version__}}
    inputs=[Path(__file__),ROOT/'scripts/standardization/dust_identifiability.py',ROOT/'scripts/salt_dust_audit/snana_extinction.py',CURRENT/'MWgaldust.c',CURRENT/'MWgaldust.h',HISTORICAL/'MWgaldust.c',HISTORICAL/'MWgaldust.h',*files]
    result['input_sha256']={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    out=ROOT/'runs/physics_audit/extinction-check.json';out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='input_sha256'},indent=2))


if __name__=='__main__':main()
