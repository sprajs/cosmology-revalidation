"""Refine source-convention comparison with common, merged-knot integration.
Piecewise-linear spectral/filter definitions are explicit; no SN flux outcomes.
"""
from pathlib import Path
import csv,json
import numpy as np
from astropy.io import fits
from ruamel.yaml import YAML
from calibration_bridge import OUT,NATIVE,RELEASE,KCOR


def integral(T, S):
    """Integral lambda*T(lambda)*S(lambda);4pointGauss on merged knots.

    Exactly integrates the piecewise cubic when both tabulated inputs are linear.
    AB S uses analytic lambda^-2 and four-point Gaussian quadrature.
    """
    if S is None:x=T[:,0]
    else:x=np.unique(np.r_[T[:,0],S[(S[:,0]>T[0,0])&(S[:,0]<T[-1,0]),0]])
    u,w=np.polynomial.legendre.leggauss(4)
    mid=(x[1:]+x[:-1])/2;half=np.diff(x)/2
    lam=mid[:,None]+half[:,None]*u
    trans=np.interp(lam,T[:,0],T[:,1],left=0,right=0)
    spec=2.99792458e18/lam**2*10**(-.4*48.6) if S is None else np.interp(lam,S[:,0],S[:,1])
    return np.sum(half[:,None]*w*lam*trans*spec)


def norm(T):
    return integral(T,np.array([[T[0,0],1.],[T[-1,0],1.]]))


def shape_l1(A,B):
    # Include native subgrid structure that is invisible on the10A released grid.
    x=np.unique(np.r_[A[:,0],B[:,0]])
    d=np.interp(x,A[:,0],A[:,1],left=0,right=0)/norm(A)-np.interp(x,B[:,0],B[:,1],left=0,right=0)/norm(B)
    cross=d[:-1]*d[1:]<0
    roots=x[:-1][cross]-d[:-1][cross]*np.diff(x)[cross]/np.diff(d)[cross]
    x=np.unique(np.r_[x,roots]);u,w=np.polynomial.legendre.leggauss(2)
    mid=(x[1:]+x[:-1])/2;half=np.diff(x)/2;l=mid[:,None]+half[:,None]*u
    d=np.interp(l,A[:,0],A[:,1],left=0,right=0)/norm(A)-np.interp(l,B[:,0],B[:,1],left=0,right=0)/norm(B)
    return np.sum(half[:,None]*w*l*abs(d))/2


def run():
    cfg=YAML(typ='safe').load((NATIVE/'filters.yaml').read_text());rows=list(csv.DictReader((OUT/'calibration-bridge.csv').open()))
    for r in rows:
        T=np.loadtxt(OUT/'release-filters'/f"{r['custom_filter']}.dat")
        oldS=np.loadtxt(r['release_reference_path']);newT=np.loadtxt(r['native_filter_path'])
        new=cfg['filters'][r['native_filter']]
        if new['magsys']=='ab':newS=None
        else:
            sp=NATIVE/cfg['standards'][new['magsys']]['path']
            if sp.suffix=='.fits':
                with fits.open(sp) as h:newS=np.column_stack([h[1].data['WAVELENGTH'],h[1].data['FLUX']])
            else:newS=np.loadtxt(sp)
        oldzero=integral(T,oldS)/norm(T)*10**(.4*float(r['release_primary_mag']))
        newzero=integral(newT,newS)/norm(newT)*10**(.4*float(r['native_primary_mag']))
        new_oldcurve=integral(T,newS)/norm(T)*10**(.4*float(r['native_primary_mag']))
        r.update(exact_piecewise_linear_zero_release=oldzero,exact_piecewise_linear_zero_native=newzero,
                 exact_reference_delta_mag=2.5*np.log10(newzero/oldzero),
                 exact_same_curve_delta_mag=2.5*np.log10(new_oldcurve/oldzero),
                 photon_weight_half_L1_merged_knots=shape_l1(T,newT),
                 released_zero_simpson_minus_exact_mag=2.5*np.log10(float(r['zero_mean_release'])/oldzero),
                 native_zero_simpson_minus_exact_mag=2.5*np.log10(float(r['zero_mean_native'])/newzero))
    with (OUT/'calibration-bridge-exact.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    result={'definition':'Merged-knot piecewise-linear filter and standard;4pointGauss per interval. Source Simpson comparison retained separately.',
            'CSP_B':{k:v for k,v in rows[0].items() if k.startswith('exact') or 'simpson_minus' in k},
            'max_released_zero_integration_difference_mmag':1000*max(abs(r['released_zero_simpson_minus_exact_mag']) for r in rows),
            'max_native_zero_integration_difference_mmag':1000*max(abs(r['native_zero_simpson_minus_exact_mag']) for r in rows),
            'max_supported_released_zero_integration_difference_mmag':1000*max(abs(r['released_zero_simpson_minus_exact_mag']) for r in rows if r['scientific_excluded']=='False'),
            'status':'No changes to observed photometry or model-training assets; differences are numerical conventions, not measured calibration biases.'}
    (OUT/'reference-integral-gate.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':run()
