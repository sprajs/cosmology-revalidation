"""Amplitude profiling kernel for a separately verified native mean oracle.

No observer model is implemented here. Native h(theta) at a reference DLMAG,
the exact exported frozen inverse covariance and measurement vector are
required. The analytic result is eligible only after native multiplicative
distance closure and flat-distance-prior gates. This module does not run
native fits, refresh covariance or infer a dust posterior.
"""
import hashlib
import json
from pathlib import Path
import numpy as np


def profile_amplitude(y, inverse_covariance, h, reference_dlmag,
                      flat_dlmag_bounds=(10.0, 60.0), other_prior=0.0):
    y=np.asarray(y,dtype=float);h=np.asarray(h,dtype=float)
    w=np.asarray(inverse_covariance,dtype=float)
    if y.ndim!=1 or h.shape!=y.shape or w.shape!=(len(y),len(y)):
        raise ValueError('Mismatched fixed row identities/dimensions')
    if not all(np.isfinite(a).all() for a in [y,h,w]):
        raise ValueError('Nonfinite inputs')
    if not np.allclose(w,w.T,rtol=1e-12,atol=1e-12):
        raise ValueError('Inverse covariance not symmetric at declared tolerance')
    np.linalg.cholesky(w)
    q=float(h@w@h);b=float(h@w@y)
    if q<=0:raise ValueError('No positive amplitude information')
    a=b/q
    if a<=0:
        return {'status':'requires_native_bounded_distance_profile','reason':'unconstrained amplitude nonpositive','a_unconstrained':a,'q':q,'b':b}
    distance=float(reference_dlmag-2.5*np.log10(a))
    if not flat_dlmag_bounds[0]<=distance<=flat_dlmag_bounds[1]:
        return {'status':'requires_native_bounded_distance_profile','reason':'distance outside verified flat-prior interval','a_unconstrained':a,'DLMAG_unconstrained':distance,'q':q,'b':b}
    residual=y-a*h
    return {'status':'analytic_interior','amplitude':a,'DLMAG':distance,
            'data_quadratic':float(residual@w@residual),
            'full_frozen_objective':float(residual@w@residual+other_prior),
            'amplitude_curvature_half':q,'amplitude_score_half':float(h@w@residual)}


def self_check():
    # Algebra-only test vectors, not an alternative SN model.
    h=np.array([.15,.4,1.,1.5,1.1,.6,.25])
    y=1.3*h
    c=np.diag(np.array([.04,.08,.06,.1,.09,.07,.05])**2)+np.outer(.02*h,.02*h)
    w=np.linalg.inv(c)
    r=profile_amplitude(y,w,h,40.)
    perm=np.array([4,0,6,1,5,2,3])
    rp=profile_amplitude(y[perm],w[np.ix_(perm,perm)],h[perm],40.)
    scale=np.array([.7,1.2,.9,1.1,1.3,.8,.95])
    ru=profile_amplitude(scale*y,w/np.outer(scale,scale),scale*h,40.)
    check={'synthetic_amplitude_error':abs(r['amplitude']-1.3),
           'direct_quadratic':r['data_quadratic'],
           'permutation_distance_error':abs(rp['DLMAG']-r['DLMAG']),
           'flux_unit_transform_distance_error':abs(ru['DLMAG']-r['DLMAG']),
           'negative_amplitude_rejected':profile_amplitude(-y,w,h,40.)['status']=='requires_native_bounded_distance_profile',
           'outside_flat_prior_rejected':profile_amplitude(1e-15*y,w,h,40.)['status']=='requires_native_bounded_distance_profile',
           'scope':'Algebra-only kernel check; native oracle and full state closure remain required',
           'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    assert check['synthetic_amplitude_error']<1e-14
    assert check['direct_quadratic']<1e-24
    assert check['permutation_distance_error']<1e-13
    assert check['flux_unit_transform_distance_error']<1e-13
    assert check['negative_amplitude_rejected'] and check['outside_flat_prior_rejected']
    Path(__file__).with_name('profile-kernel-check.json').write_text(json.dumps(check,indent=2)+'\n')
    print(json.dumps(check,indent=2))


if __name__=='__main__':self_check()
