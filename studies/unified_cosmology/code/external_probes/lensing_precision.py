"""Audit the finite precision of the public native Planck lensing text export.

The original clik binary stores differently normalized bandpowers. Comparing
numbers requires applying the same invertible diagonal row scaling to data,
window functions, and covariance. No installed data are changed by this check.
"""
import numpy as np


def audit(native, official, dls, calibration, official_loglike):
    w = np.asarray(native.bins.binning_matrix[0])
    raw = np.asarray(official.bins[:, native.pcl_lmin:native.pcl_lmax+1])
    scale = np.sum(w*raw,axis=1)/np.sum(raw*raw,axis=1)
    raw_hat = np.asarray(official.pp_hat)*scale
    hat = native.bandpowers.ravel()
    # Every released text bandpower is rounded to a grid of 1e-11.
    assert np.max(abs(hat/1e-11-np.round(hat/1e-11)))<1e-8
    rounding = hat-raw_hat
    assert np.max(abs(rounding))<=5.0001e-12
    native.get_theory_map_cls(dls, {'A_planck':calibration})
    pred=native.get_binned_map_cls(native.map_cls).ravel()
    r=hat-pred
    r_unrounded=raw_hat-pred
    ln_native=float(-.5*r@native.covinv@r)
    ln_unrounded=float(-.5*r_unrounded@native.covinv@r_unrounded)
    # Exact quadratic identity isolates the contribution of rounded bandpowers.
    rounding_shift=float(-r_unrounded@native.covinv@rounding
                         -.5*rounding@native.covinv@rounding)
    assert abs((ln_native-ln_unrounded)-rounding_shift)<1e-12
    remaining=float(ln_unrounded-official_loglike)
    # Derive a rigorous floating-export error bound, including rounded windows,
    # linear responses and inverse covariance (rather than relax a tolerance).
    ell=np.arange(native.pcl_lmax+1,dtype=float)
    dp=np.asarray(dls['pp'][:len(ell)])
    cmb=np.concatenate([np.asarray(dls[k][:len(ell)])/calibration**2
                        for k in ['tt','ee','te']])
    raw_pred=scale*(official.bins@dp-official.cor0+
                   official.cors@np.concatenate([dp,cmb]))
    raw_inv=np.asarray(official.siginv)/scale[:,None]/scale[None,:]
    raw_r=raw_hat-raw_pred
    ln_direct=float(-.5*raw_r@raw_inv@raw_r)
    assert abs(ln_direct-official_loglike)<1e-10
    dr=r_unrounded-raw_r
    di=native.covinv-raw_inv
    bound=float(abs(raw_r@raw_inv@dr)+.5*abs(dr@raw_inv@dr)
                +.5*abs(r_unrounded@di@r_unrounded)+1e-12)
    assert abs(remaining)<=bound
    return {'row_scale':scale.tolist(),'A_planck':calibration,
        'native_minus_official_loglike':ln_native-official_loglike,
        'max_bandpower_rounding':float(np.max(abs(rounding))),
        'text_bandpower_half_grid':5e-12,
        'rounding_loglike_shift':rounding_shift,
        'remaining_window_covariance_rounding_loglike':remaining,
        'derived_remaining_error_bound':bound,
        'native_chi2_reconstruction_difference':float(ln_native-native.log_likelihood(dls,A_planck=calibration)),
        'scope':'Checks public text-export precision; original likelihood and data remain unchanged.'}
