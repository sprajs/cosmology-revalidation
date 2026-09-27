"""Host likelihood pieces and explicit luminosity-evolution identification.

These functions do not turn posterior medians into observed physical ages and
do not export a brightness prior independent of overlapping supernova data.
"""
from __future__ import annotations
import numpy as np
from scipy.linalg import solve_triangular
from scipy.special import logsumexp


def host_flux_loglike(observed,predicted,covariance):
    """Gaussian flux likelihood on supplied valid bands, including negatives.

    The caller must label a diagonal catalogue-error covariance as a working
    assumption when aperture/deblending/calibration covariance is unavailable.
    Native error floors and MW corrections must not be applied a second time.
    """
    y,m,c=map(np.asarray,(observed,predicted,covariance))
    if y.ndim!=1 or m.shape!=y.shape or c.shape!=(len(y),len(y)):
        raise ValueError('One host vector and its full covariance are required')
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(m)) or not np.all(np.isfinite(c)):
        raise ValueError('Select valid bands before constructing the likelihood')
    chol=np.linalg.cholesky(c)
    delta=solve_triangular(chol,y-m,lower=True)
    return float(-.5*(delta@delta+2*np.log(np.diag(chol)).sum()+len(y)*np.log(2*np.pi)))


def conditional_brightness_loglike(residual,variance,joint_host_draws,coefficients,*,
                                  host_model_unchanged=False,log_population_ratio=None,
                                  shared_systematic_shift=None):
    """Integrate each brightness likelihood over same-row host posterior draws.

    This is p(SN data | observed host flux, FIXED author host model), not a
    prior-free host likelihood. A changed population needs log(p_new/p_old)
    at every original draw, including all selection/support factors. Ratios
    are *not* normalized per object: their evidence contribution is retained.
    Rows must be distinct physical hosts. Shared SN covariance must be carried
    by explicit conditioned systematic shifts, not discarded or diagonalized.
    Selection normalization for the SN sample is not supplied by this function.
    """
    r,v,h,b=map(np.asarray,(residual,variance,joint_host_draws,coefficients))
    if v.ndim!=1 or r.ndim!=1 or v.shape!=r.shape or h.ndim!=3 or h.shape[0]!=len(r) or h.shape[2]!=len(b):
        raise ValueError('Expected residual[N], variance[N], host_draws[N,S,P], coefficients[P]')
    if np.any(v<=0) or not all(np.all(np.isfinite(x)) for x in [r,v,h,b]):raise ValueError('Invalid numerical inputs')
    if not host_model_unchanged and log_population_ratio is None:
        raise ValueError('Changed host priors require the full joint population density ratio')
    if shared_systematic_shift is not None:
        systematic=np.asarray(shared_systematic_shift)
        if systematic.shape!=r.shape or not np.all(np.isfinite(systematic)):
            raise ValueError('Shared systematic shift must have the residual vector shape')
        r=r-systematic
    ratio=np.zeros(h.shape[:2]) if log_population_ratio is None else np.asarray(log_population_ratio)
    if ratio.shape!=h.shape[:2] or np.any(np.isnan(ratio)) or np.any(np.isposinf(ratio)):
        raise ValueError('Population ratio must be finite or negative infinity on original posterior draws')
    shift=np.einsum('nsp,p->ns',h,b)
    terms=-.5*((r[:,None]-shift)**2/v[:,None]+np.log(2*np.pi*v[:,None]))+ratio
    per_object=logsumexp(terms,axis=1)-np.log(h.shape[1])
    normalization=logsumexp(ratio,axis=1)
    weights=np.zeros_like(ratio,dtype=float)
    supported=np.isfinite(normalization)
    weights[supported]=np.exp(ratio[supported]-normalization[supported,None])
    ess=np.zeros(len(r))
    ess[supported]=1/np.sum(weights[supported]**2,axis=1)
    return dict(loglike=float(per_object.sum()),per_object=per_object,
                population_ratio_ess=ess)


def profile_linear_gaussian(y,covariance,design):
    """Profile linear means with complete supplied covariance and rank audit."""
    y,c,x=map(np.asarray,(y,covariance,design))
    chol=np.linalg.cholesky(c)
    yw=solve_triangular(chol,y,lower=True);xw=solve_triangular(chol,x,lower=True)
    coefficient,_,rank,singular=np.linalg.lstsq(xw,yw,rcond=None)
    residual=yw-xw@coefficient
    return dict(chi2=float(residual@residual),coefficient=coefficient,rank=int(rank),singular_values=singular,
                logdet_covariance=float(2*np.log(np.diag(chol)).sum()))


def within_between_design(z,host_value,edges):
    """Separate within-bin association from unrestricted bin mean brightness."""
    z,a=np.asarray(z),np.asarray(host_value)
    bins=np.digitize(z,edges)-1
    if np.any((bins<0)|(bins>=len(edges)-1)):raise ValueError('Redshift outside declared bins')
    active=np.unique(bins);intercept=np.column_stack([bins==i for i in active]).astype(float)
    means=np.array([np.mean(a[bins==i]) for i in active])
    within=a-intercept@means
    return dict(within=within,between=intercept,bin_host_mean=means,active_bins=active)


def grey_compensator(mu_trial,mu_reference,reference_weights):
    """Return B and intercept shift exactly cancelling a trial distance change.

    mu_trial + B + intercept_shift == mu_reference. B has zero reference mean.
    Host-only photometry does not independently rule out this grey direction.
    """
    difference=np.asarray(mu_reference)-np.asarray(mu_trial)
    weights=np.asarray(reference_weights,dtype=float)
    if weights.shape!=difference.shape or np.any(weights<0) or weights.sum()<=0:raise ValueError('Invalid reference weights')
    mean=float(np.average(difference,weights=weights))
    return difference-mean,mean
