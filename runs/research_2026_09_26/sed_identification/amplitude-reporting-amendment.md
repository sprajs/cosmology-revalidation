# Reporting amendment, 2026-09-26 01:28 UTC

The initial geometry ran before this amendment and is retained with its source
and manifest. No observed-flux or residual arrays were read. No basis, phase
form, drift form, cohort, coefficient target, covariance or derivative choice
is changed. The numerical-rank metadata incorrectly reported the number of
matrix rows rather than retained columns; fix that field without changing the
SVD or calculation.

The already specified .01/.02/.05/.10-mag amplitude scales are additionally
reported as hard full-grid SED RMS limits, solving the monotonic Lagrange
multiplier equation. This makes the amplitude meaning unambiguous, whereas
the original penalty scales are not hard upper bounds. Report the actual
SED RMS on observed support using equal accepted-epoch weights and normalized
positive photon-throughput weights lambda*T at each redshift. Also retain
full-grid RMS and maximum; neither metric is an empirical population prior.

Perform a finite exp(-kappa*Delta_m) broadband reintegration for the broad-drift
secondary frozen-vector matches at all four declared amplitude limits, using
both discovery-design and validation-design geometry coefficients. Hold the
SALT coordinates, covariance and nuisance projection fixed. This tests SED
linearization only and is not a nonlinear SALT refit or mechanism fit. Report
negative native spectral support separately and compare native unscaled
spectral derivatives with the official-mean fractional mapping.
