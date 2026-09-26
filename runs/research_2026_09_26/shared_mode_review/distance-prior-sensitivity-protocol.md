# Frozen contrast pushforward across the declared priors

This retrospective arithmetic check is frozen before computing these contrast
sensitivities. Read the completed root12-mode response arrays and their exact
saved high255-minus-low255 validation membership. Rebuild that same contrast
only as an arithmetic check; do not choose new objects, bins, redshifts or
weights, and do not fit light curves or a cosmology model.

Root responses use prior-whitened systematics columns including COLORLAW
amplitude0.3 and observer columns including0.02. For the separately frozen
paper-weight sensitivity, multiply response column11 (zero-based) by
sqrt(1/3)/0.3, exactly as in both discovery/validation designs used to derive
that sensitivity's saved posterior. Leave other systematics responses fixed.
For the already declared isotropic observer prior, replace its three response
columns by response_last3/0.02 @ chol(Cc), where
Cc=2*0.02^2*(B^T B)^(-1) and B maps three contrasts to zero-sum griz.

Report the same contrast under the saved discovery posteriors of null,
inherited-observer and isotropic-observer models, using both original and
paper-weight COLORLAW conventions. These are exactly the previously evaluated
model/prior choices, not new fit families. For a transformed contrast row h,
compute conditional mean h*m and standard deviation sqrt(h*V*h^T). Check
these independently by transforming the full1020-object response matrix and
forming its saved high-minus-low membership contrast.

The output is a finite-mode contribution to a local compensation contrast at
fixed alpha/beta, nominal coordinates/covariance and published selection. It
is not the contrast's total statistical uncertainty, a measured distance bias,
an identified correction, a full-systematic bound, or a cosmology posterior.
Preserve every original response and result. No probability or significance is
assigned to these contrast means.
