# Retrospective flux-residual transfer protocol

Frozen before computing the new scores on 26 September 2026. The original
12-object independent-flux aggregate held-out results have already been read;
this is not pristine preregistration. The selected pilot, full-data covariance
and published accepted mask also have known selection/nuisance leakage.

Use every original `epoch` SALT held-out fit (12 objects), with its unchanged
hashed-night train/test rows. Reconstruct its full conditional predictive
covariance and compare archived means/variances. Refuse failed/boundary fits.

Compute first-order training-refit-adjusted flux responses to two three-column
families. Observer family: g-r, i-r, z-r magnitude perturbations. Rest/phase
family: Gaussian spectral features centered at 4300 and 6000 Angstrom, each
sigma 600 Angstrom, and SALT colour-law times tanh(restphase/20 days). Each is
propagated through actual SED and passbands; spectral features are empirical
diagnostics, not identified physical laws. Coefficients have common zero-mean
Gaussian ridge scale 0.02 mag, with 0.01 and 0.05 mag sensitivity runs. These
are regularization choices, not measured calibration or population priors.

Whiten each held-out residual by its full predictive covariance. Leave out
one whole object when fitting the common feature coefficients. Score its
held-out residual using Gaussian predictive integration of the coefficients,
including their uncertainty. Compare each family's summed leave-object-out
log predictive density to the baseline standard normal. Preserve each
object's contribution, rank and response norm. Predictions share training
objects, so folds are not independent evidence.

Use 10000 standard-Gaussian residual simulations under the same fixed design,
same folds, covariance, templates, and all six family/scale candidates. Report
the upper-tail probability of the maximum score gain, which accounts for
the declared family and regularization search. A negative gain is retained.
The calibration is exact only for this fixed local predictive Gaussian
residual design; it does not regenerate clipping, selection, covariance
estimation, empirical SALT training or nonlinear parameter fits. Add explicit
one-template injections to report conditional detection power.

No distance, host-age or cosmology fit enters. Any preferred template remains
an observed residual pattern conditional on this pilot. Full covariance,
noise and redshift/field selection calibration would be needed for physical
inference. A null with poor injection power is nonidentification, not a dust
or population bound.
