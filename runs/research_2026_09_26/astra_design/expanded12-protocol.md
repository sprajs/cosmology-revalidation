# Expanded native12-mode extinction/calibration sensitivity

Freeze after observing the nine-mode result and before evaluating these added
variants. Preserve every prior result. This is a labelled post-result sensitivity
to existing released alternatives, not a newly pristine hypothesis test.

Keep the nine coupled calibration/SALT modes with amplitude scales .3. Add the
exact fitopts.yml alternatives CALSPEC (scale1, MAGOBS_SHIFT_ZP_PARAMS
0 .00714 0), MWEBV (scale1, MWEBV_SCALE .95, MWEBV_SHIFT0) and COLORLAW
(scale.3, OPT_MWCOLORLAW89, RV_MWCOLORLAW3.0). Use nominal MODEL000, fixed
nominal SALT parameters, same objects/epochs and native implementation. Evaluate
all three on discovery43 and validation1020. None is selected by score.

CALSPEC changes the observed-data convention. From snana.F90
MAGOBS_SHIFT_USRFUN, delta_m_b=.00714*lambda_mean_b/10000, with lambda_mean
in Angstrom from sntools_calib.c transmission-weighted mean. CALIB_FLUX code
multiplies both flux and quoted error by a_b=10^(-.4*delta_m_b), in single
precision. Compute the scale independently from the original KCOR bandpasses,
freeze it, verify exported data and quoted errors against this known mapping,
and divide its exported model prediction by a_b before forming a response in
nominal units. Record float32 rounding discrepancies explicitly; do not infer
an arbitrary mapping from residuals. Require unchanged SALT coordinates, z,
epochs and bands. Nominal covariance and nuisance projector remain fixed in
nominal flux units. Variant covariance is not added or silently substituted.

MWEBV deliberately changes already-processed E(B-V) to .95 times its nominal
value (no second .86 factor). COLORLAW changes its named extinction law and
RV only; observed flux/error units remain nominal for these two variants.
Gate those intended input changes separately from invariants.

Repeat the discovery-conditioned shared latent calculation with12 calibration/
extinction coordinates: H_sys and H_sys+observer, using the unchanged .02 prior
in original contrast coordinates. Repeat the already declared equal-total-power
isotropic zero-sum griz prior sensitivity. Condition both models jointly on
discovery and integrate the common posterior once across validation, retaining
all cross blocks. Also retain the centred fixed-original-direction diagnostic.
No new field-sign independence assumption, amplitude selection or validation fit.

This is still a finite-variation Gaussian second-moment approximation with local
nominal nuisance projection. Other released population, selection, velocity,
and noise modes, full calibration/dust priors, nonlinear variant parameter
transport, SALT retraining likelihood and BBC are not regenerated. The result
cannot certify a complete systematic budget or identify a physical correction.
