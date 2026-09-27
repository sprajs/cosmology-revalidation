# Expansion and residual luminosity at native accuracy two

This consumer is prepared and tested, but no observational accuracy-two history has been calculated. It accepts only a freshly qualified, distinct accuracy-two posterior. It keeps every original selected slot and its untrimmed native/proposal weight. An accuracy-one report cannot be substituted as its input.

One CAMB background per slot supplies both the expansion history and the supernova distance prediction. The background uses the actual finalized native multipole settings and the same three accuracy boosts of 2. It applies the previously validated thermal initialization adjustment and compares the rebuilt expansion rates, deceleration and jerk, matter density and sound horizon against that slot's native accuracy-two record. The original derivative, Friedmann and cosmic-age quadrature checks retain their thresholds. No new CMB spectra are calculated.

Cobaya increases the requested `lmax` when necessary and consumes `halofit_version` separately during parameter initialization. The consumer therefore uses the recorded finalized `lmax`, restores the explicitly declared nonlinear prescription if it was removed from the recorded argument dictionary, and checks the actual native `Params.max_l`. It rejects other changes to declared settings. This preserves the parent configuration; the nonlinear prescription itself is not evaluated by this background-only calculation.

The luminosity convention is positive magnitude drift for dimmer standardized supernovae at a fixed distance. The full released supernova covariance, heliocentric redshift factor, common-magnitude integral and likelihood normalization remain unchanged. The background must reproduce the native accuracy-two supernova log likelihood within $10^{-5}$ and its quadratic residual within $2\times10^{-5}$.

For a smooth luminosity sensitivity, the four spline coefficients were integrated out of the cosmological likelihood. Their conditional distribution is recovered as

$$V=(\sigma^{-2}I+F^TF)^{-1},\qquad E[c\mid\theta,d]=-VF^Tr,$$

where $F$ and $r$ use the same full-covariance whitening and projection of the unconstrained magnitude offset. The marginal history includes both $BV B^T$ and the covariance of the conditional means between cosmological points. Its covariance with $q_0$ uses only the covariance of those conditional means with the native $q_0$: conditional coefficient noise cannot correlate with a quantity fixed at the same cosmology.

The baseline imposes zero luminosity drift. The linear sensitivity uses the sampled coefficient directly. The smooth sensitivities retain the original 0.1 or 0.3 mag coefficient priors and their conditional uncertainty. All impose $B(0)=0$; none identifies stellar age as the physical cause of luminosity evolution.

Reported intervals would be pointwise credible bands conditional on the parent model and priors, not simultaneous bands. Acceleration fractions refer only to the fixed past-redshift grid. The cosmic age is inferred within the assumed CMB cosmology, rather than measured by an independent host-galaxy clock. A failed numerical or weighted-stability check withholds both admitted histories.

The synthetic validation uses analytic matter–radiation–Lambda backgrounds, independent full-intercept GLS calculations, a deterministic conditional-noise construction for total covariance and covariance with $q_0$, and a complete 2,000-slot cache fixture. The fresh observational-parent boundary is explicitly mocked in that fixture. Native CAMB and likelihood construction calls are forbidden. Changed cache bytes, wrong numerical accuracy, altered weights, alternate supernova likelihoods, failed native closure and retries of an existing cache are rejected.

After a full accuracy-two posterior qualifies, the reproduction order is:

```bash
# PY is the pinned modern environment; set all thread counts to one.
$PY studies/unified_cosmology/code/inference/native_accuracy_history_validate.py
$PY studies/unified_cosmology/code/inference/native_accuracy_history.py prepare \
  --work NATIVE2_WORK --summary NATIVE2_REPORT.json --cache NEW_HISTORY_CACHE
$PY studies/unified_cosmology/code/inference/native_accuracy_history.py execute \
  --cache NEW_HISTORY_CACHE --output NEW_HISTORY_REPORT.json
$PY studies/unified_cosmology/code/inference/native_accuracy_history.py report \
  --cache NEW_HISTORY_CACHE --output NEW_HISTORY_REPORT.json
```

Use `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1`. The final command verifies a completed cache without recalculating backgrounds. Claims, failures, record seals and the final hash manifest are preserved; an interrupted attempt is not retried in its existing directory.
