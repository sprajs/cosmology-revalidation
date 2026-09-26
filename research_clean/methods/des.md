# DES calibrated fluxes, predictors and shared calibration

## Independent flux fitting

`des-flux` operates on 12 exported fixed objectives: calibrated observations, redshifts, Milky Way E(B−V), selected epochs, covariance, starting parameters and the actual native peak-time prior centre. These are explicitly derived inputs. The original detector-image reduction, native covariance construction and survey selection are not rerun here.

The retained engine implements the photon-count integral using the frozen DES SALT3 surfaces, colour law, observer-frame F99 foreground extinction, DES throughput and AB primary spectrum. It has two interpolation/integration paths: direct SciPy surfaces and sncosmo. Throughput outside model support, negative model variance and invalid phases cause failures instead of extrapolated scientific claims.

The four fitted parameters are log amplitude, width, colour and peak date. The objective uses the frozen correlated flux covariance and the exported 10-day peak prior. Bound checks, two starts, finite-difference Hessians and a half-step Hessian comparison accompany each fit. A Hessian must be positive definite for the local covariance to be used.

Default output: 24 full fits across 12 objects and two interpolators, plus 36 deterministic held-out-epoch fits across SALT colour, F99 colour and a phase-dependent colour family. Epoch folds use the original CID/date hash; they are not selected using residuals. Conditional held-out means include the measured covariance between training and test epochs. Laplace prediction adds parameter uncertainty through the Jacobian.

```bash
uv run --frozen python research.py run des-flux --name des-flux
```

Expect all 24 full fits and 36 held-out fits to be numerically valid with this bundle. Maximum absolute difference from the native reference is approximately `0.00112 mag` in fitted amplitude and `0.096 day` in peak date. Different empirical colour families sharing the same observations and model assets do not constitute independent physical measurements of dust.

## Selected-sample prediction

`des-predictors` uses locally derived pre-BBC `(mB,x1,c)` summaries and full 3×3 covariance, with an exact frozen 1,063-object high-purity DES cohort. Fold zero contains 213 held-out objects; 850 are used for fitting. Its flexible distance-reference offsets are nuisance functions, not an inferred cosmological history.

Candidate predictors are width, colour, host-mass mixture and specified interactions/nonlinear extensions. The observed covariance projects through `v=(1,-dmean/dx1,-dmean/dc)` and contributes `v^T C v` to the magnitude variance. Two host branches are mixed using the fixed host-mass probability. Gaussian and four-degree-of-freedom Student-t residual laws are supported; the Student-t scale is chosen so its variance equals the stated variance.

Default: Tripp width+colour, coasting reference, Gaussian residuals, seed 2026092160, four sequential NUTS chains, 1,200 warmup and 2,000 retained draws per chain. The kernel retains the algebraically equivalent projection form that avoids the observed compiled JAX indexing/negation problem.

Scoring requires all R-hat values ≤1.01, effective sample sizes ≥400, zero divergences, no maximum-depth transitions, mean acceptance ≥0.6 and mobile finite chains. Failed chains are saved and the command stops before held-out scoring. The clean default replay passed, with held-out log-score sum about `91.54`; stochastic output can vary with numerical platform.

```bash
uv run --frozen --extra predictors python research.py run des-predictors --name tripp
uv run --frozen --extra predictors python research.py run des-predictors --name colour --config configs/des-colour.json
```

Compare different models on identical held-out physical objects; aggregate paired score differences by object. A good selected-sample predictor is not a selection-normalized parent-population model and does not determine cosmology or a physical dust law.

## Shared calibration modes

`calibration` starts from bundled sufficient statistics for 43 discovery and 1,020 validation SNe: projected Gram matrices `F`, projected residual vectors `u`, and an explicitly defined high-minus-low-redshift distance-response contrast `h`. These products are derived from previous flux fits and fixed released systematic modes. Their upstream projection is not reconstructed from raw pixels in this edition.

For unit-normal mode amplitudes, posterior covariance is `(I+F)^-1`, mean is `(I+F)^-1 u`, and response variance is `h^T(I+F)^-1 h`. The code evaluates 12 released calibration modes alone, plus two specified observer-filter prior conventions. It retains the joint determinant and quadratic terms when computing the validation evidence increment; residual constants common to those comparisons cancel.

Expect conditional response SDs of about `0.00942`, `0.01207` and `0.01236 mag` respectively after all design information is used. The monotonic decrease from prior to discovery to combined variance is checked. Unknown modes, arbitrary grey evolution, survey selection and population changes are outside this finite Gaussian model. These conditional means/variances cannot be promoted to a measured calibration correction.
