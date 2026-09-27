# Sensitivity to unmeasured supernova luminosity evolution

The current provisional baseline sample does **not** provide enough importance-weight overlap to establish alternative cosmological constraints under the three declared luminosity-evolution priors. This is a sampling limitation, not evidence that luminosity evolution is present or absent.

The diagnostic uses 2,000 iterations selected across four frozen baseline-chain snapshots, preserving Metropolis multiplicities. It calculates exact CAMB background distances, then integrates the additional luminosity parameters analytically. The provisional calculation has no exact CMB correction weights and is not a qualified posterior analysis.

| Declared luminosity prior | Raw weight ESS / 2,000 | Pareto k | Largest normalized weight |
|---|---:|---:|---:|
| ε log(1+z)/log(2), ε uniform [−0.5, 0.5] mag | 254.5 | 0.674 | 4.57% |
| Four smooth coefficients, Gaussian σ=0.1 mag | 124.7 | 0.698 | 6.57% |
| Four smooth coefficients, Gaussian σ=0.3 mag | 159.8 | 0.707 | 5.61% |

All three fail the registered weight-ESS requirement of 400. The broader smooth prior also fails k<0.7. Several individual chains have low weight ESS, and weighted chain-mean/batch-error checks fail. Raw weight ESS does not account for chain autocorrelation. These observations support sampling the alternative targets directly rather than reporting intervals from the current reweighting. Final exact CMB corrections could change these diagnostics; finite overlap checks cannot establish coverage of unseen modes.

For computational planning, six of the 2,000 observed baseline points lie outside the spectral interpolation envelope of 4; the maximum absolute whitened coordinate is 5.659. Their alternative-weight fractions are 0.215%, 0.237% and 0.242%, respectively. The observed reweighted centres therefore remain mostly inside the existing domain. These fractions do not constrain unseen alternative tails or the rate at which new proposal steps will require exact spectra.

The smooth basis is the existing natural cubic spline through z=0, 0.1, 0.4, 0.8 and 1.3, fixed to zero at z=0, with four independent Gaussian coefficients. All models retain a separately integrated free constant magnitude. Positive ε adds a positive distance-modulus term and corresponds to dimmer standardized supernovae at z=1. These are explicit priors for **unmeasured grey luminosity changes**, not measured age corrections, dust estimates or an empirical resolution of the host-age dispute.

## Calculation and validation

Let r be geometric predicted distance modulus minus observed distance modulus, C the released full covariance, and P the corresponding precision after projecting out the common magnitude. For the linear mode f=log(1+z)/log(2), define a=fᵀPf and b=fᵀPr. The SN likelihood ratio is the normalized bounded integral

\[
R_{\rm linear}=\int_{-0.5}^{0.5}\exp[-b\epsilon-a\epsilon^2/2]\,d\epsilon.
\]

The prior width is exactly one magnitude. The implementation evaluates this as a stable truncated-normal integral, with conditional mean −b/a. For projected, whitened spline matrix Z and residual q, each Gaussian-prior ratio is

\[
\log R_\sigma=\tfrac12v^T(\sigma^{-2}I+Z^TZ)^{-1}v
 -\tfrac12\log\det(I+\sigma^2Z^TZ),\qquad v=Z^Tq.
\]

The determinant term is retained. For final inference, each ratio multiplies the **original untrimmed exact/proposal importance weight**. The other likelihood factors remain unchanged. Pareto smoothing is used only to diagnose the weights, never to replace them in reported summaries.

Independent checks cover 24 correlated-covariance examples, adaptive numerical integration of the bounded coefficient, direct dense GLS with augmented covariance, coefficient signs, constant-mode cancellation, remote normal tails and large common magnitude offsets. Maximum toy-case log errors are below 5×10⁻¹³. On the actual 1,820-supernova covariance, the two smooth ratios agree with direct augmented-covariance calculations within 1×10⁻¹¹. Exact-background SN likelihoods agree with seven already native-evaluated cosmological points within 6×10⁻¹³. Known Gaussian importance-weight examples also test the overlap gates, including uniform and deliberately concentrated weights. A separately authored 12-case review found no sign, projection or normalization defect.

The provisional chain's printed SN likelihoods agree with the new calculation within 3.8×10⁻⁶, consistent with printed chain precision. Every attempted point is retained; there were no evaluation failures and no native CMB spectrum calculations in this lane.

## Qualification and reproduction

The [design](../code/inference/luminosity-sensitivity-design.json) fixes the priors and gates before the diagnostic. Each alternative needs passed parent-chain convergence and exact-correction checks, at least 2,000 exact points, weight ESS≥400, Pareto k<0.7, and the specified independent-chain and weighted batch-stability checks. Provisional mode can never set `posterior_qualified=true`. Machine-readable weighted means and quantiles are retained for diagnosing support; they must not be promoted to cosmological intervals while qualification fails.

The [results](../results/inference/luminosity-sensitivity.json) and [numerical validation](../results/inference/luminosity-sensitivity-validation.json) preserve source and input identities. Point selections, covariance-sized arrays and per-point calculations remain in ignored `.work` storage. Chain parameter names come from the frozen per-rank manifest and chain headers, not the shared sampler YAML file.

With the recorded modern environment and normalized SN data prepared, run the validation from the repository root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/luminosity_sensitivity_validate.py
```

`luminosity_sensitivity.py --provisional-folder CHAIN_FOLDER --output NEW_WORK_DIRECTORY` captures a new unqualified diagnostic. `--frozen-selection SAVED_SELECTION_JSON` replays exactly the saved point cohort. Final use replaces that input with `--exact-folder EXACT_CORRECTION_DIRECTORY --exact-summary SUMMARY_JSON --parent-diagnostics DIAGNOSTICS_JSON`. Both exact inputs are mandatory and their selection/source identities are checked. `--summary FILE` selects the compact output path; the default is the public sensitivity result. The output work directory must be new and beneath `.work`. None of these commands changes an active cosmological target or chain.
