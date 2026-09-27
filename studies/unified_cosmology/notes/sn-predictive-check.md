# Conditional marginal-Gaussian quadratic check

This check asks whether the aggregate supernova residuals are unusually narrow
or broad relative to the released Gaussian distance model at the cosmological
parameters supported by the joint fit. No observational result is reported here:
the calculation accepts only a native-corrected posterior that passes all
registered qualification gates.

For N supernova distances, let Σ be the full released covariance, augmented by
σ²BBᵀ when a smooth Gaussian luminosity prior is present. At fixed cosmological
and sampled nuisance parameters, removing the one free global magnitude leaves
a replicated quadratic with distribution χ²(N−1). Properly integrated Gaussian
luminosity modes enlarge Σ; they do not remove additional degrees of freedom.
Cosmological parameters are held fixed for each conditional replicate, so their
number is not subtracted either.

Each native record already contains the observed quadratic `sn_chi2`. Before
using it, the check verifies the complete normalized supernova log likelihood:

$$
\log L_{\rm SN}=-\frac12\left[
 Q+\log\det\Sigma+\log(\mathbf1^T\Sigma^{-1}\mathbf1)
 +(N-1)\log(2\pi)\right].
$$

This density accounting uses the same improper flat global-magnitude measure
as the fitted target. It is not a proper normalization over the full data space
or an absolute evidence normalization.

The calculation retains the full covariance and verifies the exact released
data hash, luminosity configuration and first-party likelihood source. It makes
no new background or CMB calculations and introduces no fitted error multiplier.

For each observed quadratic Q, the reported lower fraction is
Pr(Qrep ≤ Q), and the upper fraction is Pr(Qrep ≥ Q), under χ²(N−1). Their
averages use the original untrimmed native/proposal importance weights. A small
lower fraction indicates a quadratic smaller than expected from that conditional
Gaussian model; a small upper fraction indicates a larger quadratic. The report
also retains the observed quadratic distribution, Q/(N−1), separate-chain
averages and chronological batch Monte Carlo uncertainty.

These averages use the same observations that selected the posterior parameters.
They are **not calibrated frequentist p-values**, evidence ratios or a proof of
covariance accuracy. An aggregate statistic can also miss redshift-dependent,
population-dependent or directional errors. The released selection corrections,
classification treatment, distance Gaussian approximation and covariance remain
assumptions. A discrepancy would motivate investigation, not automatic error
rescaling or a physical host-age interpretation.

For the baseline and explicit linear-luminosity model, replication represents
projected distance noise conditional on the sampled parameters, including ε
where it is sampled. For the smooth models it is **mixed replication**: a fresh
Gaussian luminosity-coefficient vector is drawn from its declared prior for
each entire replicated dataset. This is not prediction of the same universe
or the same hosts conditioned on a shared inferred luminosity curve. The smooth
curve's conditional posterior is a separate calculation in the luminosity-history
analysis.

Independent synthetic validation uses 300,000 correlated Gaussian replicates
across two sample sizes and three luminosity-prior widths. It checks the
N−1 mean, variance and nominal conditional tails; independent dense GLS and QR
projections; and invariance to a 53-mag global offset. Direct execution of the
released SN likelihood agrees to 1.25×10⁻¹³ in log density. A separate synthetic
consumer fixture verifies unchanged replays and rejection of altered cache
values, covariance data, luminosity settings and native density accounting.
Unqualified inputs are rejected before data access. These tests validate the
conditional calculation, not the calibration of posterior-averaged tail values
or adequacy of actual observations.

The frozen assumptions are in
[`sn-predictive-design.json`](../code/inference/sn-predictive-design.json), and
the numerical checks are in
[`sn-predictive-validation.json`](../results/inference/sn-predictive-validation.json).

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/sn_predictive_validate.py

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/sn_predictive_check.py \
  --chain-folder QUALIFIED_CHAIN_FOLDER \
  --correction-summary QUALIFIED_CORRECTION_SUMMARY \
  --cache .work/unified-cosmology/inference/sn-predictive/TARGET \
  --output studies/unified_cosmology/results/inference/sn-predictive-TARGET.json
```

The compact result links to hash-bound input lineage and every per-point
quadratic, tail fraction and original weight in the ignored cache. Failed
density checks withhold all posterior tail summaries. No qualified observational
run has been substituted with a synthetic or provisional chain.
