# Native-verified cosmological fit candidates

For the declared contemporary CMB + DESI DR2 BAO + Dovekie distance likelihood, the best tested CPL point has a penalized objective **14.318 lower** than the best tested flat ΛCDM point. Both points have been evaluated using native CAMB spectra. These are bounded-search candidates, **not certified maximum-likelihood or profile optima**, and the difference is not converted into a significance, posterior probability or Bayes factor.

| Native calculation | ΛCDM candidate | CPL candidate |
|---|---:|---:|
| Penalized objective, S | 512.505041 | 498.186692 |
| −2 × sum of component log likelihoods | 511.091333 | 497.291479 |
| Gaussian calibration penalty | 1.413708 | 0.895213 |
| Ωm | 0.303871 | 0.312561 |
| H0, km s⁻¹ Mpc⁻¹ | 68.144920 | 67.453176 |
| w0 | −1, fixed | −0.811825 |
| wa | 0, fixed | −0.688852 |
| q(0) | −0.544037 | −0.336984 |

These parameter values describe individual fitted points; the table contains no uncertainty estimates. The candidate difference separates into **13.799854 from the likelihood** and **0.518495 from the calibration penalty**. Re-evaluating the ΛCDM candidate in the CPL implementation with w0 = −1 and wa = 0 gives exactly the same native objective, checking the model embedding and nuisance treatment.

## Objective and scientific scope

The calculation minimizes

\[
S(\theta)=-2\sum_k\log L_k(\theta)
 +\sum_j[(a_j-1)/\sigma_j]^2.
\]

The four calibration parameters are A_planck, P_act, Tcal and Ecal, with Gaussian standard deviations 0.0025, 0.003, 0.0036 and 0.0095. Each penalty is applied once. Uniform-prior support is enforced, but prior-density normalization constants are excluded from S. Across all valid candidate-search evaluations, S + 2 log posterior is constant within each model to less than 4.6 × 10⁻¹³. Thus the objective is a precisely defined **penalized likelihood**, not an unpenalized likelihood maximum or a comparison of posterior density heights across dimensions.

The component likelihoods retain their native normalization conventions. In particular, the supernova likelihood includes the fixed covariance normalization and an analytic integral over an unconstrained magnitude offset. Consequently, the absolute S values are not a global goodness-of-fit χ² and must not be divided by a nominal number of degrees of freedom. Differences use exactly the same data and component conventions.

The target is the frozen `modern_fast.configuration` with `calibration='official_planck'`, Dovekie total covariance, and no supernova luminosity evolution. The sound horizon and background geometry come from the same CAMB cosmology as the spectra. The ΛCDM configuration already supplies w0 = −1 and wa = 0 to both the surrogate and native theory; no caller-side interface repair was needed. Flat geometry, the stated neutrino prescription, CAMB's CPL domain and all likelihood settings remain those of the [shared inference](joint-inference.md). This is conditional on the released supernova corrections and the [cross-probe dependence assumptions](probe-dependence.md); it does not resolve an unknown age-dependent luminosity correction or recover a complete raw-data survey likelihood.

## Search accuracy and limits

Each model has two starts in coordinates scaled by the released author-chain covariance. The author means and covariance supply only numerical initialization and scaling, not a data factor or scientific prior. The second start is a fixed seeded displacement. A computational box of ±4 in Cholesky coordinates bounds the search; all reported candidates are well inside it, with maximum absolute coordinates 1.257 for ΛCDM and 0.713 for CPL.

All four L-BFGS-B runs exhausted their **700-call budget** without certifying convergence. There were no invalid prior/theory evaluations and no exact-spectrum fallbacks. The two native candidate objectives differ by 0.057635 for ΛCDM and 0.056275 for CPL. This agreement is a useful local check, not a bound on further possible improvement or on unvisited modes. Because either model could improve further, the candidate difference is neither an upper nor a lower bound on the true profile improvement. Native gradients and stationarity were not measured.

The cubic spectral approximation is used only to propose candidate positions. At the two ΛCDM candidates, native S exceeds the approximation by 0.672033 and 0.810550; native evaluation **reverses their ranking**. At the two CPL candidates, the corresponding differences are 0.006755 and 0.004951. The reported comparison therefore selects each best candidate using native values. It does not substitute surrogate maxima for exact maxima or infer posterior reliability from these seven checks.

Seven native CAMB points were evaluated: the author mean and two candidates for each model, plus the ΛCDM embedding in CPL. The proposal used 2,800 optimizer calls in total; theory caching reduced these to 908 ΛCDM and 960 CPL surrogate-spectrum computations. Both native theory and the separately validated exact lensing matrix reassociation remain unchanged. The latter is an algebraic speed optimization, distinct from the spectral approximation.

## Reproduction and records

Prepare the pinned survey inputs, modern external-probe environment, author proposal files and `cubic-0509.npz` using the [inference instructions](../code/inference/README.md), then run from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/modern_profile.py
```

The [machine-readable result](../results/inference/modern-profile.json) preserves all starts, budgets, candidate parameter vectors, native component log likelihoods, calibration penalties, source and input identities, installed versions, and the exact model-embedding check. The nine active theory/likelihood source files are hashed before and after execution and were unchanged. Detailed search evaluations and the execution design are reproducible ignored artifacts under `.work/unified-cosmology/inference/modern-profile/`. Regenerated timestamps and runtimes need not reproduce their byte hashes; the numerical comparisons and target identities are the relevant checks.
