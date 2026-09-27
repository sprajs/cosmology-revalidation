# Bounded native stationarity test

The native refinement finds slightly better tested points, but **does not establish a stationary local fit for either cosmological model**. The best visited penalized objectives are 512.401055 for ΛCDM and 498.151420 for CPL, a difference of 14.249635. This remains a comparison of finite tested points, not a certified profile improvement, significance or evidence ratio. The [original candidate calculation](modern-profile.md), including its earlier difference of 14.318349 and exhausted search budgets, is preserved unchanged.

| Native diagnostic | ΛCDM | CPL |
|---|---:|---:|
| Best objective among all tested points | 512.401055 | 498.151420 |
| Final point used for gradient checks | 512.442534 | 498.185816 |
| Largest absolute scaled cosmology gradient | 0.510482 | 0.283096 |
| Gradient tolerance | 0.100000 | 0.100000 |
| Largest-gradient direction, step 0.08 | −0.510482 | +0.283096 |
| Same direction, step 0.04 | −0.474130 | +0.058875 |
| Step-halving discrepancy | 0.036351 | 0.224221 |
| Step-halving tolerance | 0.050000 | 0.050000 |
| Failed fixed-spectrum nuisance-gradient checks | 7/35 | 4/45 |
| Native spectrum calculations | 35 | 45 |

The coordinate directions are columns of the author-covariance Cholesky factor, not unit changes in individual physical parameters. In both models, the largest-gradient column is labelled by τ. ΛCDM's derivative remains clearly nonzero at half the step. CPL's derivative changes substantially with step size; that check cannot distinguish unresolved optimization and finite-difference sensitivity well enough to certify stationarity. The tested nuisance gradients also fail the declared criterion at some points. None of these numerical failures is evidence for or against a cosmological model.

The best tested CPL point has w0 = −0.810295 and wa = −0.695276. These are coordinates of one finite evaluation, with no associated uncertainty estimate. The best visited points differ from the final points where complete native gradient checks were performed. Their objective difference is neither an upper nor a lower bound on the improvement at the unknown true profile optima.

## Calculation and unchanged target

The objective and data are exactly those of the [original penalized comparison](modern-profile.md): the contemporary CMB configuration, DESI DR2 BAO and Dovekie total distance covariance, the four calibration penalties applied once, and no supernova luminosity evolution. No active theory, likelihood, covariance, prior or numerical-library source was edited. The original source hashes were checked before and after the calculation.

Each model's surrogate search starts at its best previously native-evaluated point. L-BFGS-B reports relative-function convergence after 888 ΛCDM and 336 CPL calls, below the new 8,000-call limit. Those optimizer messages do not certify a small gradient. The native checks in the table take precedence.

At every native cosmological point, five nuisance parameters—A_planck, P_act, Tcal, Ecal and A_fg—are optimized using cached spectra. A local linear Taylor approximation to the **profiled native-minus-surrogate objective** guides two further searches within boxes of radius 0.35 and 0.25 in the scaled cosmology coordinates. The initial correction uses forward differences at step 0.08; both trial points receive complete central native-gradient checks at that step, followed by a half-step check in the final largest-gradient direction. Nuisance gradients are measured independently of the optimizer status. The stated gates require maximum absolute cosmology gradient ≤ 0.10, nuisance gradient ≤ 0.05, a step-halving difference ≤ 0.05, and an interior final trial.

These are bounded local search boxes, not an exact trust-region optimizer with a guaranteed decrease at every native trial. The first CPL trial increases the native objective; it is retained, not substituted for a best-fit result. Both final points are interior to their search boxes. Neither passes all numerical gates. The calculation stops at the declared **80 native spectrum evaluations with one worker**; no additional native evaluations were used to rescue a failed check. Calibration-only likelihood evaluations reuse spectra, and the surrogate still computes exact background distances; those cheaper operations are recorded separately.

This result establishes the limits of this bounded refinement. Better optimization or a demonstrated treatment of numerical derivative sensitivity would be needed before claiming a locally stationary profile. Posterior inference has separate convergence and exact-weight validation requirements, and cannot be qualified or disqualified solely by this optimizer's outcome.

## Reproduction and retained evidence

After preparing the same frozen environment, likelihood inputs, author proposal covariance and cubic surrogate as the original comparison, run from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/modern_profile_refine.py
```

The [result record](../results/inference/modern-profile-refinement.json) contains both models' source/configuration identities, every native trial and gradient, optimizer statuses, failed gates, nuisance fits, actual native call counts and the original-result hash. Complete likelihood evaluations, including unsuccessful nuisance iterations, are retained in ignored `.work/unified-cosmology/inference/modern-profile-refine/`. Direct reconstruction verifies the recorded calibration penalties, native objective values and central finite differences. The code and evidence preserve the unsuccessful stationarity result rather than reporting optimizer success as scientific convergence.
