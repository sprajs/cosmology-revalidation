# Evolving dark energy in the shared cosmology

The shared Dovekie supernova, DESI DR2 BAO and Planck–ACT–SPT CMB observations give the following conditional CPL inference. The dark-energy equation of state is $w(z)=w_0+w_a z/(1+z)$, and standardized supernova brightness is held constant apart from the corrections already in the released distances.

| Quantity | Posterior mean ± standard deviation | Equal-tail 95% interval |
|---|---:|---:|
| H₀ [km/s/Mpc] | 67.460 ± 0.565 | [66.395, 68.548] |
| Matter fraction Ωm | 0.31272 ± 0.00549 | [0.30237, 0.32323] |
| w₀ | −0.8059 ± 0.0555 | [−0.9126, −0.6983] |
| wₐ | −0.7222 ± 0.2109 | [−1.1481, −0.3222] |
| Present deceleration parameter q₀ | −0.3310 ± 0.0614 | [−0.4490, −0.2121] |
| Present jerk j₀ | −0.2173 ± 0.3119 | [−0.8007, +0.4075] |
| Sound horizon rᵈ [Mpc] | 147.217 ± 0.214 | [146.815, 147.640] |

These values pass the sampling and native importance-weight checks **at the declared numerical accuracy**. The completed higher-accuracy screen detects numerical sensitivity across 32 fixed posterior points: the centered total log-likelihood difference has RMS **0.18408** and maximum absolute value **0.39972**, above the declared tolerances of 0.05 and 0.10. Correcting an auxiliary background-comparison mismatch does not remove that likelihood variation. The saved-record review verifies the arithmetic and provenance without performing new physical calculations. These posterior errors are therefore provisional, and the shift in cosmological parameters at higher accuracy remains unmeasured. [Measurement and complete input bindings](../results/inference/modern-cpl-independence-measurement.json); [reviewed precision evidence](../results/inference/native-precision-reviewed-modern-cpl-independence.json).

A separate source audit identifies an upstream CAMB spline-indexing defect and prepares matched control and repaired builds. No physical comparison of those builds was completed; the defect has not been shown to cause this accuracy variation or any cosmological shift. The [source and build note](camb-spline-repair.md) and [measurement gaps](../../../docs/measurement-gaps.md) preserve that distinction. All values above retain the original likelihood target and weights.

The inferred expansion is accelerating today: the entire quoted q₀ interval is negative. Its acceleration need not be increasing with time. The jerk interval permits both signs, and about **75.6%** of the sampled posterior has both q₀ < 0 and j₀ < 0. The latter means positive scale-factor acceleration that is decreasing with time; it does not mean a shrinking universe or presently decelerating expansion. The numerical Monte Carlo uncertainty in this fraction is approximately 0.011 across the declared block partitions. It is a model-conditioned posterior fraction, not a frequentist significance or a model probability.

The dark-energy parameters closely agree with the published Dovekie combination, which reports $w_0=-0.803\pm0.054$, $w_a=-0.72\pm0.21$ and $H_0=67.47\pm0.55$. The paper uses posterior medians with integrated 68.27% intervals; the first column above uses means and standard deviations. These are separate numerical analyses of substantially shared observations, not independent experimental confirmation. [Dovekie, version 3, Table 10](https://arxiv.org/html/2511.07517v3#S10.T10).

The narrower ΛCDM inference and this CPL inference answer different conditional questions. The fact that a CPL marginal interval excludes one ΛCDM coordinate does not by itself measure the odds between the models. No model evidence, calibrated likelihood-ratio significance, or empirical age-dependent luminosity correction is inferred here. In particular, ΛCDM constrains jerk close to one by its form, whereas CPL allows its sign to be measured conditionally.

Both inferences assume spatially flat general relativity, the stated fixed neutrino sector, standard recombination, a power-law primordial spectrum, released supernova corrections/covariance and the declared CMB likelihood factorization. The additional CPL implementation domain is $w_0+w_a\leq0$. No selected posterior point is within 0.05 of that boundary, but this does not test the excluded region. The separate brightness-drift alternatives are required before interpreting the inferred dark-energy evolution as robust to supernova luminosity assumptions.

Four independent proposal chains supply 2,000 explicitly retained selected slots. Native evaluation gives raw importance ESS **1965.3**, Pareto-tail diagnostic **k = 0.2773**, and all declared chain and batch checks pass. Repeated stratified-selection slots are retained with their original weights. These finite diagnostics do not guarantee discovery of remote unvisited modes. [Native correction result](../results/inference/modern-cpl-independence-exact-correction.json).
