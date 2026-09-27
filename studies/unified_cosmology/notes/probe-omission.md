# How much each probe contributes

This helper defines three conditional comparisons: remove the released
supernova likelihood, remove DESI BAO, or remove both current lensing blocks
while keeping primary CMB, BAO and supernovae. It has been tested synthetically;
**no observational omission result is reported yet**. A production calculation
first requires an unchanged parent that passes every chain and native-CAMB
importance-correction qualification gate.

The exact component names and allowed configurations are fixed in
[the design](../code/inference/probe-omission-design.json). The lensing comparison
removes `act_dr6_lenslike.ACTDR6LensLike` with its joint ACT+Planck v1.2 covariance
and `SPT2023_lensing`, the Pan et al. (2023) reconstruction of 2018 SPT data,
using its `Lens_and_CMB` response convention.
That latter component contains lensing band powers with primary-CMB response
corrections; it is not a second primary-CMB dataset. This helper does not replace
it with the newer SPT-MUSE product or choose arbitrary likelihood subsets.

For a stored point, let the original correction weight be
$\log w=\log p_{\rm native,parent}-\log p_{\rm proposal}$. If $\mathcal R$ denotes
the declared removed components, the new weight is exactly

$$
\log w_{\rm omit}=\log w-\sum_{r\in\mathcal R}\log L_{r,\rm native}.
$$

The calculation uses the stored **native component log likelihoods**, including
their original normalization constants. It checks native/proposal prior
accounting and the equivalent target-minus-proposal identity. No CMB spectra,
background distances or new chains are calculated. Hardware/version queries
may occur when checking a GPU parent's identity; no GPU arithmetic is needed.

Every original cosmological and nuisance prior and sampled dimension is retained.
In an explicit linear-luminosity parent, removing supernovae leaves ε as an
unused coordinate with its original normalized uniform prior. Removing lensing
likewise retains the unused foreground amplitude `A_fg` and its normalized
prior. The helper does not condition these coordinates on fixed values or
remove their probability factors. They remain in all weight-stability diagnostics
but are omitted from published parameter intervals, covariance and mean changes;
a recovered prior is not a measured constraint. Smooth luminosity coefficients were already
integrated into the supernova factor and are not newly sampled. The removed
supernova magnitude integral used an improper flat-M measure; its arbitrary
constant cancels in normalized posterior weights, but does not authorize an
evidence ratio.

Each alternative must independently pass the existing requirements: at least
2,000 finite points, four qualified parent chains, raw weight ESS at least 400,
Pareto k below 0.7, and all per-chain weight, chain-mean and chronological batch
checks. The weights are never clipped or Pareto-smoothed. A failed comparison
retains diagnostic weights and failure reasons but returns no admitted posterior,
acceleration fractions, covariance or mean-shift measurement. A failure means
that new target sampling or another validated bridge is needed; it does not
show that the omitted probe is necessary for acceleration.

Passing comparisons report conditional parameter intervals, the sign fractions
of q and jerk, their covariance, and paired changes in parameter means relative
to the same parent cohort. These are not independent-fit significance tests.
Finite overlap diagnostics cannot certify distant modes the parent never
visited. Original physical support also remains: flat GR, fixed neutrino
physics and the analytic CAMB CPL restriction $w_0+w_a\leq0$ where applicable.
Removing a probe does not investigate the excluded theory domain.

This tests how the included observations affect inference **within the declared
factorized model**. It does not estimate the missing cross-probe covariance,
calibrate a discrepancy p-value, resolve a physical age correction or validate
the upstream supernova reduction. The [dependence review](probe-dependence.md)
states the relevant shared-information assumptions.

Validation covers 144 CPU/fast-CPU/GPU configuration comparisons, three
16,000-point known-Gaussian removal experiments, exact null removal, additive
component-normalization cancellation, pairwise lensing accounting, deliberately
concentrated weights and too-small cohorts. A clearly labelled 2,000-point
synthetic consumer fixture mocks parent qualification and configuration construction
solely to test the cache and lineage interface; production contains no qualification bypass. Changed
parent bytes, covariance-data bytes, cache payloads and target identities are
rejected. Native spectra, backgrounds and likelihood initialization are guarded
against during validation. [Validation record](../results/inference/probe-omission-validation.json).

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/probe_omission_validate.py

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/probe_omission.py \
  --chain-folder QUALIFIED_CHAIN_FOLDER \
  --correction-summary QUALIFIED_CORRECTION_SUMMARY \
  --omit sn \
  --cache .work/unified-cosmology/inference/probe-omission/PARENT-without-sn \
  --output studies/unified_cosmology/results/inference/probe-omission-PARENT-without-sn.json
```

Use `--omit bao` or `--omit lensing` with a distinct cache/output for the other
two declared targets. The compact report binds the qualified parent, target,
source, environment and every cached per-point component calculation. A replay
checks existing bytes instead of silently replacing a different calculation.
