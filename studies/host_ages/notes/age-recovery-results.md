# What standardization and cosmology absorb in controlled age tests

These experiments isolate specific responses of the measurement. They do not regenerate a complete survey, train a new SALT surface, or measure a physical age correction. The [cosmology design](../code/age_recovery/cosmology-protocol.json) was frozen before its new outcomes. The [flux design](../code/age_recovery/protocol.json) records its initial choices and an explicitly subsequent physical-dust extension.

## Actual cosmological fitting on the observed age sample

The injection uses the 196 unique G11-first SDSS matches, their actual redshifts, width, colour, host masses and full Pantheon+ STAT+SYS covariance. Injected ages are the observed age medians **treated as known coordinates**. This tests recovery conditional on those coordinates; it does not assert that they are the true ages. Every amplitude receives the same 500 multivariate Gaussian noise realizations. The fit profiles linear nuisance coefficients while optimizing the actual nonlinear flat-ΛCDM luminosity distance, with a free intercept and 0.001 < Ωₘ < 0.999.

| Fitted standardization before testing its residuals | Retained fraction of a small injected age slope | Mean residual slope for injected −0.030 | Joint age-fit recovery | Null rejection fraction |
|---|---:|---:|---:|---:|
| Intercept and cosmology | 98.77% | −0.02969 | −0.03005 | 4.4% |
| Also width, colour and host-mass step | 64.03% | −0.01927 | −0.03010 | 5.6% |

Slopes are mag/Gyr. Null rejection uses nominal 95% local Gaussian intervals. These cover the known injected age slope in 95.6% and 94.4% of draws, respectively; Monte Carlo errors on these proportions are about one percentage point. Detection power at −0.010 is 60.8% and 42.6%; at −0.030 it is 100% in this known-age experiment. This is **optimistic design sensitivity**, not power after uncertain host-age inference. Ωₘ boundary rates are retained in the result record; the standard fit reaches its boundary in up to 3.2% of the tested draws and the joint fit in up to 0.4%.

The directly fitted cosmology is not equivalent to forcing the residual mean to zero in arbitrary redshift bins. On this measured design, the separate covariance-weighted bin-removal diagnostic retains 96.36% of the slope. The older approximately 4.1% cosmology-projection loss used plotting-error weights: the full covariance instead gives 1.23%. These statements concern distinct weighting rules, not a changed physical age distribution.

A small loss of residual slope can coexist with a substantial cosmological shift. Under the declared Ωₘ = 0.3 truth and −0.030 injection, the mean standard fit gives Ωₘ ≈ 0.214, or ≈ 0.209 when it also fits width, colour and host mass. Jointly including the injected age coordinate restores the correct cosmology statistically. Conversely, appending a full −0.030 template to a **zero-age-effect** population gives mean Ωₘ ≈ 0.398. An imposed template can therefore generate a spurious shift. This conditional counterexample establishes neither which population is real nor the size of its residual correction.

The simulation does not refit host ages, regenerate BBC corrections, change event eligibility or impose a physical dust population. It deliberately separates the consequence of cosmological/nuisance fitting from those unresolved mechanisms. [Compact numerical record](../results/age_recovery/cosmology-recovery.json).

## Perturbations before light-curve fitting and synthetic detection

We generated calibrated fluxes on all twelve existing DES objective cadences. Each intervention receives the same twelve noise draws at native brightness and at one-quarter brightness, giving **1,728 event–intervention–noise combinations**. A declared synthetic detection cut is recalculated for each combination: two measurements above 5σ in each of two bands. Every exclusion, fit and failure is retained. All selected fits converge, no fitted parameter reaches a boundary, and two optimizer starts differ by less than 5.3 × 10⁻⁹ in χ².

The independent SciPy photon integral and sncosmo inference share the same trained SALT3 surface. Their agreement is a numerical check, not independent training physics. Independent review caught an earlier shared-quadrature path; the final campaign was rerun using native `Model.bandflux` for inference. An additional 24-vector integration check differs by at most 0.00645 of a per-epoch flux uncertainty. Separate adaptive cosmology integration reproduces eight fits within 9.1 × 10⁻⁸ in Ωₘ. [Independent checks](../results/age_recovery/independent-checks.json). The F99 and phase-dependent interventions change the generated spectra relative to the inference model. Covariance and the prior on peak time remain frozen. These cadences contain previously accepted epochs; the experiment cannot reproduce rejected epochs, DES discovery, classification or follow-up selection.

At native brightness, on the qualified paired intersection:

| Intervention | Mean change in fitted standardized magnitude | Selection compared with 140 baseline detections |
|---|---:|---:|
| Grey dimming +0.090 mag | +0.09118 mag | 139 detections |
| Grey brightening −0.090 mag | −0.09094 mag | 142 detections |
| Phase-dependent colour amplitude +0.030 | −0.01301 mag | 140 detections |
| Full F99 screen, E(B−V) = 0.030, Rᵥ = 3.1 | +0.03919 mag | 138 detections |

The standardized magnitude here uses fixed α = 0.148 and β = 3.112 after refitting the four light-curve parameters. The screen gives fitted dimming of 0.12183 mag and colour change 0.02651; the fixed colour correction removes much, but not all, of the dimming. These are response measurements under specified coefficients, **not evidence that the survey's fitted dust treatment leaves a 0.039-mag bias**. It would need its population refit and selection correction regenerated.

A separate F99 colour basis subtracts its B-band extinction by construction. Its approximately −0.087-mag standardized response must not be interpreted as a physical dust screen or dust overcorrection. Independent review prompted the explicitly labelled full-screen extension above, restoring the grey extinction term. At quarter brightness only 23/144 baseline event–noise pairs pass the cut; the physical screen leaves 18, showing why keeping the original sample fixed would miss part of the response. The paired intersections are smaller and are reported separately from selection changes. Descriptive redshift-bin responses and their off-diagonal Monte Carlo covariance are retained in a [separate record](../results/age_recovery/flux-redshift-response.json); these fixed-cadence responses do not estimate the full-survey bias function. [Flux record](../results/age_recovery/flux-recovery.json).

## Required precision and the full-survey gate

Before the new recovery outcomes, the target was set to keep a declared smooth distance-bias shape below **0.1 statistical standard deviations** in the cosmological quantities of interest. For `g(z) = z/(1+z)`, normalized to unit contrast from z = 0.05 to 1, the local flat-ΛCDM SN likelihood allows a contrast of only **1.94 mmag**. For SN+BAO CPL the corresponding limits are **2.77 mmag for q₀, 2.61 mmag for w₀ and 6.26 mmag for wₐ**. These are local Fisher precision budgets at the stated fiducial parameters; they are not shape-independent bounds or empirically measured age corrections.

A fresh public check finds that the DES-SN5YR and Dovekie repository heads are unchanged from the archived versions. All **250** previously unresolved Dovekie simulation objects still return HTTP 404 from their pinned public media endpoints. This is a specific access failure, not proof that no other copy exists. The public scaffold also lacks verified execution provenance for the full revised analysis. Historical classifier reproduction fails its published-equivalence gate. These gaps prevent claiming complete selection, retraining and BBC recovery. [Current endpoint record](../results/age_recovery/survey-availability.json).

## Reproduction

From the repository root, after restoring the frozen main inputs described in [the data guide](../../../docs/workflows.md):

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/age_recovery/cosmology_recovery.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/age_recovery/flux_recovery.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/age_recovery/flux_response.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/age_recovery/independent_checks.py
.venv/bin/python studies/host_ages/code/age_recovery/survey_availability.py --archive /path/to/original-archive
```

Large draw tables and event ledgers are regenerated under `.work/age-recovery/`; code, compact findings, input identities and prospective designs are versioned. Public endpoint status naturally depends on the time of a rerun.
