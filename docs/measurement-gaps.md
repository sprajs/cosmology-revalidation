# What the evidence establishes, and what remains unconfirmed

**Evidence closeout · 27 September 2026**

The current evidence closely recovers the published cosmological calculations and favors present acceleration within the tested models. It does not establish an additional age-dependent supernova correction, settle whether standard corrections absorb all age effects, or determine whether the positive acceleration is increasing or decreasing. The local distance calibration remains discrepant with the tested joint ΛCDM prediction. None of these findings identifies a new physical cause.

This document records the limits of the existing evidence. Research execution has stopped at the user's request; no new analysis or retry is scheduled. The [main manuscript](../README.md) gives the measured results. Missing inputs, incomplete numerical calculations and unresolved physical identification are distinguished below because they prevent different conclusions.

## 1. How final are the combined cosmological numbers?

The completed CMB + DESI DR2 BAO + Dovekie calculation gives **H₀ = 68.164 ± 0.251 km/s/Mpc** in flat ΛCDM and **67.460 ± 0.565** in flat CPL. The CPL value closely matches the published **67.47 ± 0.55**. Its present deceleration parameter is **q₀ = −0.331 ± 0.061**, and its jerk has a 95% interval **[−0.801, +0.408]**. Negative q means acceleration; the two possible jerk signs mean that positive scale-factor acceleration may be increasing or decreasing. These are independent calculations using shared observations, not independent experimental confirmations.

Sampling and importance-weight checks pass for both completed posterior samples at their original numerical settings. However, increasing the CMB integration accuracy changes the likelihood across 32 fixed points per model. After removing an auxiliary background-comparison error, the centered log-likelihood RMS changes remain **0.215 in ΛCDM** and **0.184 in CPL**, above the declared 0.05 tolerance. Maximum centered changes are **0.411** and **0.400**. This is real numerical sensitivity in the tested implementations.

**Unconfirmed:** how much the parameter means and intervals change after numerical convergence. These diagnostics are not parameter shifts, and they cannot be converted directly into an extra H₀ error bar. No complete posterior at the higher accuracy was obtained. The quoted tight joint uncertainties are provisional; the numerical issue is not evidence that acceleration disappears.

A separate source audit finds the spline-indexing defect fixed by [CAMB commit 56a95f78](https://github.com/cmbant/CAMB/commit/56a95f78fdd72709c5af6b668141ec0c617777c3) in the release used here. Package provenance, the indexing mechanism and matched control/repaired builds are verified. **Its effect on our spectra, likelihood or cosmological parameters is unmeasured.** It has not been demonstrated to explain the accuracy sensitivity. Distinguishing the repair from compiler effects requires the prepared three-way comparison of the original wheel, unpatched local build and repaired local build, followed by an adequate posterior calculation if warranted.

Evidence: [ΛCDM measurement](../studies/unified_cosmology/notes/joint-lcdm-calibration-results.md), [CPL measurement](../studies/unified_cosmology/notes/joint-cpl-results.md), [reviewed ΛCDM precision](../studies/unified_cosmology/results/inference/native-precision-reviewed-modern-lcdm-independence.json), [reviewed CPL precision](../studies/unified_cosmology/results/inference/native-precision-reviewed-modern-cpl-independence.json), [CAMB source and build audit](../studies/unified_cosmology/notes/camb-spline-repair.md).

## 2. Is an age correction additional, or already included?

There is evidence that host observations contain stellar-population information. Across 609 hosts, photometric age and the spectroscopic 4000 Å break have rank correlation **0.659**. That does not make photometric age an independently measured progenitor age.

The inferred brightness relation also changes with correction accounting. For 196 matched objects with ages treated as exact, the full-covariance slope is **−0.00495 ± 0.00456 mag/Gyr** using released corrected residuals, and **−0.01294 ± 0.00456** after reversing the exported bias correction. Allowing Gaussian age uncertainty gives a nominal 95% profile interval **[−0.056, +0.008] mag/Gyr**, including both zero and the proposed −0.030 scale. A separate 165-object spectroscopic brightness sample gives **+0.00198 ± 0.01993 mag/Gyr**, without demonstrated held-out predictive improvement.

A controlled survey model shows why a residual age slope cannot settle the dispute: a substantial residual relation can survive even when much of the mean brightness shift with redshift is absorbed by the standard bias treatment. The simulations do not yet have adequate coverage and observational support to measure that absorption in the real survey.

**Unconfirmed:** the additional redshift-dependent distance correction, if any, after the full existing standardization, dust, classification and selection treatment. Neither “an extra correction is required” nor “the standard correction already removes every age effect” has been established. Reversing an exported correction is an accounting experiment, not reconstruction of the counterfactual survey.

Resolving this needs observationally constrained joint age, dust and supernova population distributions, their measurement uncertainties, and a selection model that predicts held-out observables and recovers injected effects across the relevant population. The required target is the remaining mean distance bias with redshift, not just a residual-versus-age regression coefficient.

Evidence: [age and bias results](../README.md#4-host-ages-bias-corrections-and-dust), [host observations](../studies/unified_cosmology/notes/host-likelihood.md), [survey and simulation limits](../studies/unified_cosmology/notes/survey-selection.md).

## 3. Do the galaxy observations provide independent absolute ages?

The recovered galaxy data constrain spectral features and population differences, but dust, metallicity, star-formation history and flux response remain coupled. In the declared flexible stellar-mixture family, the five-band spectral constraints for **54 of 55** matched DESI hosts admit formed-mass ages both below 1 Gyr and above 10 Gyr. This is a limit of those compressed measurements and that model family, not a proof that full spectra contain no further age information.

The 1,088 distinct OzDES host spectra lack the relative flux calibration needed to interpret their raw count ratios as calibrated ages. Simple repeat-spectrum models also fail for many objects. Some supplied photometric age posteriors already include informative dust/metallicity priors and a cosmological age ceiling; treating them as independent age likelihoods would reuse those assumptions.

**Unconfirmed:** a sufficiently precise, independent progenitor-age distribution across redshift, and its relationship to global or local host age. Calibrated full spectral response, flexible stellar-population models, local host information where available, and explicit prior removal or reconstruction would be needed. Imposing a cosmological clock cannot independently validate that same cosmology.

Evidence: [calibrated host measurements](../studies/unified_cosmology/notes/calibrated-hosts.md), [physical age constraints](../studies/unified_cosmology/notes/calibrated-host-physics.md).

## 4. Would freely evolving supernova brightness change the conclusion?

The complete shared-probe results hold standardized brightness constant apart from the corrections already in the released distances. An externally imposed age template can shift a separate distance-fit central value to **q₀ = +0.063**, but its 95% interval **[−0.081, +0.209]** crosses zero. Its amplitude and physical transfer to the observed population have not been measured.

The new joint linear and smooth brightness-evolution calculations were not completed. Their partial native likelihood records are preserved, but no qualified posterior is reported from them.

**Unconfirmed:** robustness of the full CMB–BAO–supernova result to the specified free brightness histories, and still more to a physically inferred age correction. Completing those statistical alternatives would measure conditional degeneracies; observational identification from sections 2–3 would still be needed to call any fitted drift an age correction. The present gaps do not imply that all theories are equally probable.

Evidence: [luminosity assumptions](../studies/unified_cosmology/notes/luminosity-sensitivity.md), [stopped calculation record](../studies/unified_cosmology/results/inference/research-closeout.json).

## 5. What is unresolved in the local calibration discrepancy?

The released distance-ladder matrix gives **H₀ = 73.043 ± 1.007 km/s/Mpc**, recovering its published baseline. The separate calibrated Pantheon+SH0ES cosmology gives **73.550 ± 1.017**. These overlap other supernova samples and are not independent priors to multiply into them.

With the full released covariance, the ΛCDM CMB–BAO–noncalibrator fit poorly predicts the 77 withheld calibrator rows: the declared common-offset statistic has predictive equal-tail area **2.56 × 10⁻⁸**, with approximately **7% relative Monte Carlo error**. This is a conditional discrepancy. It is not a frequentist significance, model odds or a diagnosis of the cause; the CMB numerical limits also apply.

Two gaps remain. First, four hosts with distinct calibrating supernovae have shared released covariance **1.97–2.00 times** the reconstructed Cepheid-only host variance. Shared Cepheid nuisance marginalization and a common supernova brightness do not explain the difference, but the exact released covariance recipe has not been recovered. An additional host-shared term is possible. We have neither declared double-counting nor halved that covariance.

Second, the fully calibrated CPL replacement has effective reweighted sample size only **3.74/2,000**, so its posterior is withheld. Separately, the CPL replacement using only noncalibrator supernovae passes its posterior checks, but its prediction for the withheld calibration has insufficient tail and predictive-density precision. Those predictive quantities are withheld. A separate proposal for the fully calibrated target was prepared, but no production posterior was completed.

**Unconfirmed:** the physical source of the ΛCDM calibration discrepancy, the exact covariance construction, and the reliable calibrated CPL comparison. The missing covariance recipe and adequate direct sampling address different parts of this gap.

Evidence: [ladder reconstruction](../studies/unified_cosmology/notes/distance-ladder.md), [covariance comparison](../studies/unified_cosmology/notes/calibration-covariance-estimands.md), [joint calibration results](../studies/unified_cosmology/notes/joint-lcdm-calibration-results.md), [CPL sampling closeout](../studies/unified_cosmology/notes/anchored-proposal-refinement.md).

## 6. What prevents a complete measurement from raw observations?

Reproducing released distances and likelihoods does not reconstruct their entire calibration, training and survey selection. Rejected or unreported epochs, historical covariance transformations, selection efficiencies, contaminant populations and shared calibration terms remain incomplete. Reproducing all **17,733** released classifier probabilities resolves an interface discrepancy; it does not validate the selection model.

The detector and light-curve investigations identify specific failures, not a universal correction. One peak-time Hessian understates its profile-supported uncertainty; one pixel dominates an HST dark-image difference. Infrared model comparisons fail the declared sampling-precision checks. None justifies rescaling every supernova error, choosing a new dust law or claiming an observed infrared model preference.

The combined inference also retains physical assumptions: flat general relativity, a fixed neutrino sector, standard recombination, a power-law primordial spectrum, the stated dark-energy domain and declared treatment of dependence between CMB likelihoods. A joint correlated-lensing calculation was prepared and validated on test inputs, but no qualified new posterior establishes its impact. Failed probe-omission reweighting is inadequate sampling, not evidence that the omitted probe is inconsistent.

The July 2026 Lyα replacement yields only modest shifts within its tested approximation, but uses a rounded correlated Gaussian distance summary rather than the full unrounded author likelihood. These are specified limitations of the available measurement, not grounds for silently inventing unavailable inputs.

Evidence: [dependence between probes](../studies/unified_cosmology/notes/probe-dependence.md), [Lyα approximation and results](../studies/unified_cosmology/notes/lya-fullshape-results.md), [measurement investigations](../studies/README.md).

## 7. Exact stopping state

The main unfinished calculation trees stopped on **27 September 2026 at 19:47:18 UTC**. The separate calibrated-CPL pilot was also interrupted during closeout. Completed records, original failures and source hashes remain preserved; no partial cohort has been promoted to a final measurement.

| Calculation | Evidence retained | Scientific status |
|---|---|---|
| Fixed-brightness ΛCDM and CPL | 2,000 native evaluations per model | Sampling/weight checks pass at original accuracy; numerical convergence unresolved |
| Accuracy 1 versus 2 | 32 points per model | Both detect likelihood sensitivity |
| ΛCDM accuracy 3 | 4 completed of 32 planned; 7 claimed | Incomplete; interrupted attempts remain unknown |
| Linear brightness evolution | 1,822 of 2,000 native record files | No qualified complete posterior |
| Smooth brightness evolution | 1,501 of 2,000 native record files | No qualified complete posterior |
| Calibrated-CPL proposal pilot | 43 of 400 candidate records; 24 finite, 19 nonfinite | Incomplete; no production chains or calibrated posterior |
| CAMB spline comparison | Verified sources and matched builds; 0 of 12 physical comparisons | Repair impact unmeasured; execution prototype unvalidated |

The [main stop receipt](../studies/unified_cosmology/results/inference/research-closeout.json), [calibrated-CPL receipt](../studies/unified_cosmology/results/inference/anchored-refinement-closeout.json) and [CAMB status](../studies/unified_cosmology/results/inference/camb-spline-executor-status.json) retain the detailed boundaries. The incomplete calculations are resumable evidence, not promises of continuing work.
