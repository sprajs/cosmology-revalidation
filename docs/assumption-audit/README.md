# Cross-model audit of extinction, dust, ages and correlated uncertainties

21 September 2026. This follow-up examines the current scientific worktree, primary release files, original source implementations and new numerical checks. It preserves the earlier frozen investigation and the separately ongoing SALT/SNM correction work. The detailed evidence is divided into [host ages/Pantheon](pantheon.md), [DES/Dovekie](des.md), [foreground maps and covariance](foreground-and-covariance.md), and [BAO/CMB](bao.md).

For continuation, read the [saved agent handoff](HANDOFF.md). The [final verification](../../runs/assumption_audit/final-verification.json) and [manifest](../../runs/assumption_audit/final-manifest.json) bind the completed evidence, scripts and reports; scientific limits remain explicit below.

**The audit finds a confirmed inconsistency in a historical host-age pipeline, incompatible/overlapping covariance products, and several important restrictions on the quoted uncertainties. It does not establish a replacement cosmology or prove that every later age catalogue shares the historical error.** Some suspected issues fail the numerical tests or have small effects, and those results are retained.

## Findings ranked by the evidence and the work they require

| Finding | Evidence/status | Practical consequence |
|---|---|---|
| **Historical host ages use a different SFH from the fitted spectral model.** | Confirmed in Rose et al. 2019's pinned MC-Age/FSPS source and original archived posterior age columns. Independently reviewed twice. | Audit the age conversion and its descendants before treating those age estimates as a calibrated axis for luminosity evolution. |
| **The original retained age chains include draws outside their stated prior.** | All 105,059,947 archived draws checked: 1,288,021 outside support, affecting 86 of 103 hosts; some have negative SFH timescales. | Recalculating ages cannot turn these into valid posterior samples. The age-conversion comparison explicitly conditions on valid draws and reports exclusions separately. |
| **Host attenuation described in the paper differs from the executed model.** | Young stars receive two dust components in the code, with different power-law indices; the printed R19 formula is different. | The age fit's dust model must be taken from executed inputs. A manuscript formula alone does not reproduce it. No SN luminosity bias is inferred from this discrepancy alone. |
| **Pantheon+ individual dust-systematic matrices have a different statistical baseline from final STATONLY.** | Two direct positive-semidefiniteness counterexamples plus independent reconstructions of the common background. | Subtracting them as if they shared one baseline can misattribute statistical changes to dust. Main local cosmology uses the final total and avoids this trap. |
| **DES single-systematic products are not a disjoint uncertainty budget.** | `CAL_SALT3` already contains `CALSPEC`; subtracting the duplicate closes the total at float32 precision. | Summing every supplied component would double count a calibration contribution. Main local use of the published total avoids it. |
| **DES VPEC is not a valid independent additive covariance as supplied.** | A two-object principal submatrix and 47 negative whitened eigenvalues well beyond float32 noise establish this; the total is positive definite. | Resolve the generation/baseline/precision provenance. An illustrative PSD completion moves Ωm by only 0.000176, so this test does not explain a major cosmology shift. |
| **Shared passbands are assigned independent wavelength perturbations.** | PS1MD/Foundation use common named filter transmissions but receive separate draws in source and actual released calibration inputs. | Test common hardware/filter uncertainty jointly with genuinely survey-specific components. Physical impact remains uncomputed. |
| **A global dust scale does not cover the observed foreground-map alternative.** | Actual CSFD/SFD samples at Pantheon+ positions leave 99.1% of map-difference variation unexplained by a constant plus one global reddening scale. | Propagate spatial foreground alternatives coherently through SN and host photometry. This percentage is not a missing fraction of distance or cosmology variance. |
| **Host-age parameters are coupled; shared age–distance errors remain unpropagated.** | Fixed dust geometry, informative priors, unchanged foreground-corrected photometric errors, diagonal host likelihood, and summary-age regression inputs are verified in code. The age–distance cross-covariance has not been estimated. | Use the joint host posterior and shared SN/host nuisance responses. Correct SN covariance alone cannot supply missing uncertainty on the age axis. |
| **A shared age slope cannot be replaced by diagonal per-SN errors.** | Within the fixed C14 age-template scenario, exact Gaussian marginalization and posterior reweighting give q0 uncertainty 0.0967 versus 0.0743, P(q0<0) 0.330 versus 0.198. | The current shared-slope implementation is correct; the diagonal counterfactual would lose nearly all of that uncertainty's effect. Other age-model uncertainty is still conditional. |
| **BAO-only q0 is materially prior dependent.** | Extending only wa's lower bound from −3 to −10 changes mean q0 from about +0.035 to +0.30 in two independent ensembles. | BAO does not directly measure today's acceleration sign. Its lowest effective redshift is 0.295 and q0 comes from a model extrapolation. |
| **The adopted DES age systematic spans one model variation.** | Released W22 covariance is rank one; the methods specify an SFH/DTD-driven stretch-population alternative. | A label such as “Model SN age” does not establish uncertainty coverage for general residual luminosity-age evolution. Nor does it authorize applying the full historical slope again. |

The reports also record a nonfunctional shipped DES likelihood adapter, a factor-ten README/paper discrepancy in alpha's uncertainty, a C25 percentile/“one sigma” discrepancy, and Gaussian age summaries with appreciable probability outside physically allowed ages. These are separately classified; they are not silently combined into a claimed cosmological bias.

## Why the historical age issue is different from an ordinary model choice

The fit passes a slope `s` to FSPS whose late-time star-formation rate is proportional to `K*(1+s*u)`. The age postprocessor instead integrates `K+s*u`, where `K` is the rate at the transition and `u` is elapsed time. Their early histories agree after a common normalization, but the late histories generally do not. Normalization and time-unit conversion cannot reconcile them.

A parameter set inside the original priors produces an age of 1.41484 Gyr in the postprocessor and 10.68019 Gyr for the SFH used by the spectral fit. That extreme is a **parameter-level demonstration**, not a catalogue-wide bias. The original observational SN5916 posterior provides a more direct check: its archived median age is 6.04072 Gyr, versus 7.49920 Gyr when each of the same fitted SFHs is integrated consistently. Photometry and posterior parameters were not refitted.

The complete archive contains 103 hosts and 105,059,947 draws. Of these, 1,288,021 (1.226%) violate the executed code's prior bounds, affecting 86 hosts. On the explicitly valid subset, the median change across host posterior medians is **+0.41962 Gyr**, the mean is **+0.88175 Gyr**, and **29 of 103 hosts change by more than 1 Gyr** in absolute value. Five cross a 4-Gyr median-age threshold. Removing invalid draws alone changes a host median by at most 0.07256 Gyr, so it does not explain the conversion mismatch. These transformed posteriors have not been validated as new physical age estimates.

The full audit, including its 103-original/102-revised sample crosswalk and rare numerical-integration outliers, is documented in [the age report](pantheon.md), [the per-host results](../../runs/assumption_audit/pantheon/r19-global-validated-age-audit.csv), and [the summary](../../runs/assumption_audit/pantheon/r19-global-validated-age-summary.json). A separate [native-normalization quadrature check](../../runs/assumption_audit/age-independent-crosscheck.json) verifies 24 observed posterior draws without using the analytic moment decomposition for the reference calculation.

An [independent final review](../../runs/assumption_audit/des/r19-prior-peer-review.json) verified aggregate counts, original-code prior boundaries, the short archive member and additional native-FSPS quadrature. Here “valid” means within the stated prior; it does not establish convergence or restore correct posterior weights. Filtering cannot validate the surviving chain.

**The modern propagation boundary matters:** Chung et al. 2025 describes reusing R19's method/settings with a newer FSPS, but its exact revised SFH posterior arrays and executed age-conversion code are not available in this record. The evidence establishes the error in the original released pipeline and outputs. It does not establish the identical error in the revised C25 ages used by the present matched-Pantheon regressions. Those revised ages need a provenance check, not an invented scalar repair. Both sides of the age-versus-standardization argument can depend on the same uncertain age axis.

## What has and has not been correlated

Three different concepts must be separated:

1. **Correlated observational errors:** repeat observations of the same explosion, shared photometric calibration, dust maps, coherent velocities, and joint host-SED errors. Use a covariance or shared nuisance parameters for the same underlying uncertainty; avoid adding that uncertainty twice.
2. **Degenerate model responses:** dust, intrinsic colour, metallicity, age, selection and cosmology can produce similar observed trends. Similar distance-response vectors do not by themselves mean their prior errors have measured correlation.
3. **Reused evidence:** Pantheon+/DES share events, calibration and model training; historical age studies reuse SNe and host photometry. Different reductions, papers or correction hypotheses are not independent confirmations.

For a residual model `y=b*A_true+e` and measured age `A_hat=A_true+u`, measurement error in `y-b*A_hat` has covariance

`C_e + b² C_u − b(C_eu + C_ue)`.

This simple identity shows what is missing when age and SN errors share a dust/calibration driver. It does not estimate the cross term. The actual latent-age analysis additionally needs its population prior, host posterior/SFH information and selection normalization. Inserting an arbitrary correlation coefficient would not resolve the missing measurement response.

The existing shared-slope fit treats its normal prior as independent of the SN likelihood and holds the evolution template fixed. Shared-slope propagation is correct under those conditions. The prior's independence, SFH/DTD choices, cosmological clock, survey-selected age distribution, mean-versus-median convention and correlation with already applied corrections remain distinct assumptions. They cannot all be represented by one ±0.004 coefficient error.

## Tests that did not support a large new bias

The BAO vector's observable order, combined-tracer bookkeeping, covariance sign/scale and full likelihood all pass independent checks. In particular the reversed Ly-alpha DH/DM ordering is handled correctly. Adding the **published-scale** cross-bin theory correlations changes joint q0 by only about 0.0001–0.0004. The original within-bin correlations matter and are already retained. This test does not bound an arbitrary foreground-field error.

DESI uses dust corrections upstream in target selection and clustering weights; it is not intrinsically dust independent. Primary DESI validation nevertheless finds the acoustic peak substantially more stable than broadband clustering under the examined imaging/dust choices. No dust-induced BAO peak error is demonstrated here. The SN+BAO likelihood assumes zero error cross-covariance between probes; the compressed products do not supply enough responses to measure a common foreground contribution.

The CMB-reference chains include foreground nuisance parameters and an ACT–Planck joint lensing covariance. They are not foreground-free direct observations of q0. The local low-redshift radiation approximation changes derived q0 only about 0.0001. These chains were inspected and transformed, not independently reconstructed from maps.

The original main fits correctly preserve full release covariance, profile the SN intercept, keep an uncalibrated BAO scale free, and label extra age-template fits as conditional. The new covariance traps above concern ancillary component files or interpretation; they are not evidence that those main likelihood operations were implemented incorrectly.

## Coverage across the current model families

| Current family/data path | Audit coverage and remaining condition |
|---|---|
| Historical G11/R19 and revised C25 host ages | Original source priors/extinction, posterior conversion, sample overlap and error conventions inspected; original R19 chains reprocessed. Revised C25 executable/posteriors remain a specific missing dependency. |
| Pantheon+ distances and matched-age regressions | Final versus grouped covariance, repeat observations, foreground-map structure, joint-error assumptions and shared-slope propagation checked. No corrected age–SN joint posterior is claimed. |
| Original DES and Dovekie | Matrix formats, nesting/positivity, dust and age variation coverage, host threshold errors, calibration correlations and release adapter tested. Raw SALT/SNM correction remains with the separate task. |
| Cosmic SFH×DTD evolution templates | Current code and earlier independent integrations inspected: fixed clocks, cosmic-volume versus selected-sample weighting, response moments and transported slope are conditional assumptions. Numerical integration agreement is not an astrophysical uncertainty estimate. |
| CPL, LCDM and finite-resolution kinematic fits | Flat geometry/propagation, free intercept, fixed covariance/template and finite function-family restrictions remain. Shared correction uncertainty and CPL prior support were newly tested. q0, bin-averaged q and any past acceleration are different estimands. |
| Directional C1/C2 | Existing independently checked coordinate/frame/covariance record inspected; a fixed age curve, homogeneous population model, cubic expansion, scale boundary and lack of calibrated dipole significance remain. This audit does not redo the separate latent-light-curve correction or claim a new directional posterior. |
| DESI BAO and SN combinations | Fresh official inputs compared byte-for-byte, independent distance likelihood, within/cross-bin covariance, prior support, dust-selection validation and zero cross-probe error covariance examined. |
| DESI CMB-reference chains | Executed foreground/neutrino/geometry/lensing configuration, chain weights and derived-q approximation inspected. No new CMB-map or corrected CMB likelihood run. |

## Concrete next actions supported by this audit

- Trace the revised C25 age-conversion code and full host posterior arrays against the historical SFH mismatch. Use the accompanying parameter and archived-chain reproductions; do not ask only whether “the same method” was used.
- Compare host ages under the executed versus printed attenuation prescriptions, retaining joint metallicity/dust/SFH information and a shared foreground nuisance. Propagate those changes to the same SNe and selection model before estimating a replacement age–luminosity slope.
- Obtain the exact statistical baseline behind Pantheon+ grouping files and the DES VPEC generation ledger/full-precision component. Treat nested calibration groupings explicitly in any reconstructed error budget.
- In the separate lower-level SN pipeline, test spatial foreground alternatives and common PS1/Foundation passband components coherently. Track changes to both the distances and their covariance, and avoid imposing an independently derived correction twice.
- Continue using BAO with its full within-bin covariance and explicit prior/geometry assumptions. A dust-field joint-probe test needs both estimators' responses; ordinary sums of released chi-squares cannot answer it.

These steps are ordered by identifiable evidence gaps, not by which change would favour acceleration or deceleration. No authors were contacted and no material was published. The linked scripts, local input hashes, validation records and source citations allow the checks to be repeated; the original observations and previous results remain intact.
