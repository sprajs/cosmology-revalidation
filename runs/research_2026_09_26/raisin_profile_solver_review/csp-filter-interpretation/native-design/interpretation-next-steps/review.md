# What the stable filter response identifies, and the next decisive computation

This is a post-result scientific review of already saved artifacts, not a new result freeze. No native fit, synthetic score, observed dust posterior, calibration change or cosmology calculation was run. It preserves the source/measurement/model distinctions in the linked audits.

## What is established and what cannot be bounded

The independently verified full42 experiment identifies the **finite processing response of these released measurements**, under the fixed SNooPy B18 training realization, fixed headers, existing RC1/RC2 KCOR operators and a numerically stable native covariance iteration. Its mean raw distance response is +0.483432 mmag over42; the largest individual absolute response is5.489918 mmag. Ten exact unchanged controls and iteration/start closure establish that this particular finite response is not a failed optimizer artifact. They do not establish a physical upper bound on a passband or training error.

Let d(O,θ;y) denote this numerical distance recipe for operator O, trained model θ and data y, and let b(O,θ,P,S) denote its simulation bias correction under population P and selection S. At fixed physical membership, the changed corrected distance has the identity

Δ[d−b] = {d(O1,θ0;y)−d(O0,θ0;y)}
         + {d(O1,θ1;y)−d(O1,θ0;y)}
         − {b(O1,θ1,P1,S1)−b(O0,θ0,P0,S0)}.

The first bracket is the quantity measured here, with O1 specifically **RC2**, not a verified physical WIRC operator. The second includes the response of empirical training to the calibration convention. The third includes population and selection propagation, and is not known from raw fits. Physical WIRC minus RC2 would supply another operator term, and recalibrated photometric coordinates would require a consistently transformed y/error model. Changes in cohort selection or final weights need additional terms; the identity above holds for a fixed set of objects.

Without an independently justified admissible set for these unmeasured terms, the raw response is neither an upper nor a lower bound on the total correction, and its sign need not be retained. This is not evidence that the total correction is large. No training uncertainty magnitude has been measured here. A common grey training normalization alone cancels a common free intercept/high-minus-low contrast; chromatic, phase-dependent and sample-population training responses generally do not. Allowing arbitrary grey redshift evolution is a different identification problem, not established by this filter test.

A valid bound would require, for example, a source-linked retraining response plus an independently calibrated spectral/passband discrepancy set, and a generator/selection/bias response evaluated over that declared set. The envelope of those counterfactuals would be conditional on its specified set. Treating ±the measured RC1→RC2 shift as a prior width, summing unrelated sensitivity magnitudes, or calling a phase-contrast envelope a distance bound would supply missing assumptions without measuring them.

The baseline's poor native fit quadratics matter:105.490/14 and64.888/38 in the original affected/control pilot are model-fit limitations, not grounds for selecting another operator by smaller Q. Model-dependent covariance changes between arms. Numerical closure does not validate the latent template/population, and no normalized alternative likelihood was fitted.

## Priority and stop conditions

| Next path | What it can add | Required gate | Recommendation |
|---|---|---|---|
| Coherent archived-simulation timing pilot | Paired NIR amplitude response to an actual jointly fitted peak on the same reconstructed NIR measurements | Historical baseline, printed-input rounding, masks, stable amplitude solution and copy closure; peak/noise provenance retained | Highest priority for the next bounded native computation |
| Physical WIRC standard/passband bridge | A source-defined physical observation operator, separable from RC2 approximation | Exact primary reference realization and zero point, independent photon integration, nominal KCOR reproduction | Continue source/calibration closure now; native intervention only after this gate |
| Dust/population generative recovery | Sensitivity to distance/trained-colour assumptions under an explicit valid measurement/selection model | Signed cadence or correct censor/selection likelihood, complete noise/errors, model calibration and training overlap | Continue limited model-only algebra/recovery design; defer large observed population fits |

This ordering reflects current feasibility and identification value, not the size of a selected effect. The tiny raw RC2 response reduces the value of repeating that same fixed-model label perturbation. It does not diminish the need to establish instrument identity before calling a calibration exact.

### Timing: useful next computation, bounded interpretation

The newly recovered November2021 NIR/joint FITRES and LCPLOT form a coherent source-selected pair. Its eight-object **v11_04d baseline gate now passes**, unlike the preserved failed comparison of2022 FITRES to2021 LCPLOT. That makes the original first-eight pilot the most direct next engineering/scientific step after its still separate rounding screen. Keep all500 archive rows in the provenance ledger; only the already fixed first eight enter the pilot response, with no substitution based on outcome.

After those gates, compare the NIR distance obtained from near-truth archived fixed peak with the NIR distance obtained from that object's coherent archived joint optical+NIR fitted peak, holding reconstructed NIR measurements fixed. Keep the baseline and timing masks and any failures explicit, apply the same tested numerical convergence criterion to both, and report per-object and paired mean/dispersion responses. This is a conditional response screen, not a population bias correction. The first eight are an engineering subset, not a random or representative cosmology sample.

The joint peak uses NIR measurements, so its error and NIR amplitude noise are correlated. Matching IDs and truth/initial fields does not independently prove epoch-level cross-branch noise identity; original optical photometry is absent. The joint fit also starts near truth and may have iteration-recentered peak priors. It is a source-supported simulated alternative, not the identified creator or error distribution of observed headers.

Locally, for ε=estimated peak−true peak,

D(y,t+ε)−D(y,t) ≈ g(y)ε + ½h(y)ε².

A zero mean ε does not imply zero mean distance response: curvature and Cov(g(y),ε) can contribute. An independent Gaussian peak smear with the same marginal SD omits this correlation. Conversely, multiplying an observed finite derivative by the archived peak RMS is not a bias estimate. A later paired generative test needs complete optical+NIR draws, the actual peak-estimation and fit-selection operations, and common latent/noise draws across branches. Existing OIR empirical scatter may already contain timing-related variation; adding a new independent timing variance without checking its definition can double allocate scatter. Re-estimate the matched simulation bias and uncertainty before drawing a cosmological consequence.

### WIRC: finish the physical coordinate bridge before adding an operator

The measured WIRC throughput and a published count zero point exist. However, applying the pinned physical zero-point formula to BD17 gives RC2 magnitude8.420644, versus released8.4192: a1.444mmag reference mismatch that is not resolved by numerical integration. The same formula gives WIRC−RC2 BD17+1.647mmag. These are differences between documented conventions, not measured errors in the SN photometry. The source-derived P98 shared stellar transformation also does not make arbitrary SN spectra identical under WIRC and RC2.

The next decisive read-only step is the exact Vega/BD17/reference-spectrum and natural-system realization audit now assigned to Sol. If that closes, preserve a new distinct WIRC band, build a versioned KCOR with source-defined MAGREF, and require independent standard-star integration and unchanged nominal-band controls before any affected SN response. Reusing RC2's8.4192 with WIRC throughput alone would introduce an unmeasured calibration assumption.

Existing CSP-II measured spectra provide useful same-spectrum and phase-difference operator checks: grey normalization cancels. They do not identify the archived training model's response or an absolute distance correction. Object-level held-out spectra and actual dual-instrument measurements with calibration/telluric covariance are stronger calibration tests than adopting a larger learned template without an independent training ledger. A physical operator plus fixed old templates remains a diagnostic until training and the matched bias simulation are also addressed.

### Dust: identification first, no generic substitute

The public M20 broadband operator, proper one-dimensional grey-amplitude integral and model-only derivative gates are executable. They show that removing the external distance prior does not remove the trained intrinsic-colour prior, and that free intrinsic residual coordinates can closely mimic dust locally. Broad amplitude makes grey luminosity versus distance unresolved; it does not deliver cosmology-independent dust truth.

The released RAISIN optical data omit negative-flux rows, and precursor quoted errors are not always identical. An ordinary Gaussian likelihood on the retained positive rows is therefore not the validated measurement model. A defensible next recovery would use complete signed cadence/noise or a correct known-cadence selection/censor likelihood, the same calibrated forward operator, and fixed truth/prior scenarios chosen before outcomes. It should compare external-distance and proper broad-distance arms under the **same** training and noise model. No generic dust law/template replacement or observed posterior is justified by the current calibration result.

A small model-only proper-amplitude quadrature/recovery check can advance numerical readiness. A large population HMC run before signed-flux, selection and independent training/colour support would add precision conditional on unresolved assumptions rather than resolve the user's luminosity/dust concern.

## Evidence inspected

- `native-design/stable-filter-response/{review.md,result.json,manifest.json}` and root's independent `csp_native_filter_response/root-state-review/stable-response/result.json`.
- `docs/research-2026-09-26/csp-filter-bias-propagation.md`: generator/SIMLIB, optical/NIR selection, signed bias subtraction and unknown execution/training linkage.
- `docs/research-2026-09-26/csp-wirc-operator-feasibility.md`: physical throughput, standard-star integrals and unresolved reference realization.
- `docs/research-2026-09-26/csp-spectral-identification.md`: measured-spectrum constraints and independence/calibration limits.
- `docs/research-2026-09-26/raisin-timing-asset-recovery.md` and `astra_design/raisin_timing_assets/baseline_2021/baseline-gate.json`: recovered historical source pair and passing latest baseline.
- `docs/research-2026-09-26/raisin-nir-timing-sensitivity.md`: conditional observed header responses, not measured timing errors.
- `docs/research-2026-09-26/bayesn-distance-identification.md`: proper amplitude measure, inherited training, forward/Fisher checks and blocked observed likelihood.

This note's manifest freezes the exact reviewed artifacts. Later timing/standard-star results should be treated as new evidence, not silently retroactively included in this assessment.
