# Cosmology Revalidation

### Measurements, corrections and the evidence for cosmic acceleration

**Working manuscript · 26 September 2026**

[Methods](docs/methods/README.md) · [Research studies](studies/README.md) · [Data and code](docs/workflows.md)

## Abstract

The released Pantheon+ supernova distances give **Ωₘ = 0.3323 ± 0.0182** in a spatially flat ΛCDM model, closely matching the published value of 0.334 ± 0.018. The corresponding present deceleration parameter is **q₀ = −0.5016 ± 0.0273**: these distances favour accelerating expansion under that model. Combining supernovae with baryon acoustic oscillations in a more flexible dark-energy model also favours acceleration. A separate inequality test of BAO alone is inconclusive.

Host-age trends depend strongly on which brightness corrections are included. In 196 matched objects, the full-covariance residual–age slope changes from **−0.00495 ± 0.00456** to **−0.01294 ± 0.00456 mag/Gyr** when the exported bias correction is reversed. An assumed population-age template can shift the inferred q₀ to **+0.063**, but its 95% interval crosses zero and the template is not an empirically established correction. Optical and infrared distances show a related dependence on correction accounting: their high-minus-low-redshift contrast changes from **+0.0012 to +0.0754 mag** when the exported mass and bias terms are reversed.

The measurement investigations reveal specific limitations. Light-curve mean fluxes agree closely between independent calculations, but some uncertainty prescriptions disagree. A native local peak-time Hessian understates the profile-supported uncertainty by a factor of **4.55** in one supernova. An apparent infrared timing precision largely repeats the supplied initializer. In an HST dark-image comparison, one pixel contributes **82.9%** of the squared difference, making a universal detector-error correction unjustified. These findings identify weaknesses in particular estimators and interpretations; they do not establish a new cosmological correction or a failure of acceleration. The calculations are independent reanalyses of shared published measurements, not independent observations.

## 1. What are we trying to measure?

A Type Ia supernova is not a perfectly standard light bulb. Its measured brightness depends on distance, its intrinsic properties, dust, the observing system and which events enter the sample. A useful accounting equation is

$$m_{\mathrm{observed}} = M_{\mathrm{intrinsic}} + \mu(z) + \Delta_{\mathrm{dust}} + \Delta_{\mathrm{calibration}} + \epsilon.$$

Here magnitudes increase as objects become fainter; μ is the distance modulus, z is redshift, and ε represents measurement variation. Selection changes which realizations of this equation are observed. Published supernova distances have already fitted or corrected several of these terms. Adding a further correction requires checking what has already been removed.

Our expansion diagnostic is the **deceleration parameter**, q. Negative q means accelerating expansion. A model can infer today's q₀ by extrapolation, while a flexible reconstruction may constrain only an average over a finite redshift interval. Those are different measurements. Supernovae with a freely fitted absolute-brightness offset constrain relative distances; they do not separately measure the absolute expansion rate H₀.

Baryon acoustic oscillations (**BAO**) supply another distance measurement, expressed relative to a sound-horizon ruler. The present joint fit leaves the ruler normalization free. It contains **no cosmic microwave background (CMB) likelihood**.

## 2. Measurements and assumptions

We compare three levels of evidence: released distances; calibrated optical and infrared fluxes; and selected detector images and raw reads. A released distance already contains model fits and corrections. A flux fit tests more of that construction, while a detector comparison tests only the exposures and measurement weights actually examined.

The distance analysis uses Pantheon+ and DESI DR2 BAO. The measurement studies use the frozen DES-SN5YR release, RAISIN and CSP photometry, published host-age tables, and selected HST observations. The DES data are not the later Dovekie reduction, and the current HST calibration references are not identical to those used for the historical RAISIN distances. The [source catalogue](docs/sources.md) specifies the releases.

Unless stated otherwise, uncertainties attached to cosmological parameters are posterior standard deviations; regression errors are conditional standard errors; and intervals marked 95% use the method specified alongside them. Simulated examples test an explicit generator. They are not measurements of the corresponding bias in a survey.

## 3. Expansion recovered from released distances

### 3.1 A direct check of the simplest fit

The redshift cut z > 0.01 leaves **1,590 distance rows representing 1,473 distinct supernova identifiers**. Multiple rows for an object remain in the supplied covariance; treating them as independent objects would be incorrect.

For spatially flat ΛCDM, we calculate luminosity distances, retain the full supplied covariance, and fit out the unknown absolute-brightness offset. Direct integration of the matter-density posterior gives

$$\Omega_m = 0.33226 \pm 0.01821,\qquad q_0 = -0.50161 \pm 0.02732.$$

The quoted uncertainties are posterior standard deviations conditional on this model, covariance and prior. Direct integration and independent likelihood calculations agree to numerical precision. [Calculation record](validation/reports/core.json).

The published Pantheon+ supernova-only value is Ωₘ = 0.334 ± 0.018. Our mean differs by about **0.10 of that quoted standard deviation**. This is descriptive agreement between closely related analyses, not an independent tension test. [Brout et al., Pantheon+ cosmological constraints](https://arxiv.org/abs/2202.04077).

We also generated 200 Gaussian datasets with Ωₘ = 0.33 and the frozen covariance. Mean recovery bias was **+0.00036**, with Monte Carlo standard error **0.00122**; 138/200 nominal 68% intervals covered the injected value. This tests conditional recovery at one truth. It does not test whether the actual covariance or survey selection model is correct.

A remaining diagnostic deserves attention: χ² = **1402.92** for approximately **1588 degrees of freedom**, or χ²/dof = **0.883**. Under a fixed Gaussian covariance and a locally linear mean model, the lower-tail probability is about 0.00033. The residuals are smaller than that simple reference expects. Covariance construction, fitting and selection need investigation before attributing this to any physical correction. We have **not** rescaled the errors to force χ²/dof to one.

### 3.2 How much does the answer depend on the model?

“CPL” below means a dark-energy model with equation of state w(z) = w₀ + wₐz/(1+z). The flexible model fits five redshift bins of q with its stated smoothing prior. The final row adds a fixed age-based distance template; it is a sensitivity experiment.

| Data and assumptions | Mean expansion diagnostic | 95% posterior interval |
|---|---:|---:|
| Supernovae, flat ΛCDM | q₀ = **−0.502** | [−0.554, −0.447] |
| Supernovae + BAO, flat CPL | q₀ = **−0.436** | [−0.590, −0.276] |
| Supernovae, flexible expansion | q over 0 < z < 0.1 = **−0.326** | [−0.588, −0.064] |
| Supernovae + BAO, CPL with imposed age template | q₀ = **+0.063** | [−0.081, +0.209] |

![Conditional expansion estimates under four sets of assumptions; the imposed correction changes the inference, while its interval spans zero.](docs/figures/expansion.png)

*Figure 1. Points are posterior means and lines are 95% intervals. The third row concerns a finite redshift bin, not an independently measured instantaneous q₀. The fourth row conditions on an assumed correction without propagating uncertainty in its physical construction.*

The fixed-correction case has about **20%** of sampled posterior mass at q₀ < 0, and its 95% interval spans both signs. It does not establish present deceleration at 95% credibility.

The models use different data or correction assumptions, so their raw χ² values are not a ready-made ranking. The joint analysis also cannot reproduce published **BAO + CMB + supernova** claims without the CMB contribution. [DESI DR2 cosmological analysis](https://arxiv.org/abs/2503.14738). Exact settings and results: [joint CPL](results/alternatives/joint-cpl/run.json), [flexible expansion](results/alternatives/flexible-expansion/run.json), [fixed-template experiment](results/alternatives/age-template/run.json).

### 3.3 A less model-specific BAO check

A separate calculation asks whether the anisotropic BAO measurements can satisfy a set of inequalities implied by flat, nonaccelerating expansion. The distance from the allowed set, measured with the released covariance, is **11.467067**. The conservative simulated cone-tail fraction is **0.4166**, with binomial 95% interval approximately [0.403, 0.430].

This particular test does not reject its composite null. It also does not establish nonacceleration: it has different assumptions, information and power from the parametric joint fit. The tail fraction is neither the probability that the universe decelerates nor a posterior for q₀. The test uses 12 anisotropic entries at six redshifts; it excludes the isotropic BAO entry. [BAO result](results/baseline/bao-shape/summary.json).

## 4. Host ages, bias corrections and dust

### 4.1 A discrepancy requires matching the quantity first

The age supplement contains 199 G11 and 102 R19 entries. Matching and the stated redshift/sample rules leave **196 unique objects** in the default diagnostic. Objects appearing in both catalogues use G11 ages in this comparison; the alternative catalogue choice is treated separately.

With age treated as exact, the relation between Pantheon+ residual and age is weak after the released corrections. Reversing the exported bias term changes the slope:

| Residual definition | Diagonal weighted slope | Full-covariance slope |
|---|---:|---:|
| Released corrected residual | −0.00417 ± 0.00642 | −0.00495 ± 0.00456 |
| Exported bias correction reversed | −0.01327 ± 0.00642 | −0.01294 ± 0.00456 |

Slopes are in **magnitudes per billion years**; uncertainties here are conditional standard errors. Neither row includes uncertainty in the ages themselves. The change in the diagonal slope, **−0.00910 mag/Gyr**, is exactly the contribution of the reversed bias term under the same regression weights. It is an accounting identity, not proof that the original correction removed a genuine age effect.

The published age analysis uses older residuals associated with its G11/R19 samples, a redshift adjustment for G11, and an errors-in-variables regression. Our corrected Pantheon+ residual is a different outcome. Regressing the literal author residuals on the matched sample already gives a steeper diagnostic slope, approximately **−0.0274 mag/Gyr**, even before recreating the authors' estimator. Comparing −0.0042 directly with a published age slope would therefore conflate changed residuals, errors and methods. [Chung et al., source analysis](https://arxiv.org/abs/2411.05299).

We also checked a Gaussian latent-age likelihood using a separate joint block covariance calculation. Its numerical likelihood is correct, but one fit reaches the intrinsic-scatter boundary. Such an optimum is not a calibrated detection significance. Exact reproduction of the published regression still requires the original age probability distributions and full fitting choices. [Age audit](validation/reports/core.json).

### 4.2 The uncertainty columns are not interchangeable

In the selected age sample, the median quoted distance-modulus error is **0.196 mag**, while the median square root of the full covariance diagonal is **0.139 mag**. We traced this difference back to the original author release: relevant table columns match literally, and the full covariance has the same SHA-256 hash. The extraction did not introduce the difference.

The release distinguishes plotting errors from the covariance required for cosmological fitting. We retain the full covariance in the cosmology likelihood and report diagonal-error age regressions as separate diagnostics. An unexplained difference between uncertainty representations is worth documenting, but substituting one for the other would silently change the analysis. [Release comparison](validation/reports/release-comparison.json); [Pantheon+ release instructions](https://github.com/PantheonPlusSH0ES/DataRelease/tree/7fc6805/Pantheon%2B_Data/4_DISTANCES_AND_COVAR).

### 4.3 A population template is not an observed correction

The population calculation combines an assumed star-formation history with a distribution of delays between star formation and supernova explosion. These describe the ages of potential progenitors before survey selection, not the measured ages of the observed host galaxies.

For the default delay model, the median delay changes by about **5.44 Gyr** between z = 0 and z = 1. Multiplying this by an imposed **0.030 mag/Gyr** slope produces a correction magnitude of **0.163 mag** at z = 1. A different supplied delay model gives about **0.074 mag**. Neither number is a measured correction for the selected survey population.

This distinction matters because host-population age, progenitor delay and supernova luminosity are not the same latent quantity. Dust, selection, redshift and already-applied corrections can change the mapping between them. The large cosmology shift in Figure 1 establishes sensitivity to a template; it does not identify the template as the right physical model.

Dust geometry and intrinsic colour also remain degenerate in the constructed models: different mixtures can produce similar observed colour–brightness relations. A separate physical check finds that extrapolating the tested extinction prescriptions to sufficiently low Rᵥ can yield **negative extinction**. A passive absorbing screen cannot brighten a source in that way. This is a failure of the extrapolated model domain, not a measurement of the resulting cosmology bias; clipping the extinction would define another model whose population and selection effects would need to be assessed. [Dust findings](studies/dust/README.md). [Population outputs](results/baseline/populations/summary.json); [dust outputs](results/baseline/dust/summary.json).

## 5. DES: from calibrated flux to prediction and calibration

Two independent fitting paths for twelve DES supernovae give **24 full light-curve fits**. They closely reproduce the calibrated mean fluxes; finer wavelength integration and alternative fit starts change the answers only slightly. Differences from the frozen native reference reach about **0.0011 mag** in fitted brightness and **0.095 day** in peak time. These are bounded implementation comparisons on the same accepted epochs, not evidence that the original detector reduction was wrong.

Agreement in mean flux does not imply agreement in uncertainty. For **CID1896213**, the native local Hessian gives a peak-time uncertainty of **0.124 day**, while an independent smooth Hessian and a profile of the likelihood give **0.564 day**. The native MINOS interval and published scalar uncertainty also agree with the larger value. The discrepancy therefore concerns that local Hessian estimate, not every published timing error. Replacing the effective SALT colour law with a fixed F99-shaped alternative does **not** give a decisive held-out predictive improvement in this twelve-object sample. [Light-curve findings](studies/light_curve_fitting/README.md).

The predictor experiment uses **850 training and 213 held-out supernovae**. Adding light-curve width to colour improves the summed held-out log predictive score by **40.10 nats**. Higher score means the model assigned more probability to the held-out observations. Both paired-object resampling and a field-level bootstrap keep the improvement positive. The fitted colour coefficient is **2.473 ± 0.078**, conditional on this selected sample and predictor model.

This supports width as a useful predictor in this selected sample. It does not identify a universal physical luminosity law or validate a survey selection correction. Only ten observing fields are available for the field bootstrap. [DES numerical and statistical audit](docs/validation/des.md).

Allowing shared calibration patterns changes both the estimated distance contrast and its uncertainty. The conditional uncertainty in the chosen distance contrast is **0.00942 mag** with systematics-only modes, **0.01207 mag** with inherited observer modes, and **0.01236 mag** with an isotropic observer prior. The corresponding means change materially too. Thus an apparently precise calibration result remains dependent on which residual patterns the model permits. We cannot uniquely assign every fitted mode to an instrumental calibration error or transport it directly into cosmology.

The inputs here are the frozen **DES-SN5YR** assets. They do not constitute a reconstruction of the later DES-Dovekie reanalysis. Comparing this diagnostic with a later published cosmology result requires matching the calibration release and full likelihood first. [DES-SN5YR paper](https://arxiv.org/abs/2401.02929); [DES-Dovekie reanalysis](https://arxiv.org/abs/2511.07517v3).

## 6. Infrared distances, timing and measurement selection

### 6.1 Paired distance accounting

The RAISIN comparison contains the same **79 supernovae**, with 42 at low redshift and 37 at high redshift, across the optical, near-infrared and combined releases. For optical minus infrared distance, the high-minus-low mean difference is **+0.00124 mag** as released. Reversing both exported mass and bias corrections changes it to **+0.07541 mag**.

The correction accounting uses the signs in the released definitions: undoing the bias term adds the exported bias correction; undoing the mass term subtracts the exported mass correction. These operations retain the original selected sample. They do not redo the light-curve fit or selection.

A paired whole-object bootstrap gives a 95% interval of approximately **[−0.088, +0.085] mag** for the released optical contrast. The broad interval prevents interpreting its small central value as proof that all branches agree physically. Common calibration uncertainty and unknown cross-branch covariance remain outside this simple resampling calculation.

A second question asks whether published systematic covariance can be rebuilt from the available raw systematic distance-shift vectors. The relative reconstruction residual remains **13.45%** for infrared, **1.44%** for optical and **0.74%** for the combined branch. Allowing signed rather than nonnegative coefficients barely changes those residuals. This is a real reconstruction gap in the available representation. Missing historical vector transformations or covariance operations must be resolved before declaring the published matrix erroneous. [Infrared audit](docs/validation/infrared.md).

### 6.2 Timing and filter provenance

Of 30,000 infrared simulation rows, 29,995 join to the corresponding combined-fit simulation. **Every matched infrared peak time equals its supplied initializer.** The apparent 0.010-day timing residual therefore does not establish independently recovered timing precision. The combined-fit timing residual has a standard deviation of about **0.549 day**. Selection also changes simulated distance residual means, demonstrating why selected and unselected simulation summaries cannot be interchanged.

The CSP lineage audit recovers **5,491** raw/intermediate photometry rows with the documented filter relabelling. Three listed objects lack rows in the source archive; all three lie outside the 79-object cosmology cohort. They are coverage gaps, not demonstrated converter failures.

For the tested spectra, changing the J-band filter representation produces phase-dependent differences reaching about **0.069 mag**. This is a spectral sensitivity result. It is not yet a correction to a fitted supernova distance.

### 6.3 Why signed measurements matter

In a controlled simulation with true amplitude **0.4**, the three methods recover:

| Treatment of measurements | Mean recovered amplitude, 128 simulations |
|---|---:|
| Keep all signed measurements | **0.4069** |
| Keep positive values and ignore the selection | **1.2023** |
| Model censoring with the omitted-epoch schedule known | **0.3960** |

The correctly censored and fully signed methods differ by −0.0110, with paired Monte Carlo standard error 0.0089. The naive positive-only estimator has a large upward bias under this generator. All released positive fluxes alone cannot establish which censoring mechanism a real survey used; the omitted observations and their selection probabilities matter.

The real signed-baseline check covers ten selected objects. Its residual χ²/dof is **1.056** with a 180-day exclusion around peak, and **1.029** with a 365-day exclusion. Whole-object bootstrap intervals include one. Most of the modest excess is associated with one object. These nested cuts and repeated measurements are not independent confirmations of a universal noise correction. [Selection and baseline audit](validation/reports/raisin.json).

## 7. HST: a measurement depends on its spatial and temporal weights

The science-image experiment uses eight current calibrated images, four in each of two visits. One pair per visit defines the source masks; the other pair supplies a signed repeat difference. There are **170 retained aperture sites**, but only **one held-out temporal repeat pair per visit**. More apertures do not create more independent detector realizations.

The squared repeat difference divided by the propagated diagonal error variance is **0.708** for the search visit and **0.729** for the template visit. The independently propagated flat-reference contribution is only about **0.005%** of pair variance, far too small to explain a 27–29% deficit. Unknown spatial covariance, shared terms and the sampling of detector states remain open explanations. The strict all-pixels-valid mask has no support; these results concern the documented masked-coverage sample.

The dark-image result appears contradictory until we inspect what is being averaged. Over 993,750 eligible pixels the variance ratio is **4.202**, while 247 fixed background-subtracted apertures give **0.547**. One pixel contributes **82.9%** of the squared pixel difference and has zero weight in every aperture. Removing fixed detector blocks makes the pixel ratio range from **0.720 to 4.216**. This sensitivity is evidence against summarizing the pair with a universal error multiplier.

![Science and dark variance ratios, with ranges from deleting fixed spatial groups. The dark-pixel result has far greater spatial sensitivity.](docs/figures/spatial-influence.png)

*Figure 2. Dots use the full eligible samples. Lines show post-hoc deletion of fixed spatial groups, not confidence intervals. Pixel and aperture statistics apply different spatial weights. None of the deletions replaces the primary estimate.*

The raw-read comparison covers fourteen exposures, seven disjoint primary pairs and one overlapping sensitivity pair, with **11,856** signed aperture slopes. Including versus discarding temporal cross terms changes the relevant moment contraction by factors ranging from **0.668 to 2.314**. Covariance can change the direction, as well as the size, of an error in a diagonal approximation.

Those raw moments mix detector state, events and noise; they do not isolate electronic read noise. Eleven exposures have indeterminate exposure flags, so nominal read times alone do not establish correct physical timing. Current science headers also identify a 2026 nonlinearity reference, which does not reproduce the historical RAISIN reduction. The dark slope comparison conditions on two fixed CALWF3 output images; it does not independently establish the correctness of the entire detector reduction. [HST audit and source discussion](docs/validation/hst.md); [STScI calibration handbook](https://hst-docs.stsci.edu/wfc3dhb/chapter-3-wfc3-data-calibration/3-3-ir-data-calibration-steps); [2026 nonlinearity update](https://arxiv.org/abs/2602.12110).

## 8. Agreement, disagreement and unresolved interpretation

The released-distance result is reassuringly close to the published Pantheon+ cosmology. It says that the stated distance likelihood supports the stated flat-ΛCDM answer. It does not independently verify the survey reductions, population models or covariance construction that supplied those distances.

The host-age comparison is less settled. We find a weak trend in the corrected Pantheon+ residuals and a stronger one when the exported bias term is removed. Chung et al. find a stronger relationship using different residuals and age-error modelling. The difference is real between these calculations, but it is not yet a controlled disagreement with their estimator. Nor does the accounting change prove that a survey bias correction has removed an astrophysical age effect. That requires a joint explanation of the measured ages, luminosity population and selection.

Some conclusions are more direct. The peak-time Hessian and likelihood profile disagree for a specific object; the profile is supported by an independent fit and the native MINOS result. The infrared simulations do not demonstrate independent 0.010-day timing recovery because the fitted values reproduce the supplied initializers. The dark-image pixel statistic is too concentrated to support a detector-wide variance multiplier. These are identifiable failures of particular uncertainty summaries or interpretations.

Other attempts remain inconclusive. The available systematic shifts do not fully reconstruct the infrared covariance, but the missing historical transformations prevent identifying the source of the gap. The broader BAO inequality test does not reject nonacceleration under its assumptions. Alternative colour laws do not establish a preferred physical dust correction in the small DES pilot. Classifier reconstruction and several native infrared fits have not met the accuracy or optimizer-stability requirements needed for reliable selection inference. Failed or incomplete identification is part of the scientific result, not evidence for a particular alternative cosmology. [Supporting investigations](studies/README.md).

## 9. Conclusions

**Acceleration is recovered from the released distances under the tested uncorrected distance models.** The flat-ΛCDM matter density agrees with the published result, and a supernova-plus-BAO fit with evolving dark energy also gives negative q₀. The broader BAO-only test supplies no corresponding rejection of its nonaccelerating class.

**Corrections can matter enough to change the cosmological interpretation, but their physical identification remains incomplete.** The assumed age template moves the central q₀ estimate across zero, yet neither its amplitude nor its transfer to the selected supernova population has been established. We have not measured a correction that overturns acceleration, isolated a new dust law, or obtained a selection-complete unified cosmology measurement.

**The main unresolved issue is connecting the measurement levels without losing their shared uncertainties.** Detector repeats constrain particular spatial and temporal combinations. Light curves add spectral, timing and calibration assumptions. Population and selection models turn those fits into distances. A defensible joint measurement must propagate those dependencies rather than treating each conditional result as an independently measured correction. The [experimental programme](docs/experimental-plan.md) specifies the measurements needed to resolve these ambiguities.

## Data and code availability

The [reproduction guide](docs/workflows.md) gives installation and execution instructions; the [methods](docs/methods/README.md) describe the statistical assumptions. Original acquisition, extraction, fitting and diagnostic code is organized by topic in [studies](studies/README.md), including investigations with negative or unresolved results. Downloadable inputs and third-party implementations are identified by their [sources and hashes](docs/sources.md), rather than bundled as research code. Compact reference summaries and manuscript figures are retained; bulk data, chains and intermediate products are acquired or regenerated locally. The distinction between the validated manuscript calculations and additional research requiring upstream inputs is explicit in the study guide.
