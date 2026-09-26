# Cosmology Revalidation

### A reproducible examination of measurements, corrections and expansion history

**Working manuscript · 26 September 2026**

[Run the workflows](docs/workflows.md) · [Methods](docs/methods/README.md) · [Results](results/README.md) · [Validation](validation/README.md)

## Abstract

We re-executed the complete reproducible workflow set: **19 workflows, five supplied alternative configurations and five additional random-seed runs**. All 391 input files passed their recorded hashes; all 19 default result summaries exactly matched the earlier edition. Separate numerical implementations checked the main distance likelihood, age regressions, population integrals, light-curve fits, calibration calculations, infrared comparisons and detector statistics. These checks found no new numerical defect in the audited calculations.

The released Pantheon+ distances give a matter-density fraction **Ωₘ = 0.3323 ± 0.0182**, consistent with the published flat-ΛCDM result. Acceleration remains favoured in the tested released-distance models. An imposed population-age correction substantially changes the inferred present expansion, but its amplitude and transfer to the observed supernova population are not established by these checks. Detector and selection investigations identify concrete weaknesses in how simple summaries can be interpreted: one pixel dominates a dark-image statistic, a simulation's near-perfect timing largely repeats its initializer, and discarding negative measurements creates a large bias in a controlled example.

This is a numerical revalidation of the **current research workflows**, with new independent diagnostics. It is not yet a complete reconstruction of the surveys or a unified cosmology measurement. The sections below separate observations, fitted quantities, assumptions and unresolved explanations.

## 1. What are we trying to measure?

A Type Ia supernova is not a perfectly standard light bulb. Its measured brightness depends on distance, its intrinsic properties, dust, the observing system and which events enter the sample. A useful accounting equation is

$$m_{\mathrm{observed}} = M_{\mathrm{intrinsic}} + \mu(z) + \Delta_{\mathrm{dust}} + \Delta_{\mathrm{calibration}} + \epsilon.$$

Here magnitudes increase as objects become fainter; μ is the distance modulus, z is redshift, and ε represents measurement variation. Selection changes which realizations of this equation are observed. Published supernova distances have already fitted or corrected several of these terms. Adding a further correction requires checking what has already been removed.

Our expansion diagnostic is the **deceleration parameter**, q. Negative q means accelerating expansion. A model can infer today's q₀ by extrapolation, while a flexible reconstruction may constrain only an average over a finite redshift interval. Those are different measurements. Supernovae with a freely fitted absolute-brightness offset constrain relative distances; they do not separately measure the absolute expansion rate H₀.

Baryon acoustic oscillations (**BAO**) supply another distance measurement, expressed relative to a sound-horizon ruler. The present joint fit leaves the ruler normalization free. It contains **no cosmic microwave background (CMB) likelihood**.

## 2. Data and validation strategy

The repository freezes the inputs and numerical settings needed for repeatable calculations. Some inputs are measurements, some are products of earlier fits, and some are simulations. Their roles must remain visible.

| Research component | Starting evidence | What this revalidation checks |
|---|---|---|
| Expansion and BAO | Released corrected distances and covariance matrices | Likelihood algebra, direct integration, shape constraints and sampler stability |
| Host age and populations | Published age summaries; specified star-formation and delay models | Extraction, matching, regression, numerical convolution and correction sensitivity |
| DES light curves | Calibrated fluxes, SALT3 model assets and frozen accepted epochs | Flux integration, fits, held-out prediction and shared calibration inference |
| RAISIN and CSP infrared work | Released distances, photometry, filter curves and archived simulations | Object pairing, correction signs, covariance reconstruction, timing and passband integrals |
| HST detector work | Current calibrated images, raw dark reads and reference files | Signed repeat extraction, geometric weights, reference-error propagation and spatial influence |
| Signed selection | Released signed photometry and an explicit synthetic generator | Baseline diagnostics and the consequences of dropping negative measurements |

We used three levels of validation. **Replay** establishes that the recorded calculation can be repeated. **Independent recomputation** changes the numerical route without importing the tested kernels. **Scientific comparison** checks whether the data selection, statistical model and measured quantity actually match the published claim. Passing the first two does not automatically establish the third.

The validated scientific inputs and published result files retain their original bytes. The repository has since been reorganized around the executable workflows; [provenance](provenance/README.md) records the original snapshot and the new locations. Detailed coverage and commands are in the [validation guide](validation/README.md); machine-readable evidence is in the [manifest](provenance/history/revalidation.json). No single “accuracy percentage” is assigned across these very different checks.

## 3. Expansion recovered from released distances

### 3.1 A direct check of the simplest fit

The redshift cut z > 0.01 leaves **1,590 distance rows representing 1,473 distinct supernova identifiers**. Multiple rows for an object remain in the supplied covariance; treating them as independent objects would be incorrect.

For spatially flat ΛCDM, we independently calculated luminosity distances with Astropy, used the full covariance, fitted out the unknown brightness offset, and integrated a dense one-dimensional posterior grid. This gives

$$\Omega_m = 0.33226 \pm 0.01821,\qquad q_0 = -0.50161 \pm 0.02732.$$

The quoted uncertainties are posterior standard deviations conditional on this model, covariance and prior. The best-fit matter density agrees with the workflow to **1.2 × 10⁻⁸**, and χ² agrees to **3.1 × 10⁻⁹**. Halving the grid density does not materially change the answer. A second sampler seed also agrees with the integrated posterior. [Calculation record](validation/reports/core.json).

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

All four models passed their within-ensemble autocorrelation checks and were rerun with another seed. That supports numerical stability of these summaries; it does not establish model adequacy. The fixed-correction case has about **20%** of sampled posterior mass at q₀ < 0, and its 95% interval spans both signs. It does not establish present deceleration at 95% credibility.

The models use different data or correction assumptions, so their raw χ² values are not a ready-made ranking. The joint analysis also cannot reproduce published **BAO + CMB + supernova** claims without the CMB contribution. [DESI DR2 cosmological analysis](https://arxiv.org/abs/2503.14738). Exact settings and results: [joint CPL](results/alternatives/joint-cpl/run.json), [flexible expansion](results/alternatives/flexible-expansion/run.json), [fixed-template experiment](results/alternatives/age-template/run.json).

### 3.3 A less model-specific BAO check

A separate calculation asks whether the anisotropic BAO measurements can satisfy a set of inequalities implied by flat, nonaccelerating expansion. An independent constrained optimizer reproduces the projection statistic **11.467067**. The conservative simulated cone-tail fraction is **0.4166**, with binomial 95% interval approximately [0.403, 0.430].

This particular test does not reject its composite null. It also does not establish nonacceleration: it has different assumptions, information and power from the parametric joint fit. The tail fraction is neither the probability that the universe decelerates nor a posterior for q₀. The test uses 12 anisotropic entries at six redshifts; it excludes the isotropic BAO entry. [BAO result](results/baseline/bao-shape/summary.json).

## 4. Host ages, bias corrections and dust

### 4.1 A discrepancy requires matching the quantity first

The age supplement contains 199 G11 and 102 R19 entries. Matching and the stated redshift/sample rules leave **196 unique objects** in the default diagnostic. The supplied alternative prioritizing R19 ages was also rerun.

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

The population calculation combines an assumed star-formation history with a distribution of delays between star formation and supernova explosion. Independent integration reproduces the mean and median delays for all three supplied delay models at z = 0, 1 and 2.5, with differences below **0.00003 Gyr**.

For the default delay model, the median delay changes by about **5.44 Gyr** between z = 0 and z = 1. Multiplying this by an imposed **0.030 mag/Gyr** slope produces a correction magnitude of **0.163 mag** at z = 1. A different supplied delay model gives about **0.074 mag**. Neither number is a measured correction for the selected survey population.

This distinction matters because host-population age, progenitor delay and supernova luminosity are not the same latent quantity. Dust, selection, redshift and already-applied corrections can change the mapping between them. The large cosmology shift in Figure 1 establishes sensitivity to a template; it does not identify the template as the right physical model.

The dust workflow checks absorption geometry and a latent-colour degeneracy. Independent depth integration reproduces the mixed-source transmission formula. These constructed examples explain how different dust and intrinsic-colour assumptions can resemble one another; they do not estimate a new empirical dust law. [Population outputs](results/baseline/populations/summary.json); [dust outputs](results/baseline/dust/summary.json).

## 5. DES: from calibrated flux to prediction and calibration

The flux workflow refits **24 full light curves** and evaluates **36 conditional held-out predictions**. Independent reconstruction agrees with the fit objectives and predictive densities to better than 5 × 10⁻¹⁴. Finer wavelength integration, alternative fit starts and Hessian step sizes give small changes. Differences from the frozen native reference reach about **0.0011 mag** in fitted brightness and **0.095 day** in peak time. These are bounded implementation comparisons on the same accepted epochs, not evidence that the original detector reduction was wrong.

The predictor experiment uses **850 training and 213 held-out supernovae**. Adding light-curve width to colour improves the summed held-out log predictive score by **40.10 nats**. Higher score means the model assigned more probability to the held-out observations. Both paired-object resampling and a field-level bootstrap keep the improvement positive. A second seed changes the default score by −0.020 nat, within its approximate Monte Carlo uncertainty. Four-chain rank-normalized convergence checks pass, with maximum split R-hat about **1.0023**.

This supports width as a useful predictor in this selected sample. It does not identify a universal physical luminosity law or validate a survey selection correction. Only ten observing fields are available for the field bootstrap. [DES numerical and statistical audit](docs/validation/des.md).

Shared calibration inference was independently rebuilt from its design matrices. The conditional uncertainty in the chosen distance contrast is **0.00942 mag** with systematics-only modes, **0.01207 mag** with inherited observer modes, and **0.01236 mag** with an isotropic observer prior. The corresponding means change materially too. Thus an apparently precise calibration result remains dependent on which residual patterns the model permits. We cannot uniquely assign every fitted mode to an instrumental calibration error or transport it directly into cosmology.

The inputs here are the frozen **DES-SN5YR** assets. They do not constitute a reconstruction of the later DES-Dovekie reanalysis. Comparing this diagnostic with a later published cosmology result requires matching the calibration release and full likelihood first. [DES-SN5YR paper](https://arxiv.org/abs/2401.02929); [DES-Dovekie reanalysis](https://arxiv.org/abs/2511.07517v3).

## 6. Infrared distances, timing and measurement selection

### 6.1 Paired distance accounting

The RAISIN comparison contains the same **79 supernovae**, with 42 at low redshift and 37 at high redshift, across the optical, near-infrared and combined releases. For optical minus infrared distance, the high-minus-low mean difference is **+0.00124 mag** as released. Reversing both exported mass and bias corrections changes it to **+0.07541 mag**.

Independent accounting verifies the correction signs: undoing the bias term adds the exported bias correction; undoing the mass term subtracts the exported mass correction. These operations retain the original selected sample. They do not redo the light-curve fit or selection.

A paired whole-object bootstrap gives a 95% interval of approximately **[−0.088, +0.085] mag** for the released optical contrast. The broad interval prevents interpreting its small central value as proof that all branches agree physically. Common calibration uncertainty and unknown cross-branch covariance remain outside this simple resampling calculation.

A second question asks whether published systematic covariance can be rebuilt from the available raw systematic distance-shift vectors. The relative reconstruction residual remains **13.45%** for infrared, **1.44%** for optical and **0.74%** for the combined branch. Allowing signed rather than nonnegative coefficients barely changes those residuals. This is a real reconstruction gap in the available representation. Missing historical vector transformations or covariance operations must be resolved before declaring the published matrix erroneous. [Infrared audit](docs/validation/infrared.md).

### 6.2 Timing and filter provenance

Of 30,000 infrared simulation rows, 29,995 join to the corresponding combined-fit simulation. **Every matched infrared peak time equals its supplied initializer.** The apparent 0.010-day timing residual therefore does not establish independently recovered timing precision. The combined-fit timing residual has a standard deviation of about **0.549 day**. Selection also changes simulated distance residual means, demonstrating why selected and unselected simulation summaries cannot be interchanged.

The CSP lineage audit recovers **5,491** raw/intermediate photometry rows with the documented filter relabelling. Three listed objects lack rows in the source archive; all three lie outside the 79-object cosmology cohort. They are coverage gaps, not demonstrated converter failures.

Independent analytic integration reproduces **140 passband calculations** to better than 3 × 10⁻¹⁵ mag. For the tested spectra, changing the J-band filter representation produces phase-dependent differences reaching about **0.069 mag**. This is a spectral sensitivity result. It is not yet a correction to a fitted supernova distance.

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

Independent integration of aperture geometry and re-extraction from FITS reproduce the reported measurements. The raw-read audit checks fourteen exposures, seven disjoint primary pairs and one overlapping sensitivity pair. All **11,856** signed aperture slopes agree to numerical precision. Including versus discarding temporal cross terms changes the relevant moment contraction by factors ranging from **0.668 to 2.314**. Covariance can change the direction, as well as the size, of an error in a diagonal approximation.

Those raw moments mix detector state, events and noise; they do not isolate electronic read noise. Eleven exposures have indeterminate exposure flags, so nominal read times alone do not establish correct physical timing. Current science headers also identify a 2026 nonlinearity reference, which does not reproduce the historical RAISIN reduction. Native CALWF3 processing itself was not rebuilt in this pass: its two dark output images are frozen starting inputs. [HST audit and source discussion](docs/validation/hst.md); [STScI calibration handbook](https://hst-docs.stsci.edu/wfc3dhb/chapter-3-wfc3-data-calibration/3-3-ir-data-calibration-steps); [2026 nonlinearity update](https://arxiv.org/abs/2602.12110).

## 8. What has been resolved, and what remains open?

| Question | Current assessment | What would change that assessment? |
|---|---|---|
| Do the calculations repeat? | Yes: all 19 default summaries match exactly; all supplied alternatives complete. | A changed input, code path or environment must receive its own validation record. |
| Is the simple released-distance cosmology numerically correct? | Independent integration, likelihood checks and conditional simulations agree. | A failure under upstream covariance or selection reconstruction would require revisiting the inference. |
| Does the age discrepancy prove either analysis wrong? | No: residual definitions, uncertainties and estimators differ. | Matched samples, original age distributions and an explicit correction ledger in a joint likelihood. |
| Is the infrared systematic covariance fully reconstructed? | No: a documented residual remains with the available vectors. | Historical vector operations and complete covariance-building instructions. |
| Is the dark-image variance anomaly a universal detector correction? | No: a single pixel dominates one statistic; apertures measure a different quantity. | Independent repeats that predict held-out measurements with their actual weights and covariance. |
| Has a physical correction overturned acceleration? | No such correction is identified here. A fixed template changes the answer conditionally. | A correction learned independently of cosmology, passing selection and transfer tests with uncertainty propagated. |

Several useful statistical methods strengthen this work: direct integration checks a sampler without trusting its draws; paired resampling compares models on the same objects; field and spatial grouping expose dependence; fixed-group deletion reveals concentration; constrained projection tests a broader class of expansion histories; and shared-mode inference makes calibration degeneracies explicit. These are established methods applied to this problem, not claims of newly invented statistics.

The new contribution of this pass is the complete fresh execution, independent audit implementations, recovery and influence diagnostics, and a clearer reconciliation of apparent discrepancies. The underlying research findings often predate this manuscript. Successful revalidation does not turn them into newly discovered physical effects.

## 9. Towards a unified measurement

The next useful experiment should constrain an ambiguity that currently prevents the pieces from being combined. The [full experimental programme](docs/experimental-plan.md) gives the detailed sequence. Its immediate priorities are:

1. **Predict independent detector repeats.** Estimate temporal and spatial covariance from timing-verified training exposures, then predict unseen aperture contrasts, including shared-template uncertainty. Compare reference versions in controlled native reductions.
2. **Recover injected sources through selection.** Inject signed fluxes before detection and fitting; retain rejected objects and epochs. Require recovery of brightness, timing and uncertainty coverage across brightness, colour, host properties and redshift.
3. **Fit populations and corrections together.** Use age probability distributions and shared latent dust, luminosity and calibration parameters. A proposed extra correction must account for corrections already applied and predict observations outside its training sample.
4. **Connect to distance and geometry.** Propagate shared calibration and selection uncertainty into the distance likelihood. Compare supernovae and BAO first; add an explicit, validated CMB likelihood only when making CMB-dependent claims.
5. **Test the complete inference before interpreting it.** Use repeated synthetic datasets, alternative plausible population generators, independent random seeds and survey/field holdouts. Record sensitivity to priors and missing data. Keep cosmology summaries blinded during upstream model choices where practical.

A measurement should advance through these steps only when its predictions survive the relevant held-out check. A more elaborate fit alone does not remove an identification problem.

## 10. Reproduce the results

The [workflow guide](docs/workflows.md) describes the executable research package. Its [source guide](docs/sources.md) records input provenance and acquisition limits. **A Git-only clone does not include all bulk scientific inputs or the new full sampler arrays.** They remain locally hashed; public inputs can be restored where acquisition metadata supports it, and derived bundles require their recorded upstream preparation. Missing inputs must not be replaced with invented data.

With the complete data bundle present:

```bash
uv sync --frozen --extra predictors
.venv/bin/python research.py verify

# Use a new name: existing workflow outputs are never overwritten.
.venv/bin/python validation/run.py --name my-revalidation
```

The last command performs all 29 runs and the four independent audit programs. It refreshes audit reports, so preserve the dated Git version before using it for another campaign. The figures in this manuscript are generated by `validation/figures.py` from the **published reference results**. See the [validation guide](validation/README.md) for partial reruns, coverage and interpretation.

| Read further | Purpose |
|---|---|
| [Workflow guide](docs/workflows.md) | Installation, commands and input requirements |
| [Methods](docs/methods/README.md) | Assumptions and outputs, ordered from measurements to cosmology |
| [Results](results/README.md) | Baseline analyses, alternative models and robustness runs |
| [Validation](validation/README.md) | Independent checks and reproducibility evidence |
| [Physical foundations](docs/physics/README.md) | Equations and their interpretation limits |
| [Experimental plan](docs/experimental-plan.md) | Next experiments towards a unified measurement |
| [Sources and provenance](docs/sources.md) | Input attribution, hashes and historical source records |

This repository is a working research record. Earlier exploratory work is preserved in Git history; [the provenance guide](provenance/README.md) explains how to retrieve it. Third-party material retains its original attribution and licensing; see the [license catalogue](provenance/licenses.json). No authors were contacted during this revalidation.
