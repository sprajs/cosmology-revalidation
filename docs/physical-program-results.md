# What the new galaxy observations establish about an age correction

**27 September 2026.** Host ages trace real stellar-population differences, but the present observations do not yet determine a unique additional supernova brightness correction. Independent spectra validate age ordering; spatially resolved spectra reveal a measurable difference between galaxy centres and supernova sites; flexible stellar-population fits expose how weakly broad-band light determines absolute age without restrictions on star-formation history. The central issue is now the physical mapping from these observations to residual supernova luminosity, including dust and survey selection.

The original [reconciliation and recovery results](age-correction-results.md) remain the starting point. This extension recovers actual spectra, local environments, high-redshift photometry and infrared images and tests them against physical models. It does not replace missing author likelihoods with unlabelled Gaussian ages, or change the baseline cosmology by appending an unvalidated age template.

## 1. The age information is real, but the absolute scale is not certified

On **609 distinct TITAN/ZTF hosts**, the correlation between published photometric age and the observed narrow 4000 Å break is **0.659**, with a physical-host bootstrap 95% interval **[0.605, 0.706]**. The relation survives control for host mass and redshift, and independent integration of the actual spectra gives essentially the same answer. Hδ absorption provides a complementary association in the opposite direction.

Those measurements support a stellar-population ordering. They do not separately identify age, metallicity, recent star formation and attenuation, nor do they measure a supernova progenitor's delay time. On the same 30 G11 and 42 R19 hosts, original and updated age estimates both correlate with spectroscopy; the paired comparisons do not establish that the revised ages are more accurate.

In the **165 common supernovae** with usable spectra and the declared ZTF brightness selection, the incremental age coefficient is **+0.00198 ± 0.01993 mag/Gyr** after the stated host and light-curve controls. Including spectral indices gives **+0.00515 ± 0.02044**. Neither addition improves held-out brightness prediction. These are conditional standard errors in a small selected sample. The brightness contrast is internal to the release and assumes its documented common-zero-point blinding description; the exact x0 transformation remains unrecovered, so this is not a certified unblinded luminosity test. The result does not prove that age physics is absent: some controls may absorb its effects, and the measured central population is often not the explosion environment. [Spectroscopic analysis](../studies/host_ages/notes/galaxy-validation-results.md).

## 2. Central spectra do not fully describe the supernova environment

Among the 609 hosts, **411 supernovae lie outside the three-arcsec fibre's radius**. MaNGA maps provide a direct spatial check. In **34 distinct ZTF galaxies**, Dn4000 at the supernova position minus the central value has median **−0.1234**; HδA changes by **+1.434 Å**. Both changes indicate a systematically younger-looking population at the explosion site under ordinary stellar-population interpretation. The reported spreads across objects are not confidence intervals on a universal correction.

Global photometric age still correlates with local spectral features in the smaller matched set. Thus “global ages contain no information” and “global ages exactly measure progenitor age” are both stronger than these observations support. Seeing, galaxy structure and stellar-population modelling must be included explicitly.

![Independent age ordering, central-versus-local spectral features, and broad physical age compatibility regions.](figures/host-age-physics.png)

*Panel A uses 609 distinct hosts; its interval is a bootstrap 95% interval for rank correlation. Panel B shows released spectral-index measurement errors for 34 spatial pairs. Panel C shows the median and 5th–95th percentiles of age-region widths among statistically compatible fits, not uncertainty on each sample's median. The samples have different redshifts and age ceilings.*

## 3. Allowing flexible star formation changes what a precise age means

We fit signed galaxy fluxes with nonnegative mixtures of stellar populations. A convex optimization finds the youngest and oldest formed-mass mean ages compatible with the observed flux errors, over a declared dust/metallicity family. This avoids imposing one narrow parametric star-formation history.

| Observed aperture sample | Eligible | Compatible with statistical errors alone | Median compatible age-range width |
|---|---:|---:|---:|
| SDSS central fibres | 788 | 768 | 12.15 Gyr |
| Roman local SDSS-family apertures | 485 | 430 | 10.93 Gyr |
| DES eight-band hosts | 331 | 124 | 7.10 Gyr |

These are conditional confidence regions for a specified finite physical family, not Bayesian age posteriors or 1,604 independent galaxies. The large widths show that broad-band light alone does not tightly identify an arbitrary star-formation history. Stronger age constraints can be obtained by adding physical information or priors, but their contribution must remain visible.

The model also fails an informative precision test: 191 of 698 usable fibre comparisons have disjoint photometric and measured Dn4000 regions. Adding an **assumed** 0.03-mag model/calibration term makes the mismatch much rarer and admits 307/331 DES hosts, but does not make their ages precise. That sensitivity is not an empirically calibrated error model. Grid refinement, alternative mass weighting and foreground conventions retain the broad-age conclusion in the tested subset. A different cosmological clock changes the upper age bounds, so the clock must be propagated if ages are used to infer cosmology. [Equations, physical assumptions and numerical tests](../studies/host_ages/notes/physical-ages-results.md).

## 4. High-redshift host data exist; transport is still a measurement problem

The recovered DES deep fields contain **1,906,261 galaxies**. Secure matching and quality/association requirements yield **331 selected supernova hosts with eight-band photometry**, including **182 at z ≥ 0.6**. This establishes observed high-redshift support beyond the older low-redshift age tables.

The input audit found two consequential representation problems. Roman's full-precision SDSS-family local fluxes use maggies despite a contradictory text-table unit label; the conversion is recoverable from their flux–magnitude identity. Older DES host magnitudes differ from the corrected official release by median **0.61–0.78 mag** across griz. The new deep photometry and an independent published catalogue agree much more closely with the corrected release. These findings concern host measurements; they do not by themselves establish an error in supernova distances.

The available high-redshift hosts cannot simply be treated as a random sample. Coverage differs across fields, and the tested availability model does not explain those differences well in held-out fields. Matching redshift and global colour/mass moments also fails to improve held-out prediction of local colour across SDSS and SNLS. Inverse-probability weights bring some observed means into agreement, but that does not establish that age is independent of missingness. [Data recovery, calibration and transport analysis](../studies/host_ages/notes/host-transport-results.md).

## 5. Infrared images constrain what can honestly be attributed to dust

We recovered the public SPIRE images behind the proposed optical-to-infrared extension. **265 of the 331 eight-band DES hosts** have usable measurements at 250, 350 and 500 μm. Signed negative values are retained. Typical local scatter is **6.5–7.2 mJy**, exceeding the quoted instrumental errors, with substantial cross-band correlation.

At 250 μm, 260/265 hosts have another optical galaxy within half a beam FWHM; all 265 do at the two longer wavelengths. An optical neighbour is not necessarily an infrared emitter, but the measurements do not identify which galaxy supplies each beam's light. Assigning all positive flux to the host, or treating non-detections as independent zero-flux observations, would invent dust information.

These observations support a next physical measurement: simultaneous image modelling of hosts and neighbours, with actual PSFs, correlated confusion, shorter-wavelength infrared information and explicit dust-heating assumptions. Until then, the SPIRE measurements are conditional beam constraints, not deblended host luminosities or exact reproduction of the unreleased 501-host photometric catalogue.

## 6. A physically injected age term is only partly absorbed

A grey luminosity term applied before measurement and selection is not automatically erased in the tested survey model. With physically positive fixed-Rv = 3.1 dust, the frozen nominal correction transmits **−0.02988 ± 0.00131 mag/Gyr** of an injected −0.030 slope. Independently refitting standardization and the correction leaves **−0.02120 ± 0.00224 mag/Gyr**, approximately **71%**. Both are paired responses relative to nominal closure; the baseline itself is not proven perfectly unbiased.

The experiment uses 60,000 generated attempts across nine native simulations, including 18,000 in the physically supported dust population. It propagates the intervention into noisy photons, detection, host-redshift acceptance and native light-curve fitting. Correction training and evaluation have independent seeds; the uncertainty resamples both pools. A separate review of 29,576 noiseless paired epochs verifies the injected sign and amplitude to 3.8 × 10⁻⁶ mag.

This is a conditional counterexample to automatic complete absorption, **not an observed age correction**. The simulated progenitor delays are assumptions, mock host masses have no measurement error, the trained light-curve surface is fixed, and the primary correction is a declared nearest-neighbour model. The supported redshift-bin mean shifts remain imprecise. Native redshift-only BBC completes with positive-definite conditional covariance; production-dimensional BBC fails genuine training/interpolation support at this Monte Carlo volume. The pure-Ia experiment does not establish classifier contamination or BEAMS calibration. [Experiment, native results and failed gates](../studies/host_ages/notes/survey-physics-results.md).

## 7. What remains before a unified cosmological correction

The evidence supports neither a universal claim that dust standardization automatically removes every age effect nor a claim that the full historical age template must be added. A residual slope depends on the age definition, selection, existing nuisance fits and the physical population admitted by the data.

The remaining scientific requirements are specific:

1. Model the **local** stellar population and its dust/metallicity/SFH uncertainty, validated against actual spectra and spatial information. A central-age summary is insufficient.
2. Separate stellar attenuation, nebular attenuation and supernova line-of-sight extinction. Infrared image confusion and shared photometric calibration must enter the likelihood.
3. Infer the selected high-redshift host population with supported observation and follow-up selection, preserving absent and ambiguous counterparts.
4. Test physical population changes before detection, refit standardization and correction training on independent samples, and establish interpolation support and coverage. A sparse surrogate correction or a fixed trained light-curve surface must remain labelled.
5. Only then propagate the supported residual distance bias and its correlated uncertainty into the cosmological likelihood, including the dependence of inferred stellar ages on cosmic time.

The released-distance acceleration result is unchanged. The new observations substantially sharpen which assumptions fail or remain unmeasured, but they do not yet identify a replacement distance correction or a complete selection-aware cosmology measurement.
