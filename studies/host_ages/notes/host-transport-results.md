# Observed host populations across surveys

We recovered public observed photometry for 881 nearby and distant supernova environments and eight-band measurements for 331 selected DES supernova hosts. These data replace an assumed cosmic host population with measurements of actual selected hosts. They do not, by themselves, identify stellar ages, progenitor delay times, or an additional supernova brightness correction.

The analysis also exposed two consequential catalogue problems: Roman et al.'s SDSS-family FITS fluxes are in maggies despite a nanomaggies label, and older DES host magnitudes were superseded by a documented correction. Neither finding establishes an error in the released supernova distances.

## Observations and sample construction

| Input | Observed support | What is available |
|---|---:|---|
| [Roman et al. 2018](https://arxiv.org/abs/1706.07697), [CDS catalogue](https://cdsarc.cds.unistra.fr/ftp/J/A+A/615/A68/ReadMe) | 881 objects: 396 SNLS5, 389 SDSS, 55 CfAIII, 34 CfAIV, 7 CSP; z=0.01–1.06 | Local 3-kpc fluxes/errors, local/global magnitudes, published rest-frame U−V and mass summaries |
| [DES Y3 deep fields](https://des.ncsa.illinois.edu/releases/y3a2/Y3deepfields), [Hartley et al.](https://arxiv.org/abs/2012.12824) | 1,906,261 galaxies in C3, E2, X3 | Calibrated signed DECam ugriz and VISTA/VIRCAM JHKs model fluxes/errors, foreground values and measurement flags |
| [DES-SN5YR public release](https://github.com/des-science/DES-SN5YR/tree/c9a4fcafc4cbd19bd750dee47fc76194a45c181f) | 1,623 DES events in the actual Dovekie metadata | Exact event IDs, selected-sample membership, host coordinates and corrected griz magnitudes |

The release's descriptive count is not used as a selection mask: the pinned metadata contain 1,623 `IDSURVEY=10` rows. Every one matches an exact numeric event ID in both the original candidate catalogue and the updated SMP catalogue. We match released **host centroids**, not supernova coordinates, to individual deep-catalogue galaxies.

The frozen primary rule requires a nearest counterpart within 1 arcsec, no second counterpart within 1 arcsec, unflagged optical/NIR measurements and masks, a galaxy classification, and host directional-light-radius distance below four. There are 474 unambiguous selected-event counterparts, 331 passing the quality/association rule, with 331 distinct deep-catalogue hosts. Tightening the matching radius to 0.5 arcsec gives 327 quality hosts; widening it to 2 arcsec gives 320 because the ambiguity veto removes newly competing counterparts. No flux-detection threshold is imposed. Explicit ±9.999×10^9 native missing-value placeholders are removed; real negative flux measurements remain. None of the 331 selected hosts contains these placeholders.

The primary hosts comprise 28 at 0.06≤z<0.3, 121 at 0.3≤z<0.6, and 182 at 0.6≤z<1.2. The complete candidate catalogue is heterogeneous: an unselected candidate is not automatically a rejected Type Ia supernova. This catalogue therefore does not supply the discovery/classification selection function.

## Calibration and measurement interpretation

For all positive SDSS-family local flux/magnitude pairs, `m+2.5 log10(f)=0` within 10^-9 mag. In the catalogue-declared AB system, the native flux is therefore a **maggie**, requiring multiplication by 10^9 to obtain nanomaggies. This numerical identity establishes the unit scale; it does not independently establish that small native-SDSS-to-AB zero-point offsets were already applied. The paper ties its images to SDSS field-star magnitudes, so absolute zero-point uncertainty remains a separate sensitivity. The ASCII table rounds many of these fluxes to zero; the FITS table preserves them. The reusable interface retains 48 legitimate negative local measurements. Global fluxes derived from magnitudes have only first-order Gaussian flux errors; the original measurement representation remains a magnitude likelihood.

SNLS local measurements instead satisfy a native Vega zero point of 30 within 5.1×10^-6 mag, consistent with the paper. The catalogue's additional conversion description has not been reconciled to an absolute AB calibration. We preserve native SNLS flux and block its use as calibrated AB input. The paper applies Planck Milky Way extinction during fitting but does not release per-object foreground values in this table. Local and global apertures overlap; their errors cannot be multiplied as independent likelihoods.

The DES deep catalogue uses AB zero point 30 in every band, including JHKs. Dividing its native flux by 1,000 gives AB nanomaggies. The flux/magnitude identity is accurate to 7.2×10^-15 mag. The `CALIB` and `DERED_CALIB` fluxes are both retained. Their ratio gives the released foreground coefficients A_band/E(B−V)_SFD98 = 3.9627, 3.186, 2.140, 1.569, 1.196, 0.705, 0.441 and 0.308 in ugrizJHKs. A second foreground correction must not be applied to the latter product.

Pinned public SVO response curves supply the correct DECam ugriz and VISTA/VIRCAM JHKs instrument families. They approximate the response of the deep coadds; they are not the actual exposure/CCD-weighted response. SVO's Vega reference metadata are not applied to these AB-calibrated measurements. The source model shares aperture shape across bands, and joint deblending/calibration covariance is not released.

The [SMP README](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/0_DATA/DES-SN5YR_DES/DES-SN5YR_DES.README) documents a December 2025 host-magnitude correction. For 329 eligible common hosts, median old-minus-corrected griz magnitudes are +0.779, +0.700, +0.624 and +0.606 mag. The deep-minus-corrected medians are instead +0.028, −0.024, −0.030 and −0.014 mag. An exact event-ID/legacy-name link gives 206 matches to [Wiseman et al. 2020](https://cdsarc.cds.unistra.fr/ftp/J/MNRAS/495/4040/ReadMe), with median corrected-minus-published differences below 0.0014 mag. The broad tails between deep and transient-free host measurements remain much larger than formal magnitude errors. Aperture, blending and possible transient contamination must be modeled or tested; a small median offset does not establish interchangeable likelihoods.

## What transports across the observed samples?

The Roman rest-frame U−V colours are outputs of its published SED interpolation. They are **proxies conditioned on that procedure**, not independent age measurements. Of 881 objects, 800 have usable local/global colours, errors and mass summaries. Their mean local-minus-global colour differences are +0.0307 mag at 0.06≤z<0.3, +0.0101 at 0.3≤z<0.6, and +0.0474 at 0.6≤z<1.2. Unique-supernova bootstrap intervals are respectively [0.0046,0.0565], [−0.0193,0.0356] and [0.0194,0.0777] mag. These intervals describe sampling variation, not missing joint SED covariance.

In the common 0.2≤z<0.4 range, we balanced 166 SDSS hosts to the mean redshift, global U−V and stellar mass of 61 SNLS hosts. The effective SDSS sample falls to 40.3. The transported local-minus-global offset differs from the actual SNLS offset by −0.0858 mag, with bootstrap interval [−0.2027,+0.0279]. The reverse transport has effective sample 21.9 and difference +0.0386 [−0.0112,+0.0916] mag; one of 500 bootstrap targets lies outside the sampled source convex hull and is explicitly recorded as unsupported.

An SDSS-trained linear prediction of local colour from the same global observables has held-out SNLS RMSE 0.190 mag, worse than using the target's global colour directly (0.143 mag). The reverse prediction also fails to improve materially (0.258 versus 0.256 mag). Matching three mean observables therefore does not validate exchangeability of local environments across these surveys. No colour difference is converted to Gyr or a distance correction.

## Availability of eight-band measurements

There are 701 already-selected DES events assigned to C3, E2 or X3, counting composite release fields such as C1+C3. Their eight-band quality coverage is 144/272, 66/167 and 121/262 respectively. We froze a ridge-logistic availability model using redshift, corrected host i magnitude and observed g−i, and predicted each field using only the other two.

| Held-out field | Measured availability | Predicted mean | AUC | Brier score; constant-source baseline |
|---|---:|---:|---:|---:|
| C3 | 0.529 | 0.479 | 0.515 | 0.260; 0.258 |
| E2 | 0.395 | 0.472 | 0.571 | 0.242; 0.249 |
| X3 | 0.462 | 0.473 | 0.596 | 0.242; 0.249 |

The model misses field-specific completeness. Descriptive inverse-probability weighting has effective sample 320 and moves the available mean redshift from 0.613 to 0.587, close to the full selected-field mean 0.590. It similarly moves mean i from 21.272 to 21.087 versus 21.073 in the parent. This agreement does not establish that missingness is independent of age: morphology, blending, observing depth and unmeasured stellar properties remain. These weights are not an empirical correction of the cosmic or selected age distribution.

## An independent infrared measurement

[Ramaiya et al. 2025](https://arxiv.org/abs/2509.12035) identifies a useful physical extension: optical-to-Spitzer/Herschel photometry constrains emission from dust as well as attenuated starlight. Its paper points to public survey data but does not link its assembled 501-host flux table, gold-sample mask or joint posterior samples. Its reported supernova results are therefore not replicated here.

We nevertheless recovered and measured the underlying public images. The HELP flat-file endpoint returned HTTP 402 and its VO service timed out. The [NASA IRSA HerMES archive](https://irsa.ipac.caltech.edu/data/Herschel/HerMES/) supplies the exact original CDFS-SWIRE-NEST and XMM-LSS-NEST SMAP v6.0 files identified in HELP's pinned [source-map index and processing code](https://github.com/H-E-L-P/dmu_products/tree/0e451df5da3ec721b87f56cae5efd2fef17cf30d/dmu19/dmu19_HELP-SPIRE-maps). HELP copies the original image values without rescaling; we use those original images with their native error and flag maps.

The frozen extraction samples each known optical host position, preserving signed nearest-pixel values in Jy/beam. It subtracts the mean of 128 common nearby sky positions, placed deterministically in a 60–180 arcsec annulus, and retains median subtraction as a sensitivity. The full three-band covariance of these positions measures local sky/confusion variation. It is not divided by 128: the goal is the variation of a possible contaminating beam, not an artificially precise mean. The sky samples are themselves spatially correlated, so this empirical covariance is a diagnostic rather than an exact Gaussian likelihood. A Gaussian beam calculation separately records the roughly 5% typical pixel-centering attenuation; it is an explicit response approximation.

There are 1,100 selected DES hosts within the rectangular map bounds. After native flags and local-sky coverage requirements, 1,003, 1,004 and 1,004 have usable 250, 350 and 500 μm samples. The existing eight-band sample contributes 265 hosts with all three infrared bands.

| SPIRE band | Eight-band hosts | Negative background-subtracted values | Median sky scatter | Median instrumental error | Above three times local sky scatter |
|---|---:|---:|---:|---:|---:|
| 250 μm | 265 | 93 | 6.52 mJy | 2.79 mJy | 19 |
| 350 μm | 265 | 114 | 7.18 mJy | 2.67 mJy | 16 |
| 500 μm | 265 | 131 | 6.95 mJy | 3.15 mJy | 10 |

The last column is a descriptive threshold count, not a calibrated detection significance: the sky distribution is asymmetric and confused. Median local correlations are 0.744 between 250/350 μm, 0.553 between 250/500 μm and 0.725 between 350/500 μm. Treating the three bands as independent would overstate information.

Among these 265 hosts, 260 have another unflagged optical galaxy inside half the 250 μm FWHM; all 265 do at 350 and 500 μm. This does not prove each optical neighbor emits in the infrared, but demonstrates why beam flux cannot be assigned uniquely to the supernova host. Measurements, coordinates, neighbor counts, signed sky samples and covariance are supplied for subsequent joint image modeling. We do not silently interpret a positive beam value as intrinsic host dust luminosity or replace negative values with upper limits.

A defensible next physical inference requires PSF- and confusion-aware deblending using all plausible neighboring sources, shorter-wavelength infrared measurements or an explicit dust-temperature/emissivity prior, and a model for stellar/AGN heating and attenuation geometry. Even a well-measured total dust luminosity would constrain an energy budget; it would not uniquely recover an arbitrary star-formation history or progenitor age. This branch establishes actual observing support and these limitations, rather than fitting an extra age correction to distance residuals.

## Availability and validation

The [TITAN release page](https://titan-snia.github.io/dr1.html) still displays “Stay tuned!” at the recorded check. A deeper inspection of the author's public repositories confirms the previously recovered `age-of-titans` snapshot at `9b9ed9e5faf8b2f6267c0dce8cea97704d9b95ea`. Its historical host summaries remain useful observations conditioned on the author's model; they do not provide the final 6,983-object selection or joint host SFH/dust posteriors. Unrelated galaxy-imaging files in another public repository do not establish availability of those likelihoods. Negative availability statements here are limited to the specifically inspected public sources.

Independent checks verify publisher MD5 checksums for all three deep catalogues; input, code and output SHA-256 hashes; native flux/magnitude and foreground identities; direct angular separations; entropy-balancing moments; and held-out AUC by pairwise comparisons. The infrared check independently computes the tangent-plane projection, checks exact native pixel/error/flag values, and reconstructs every eligible covariance from the saved signed sky samples. These checks establish calculation integrity, not the physical adequacy of unmeasured covariance, approximate response curves or conditional stellar models.

The explicit missing-value repair changes the full DES parent interface but leaves every flux and its order in the 331-host subset unchanged. `selected-photometry-integrity.json` records that proof and hashes the selected-only table. The raw downloads and generated observation tables stay in ignored `.work/host-transport`; the repository contains acquisition code, pinned source hashes, compact findings and interpretation.

## Reproduction and interfaces

Run from the repository root. This branch uses an isolated Python 3.12 environment and does not alter the repository's main locked environment:

```bash
uv venv --python 3.12 .work/host-transport/.venv
uv pip install --python .work/host-transport/.venv/bin/python -r studies/host_ages/code/host_transport/requirements.txt
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4
.work/host-transport/.venv/bin/python studies/host_ages/code/host_transport/acquire.py
.work/host-transport/.venv/bin/python studies/host_ages/code/host_transport/roman.py
.work/host-transport/.venv/bin/python studies/host_ages/code/host_transport/deep.py
.work/host-transport/.venv/bin/python studies/host_ages/code/host_transport/smp.py
.work/host-transport/.venv/bin/python studies/host_ages/code/host_transport/filters.py
.work/host-transport/.venv/bin/python studies/host_ages/code/host_transport/proxies.py
.work/host-transport/.venv/bin/python studies/host_ages/code/host_transport/selection.py
.work/host-transport/.venv/bin/python studies/host_ages/code/host_transport/infrared.py
.work/host-transport/.venv/bin/python studies/host_ages/code/host_transport/validate_infrared.py
.work/host-transport/.venv/bin/python studies/host_ages/code/host_transport/validate.py
```

The source downloads total several GB. Reuse is hash-checked; `acquire.py --verify-only` verifies cached inputs without retrieval. A changed upstream byte stream fails the pin rather than silently replacing an input.

The main generated interfaces are:

- `roman-photometry.csv`: object ID, redshift, aperture, native flux/magnitude/error, instrument family, calibration and foreground status; use `aperture=local`, `filter_family=SDSS`, `measurement_available=true` for directly usable AB local fluxes.
- `des-deep-photometry.csv` and `des-deep-selected-photometry.csv`: exact event and host IDs, coordinates, redshift, eight signed observed/dereddened flux/error pairs, flags and foreground; the latter preserves exactly 331×8 selected measurements.
- `filters/des-deep-responses.npz`: wavelength in Angstrom and photon throughput for each actual instrument family, explicitly approximate coadd responses.
- `des-spire-photometry.csv`, `des-spire-covariance.json` and `des-spire-sky-samples.npz`: signed beam samples, instrument errors, local background statistics, empirical covariance, centering response and optical-neighbor diagnostics. These are not deblended host likelihoods.

Machine-readable results are in [`../results/host_transport`](../results/host_transport); the separately frozen protocols and reusable programs are in [`../code/host_transport`](../code/host_transport).


The separate physical-model review is also reusable. It identified a native observer-grid quadrature error reaching 0.0253 mag in a DES colour; the replacement operator uses Gaussian quadrature over both redshifted stellar-spectrum and passband knots. Independent 0.1 Å uniform integration now agrees within 1.73×10^-7 mag over 36 instrument/dust/metallicity/redshift/foreground combinations and four ages per case. The same review identified five valid Roman measurements excluded by a positive-total-flux condition; they are retained in the revised physical analysis, with uninformative zero-compatible cases treated analytically. These corrections are documented with the physical-age results.

After building the separate physical-age library and its isolated environment, reproduce that cross-check with:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 .work/physical-ages/.venv/bin/python studies/host_ages/code/host_transport/review_physical_operator.py
.work/host-transport/.venv/bin/python studies/host_ages/code/host_transport/validate.py --include-physical-review
```

The default host-data validation remains runnable without the physical-age library. The final validation record explicitly states whether the additional physical-operator proof was included.
