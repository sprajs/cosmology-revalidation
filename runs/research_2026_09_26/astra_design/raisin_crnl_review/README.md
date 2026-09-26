# RAISIN WFC3/IR count-rate nonlinearity: physical interpretation and next gate

2026-09-26. Independent Astra source/theory review; no image pixels, observed flux scores, light-curve fits, or corrections were produced. Primary documentation is saved with URLs and SHA-256 hashes in `acquisition.json`. `algebra_check.py` is a synthetic sign/background-limit check, not a RAISIN calibration estimate.

**This is a high-priority, independently testable detector-calibration question. It is not yet an established omission in RAISIN.** The exact public author converter starts after an external magnitude catalog; its missing upstream execution prevents classification as corrected or uncorrected. Sol's separate [inclusion audit](/home/szymon/Documents/ChatGPT/supernova/docs/research-2026-09-26/raisin-crnl-inclusion.md) documents that gap. A universal −0.04-mag correction is not justified.

## Source facts and what they establish

The current [WFC3 Data Handbook7.7](https://hst-docs.stsci.edu/wfc3dhb/chapter-7-wfc3-ir-sources-of-error/7-7-count-rate-non-linearity) says CALWF3 does not apply the count-rate correction. Its older ≈0.04±0.01-mag faint-source illustration assumes a bright-standard calibration and sky-dominated target. The [official photometry workflow9.1](https://hst-docs.stsci.edu/wfc3dhb/chapter-9-wfc3-data-analysis/9-1-photometry) lists aperture/encircled-energy, zeropoint and CRNL operations separately. Its time-dependent-sensitivity discussion is also separate; the still-future-tense reference to planned2024 processing must not be mistaken for an execution record in2026.

[Riess, Narayan & Calamida2019, ISR2019-01](https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs/_documents/2019/WFC3-2019-01.pdf) Table1 and body give the combined coefficient0.0075±0.0006 **mag/dex**, with no detected wavelength dependence. The report describes the rate effect becoming background-limited for very faint sources. Its white-dwarf comparison predicts F160W from other photometry/spectroscopy, omitting F160W from the SED fit; it is not a SN-distance calibration. Use the explicit mag/dex units: the abstract's percent shorthand is not an exact magnitude-to-flux conversion. The conclusion contains “0.006” where the table/body give “0.0006”; this review uses the latter. The older [ISR2010-07](https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs/_documents/2010/WFC3-2010-07.pdf) derives its0.04-mag illustration from an approximately0.01-mag/dex coefficient over4dex. Neither number is an observed RAISIN correction.

[ISR2020-10](https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs/_documents/2020/WFC3-ISR-2020-10.pdf) derives infrared zeropoints from five standard stars and describes revised CALSPEC spectra plus flat-field normalization. Its inverse sensitivities refer to infinite aperture. It provides no explicit faint-rate reference transfer in those sections. Therefore citing the2020 CALSPEC/zeropoint update is insufficient to demonstrate CRNL absorption. This is supported by the separately prescribed CRNL step, not merely a keyword search.

The [RAISIN paper§2.2,§3.7.1,AppendixA](https://arxiv.org/html/2201.07801v2) uses drizzled/subtracted FLT images,0.4-arcsec aperture photometry, P330E aperture checks, artificial stars and public bright-standard zeropoints. It reports a0.5% CALSPEC calibration uncertainty. These operations do not by themselves establish the rate correction. Post-detector artificial-star injections test extraction/noise; without an explicit detector-response forward step, they cannot measure how the incident faint-source photons were converted to electrons. A measured bright-star PSF may carry some rate-dependent shape, but rescaling it is not a validation of the faint/bright flux transfer.

The released FITOPT002 `HST_CAL` is `0.00714 × wavelength[µm]` in the [signed-input ledger](/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_differential/nml-signed-variants.json). It has no source or sky rate argument. It is neither an explicit nominal CRNL correction nor an identified CRNL prior. Partial overlap in distance space does not establish equivalent physics or justify adding another copy of its covariance.

## Response, subtraction and sign: explicit conditional model

The following is a transparent effective power-law model, not a microscopic detector derivation or an approved pixel correction. Let `q>0` be the magnitude dimming per decreasing rate dex, and let incident photon-equivalent rate `r` have response

`g(r)=A r (r/r0)^ε`, with `ε=q/2.5`.

Thus `q=.0075 mag/dex` means `ε=.003`. For a source signal `s_p` on local photon background `b_p`, the background-subtracted aperture response is

`F_meas = Σ_p w_p [g(b_p+s_p)−g(b_p)]`.

The calibration star has its own PSF, aperture, rate and background. Normalize using its response per incident source flux. Sky/host subtraction after detection cannot undo the preceding nonlinear response. If search and template backgrounds differ, static-host subtraction can leave further terms; this requires the actual exposure/template operator. Electronics offsets and dark current are not automatically equivalent to extra incident photons: that equivalence needs detector evidence.

For `s << b`, `g(b+s)−g(b) ≈ s g′(b)` and `g′(b)=A(1+ε)(b/r0)^ε`. As a source becomes still fainter at fixed background, its incremental response tends to a constant. Continuing a correction proportional to the logarithm of its sky-subtracted flux down to zero is therefore wrong in this model. A catalog total magnitude alone does not specify its effective detector illumination. Never raise noisy negative background-subtracted flux to this power or discard it; the physical positive incident signal belongs in a forward model before measurement noise.

For one effective pixel and small `ε`, define

`k(s,b) = [(b+s)ln((b+s)/r0)−b ln(b/r0)] / (s ln10)`.

Then the uncorrected relative magnitude bias is `δm ≈ q[k(anchor)−k(target)]`. The reference unit/pivot cancels. A fainter target in a lower illumination regime has `δm>0`: correcting **the target onto the bright-anchor system** makes it brighter (`Δm=−δm`, flux multiplier `10^(0.4δm)`). Alternatively moving the bright anchor onto a faint system has the opposite signed anchor adjustment; the relative calibration is the same. “Positive/negative correction” without specifying which side and pivot is ambiguous.

`algebra-result.json` verifies the first-order derivative against the exact response to5.46e−11 and the faint-background limit to5.04e−14. Its synthetic4dex/no-sky example is0.030±0.0024mag, flux multiplier1.028016 and fixed-shape distance multiplier0.986279. These are unit/sign illustrations, **not applied RAISIN values or uncertainties**.

## When another convention can absorb the effect

A faint-anchor strategy can remove the differential effect if its independently established incident flux scale, effective pixel illumination/background, PSF/aperture and bandpass match the science regime, and the executed transfer is carried into the published data and model. A faint standard catalog calibrated from the same bright zeropoint without a rate correction inherits the problem. Similarly, an empirical zeropoint may absorb the response at one effective pivot, but a single constant cannot remove an arbitrary residual rate dependence across both bright and faint regimes.

An aperture ratio derived from bright P330E can absorb part of the wing/core response, not automatically the absolute faint/bright transfer. A CALSPEC spectral update changes the reference SED; that alone is not the detector rate operator. A shared uncertainty with zero mean does not document a missing mean correction. Independence of a new standard-star check also requires checking whether its NIR SED or calibration was itself inferred using the same CRNL prior/data.

## Decision rule for the execution audit

| Classification | Minimum evidence required |
|---|---|
| Included explicitly | Executed pixel/catalog correction or documented effective-rate zeropoint transfer, sign/pivot/units, source version and calibration references, linked numerically from images/intermediate catalog to released J/H rows. |
| Included by cancellation | Quantified equality of target and anchor effective response, including background/PSF/filter, or a bounded differential response below a predeclared precision target. A faint label alone is insufficient. |
| Omitted in a defined chain | A complete linked image-to-catalog-to-release execution uses the bright-standard scale without an intervening rate transfer, and the relevant differential illumination is established. Absence in a post-photometry converter proves only that stage. |
| Unresolved | Missing upstream catalog/execution, undocumented global offsets, or uncertain calibration pivot. This is the present RAISIN status. |

Required next inputs are the executed producer of `../phot/raisin_photometry.txt`, source commits/commands, epoch-matched ZP/EE/CRDS files, aperture source counts, local sky/host rates and template handling, and the standard-star reference observations/transfer. Current library capability, `NLINCORR=COMPLETE`, the paper's nonmention, or a covariance-group name cannot decide inclusion.

## Propagation if inclusion and response become identified

Common slope across wavelength does **not** imply equal J/H magnitude offsets: source spectra, throughput, background and calibration pivots differ. Decompose the photometric response into a J/H mean and their difference. With shape/dust/peak fixed, the mean predominantly shifts NIR distance; colour changes can shift dust/shape and optical+NIR distances when those parameters are fitted. The whole-sample gray direction is degenerate with absolute magnitude. A common HST-only shift relative to CSP is not the same direction, because instrument and redshift are associated.

For fixed corrected photometric multipliers use `F′=SF`, `C′=SCSᵀ`; a rate-dependent transformation instead needs its Jacobian and shared calibration/background uncertainties. A known slope with uncertainty `σq` and response vector `k` contributes `σq² kkᵀ`, conditional on known rates and its prior origin. Do not add this independently to an existing mode that represents the same uncertainty. Other pivot/sky/SED uncertainties and cross-survey anchors may add correlations. New means require validated refitting and selection/bias recalculation; shifting final bias-corrected distances alone does not establish a corrected likelihood.

A local cosmology response, after such a distance perturbation is justified, is `δθ=(JᵀSJ)⁻¹JᵀSδμ`, with `S=W−W1(1ᵀW1)⁻¹1ᵀW` projecting the fitted intercept. This states the assumptions (fixed metric/model, local linearization); no cosmology derivative or fit is evaluated here. A percent-level possible flux effect is not evidence for a particular cosmological shift.

## Distinguishing this from2026 archive reprocessing

[Huynh etal., ISR2025-09, updated February2026](https://arxiv.org/html/2602.12110v1) describes a pixel-based **accumulated-count** nonlinearity solution, different ramp cosmic-ray flagging, an October2025 NLINFILE and its February2026 update. Its standard-star zeropoint changes are about0.1–0.2%; this does not bound individual faint-SN processing changes or supply CRNL.

The independently acquired current first/template primary headers contain `CAL_VER=3.7.3`, `NLINFILE=a2412448i_lin.fits`, `NLINCORR=COMPLETE`, and `IMPHTTAB=8ch15233i_imp.fits`. They have different PHOTFLAM values and a0.12825-arcsec drizzle scale, unlike the paper's0.11. Those are current-product facts, not the historical RAISIN recipe. The manual[NLINCORR description](https://hst-docs.stsci.edu/wfc3dhb/chapter-3-wfc3-data-calibration/3-3-ir-data-calibration-steps) is explicitly about integrated counts. Preserve rate nonlinearity, time-dependent sensitivity and accumulated-count correction as distinct operators, and compare a reproducible historical recipe before interpreting current-product differences.

## Highest-information next experiment

First close the inexpensive catalog/execution inclusion gate. If it stays missing, do not claim historical correction recovery. A separately frozen measurement can test the differential detector response without supernova distances: calibrated faint/bright standards whose NIR predictions are independent of the tested WFC3 measurements, or repeated identical sources with independently measured different photon backgrounds. Fix count/background/PSF support and model comparisons before outcomes, retain cross-exposure covariance and temporal/persistence nuisances, and validate on held-out objects/epochs. A SN brightness-versus-redshift trend cannot identify this operator independently of luminosity, selection, dust or phase-template changes.

Only after that response and the actual photometry convention are linked is a paired, source-consistent SN reduction warranted. No physical correction, distance likelihood, or cosmological conclusion has been produced in this review.
