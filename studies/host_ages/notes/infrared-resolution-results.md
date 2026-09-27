# What angular resolution permits in the recovered infrared images

The recovered SPIRE images do not independently determine a dust luminosity for each supernova host. Actual galaxy separations and measured telescope beams show why: even freeing only the nearest optical neighbour typically multiplies the conditional host-flux standard error by **2.54, 3.39 and 4.89** at 250, 350 and 500 μm. This calculation concerns source separation under an explicit image model. It is not a fitted dust flux, a calibrated uncertainty on the sky, or evidence that every optical neighbour emits in the infrared.

## Observations and selection

The sample is the existing **265 distinct DES hosts** with accepted eight-band photometry and valid measurements in all three SPIRE images. Selection uses source associations, masks and usable local coverage, without requiring a positive infrared flux. Both the host and neighbours use the released DES deep-field coordinates. Candidate neighbours are other distinct, unflagged, galaxy-classified sources within twice the nominal beam FWHM. The median nearest separation is **4.27 arcsec**; the median candidate counts are 75, 138 and 285 at the three wavelengths.

We downloaded the official [ESA 1-arcsec Neptune beam maps](https://archives.esac.esa.int/hsa/legacy/ADP/PSF/SPIRE/SPIRE-P/) and their [calibration description](https://archives.esac.esa.int/hsa/legacy/ADP/PSF/SPIRE/SPIRE-P/README.html). A fitted elliptical core determines each beam's centre rather than assuming that the header reference pixel equals the centre. Radial averaging gives empirical FWHM values **18.42, 25.09 and 36.88 arcsec**. These are calibration-source beams, not a reconstruction of the exact exposure-weighted coadd PSF for every galaxy spectrum.

The actual recovered SMAP image pixels are **6, 8.333 and 12 arcsec**, verified from their WCS. The primary calculation samples the beam at valid pixel centres within four nominal FWHM of the host and truncates its radial support at four FWHM. This is an explicit point-sampling approximation. A separate pixel-averaging sensitivity is not asserted to reproduce the actual mapmaking response.

## Conditional information, with no fitted image flux

For an image model $y=a_h h+N a_n+S b+\epsilon$, $h$ is the host beam template, $N$ contains neighbour templates and $S$ contains a constant and two linear sky terms. The conditional errors obey $\operatorname{Cov}(\epsilon)=\sigma^2 I$. All source amplitudes are unrestricted signed parameters; there is no positivity, flux-ratio or SED prior.

After projecting out the sky, let $\tilde h$ be the host template and let $r$ be its residual after projection onto the neighbour span. The fitted amplitude has conditional variance

$$\operatorname{Var}(\hat a_h)=\frac{\sigma^2}{r^\mathsf{T}r},\qquad
\mathcal I=\frac{\operatorname{SE}(\hat a_h\mid\text{neighbours free})}
{\operatorname{SE}(\hat a_h\mid\text{host and sky only})}
=\frac{\|\tilde h\|}{\|r\|}.$$

For one neighbour this reduces to $\mathcal I=(1-\rho^2)^{-1/2}$, where $\rho$ is the correlation of the two sky-projected templates. This is **standard-error inflation**, not variance inflation. A rank-revealing SVD computes the general projection without squaring its condition number.

| Wavelength | Median nearest-pair inflation | 5th–95th percentiles | Hosts with pair inflation >2 | Median all-candidate inflation |
|---|---:|---:|---:|---:|
| 250 μm | 2.54 | 1.60–4.68 | 209/265 | 153 |
| 350 μm | 3.39 | 2.06–6.28 | 255/265 | 805 |
| 500 μm | 4.89 | 2.90–9.17 | 265/265 | 1,845 |

Percentiles describe different host configurations; they are not confidence intervals. The all-candidate column deliberately frees every included optical neighbour. Its large values are conditional on that source model and are especially sensitive to beam and pixel approximations. They must not be read as measured uncertainties or universal limits on infrared inference. Informative source/SED priors can reduce uncertainty, while making their physical assumptions consequential.

![Empirical telescope beams are wide relative to the nearest-galaxy separation; conditional pair noise inflation increases with wavelength.](../../../docs/figures/infrared-source-separation.png)

*The left panel shows measured calibration-source radial beams. The right panel is the distribution of conditional nearest-pair standard-error inflation across the 265 host positions. It contains no fitted infrared fluxes or observational confidence intervals.*

Four rotations of the native two-dimensional beam give median nearest-pair inflation ranges of **2.49–2.62**, **3.28–3.46** and **4.67–5.07**. This checks beam asymmetry for the nearest pair; it does not validate the full multi-source radial approximation.

For two continuous Gaussian beams separated by $d$, their normalized overlap is $\rho=\exp[-d^2/(4\sigma_b^2)]$. Requiring pair inflation at most two gives

$$\mathrm{FWHM}\leq d\sqrt{\frac{4\ln2}{-\ln(3/4)}}\simeq3.105d.$$

Across these separations the median pair-only resolution requirement is **13.24 arcsec**, with 5th–95th percentiles **7.02–22.75 arcsec**. It is not sufficient for separating all neighbours, a complete observing proposal, or an exposure-time calculation.

## Numerical and physical limitations

The pair projection and analytic correlation identity agree to **8.2×10⁻¹⁴** relative precision across all 795 host/band configurations. An independent normal-matrix comparison is used only where conditioning permits: 184 configurations agree within **2.5×10⁻¹⁰**. Independent validation reconstructs source selection and geometry without the analysis KD tree and compares QR, a different SVD driver and original-unit orthogonality on 12 identifier-hash-selected hosts in all three bands. White-noise injections test the declared linear estimator, not real-image coverage. The [validation record](../results/infrared_resolution/validation.json) retains all numerical checks and model-boundary sensitivities.

Small empirical beam features can supply much of the formal separation information when many nearly overlapping templates are free. Pixel averaging, radial-profile smoothing, beam support, fit radius and neighbour radius are therefore explicit sensitivities. The [independent method review](../results/infrared_resolution/method-review.json) also retains six outcome-selected stress examples, clearly distinguished from the hash-selected validation sample. No numerical success certifies the underlying pixel or source model.

In the 12 hash-selected hosts per band, independent projections agree with the primary calculation within **5.2×10⁻¹³** relative precision, and changing the SVD rank threshold from 10⁻⁶ to 10⁻¹² changes none of their primary results. However, averaging the beam over nine subpixel points increases the median all-candidate inflation by factors **1.54, 2.76 and 3.48**. Extending the measured beam support from four to five FWHM changes it by factors **1.03, 1.65 and 10.20**. These are sensitivity ratios, not calibrated corrections to an uncertainty estimate.

At 500 μm, extending the candidate radius from two to three FWHM leaves the host amplitude numerically unresolved in **11/12** configurations; reducing the fitted image radius to three FWHM does so in **8/12**. Those cases are recorded as unresolved, rather than assigning a finite precision from the inverse of a floating-point residual. Thus the declared primary matrices are numerically solvable, but the available infrared information is not robust to all reasonable source and beam boundaries. The 72,000 white-noise estimator injections validate the stated linear calculation only.

The actual images have correlated instrumental and confusion noise. Optical candidates need not be infrared emitters; infrared-only sources can be absent from this list. Source extension, the exact coadd PSF, source-spectrum-dependent beam response and foreground backgrounds remain unmodelled. Deblended fluxes would need a joint image likelihood, additional infrared constraints and justified source-population information. Converting such fluxes into an age constraint additionally requires dust-heating and star-formation assumptions. None of those missing quantities is supplied by the geometric calculation above.

The result narrows the next measurement requirement: use spatially resolved or shorter-wavelength information and a validated joint source model before treating these beam measurements as host dust luminosities. [Acquisition and reproduction commands](../code/infrared_resolution/README.md).
