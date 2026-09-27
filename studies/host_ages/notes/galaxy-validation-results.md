# Galaxy spectra as an independent check of host ages

The photometric ages contain real information about stellar populations. On 609 distinct TITAN/ZTF hosts with reliable SDSS spectra, older inferred populations have stronger 4000 Å breaks and weaker Hδ absorption. This ordering survives control for host mass and redshift, and an independent integration of the released flux spectra gives the same result. It does **not** calibrate the absolute age scale in Gyr or establish an additional supernova luminosity correction.

The central spectrum is often a poor substitute for the supernova environment. Spatially resolved MaNGA data directly show lower 4000 Å breaks and stronger Hδ at the supernova positions than in the galaxy centres. Both the photometric age information and this aperture difference must be retained in a physical model.

## Observations and sample identity

We matched released host positions within 3 arcsec and redshifts within 0.001, requiring a unique physical SDSS photometric counterpart. Repeat spectra of that counterpart were reduced to the highest signal-to-noise spectrum; possible alternative counterparts were retained and excluded rather than resolved using their ages or brightness residuals. The resulting 810 host/sample associations comprise 738 ZTF, 30 G11 and 42 R19 rows, corresponding to 795 unique spectrum keys. The public SDSS tables returned 794 spectra; one requested key, 2129–54243–58 for ZTF19aarfnto, was absent from the joint spectral-table query. All 794 returned flux spectra were downloaded. Nine ambiguous host associations were excluded.

The R19 host coordinates were recovered from the original [author data repository](https://github.com/benjaminrose/MC-Age/tree/92713be96a89da991fe53bffcc596a5c0942fc37/data). G11 and R19 overlap and are analysed separately; their counts must not be added as independent galaxies. TITAN is the pinned public 8,610-row snapshot, not a reconstruction of its final published 6,983-host selection. Exact IAU names connect TITAN to ZTF; the spectroscopy verifies the ZTF host position, while TITAN's own file contains no coordinates for independently checking its host association.

All 64-bit SDSS identifiers are stored as exact decimal strings, including nullable photometric IDs. Native model, cmodel and fibre ugriz fluxes and inverse variances are retained for forward modelling. The existing DustPedia/Pantheon crosswalk did not supply usable matched host coordinates here; this is not evidence that those galaxies lack spectroscopy.

The selected spectroscopic sample is not representative of every TITAN host. Among 2,259 accepted TITAN/ZTF associations, the 609 with reliable indices have median log host mass 10.296, compared with 9.932 among the 1,650 without them, and median ages 6.82 versus 6.55 Gyr. No population transport correction is inferred from this selection.

## Age-sensitive features and dust

The following are physical-host bootstrap 95% intervals for Spearman correlations, conditional on the observed catalogue summaries.

| Age sample | Distinct hosts | Age versus Dn4000 | Age versus emission-corrected HδA |
|---|---:|---:|---:|
| TITAN/ZTF | 609 | +0.659 [0.605, 0.706] | −0.645 [−0.694, −0.589] |
| Updated G11 | 30 | +0.775 [0.617, 0.862] | −0.759 [−0.874, −0.545] |
| Updated R19 | 42 | +0.656 [0.389, 0.844] | −0.621 [−0.795, −0.368] |

For TITAN/ZTF, partial rank correlations after host mass and redshift are +0.638 [0.580, 0.693] for Dn4000 and −0.598 [−0.654, −0.535] for HδA. This validates an ordering along an observed stellar-population sequence. Age, metallicity, recent star formation and attenuation are not individually identified by these two features.

We independently integrated mean Fν over rest-frame 3850–3950 Å and 4000–4100 Å, using fractional pixel widths, positive inverse variances and unmasked pixels. Of 794 spectra, 779 have at least 90% usable coverage in both bands. Among 598 distinct TITAN hosts with these integrations, the age–Dn4000 rank correlation is +0.658 [0.605, 0.705].

The independently integrated foreground-corrected Dn4000 exceeds the MPA-JHU catalogue value by a median 0.00581; the central 68% object distribution is 0.00118–0.01069. This is an approximate measurement check, **not an exact reproduction** of the catalogue processing. Pixel covariance, spectral model subtraction and velocity-dispersion treatment are not all reproduced. The raw band means and variances remain available so other foreground models can be applied. Our catalogue comparison uses unscaled SFD reddening and O'Donnell Rv=3.1; the separate physical-age forward model uses its explicitly stated foreground prescription.

For 254 TITAN hosts with all four BPT emission lines detected above S/N=3 and a star-forming line-ratio classification, Hα/Hβ has median 4.067. Assuming intrinsic ratio 2.86 and F99 Rv=3.1 gives median nebular Av=0.942 mag. Its rank correlation with the SED Av summary is +0.724, but it exceeds that summary by median 0.755 mag. Nebular attenuation, stellar-continuum attenuation and supernova line-of-sight extinction are different measurements; their numerical difference is not a failed dust correction. MPA-JHU emission-line extraction itself uses a stellar-continuum model. See the [SDSS measurement description](https://www.sdss4.org/dr17/spectro/galaxy_mpajhu/).

## The aperture difference is measurable

For the 609 TITAN hosts, the legacy SDSS fibre has median physical diameter 3.37 kpc, and 411 supernovae lie outside its 1.5-arcsec angular radius. These physical sizes use a stated fiducial distance relation solely to describe apertures.

We acquired MaNGA maps for 53 distinct matched galaxies and sampled the actual supernova coordinates and host centres. Nearest 0.5-arcsec spaxels required S/N>3, unmasked indices and positive inverse variance. The [published dispersion corrections](https://www.sdss4.org/dr17/manga/manga-data/working-with-manga-data/) were applied, and failed or uncovered positions were retained in the ledger.

Among 34 distinct ZTF hosts with valid paired positions:

- Dn4000 at the supernova minus Dn4000 at the centre has median **−0.1234**, with 16th–84th object percentiles **−0.2526 to −0.00884**.
- The corresponding HδA difference is **+1.434 Å**, with object percentiles **+0.161 to +2.936 Å**.
- Median supernova–centre separation is **4.27 arcsec**, compared with median r-band seeing **2.48 arcsec**.

These percentile ranges describe galaxy-to-galaxy variation, not uncertainty on the medians. Seeing mixes populations, so spaxels are not independent. On the 28 TITAN hosts with valid local spectra, the global photometric age remains correlated with local Dn4000, +0.481 [0.093, 0.752], and local HδA, −0.618 [−0.829, −0.267]. The ordering persists, while the centres generally indicate older populations than the supernova sites. Neither measurement identifies the progenitor's birth environment or delay time.

The two [FIREFLY stellar-library fits](https://www.sdss4.org/dr17/manga/manga-data/manga-firefly-value-added-catalog/) also expose model dependence. On 49 distinct galaxies with valid central ages, MaStar minus MILES mass-weighted log age has median +0.027 dex and 16th–84th percentiles −0.105 to +0.178 dex. On 48 with both central and one-effective-radius measurements, median central-minus-outer-shell log age is −0.0668 dex with MILES and +0.0581 dex with MaStar. These are the same spectra under different stellar libraries; they are not independent replications. The outer measurement is an elliptical shell, not a global integrated age.

## Incremental supernova brightness information

The same 165 physical SNe have reliable spectral indices and satisfy the earlier explicit 401-object TITAN/ZTF brightness selection. Every nested comparison uses those same 165 objects and keeps each physical host in one held-out fold. The likelihood retains released SALT covariance, nuisance-dependent measurement variance and fitted intrinsic scatter. It conditions on SED medians; the missing joint age/dust/metallicity posterior is not invented.

After nonlinear SALT width/colour terms, a soft mass step, local colour, redshift and continuous host mass/Av/metallicity, the age coefficient is **+0.00198 ± 0.01993 mag/Gyr**. Adding Dn4000 and Hδ gives **+0.00515 ± 0.02044 mag/Gyr**. These are fitted standard errors, not 95% intervals. Independent propagation of the two index errors changes the latter to +0.00496 ± 0.02045.

Held-out RMSE is 0.22448 mag for host controls, 0.22557 after adding age, 0.23261 after adding spectral indices, and 0.23499 after adding both. Age changes held-out log predictive score by −1.45 [−2.70, −0.30] without the indices and −2.45 [−4.30, −0.69] with them. Those bootstrap intervals condition on the fitted overlapping training folds. This small selected sample supplies no incremental predictive gain. It does not demonstrate that age physics is absent: the controls can be mediators, and central indices are imperfect local proxies. Calibration correlations, survey selection and high-redshift transport remain outside this fit.

## Original versus revised age estimates

A labelled post-primary-inspection comparison tested the original and updated G11/R19 age orderings against the same spectra. G11's age–Dn4000 correlation is 0.7775 with the original ages and 0.7749 with C25; the paired difference is −0.0026 [−0.154, +0.164]. For R19 the correlations are 0.6149 and 0.6563, difference +0.0414 [−0.235, +0.295]. Hδ differences are also unresolved.

Thus both age estimates trace stellar populations. These observations do not establish that the revised age scale is more accurate, or resolve the previously identified sensitivity of age–brightness slopes to the age model. The R19 CDS age-column labels are misleading; this particular check uses rank order and does not exponentiate the numerical entries.

## Numerical checks and retained limits

Source hashes, exact identifiers, coordinate/redshift cuts, duplicate handling, host folds, Fν units and band-ratio identities were checked. The actual brightness likelihood was independently reevaluated and optimized. Acquisition and analysis are reproducible from the [code directory](../code/galaxy_validation/README.md); full tables, spectra and maps remain outside Git. Compact results are in [diagnostics-summary.json](../results/galaxy_validation/diagnostics-summary.json), [manga-local-summary.json](../results/galaxy_validation/manga-local-summary.json) and [historical-age-summary.json](../results/galaxy_validation/historical-age-summary.json).

A separate read-only solver audit checked the older BAO and RAISIN uses of NNLS after an inconsistency was found in a new underdetermined stellar-fit problem. The older square BAO designs passed 2,002 actual-design cases; all 117 overdetermined RAISIN designs agreed with independently column-scaled bounded least squares. This is not a universal solver guarantee.

The audit did find a distinct RAISIN reconstruction defect: known mass-step responses cancel algebraically after subtracting the mass correction, but their floating-point remainders had been amplified by fitted weights of order 10²⁵. Removing just those known cancelled directions raises the all-systematics relative covariance residual from **0.3017 to 0.4176** in NIR, **0.2761 to 0.3003** in optical and **0.0878 to 0.1316** in the combined branch. The mass-step-only reconstruction then correctly has residual 1.0. The original results remain preserved; direct paired brightness contrasts and the BAO conclusions are unchanged. This strengthens the limitation on reconstructing a physical covariance from those exported variants. Unscaled bounded least squares also reported success for inferior solutions in three of these ill-scaled cases, so changing solvers without checking scaling and residuals would not have fixed the problem. Details are in [legacy-nnls-audit.json](../results/galaxy_validation/legacy-nnls-audit.json).
