# Stellar-age feasibility in 55 calibrated host spectra

The five measured spectral bands do **not** determine a useful formed-mass stellar age when the star-formation history, metallicity and attenuation are allowed to vary across the declared mixture family. All 55 host spectra have compatible models. For 54 hosts, the allowed set includes both an age below 1 Gyr and an age above 10 Gyr. The median interval width is 13.495 Gyr, close to the entire assumed 0.005–13.5 Gyr library range; even the narrowest interval spans 12.465 Gyr. These are ranges of models compatible with these measurements, not estimates that the galaxies truly occupy the extremes. [Compact results](../results/calibrated_host_physics/summary.json), [individual hosts](../results/calibrated_host_physics/host-age-feasibility.csv).

This calculation uses the 55 independently associated DESI DR1 host spectra described in the [calibrated-host observations](calibrated-hosts.md). Their redshifts span 0.162–0.955. Each physical host enters once. All signed flux measurements are retained, including the objects for which a displayed break or absorption index failed its signal-to-noise gate. The analysis uses five overlapping bands near the 4000 Å break and Hδ, not the full spectral likelihood; additional information elsewhere in the spectra or in broader photometry has not been ruled out.

## What was fitted

The independently compiled C3K_HR/MIST stellar library contains 70 simple stellar populations: 14 ages from 0.005 to 13.5 Gyr and five metallicities from log(Z/Zsolar) = −1.5 to +0.25. Its native resolution is about R = 3000 near 4000 Å. The former coarse C3K_LR library and its results remain unchanged. The new stellar spectra contain no nebular emission. The [build record](../results/calibrated_host_physics/stellar-library.json) identifies the source, compiler, spectral data and isolated environment; the [FSPS source documentation](https://github.com/cconroy20/fsps/blob/05b5e550ddd7ffb1e71102b75cfbfdf057c211fd/SPECTRA/C3K/readme.md) defines the spectral library.

The primary family permits nonnegative mixtures of stellar age, metallicity and attenuation, without a parametric star-formation history or shrinkage prior. Attenuation is exp[−τV(λ/5500 Å)^p], with τV = 0, 0.5, 1 or 2 and p = −0.7 or −1; the zero-attenuation duplicate is removed. Mixtures of differently attenuated populations explicitly permit some stellar mass to contribute little visible light. This finite dust family is an assumption, not a dust measurement or an energy-balance constraint. The IMF is Chabrier and the abundance pattern is solar-scaled.

For every camera, the model is sampled on the native observed wavelength grid, multiplied by the foreground law, passed through the released DESI resolution matrix, and measured with the same band weights and masks as the data. The entire camera resolution operator is applied before the camera selection. Its diagonal-storage convention and row sums are preserved; it is not silently normalized. The data are not already corrected for Milky Way extinction, so F99 attenuation with native DESI E(B−V) × 0.86 is applied to the model. The common amplitude absorbs distance and aperture normalization, without changing relative formed-mass weights among stellar columns.

The primary measurements mask ±400 km/s around all seven declared FSPS nebular-line centres within the five windows. This removes emission-prone regions and some age-sensitive stellar absorption information together. The same mask is applied to the model. Retained widths, rather than imputed line cores, determine the means and covariance. An alternative uses the original full bands with independent nonnegative amplitudes for those emission lines. The gas coefficients contribute zero stellar mass; if gas alone fits, arbitrarily small stellar components leave the entire stellar-age range possible. Neither a Case-B line ratio nor a fixed gas-to-stellar attenuation relation is imposed.

Let the signed band measurements be **y**, their full formal covariance **C**, and the matrix of model columns **A**. The set of compatible nonnegative coefficients **a** satisfies

$$
(A a-y)^T C^{-1}(A a-y)\leq\chi^2_{5,0.95}=11.07050.
$$

No fitted degrees of freedom are subtracted: this is the absolute 95% Gaussian confidence ellipsoid for the five mean fluxes, conditional on the supplied covariance. The overlapping bands share pixels, so their covariance is retained. The extremized age is

$$
\bar t=\frac{\sum_{j\in\mathrm{stars}} a_jt_j}{\sum_{j\in\mathrm{stars}}a_j}.
$$

It is weighted by stellar **formed mass**, not light, surviving mass or supernova progenitor probability. The generalized Charnes–Cooper transformation converts the extrema to second-order cone problems, with zero denominator weight for the gas columns. Original-coordinate feasible witnesses and dual gaps are checked. These sets are not Bayesian age posteriors.

## What changed under the assumptions

| Conditional model | Compatible hosts | Median age-interval width |
|---|---:|---:|
| Primary, no fitted cosmological clock | 55/55 | 13.495 Gyr |
| Original bands and free emission amplitudes | 55/55 | 13.495 Gyr |
| Wider ±800 km/s emission masks | 55/55 | 13.495 Gyr |
| Extra attenuated components, τV up to 4 | 55/55 | 13.495 Gyr |
| Imposed flat ΛCDM clock, H0 = 70 and Ωm = 0.3 | 55/55 | 9.141 Gyr |

All 55 primary upper endpoints lie within 0.01 Gyr of the library ceiling. All 55 clock-constrained upper endpoints likewise lie within 0.01 Gyr of the imposed age of the Universe at that host's redshift. The substantial narrowing in the last row therefore comes from the assumed cosmological clock, not an independent stellar-age measurement. The deliberately generous primary ceiling is not a claim that 13.5 Gyr stars physically exist at these redshifts.

The other declared alternatives vary velocity broadening from 0 to 150 or 300 km/s, apply a continuum tilt of ±2% per 100 rest-frame Å, vary the redder contributing camera by ±3%, change the foreground scale to 1.0, omit the native resolution convolution, or add an assumed 3% independent band discrepancy. All remain compatible for all 55 hosts; the broad age ambiguity persists. Those response and discrepancy choices are conditional sensitivity tests, not empirically calibrated error bounds. A finer grid of 27 ages and nine metallicities, evaluated for eight hosts selected by an ID hash before fitting, changes the primary lower endpoints by at most 0.043 Gyr and leaves the upper endpoints at the same ceiling. [All scenarios](../results/calibrated_host_physics/fits.json), [refinement and numerical checks](../results/calibrated_host_physics/validation.json).

## What this establishes, and what it does not

The recorded uncertainties in these five observed bands cannot support precise stellar formed-mass ages over this flexible family. Assigning a narrow age likelihood would require additional observations or stronger assumptions about star formation, differential attenuation, metallicity, abundance patterns and the cosmological clock. A successful fit is not evidence that those assumptions are true.

This is a fibre-aperture stellar-population calculation. It is not a global host-age measurement or a delay-time/progenitor-age likelihood. Formal diagonal-IVAR errors also omit interpixel and calibration uncertainty. Seventeen hosts span a camera boundary and 30 are multi-exposure coadds affected by the scope of the known DR1 resolution-weighting issue. C3K_HR already has finite spectral resolution; applying the native response adds smoothing rather than recovering an infinitely resolved truth spectrum. The retained response alternatives expose some consequences but do not calibrate these limitations. See the [DESI DR1 release documentation](https://data.desi.lbl.gov/doc/releases/dr1/).

No supernova brightness, age–brightness coefficient, redshift-dependent luminosity correction or cosmological parameter was fitted here. In particular, these five-band compatibility sets must not be inserted as a calibrated progenitor-age prior into the SN + BAO + CMB fit. A selected host population and a corrected supernova distance do not become independent observations by sharing an identifier.

## Numerical evidence and reproduction

The calculation completed 770 primary/sensitivity cases plus 16 refined-grid cases, with no unresolved numerical failures. All 1,572 returned endpoints pass original-coordinate constraint checks. Independent quadratic age-profile calculations on eight reduced physical sublibraries agree with all 16 cone endpoints within 1.1 × 10⁻⁶ Gyr. Forty synthetic mixtures, explicit zero-flux and gas-only tests recover the compatible age support. These test numerical correctness, not the empirical coverage of incomplete calibration errors.

An independent implementation rebuilds the native resolution response, masked fluxes and covariance for all 55 hosts and three masking choices. Direct FITS diagonal scattering agrees to 8.2 × 10⁻¹⁶ relative; masked flux/covariance reconstruction agrees to 5.5 × 10⁻¹³ and synthetic forward-band predictions to 2.7 × 10⁻¹⁴. [Independent review](../results/calibrated_host_physics/independent-review.json).

The initial resolution-offset programming error aborted before fitting. Two generic SLSQP cross-check formulations also failed numerically; their records are retained. The final independent check uses convex quadratic age profiles and reduces boundary faces exactly. No scientific acceptance gate was relaxed, and the production cone calculations did not change in response to those cross-check failures.

Use the [reproduction commands and module interface](../code/calibrated_host_physics/README.md). The [initial design](../code/calibrated_host_physics/design.json) and [implementation specification](../code/calibrated_host_physics/implementation-design.json) were frozen before the new stellar fits, after the existing spectral observations had already been inspected. This was an exploratory physical feasibility test, not a blinded or preregistered discovery test.
