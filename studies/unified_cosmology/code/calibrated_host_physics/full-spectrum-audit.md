# What additional DESI spectral information is usable?

All 55 matched hosts have at least 95% valid native coverage over rest-frame 3700–4900 Å. Their median per-pixel signal-to-noise ratios range from 0.65 to 5.02. Additional broad regions cover Hγ/G4300 for all 55, Hβ for 54, Mg b for 52 and Fe5270/5335 for 48. These are coverage counts, not stellar absorption detections or additional age measurements. No stellar-age fit was performed in this audit.

The available C3K SSP library has an approximate optical velocity resolution of 42.44 km/s (Gaussian sigma). About 96.1% of valid common-region native DESI response columns have a smaller signed second moment. Those moments are diagnostics: the released response includes negative wings, and one column has an undefined signed variance. Convolving the already broadened SSPs with DESI's response alone therefore does not give a generally valid full-resolution stellar model. The native response must not be renormalized or mistaken for a noise covariance.

## A numerical remedy that partly works

Let `R` be the released DESI response, `D` the host/foreground attenuation and redshift density factor, and `L ≈ G f` the available library spectrum. We tested broad measurements of the transformed data `W G y` against the model `W R D L`. The full induced covariance is

$$
C'=WG\,\mathrm{diag}(1/\mathrm{IVAR})\,G^T W^T.
$$

`G` explicitly convolves wavelength times flux density in log wavelength. No Gaussian replacement of `R`, deconvolution, fitted age, or minimum stellar velocity dispersion is used. The fixed masks exclude every available FSPS nebular centre within 400 km/s and erode the output support wherever the blur touches an invalid or masked input. The remaining means are signed native flux densities, not Lick indices or the earlier mixed Fν bands.

For all 70 existing SSPs and seven attenuation settings per host, the ordering difference `W (G R D − R D G) f` passes the preregistered 0.1-sigma per-band and joint-whitened gates. The worst errors are **0.02566 sigma per band and 0.04781 jointly**. The independent covariance reconstruction agrees to relative 4.75×10⁻¹⁶. Cross-band correlation reaches **0.812**; ignoring it would materially change the calculation.

The separate unresolved-line test does **not** pass for every host. An artificial line of one observed-Ångstrom equivalent width at any native pixel yields a maximum joint error of **0.12189** for target `39627682694039936`; 54 of 55 hosts pass the same 0.1 threshold. This is a numerical stress, not a measured line residual or an empirical line-strength bound. The finite resolved-SSP success does not certify arbitrary finer spectral structure or establish the true C3K line-spread function.

## The additional smoothing trial failed

After recording that outcome, we declared one fixed extra 150 km/s blur `H` for every host. This changes the data to `W H G y`, model to `W H R D L`, and propagates the corresponding full covariance. It is data degradation, not a stellar-dispersion prior. Windows and thresholds were unchanged, and masks were eroded over the complete `H G` support.

All 55 resolved-family tests still pass, with worst joint error **0.04940**. However, the unresolved-line test now passes only **52 of 55**, with maximum **0.12233**. The extra blur also removes masked absorption support: for the 51 hosts retaining the Hδ contrast, its resolved-template contrast range divided by formal noise is a median **0.349 times** its original value. The analogous Hγ/G4300 comparison is **1.073 times** for 54 hosts. These are template/operator sensitivity measures, not data-derived age information. This repair is **not adopted**, and no further smoothing was selected from these outcomes.

## The next defensible inference

A broad-feature compatibility calculation remains feasible under explicitly stated library-response assumptions. It should retain the full correlated data vector, nonnegative stellar mixtures and a union over attenuation, calibration and kinematic nuisance settings. It must report absolute model adequacy alongside formed-mass age sets. Adding observations changes the dimension of the compatibility region, so age-interval contraction cannot be assumed in advance. To measure incremental information from the same photons, the old measurement span and its covariance-orthogonal complement can be compared with separately preregistered coverage gates; the old and new bands must not be treated as independent data.

For a stronger full-spectrum claim, one of the following response remedies is needed before fitting:

- Obtain and validate an actually sharper SSP grid with compatible age support and mass conventions. The locally available C3K grid is not that grid. Public XSL models offer higher resolution but different stellar modelling and a minimum age of about 50 Myr; silently removing younger populations would impose new support restrictions.
- Model physical stellar dispersion only where it exceeds the library width, adding the *difference in quadrature* before the native response. Any excluded lower-dispersion populations must remain explicit.
- Validate a covariance-preserving data transformation against an independently known library kernel or a justified bound on unresolved spectral structure. The present finite-family audit and failed extra-blur trial do not provide that bound.

Thirty hosts have multiple exposures. DESI DR1 documents a coadded-response weighting issue, so per-exposure spectra or a corrected coaddition provide a concrete additional check. Diagonal pixel IVAR also omits calibration and resampling correlations. Camera-relative response, emission infill, nebular continuum, stellar abundance patterns and aperture differences remain relevant. Increasing the number of pixels does not solve these model deficiencies. Optical spectra also cannot exclude arbitrarily dust-hidden stellar mass without an attenuation restriction, and fibre-population ages are not supernova progenitor clocks.

See the [DESI DR1 known issues](https://data.desi.lbl.gov/doc/releases/dr1/known-issues/), [coadd data model](https://desidatamodel.readthedocs.io/en/latest/DESI_SPECTRO_REDUX/SPECPROD/healpix/SURVEY/PROGRAM/PIXGROUP/PIXNUM/coadd-SURVEY-PROGRAM-PIXNUM.html), [pinned C3K library description](https://github.com/cconroy20/fsps/blob/05b5e550ddd7ffb1e71102b75cfbfdf057c211fd/SPECTRA/C3K/readme.md), and [XSL model release](https://xsl.u-strasbg.fr/page_ssp_all.html).

## Reproduction

Run from the repository root after restoring the existing calibrated-host spectra and SSP library:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
.work/unified-cosmology/calibrated-host-physics/.venv/bin/python \
studies/unified_cosmology/code/calibrated_host_physics/spectral_resolution_audit.py

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
.work/unified-cosmology/calibrated-host-physics/.venv/bin/python \
studies/unified_cosmology/code/calibrated_host_physics/degraded_band_audit.py
```

The two pre-outcome designs and compact results bind their sources and inputs by SHA-256. Full host vectors and covariance matrices are regenerated under `.work/unified-cosmology/full-spectrum-design/`. The second audit reproduces the first operator's errors exactly before testing its single declared modification. Neither script changes the existing stellar-age results or any cosmological likelihood.
