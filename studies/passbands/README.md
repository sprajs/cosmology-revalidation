# Passbands and photometric lineage

CSP lineage recovers 5,491 raw/intermediate rows. The three listed objects absent from the source archive are outside the 79-object cosmology sample. Different J-band representations produce phase-dependent changes reaching about 0.069 mag for the tested spectra. This is a spectral sensitivity, not an identified distance correction. Native filter-response and timing studies retain their optimizer, covariance and model-domain limitations.

## Code and scope

Filter acquisition, magnitude lineage, photon integrals, spectral contrasts, calibration bridges, native filter perturbations and first-party instrumentation patches.

The [source directory](code/) retains the original calculations. See the [execution guide](../README.md#execution) before running them: historical imports and relative paths are assembled into an isolated workspace, and some studies require additional data or native software. Source preservation does not certify that an incomplete experiment succeeded.

## Supporting findings

- [CSP native covariance amplitude check](notes/csp-covariance-independent.md)
- [Direct cross-instrument calibration: metadata feasibility](notes/csp-direct-overlap-feasibility.md)
- [Original CSP DR3 products and calibration errata](notes/csp-dr3-provenance.md)
- [CSP J-filter intervention and existing RAISIN bias correction](notes/csp-filter-bias-propagation.md)
- [CSP J-filter interpretation: independent source review](notes/csp-filter-interpretation-review.md)
- [CSP native baseline: reproduction and convergence are separate](notes/csp-native-baseline-review.md)
- [Optical V labels follow the physical filter dates](notes/csp-optical-filter-identification.md)
- [CSP infrared passband identity and conditional spectral response](notes/csp-passband-lineage.md)
- [CSP 42-object peak-header provenance audit](notes/csp-peak-header-provenance.md)
- [CSP-II spectra identify passband-shape effects, not a distance correction](notes/csp-spectral-identification.md)
- [CSP J-filter processing response](notes/csp-stable-filter-response.md)
- [Physical WIRC J observer-operator feasibility](notes/csp-wirc-operator-feasibility.md)
- [CSP WIRC-J reference bridge: pinned SNooPy VegaB convention](notes/csp-wirc-operator-vega-bridge.md)
