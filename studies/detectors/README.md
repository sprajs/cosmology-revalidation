# Detector repeats and calibration

HST science-repeat variance ratios are 0.708 and 0.729 under the chosen masked aperture operator. The propagated flat-reference contribution is about 0.005%, too small to explain the deficit. A dark-pixel ratio of 4.202 and aperture ratio of 0.547 weight the detector differently: one pixel supplies 82.9% of the squared pixel signal and no aperture signal. Raw moment ratios mix state, events and noise; they do not identify a universal read-noise multiplier.

## Code and scope

MAST/CRDS acquisition, FITS header lineage, aperture geometry, raw-read extraction, first-eight-read prefix construction, native CALWF3 build/run preparation, and pixel/aperture influence checks.

The [source directory](code/) retains the original calculations. See the [execution guide](../README.md#execution) before running them: historical imports and relative paths are assembled into an isolated workspace, and some studies require additional data or native software. Source preservation does not certify that an incomplete experiment succeeded.

## Supporting findings

- [HST repeat-noise result: accounting for the quoted variance](notes/hst-ramp-variance-accounting.md)
- [Pinned CALWF3 RAW-to-FLT closure feasibility](notes/raisin-calwf3-native-feasibility.md)
- [Exact science and error replay; downstream quality flags remain separate](notes/raisin-calwf3-native-replay.md)
- [Quality and reference audit of the 14 public dark-ramp candidates](notes/raisin-dark-quality-reference-audit.md)
- [Public WFC3/IR dark-ramp metadata near the DES16E1dcx visits](notes/raisin-dark-ramp-metadata.md)
- [Four public WFC3/IR RAW ramp pilot: metadata and read-sequence closure](notes/raisin-dark-raw-pilot.md)
- [Public HST image-identification feasibility for DES16E1dcx](notes/raisin-hst-pixel-feasibility.md)
- [Current HST F160W signed-pixel repeat-noise check](notes/raisin-hst-signed-pixel-noise.md)
