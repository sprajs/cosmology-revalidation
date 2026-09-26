# Exact science and error replay; downstream quality flags remain separate

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The private CALWF3 replay reproduces **every SCI, ERR, SAMP and TIME pixel
exactly** in the two selected current archive exposures, `icxoi1bcq` and
`icxoi4hgq`. Each array contains 1,028,196 pixels. The data-quality array differs
only by additional archive bit 4096 flags. The original whole-product equality
gate therefore remains failed; the result establishes exact numerical closure
of these four CALWF3 output components, not of AstroDrizzle or the complete
archive product. The archive FLTs also contain 13 extensions versus five in
the native outputs, with downstream astrometric/headerlet metadata. The
DQ-only statement applies to the five compared image components, not to all
file bytes or header cards.

| Comparison | Search exposure | Template exposure |
|---|---:|---:|
| SCI/ERR/SAMP/TIME unequal values | 0 / 0 / 0 / 0 | 0 / 0 / 0 / 0 |
| Maximum SCI or ERR floating-point difference | 0 | 0 |
| DQ pixels differing | 3,897 | 4,262 |
| Sole XOR flag in differing DQ pixels | 4096 | 4096 |
| Pixels with native DQ=0 additionally rejected by archive DQ | 186 | 213 |
| Other DQ bits different | 0 | 0 |
| Calibration reference names, switches and units different | 0 | 0 |

The native output has no bit 4096 pixels in either exposure; the archive only
adds this bit and removes no native bits. The WFC3 handbook identifies 4096 as
the [AstroDrizzle cosmic-ray flag](https://hst-docs.stsci.edu/wfc3dhb/chapter-2-wfc3-data-structure/2-2-wfc3-file-structure).
The [AstroDrizzle instructions](https://hst-docs.stsci.edu/drizzpac/chapter-6-reprocessing-with-the-drizzlepac-package/6-3-running-astrodrizzle)
describe its insertion into input FLT quality arrays when combining an
association. This explains the processing-stage distinction, without claiming
to have rerun its rejection decisions. The published DQ arrays and all earlier
aperture results remain unchanged. Most additional flagged pixels already had
another native quality flag, which is why their effect on a DQ=0 mask is smaller
than the total differing-pixel count.

## Provenance and computation

The source is HSTCAL commit `6a1147d7bdccb7e2a7b73af276f4597c6420fbb8`,
whose 812 regular files were independently matched to their Git blob IDs.
CALWF3 reports version 3.7.3, matching both archive FLTs. Fourteen exact named
references are staged through immutable symlinks; the eight previously missing
files were acquired from official CRDS, checked for byte size and FITS structure,
and hashed. Source, CMake wheel and reference acquisition totaled 716,501,580
bytes in 31.93 seconds. No historical 2016–2017 reduction is reproduced here.

The initial configuration failed because a nested library requires Fortran.
That attempt is retained. The subsequent private build reused the existing
project Fortran compiler and did not install system packages or alter existing
environments. It used GCC with `-O2`, disabled OpenMP and required about 9.35
seconds for configuration and compilation. Generic Git-version discovery was
disabled to avoid attributing the unrelated enclosing supernova repository to
HSTCAL; its generic provenance string is `unknown`, while the independently
verified source commit and unchanged WFC3 version header remain recorded.
The binary SHA-256 is
`230dc5e2c760217e7924f0c11153278dd43e10c22fa73828cad8aff0a9c97d00`.

Two sequential calibration calls consumed 4.04 seconds. A protocol written
before execution required complete array equality, retained all pixels and
specified an eight-ULP diagnostic separately from exact reproduction. All
floating arrays are bit-identical; neither primary nor approximate whole-product
closure passes because DQ equality is required. No threshold was changed.

Artifacts are under
`runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/`,
including acquisition/build protocols, preserved failures, `replay-protocol.json`,
`replay-result.json` and per-exposure native logs. The independent review in `independent-review/`
verifies all 812 source blobs, 26 frozen replay inputs, 14 reference links and
the complete pixel comparisons; it confirms the component result without
promoting it to whole-product reproduction. The executor is
[run_calwf3_closure.py](../code/run_calwf3_closure.py).

## What this resolves

The exact ERR replay connects the audited ramp/flat equations to actual current
archive error arrays for these two exposures. Agreement of implementations does
not validate their physical covariance model. The dark-ramp experiment still
needs to distinguish total repeat variability, temporal/spatial covariance,
cosmic rays and calibration-state differences. Current FLT error construction
also does not establish the errors or count-rate correction in the missing
historical RAISIN aperture catalog. A luminosity or cosmological correction
cannot be obtained merely by rescaling these reproduced errors.
