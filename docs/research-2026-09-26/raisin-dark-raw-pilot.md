# Four public WFC3/IR RAW ramp pilot: metadata and read-sequence closure

This pilot acquired only the two frozen nearest verified BLANK-filter DARK roots, `idbx43p7q` (search visit) and `idp247tnq` (template visit), and the two frozen held-out F160W science roots, `icxoi1bcq` and `icxoi4hgq`. Their exact public [MAST Download](https://mast.stsci.edu/api/v0/_services.html) URIs and HEAD lengths were verified before the GETs and frozen in `runs/research_2026_09_26/raisin_hst_pixel_feasibility/dark_ramp_raw_pilot/protocol.json` (SHA-256 `c12defe7f5d5421ba8942d5f0d937de2e4d1af606aadd261f03d869e57382d24`). The four RAW files total 118,897,920 bytes, within the 140-MB/200-s cap; acquisition took 8.82 s with one attempt per file. Exact file SHA-256 values, MAST product sizes, response records, and byte-level provenance are in `acquisition-results.json`. No other darks were acquired.

All four FITS files parse without structural warnings. The two 16-read dark RAWs contain 81 HDUs each; the two eight-read science RAWs contain 41 each. Every read has the expected SCI/ERR/DQ/SAMP/TIME group. The final HDU ends at the exact downloaded byte count. No HDU carries a FITS `CHECKSUM` or `DATASUM` card, so structure plus local SHA-256 and HTTP length are the available integrity evidence; there is no independent archive cryptographic digest in this pilot. Inspection accessed FITS headers and checksum metadata, **not SCI/ERR/DQ pixel values or statistics**.

The chronological first eight dark SCI read headers match the respective science eight **exactly at printed header precision** for `SAMPNUM`, `SAMPTIME`, `DELTATIM`, constant SAMP `PIXVALUE`, constant TIME `PIXVALUE`, and SCI 1024×1024 geometry. The shared times in seconds are 0, 2.932291, 52.932701, 102.933113, 152.933533, 202.933945, 252.934357, and 302.934753. `SAMPZERO=2.911756` in all four primary headers. Dark FITS extension versions run in reverse chronology (16→9 for the first eight); science versions run 8→1. Sorting by `SAMPTIME` avoids aligning the wrong groups. TIME and SAMP extensions in these RAWs are constant-header HDUs (`NAXIS=0`), not image arrays; their `PIXVALUE` metadata closes with the SCI headers. The darks have eight additional later reads through 702.938171 s, which **must not** enter a future first-eight processing or mask decision.

Root's independent header parser confirms both prefix matches, complete group
membership and common gain/overscan/nonlinearity reference names. Its results
and four-file/source hashes are in `dark_ramp_raw_pilot/root-prefix-review/`.
It does not access an HDU data array.

The instrument headers agree on IR detector, full frame, `SPARS50`, MULTIACCUM, aperture IR, `CCDAMP=ABCD`, `CCDGAIN=2.5`, `CCDTAB=iref$t2c16200i_ccd.fits`, `NLINFILE=iref$a2412448i_lin.fits`, and listed `ZOFFCORR`/`BLEVCORR`/`NLINCORR`/`DQICORR`/`CRCORR`/`UNITCORR` states. RAW `ATODGN*` and `READNSE*` are **zero placeholders**, not evidence for zero physical gain/noise; actual calibration must resolve `CCDTAB` and read processing. Search dark and science share `BPIXTAB=iref$3562029fi_bpx.fits`; template dark uses `iref$3562028ni_bpx.fits`, while template science uses `iref$3562018mi_bpx.fits`. Dark `DARKCORR` and `FLATCORR` are `OMIT` and `DARKFILE=N/A`, versus `PERFORM` and named dark reference in science. The search dark's current `CRDS_CTX` differs from its science RAW (`hst_1337` versus `hst_1339`); template pair uses `hst_1337`. Raw `CAL_VER` is blank, so current product processing equivalence does not follow from these headers.

Two additional flags constrain interpretation. The nearest template dark `idp247tnq` has RAW `EXPFLAG=INDETERMINATE`, whereas the other three say `NORMAL`; its meaning and effect have not been resolved from source, and it is not silently replaced. Its DQ groups are constant-header HDUs, while template science `icxoi4hgq` has eight full 1024×1024 DQ image HDUs, explaining that science RAW's larger byte count. No DQ values were read. A later protocol must decide how bad-pixel/mask maps and exposure flags are handled before any pixel-level covariance score.

The zero-read signal step also differs: `ZSIGCORR=OMIT` in both darks and
`PERFORM` in both science exposures. The pinned source applies a thresholded
super-zero-reference correction before subsequent nonlinearity processing.
Gain multiplication and matching read times alone do not reproduce that
measurement operator. The common overscan reference is
`iref$q911321mi_osc.fits`; its reference-pixel definition still needs matching
to the source and table values.

The [archive handbook](https://archive.stsci.edu/manuals/archive_handbook/appendix2.html)
defines `INDETERMINATE` generically as an exposure time that could not be
derived successfully from telemetry with no predicted duration available.
It also says interpretation is instrument dependent. This supports a quality
provenance check, not an automatic claim that this dark is corrupt or a
decision to replace it. Its cause in this particular exposure remains open;
matching header sample times is not independent telemetry verification.

This is **read-cadence metadata closure for one pair per visit**, not a noise validation. A future bounded test would need common reference-pixel, zero-read, nonlinearity, gain and mask handling on the frozen first eight RAW groups; it must account for changing dark current, sky illumination, persistence and detector state over the roughly 11–12 days separating the darks from science. Paired RAW differences can cancel a static mean dark pattern, but calibrated IMA/FLT dark residuals may reuse a master DARKFILE partly trained on these ramps and do not independently validate that master-dark uncertainty. No remaining candidate acquisition or covariance scoring is released by this pilot.

Machine-readable details: `runs/research_2026_09_26/raisin_hst_pixel_feasibility/dark_ramp_raw_pilot/{verified-products.json,acquisition-results.json,header-inspection.json,header-comparison.json,state-ledger.json,manifest.json}`. All header comparisons were made without an image-array read.
