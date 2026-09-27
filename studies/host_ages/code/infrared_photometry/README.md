# Infrared catalogue associations for DES hosts

This acquisition and positional analysis uses the existing 265 selected C3/X3
hosts with three valid SPIRE measurements. It retrieves official SWIRE IRAC/MIPS24
and SERVS two-band catalogue records inside 75 arcsec of each host. It does not
infer a host dust luminosity or update any stellar age.

Run from the repository root with Python, NumPy, Pandas, SciPy and PyArrow:

```sh
.work/host-transport/.venv/bin/python studies/host_ages/code/infrared_photometry/acquire.py
.work/host-transport/.venv/bin/python studies/host_ages/code/infrared_photometry/analyze.py
.work/host-transport/.venv/bin/python studies/host_ages/code/infrared_photometry/matching_controls.py
```

The first command requires the earlier host-transport acquisition and infrared
analysis products, including `des-deep-crosswalk.csv`, `des-spire-photometry.csv`
and the two DES deep-field catalogues under `.work/host-transport/`. Their exact
hashes are recorded. Downloads and full output ledgers remain ignored under
`.work/infrared-photometry/`; compact results and first-party code are retained.
Subsequent acquisition verifies pinned downloaded bytes and queries rather than
silently accepting an altered upstream response.

`candidates.csv` retains the union of native catalogue source positions and exact
string identifiers. `host-candidates.csv` records all candidates within 75 arcsec
per host. `associations.csv` has every host and its three available catalogue
families, including missing and ambiguous associations. `photometry.csv` retains
native aperture-2 measurements, sentinels, errors and quality flags alongside
explicit validity and association status. Fluxes are in microJy. The reported
values are native observer-frame measurements; this branch adds no Galactic
extinction correction, K-correction or extra aperture correction. The reported
`clean_photometry` flag is an intentionally strict all-flags-clear diagnostic,
not the catalogue's official selection or a calibrated measurement likelihood.

Missing catalogue entries are neither zero flux nor numerical upper limits.
Catalogue fluxes are point-source aperture corrected and can be inappropriate
for extended galaxies. CDFS SWIRE has conflicting aperture-radius metadata in
its native error labels; no corrected error is invented. SWIRE/SERVS agreement
is a release-consistency check with shared IRAC instrumentation and explicit
calibration to SWIRE; individual exposure overlap was not audited. This does
not independently validate the flux calibration. Bandmerged and separate MIPS24
records can identify the same detection through `detid_24`.

The frozen design and acquisition-only amendment precede matching outcomes.
Root's separately designed positional controls use shifted coordinates and
independent spherical arithmetic; they estimate conditional coincidence rates,
not posterior host probabilities or catalogue completeness.
