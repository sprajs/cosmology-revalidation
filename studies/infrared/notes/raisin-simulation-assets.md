# Source assets behind the RAISIN simulation corrections

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

Seven additional files have been acquired from the pinned author repository
`djones1040/RAISIN_cosmo@b888214a5cbae38ac0bf886488ce7734f5b77e87`.
Every file matches its Git-tree blob and a recorded SHA-256 hash. The earlier
selection review correctly described these as absent from its inspected local
files; that acquisition gap is now closed for the named DES error/efficiency
files. Historical execution identity remains open.

The [acquisition record](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_simulation_assets/acquisition.json)
contains exact upstream URLs, including the
[nominal DES input directory](https://github.com/djones1040/RAISIN_cosmo/tree/b888214a5cbae38ac0bf886488ce7734f5b77e87/sim/inputs/DES).
The [content check](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_simulation_assets/result.json)
and [executable](../code/raisin_simulation_assets.py)
record the following mappings without generating or fitting any supernovae.

| Asset | Verified numerical/source content | Interpretation |
|---|---|---|
| `DES3YR_SIM_ERRORFUDGES.DAT` | Eight band/field maps, 71 rows of host surface brightness versus error scale; tabulated scale 1.0468–5.1458 | The historical reader defaults to applying these maps to both generated scatter and reported simulated error. This is not an inferred multiplier to apply again to observed errors |
| `SEARCHEFF_PIPELINE_RAISIN.DAT` | Each griz table exactly shares 32 numerical rows with the ordinary DES file, plus a unity-efficiency endpoint at SNR=1,000,000 | Epoch detection probabilities based on injected-source recovery; detection is distinct from retaining every measured positive flux |
| `SEARCHEFF_SPEC_DES_Moller_G10_v7.DAT` | Nominal file has 99 peak-i-magnitude/efficiency rows; its 79 shared magnitude rows exactly match the flat-distribution directory version | This file models a brightness-dependent spectroscopic selection; filenames alone do not certify every historical override or HST selection |
| `OIR.J19/OIR.INFO` | Author file is byte-identical to bundled 2024 SNDATA_ROOT file, SHA-256 `f68d654809a6b5423566fc6926aa53d7837fc27805b7eec17752622490031d9e` | Seven correlated spectral-node scatter amplitudes, trained from CSP fits with empirical optical shape/extinction assumptions |
| Bundled `SEARCHEFF_PIPELINE_LOGIC.DAT` | DES requires two distinct detection epochs in any griz band; generator sets `NEWMJD_DIF=0.4` days | Available local logic, not a verified copy of the historical SNDATA_ROOT used in the published run |

The historical `sntools_fluxErrModels.c` initializes map application to both
simulation and reported data errors at lines 166–179, then calculates these
separately at lines 803–815. Additional execution overrides must be retained.
No covariance between epochs is specified by the newly acquired error map.
That absence is a property of this file, not proof that extraction errors are
independent.

The historical `sntools_genSmear.c:2792–2995` constructs the OIR spectral-node
covariance from the seven sigmas and correlation matrix, adds `1e-9` to the
diagonal before applying a covariance factor of 1.3, draws a common spectral
perturbation and interpolates it in wavelength. The tabulated correlation
matrix is positive definite, with smallest eigenvalue 0.009986. Its seven
sigmas are 0.125, 0.118, 0.130, 0.124, 0.116, 0.121 and 0.136 mag; they are
not independent per-epoch measurement errors. The node covariance is not the
final filter-integrated covariance. The historical parser reads
`COLOR_SIGMA_SCALE` but does not use it in the OIR covariance construction;
the nominal file value is one, so this observation changes no nominal number.
An initial descriptive claim that the scale was applied was caught by source
review and corrected; its source/result snapshots are preserved.

An [independent raw-file/source review](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_simulation_assets/independent-review/review.md)
confirms all seven Git blobs, 13 recorded pre-review file hashes and the
numerical/source mappings. The released `model/OIR.J22/OIR.INFO` is also
byte-identical to the author and bundled `OIR.J19/OIR.INFO` files. Different
directory labels therefore do not imply different numerical models here.

These inputs make a closer source simulation feasible, but do not close the
joint cadence/extraction/selection distribution. The existing SIMLIB inherits
observed positive-epoch membership, while fresh simulated flux noise does not
repeat that selection. Original exposure metadata for omitted rows, the
historical execution manifest, optical-to-NIR timing, population assumptions
and event/follow-up selection still require explicit treatment. No distance
correction or cosmology has been inferred from this source acquisition.

## Available simulation timing path

The separately acquired and Git-blob-verified
[simulation NIR fit file](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_simulation_assets/timing/REFAC_DES_RAISIN_NIR.nml)
fixes shape=1, AV=0 and peak time while fitting JH. The nominal DES generator
sets `GENSIGMA_SEARCH_PEAKMJD=0.01` day. In the v11_03c generator this option
disables the noisy-flux peak estimator and writes a peak equal to the
generating true peak plus a Gaussian perturbation of that width. The FITS
writer and reader map it through the `PEAKMJD` field into the fixed NIR fit.
The [source-chain record](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_simulation_assets/timing/source-chain.json)
hashes inputs and identifies exact source lines. This is an available-source
configuration result, not a reproduction of the published simulation run.

The release README describes real-data header peaks as originating in
optical+NIR SNooPy fits; the exact upstream executions and timing uncertainty
remain unlinked. A 0.01-day perturbation about generating truth is not a
reproduction of that estimator. The real and simulated header-time processes
therefore require an explicit closure test before claiming their timing
effects are already represented by the bias simulation. Differences between
the author timing table and real-data headers are not known true errors.
Some timing dispersion may also enter empirical intrinsic-scatter estimates,
so adding another variance term without checking overlap would be unjustified.

## Subsequent field-applicability audit (2026-09-26)

The map defaults above apply only when a field matches. A separate exact
v11_04d [source/configuration review](raisin-prospective-noise-review.md)
finds that none of the 120 band/field pairs in the 23-entry, 6,877-row
archived DES RAISIN SIMLIB matches the eight maps: actual fields are full
object names such as `DES16E1dcx`, while map groups contain short tokens
such as `E1`. Native fallback returns unchanged true and reported errors;
this is a no-op, not a fatal error. The complete hashed
[matching ledger](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/raisin_timing_assets/prospective_noise_review/map-matching.csv)
is retained. This does not establish every historical execution or a
measured noise deficit. The pinned SIMLIB constructor imports optical
sky/PSF/zero-point metadata from a private survey library and estimates NIR
sky from private images, rather than directly deriving them from the quoted
SN errors. Their original execution and upstream correction content remain
unlinked. Renaming fields to activate the maps would therefore be a new
intervention with a possible double-counting risk, not an established
correction to apply automatically.
