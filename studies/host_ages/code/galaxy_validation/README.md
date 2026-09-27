# Independent galaxy spectroscopy

This experiment tests SED age summaries against observed stellar absorption and emission features, checks central-versus-supernova-site apertures, and asks whether the available spectral information adds held-out brightness prediction. Read the [scientific findings](../../notes/galaxy-validation-results.md) before interpreting any age or dust coefficient.

Use the repository's locked `.venv`. Commands below run from its root, use at most three download workers and two numerical threads, and write large data only to ignored `.work/galaxy-validation/`.

```bash
.venv/bin/python studies/host_ages/code/galaxy_validation/acquire.py
.venv/bin/python studies/host_ages/code/galaxy_validation/restore_inputs.py
.venv/bin/python studies/host_ages/code/galaxy_validation/match.py --archive .work/galaxy-validation/prerequisites
.venv/bin/python studies/host_ages/code/galaxy_validation/query_sdss.py
.venv/bin/python studies/host_ages/code/galaxy_validation/spectra.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python studies/host_ages/code/galaxy_validation/diagnostics.py --archive .work/galaxy-validation/prerequisites
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python studies/host_ages/code/galaxy_validation/manga_local.py
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python studies/host_ages/code/galaxy_validation/historical_ages.py --archive .work/galaxy-validation/prerequisites
```

`restore_inputs.py` restores eight inherited inputs and verifies `prerequisite-lock.json`. An old local archive is an optional read-only cache, configurable with `--archive-cache PATH`. Without it, the public ZTF ZIP is about 1.37 GB; only its three explicit tables are extracted. G11 and the pinned ZTF master list are downloaded directly. The repository's main acquisition must already supply `data/ages/{table1,table2}.dat` and `data/distances/Pantheon+SH0ES.dat`; see [data sources](../../../../docs/methods/data.md). No generated result from a previous age experiment is needed for the primary spectrum/brightness analyses. The optional existing DustPedia crosswalk contributes only its separately reported coordinate-availability check.

`acquire.py` downloads the pinned TITAN snapshot, original R19 host coordinates and age table, MPA-JHU position catalogue, two FIREFLY catalogues and method descriptions. `query_sdss.py` fetches only exact plate–MJD–fibre matches from the public SDSS SQL endpoint, storing each query and response. `spectra.py` downloads all returned reduced spectra. `manga_local.py` downloads the matched MaNGA maps. `source-lock.json` refuses changed downloaded bytes; a later catalogue release requires an explicit new edition, rather than silently changing the input.

Primary design, detailed analysis, local-aperture extension and original-versus-updated age comparison have separate timestamped design files. The last comparison was chosen after primary results were inspected. The spectral-sample selection summary is also descriptive post-inspection analysis.

## Interfaces and definitions

- `host-spectrum-age-ledger.csv`: one accepted age-catalogue/SN association per row, exact spectrum IDs, physical host key, host/SN coordinates, redshift, SED summaries, brightness selection, native photometry, emission/absorption features and independently integrated narrow bands. The same galaxy may occur in G11 and R19; do not count these as independent observations.
- `sdss-spectral-measurements.csv`: one returned spectrum per row; native SDSS `modelFlux`, `cModelFlux` and `fiberFlux` are nanomaggies and their `*Ivar` values are inverse squared flux. These have not been changed to AB fluxes or foreground-corrected here.
- `remeasured-dn4000.csv`: raw and unscaled-SFD/O'Donnell-corrected narrow-band mean observed Fν in microJy, variances in microJy², coverage and ratio. Bands are 3850–3950 and 4000–4100 Å in the rest frame. Pixel widths, including fractional boundary pixels, weight Fν, with Fν = Fλ λ_obs²/c. Zero cross-band covariance means the supplied diagonal pixel-ivar approximation, not a measured absence of calibration/resampling covariance.
- `brightness-cohort.csv` and `brightness-heldout.csv`: the common selected cohort and all held-out predictions/folds.
- `manga-central-SN-indices.csv`: actual nearest spaxel at each specified sky position, mask/coverage failures, seeing, dispersion correction and local/central indices. Adjacent spaxels are not treated as independent.
- `firefly-paired-models.csv`: identical host/galaxy identities across the MILES and MaStar models. Catalogue age fields are log10(age/Gyr); missing sentinels are excluded before exponentiation or differencing.

**CSV consumers must preserve catalogue IDs:**

```python
pd.read_csv(path, dtype={
    "specObjID": "string", "bestObjID": "string",
    "photometric_objID": "string", "photometric_flags": "string",
})
```

Default pandas inference turns nullable 64-bit IDs into inexact floating values. The released SQL bytes and our saved files retain their exact digits; never recover an identifier by rounding a float.

The brightness comparison imports the existing first-party environmental likelihood. It reconstructs the sample from the original released light-curve/host tables and TITAN summaries, preserving the SALT covariance and slope-dependent variance. It does not infer a latent age likelihood, refit detector photometry or regenerate survey selection. Conditioning on dust, colour, mass and metallicity is an incremental prediction test and can remove mediated age information.

## Numerical validation

```bash
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python studies/host_ages/code/galaxy_validation/audit_legacy_nnls.py --archive PATH_TO_FROZEN_RESEARCH_INPUTS
OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 .venv/bin/python studies/host_ages/code/galaxy_validation/validate.py
```

The optional legacy audit reads old-layout DESI BAO and RAISIN release inputs. It redirects the current first-party sources' historical root paths in memory, without editing the sources or old outputs. Its compact result is preserved here so `validate.py` can check the current result inventory without rerunning historical work. Main validation checks locked download identities, original SQL identifier strings, host association, folds, Fν equations, rank correlations, actual likelihood evaluation and a separate optimizer. It does not validate missing calibration or selection models. All inherited source files are left intact.
