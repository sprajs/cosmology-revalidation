# Calibrated DESI spectra of observed supernova hosts

This package recovers matched DESI DR1 host spectra for the existing OzDES/Dovekie cohort. It provides signed spectral data and covariance-aware diagnostics; it does not convert a break index to progenitor age or refit cosmology. See the [results note](../../notes/calibrated-hosts.md).

From the repository root, after the survey-selection OzDES crosswalk has been regenerated:

```bash
uv venv .work/unified-cosmology/calibrated-hosts/.venv --python 3.12
uv pip install --python .work/unified-cosmology/calibrated-hosts/.venv/bin/python -r studies/unified_cosmology/code/calibrated_hosts/requirements.lock.txt
.work/unified-cosmology/calibrated-hosts/.venv/bin/python studies/unified_cosmology/code/calibrated_hosts/sources.py
OPENBLAS_NUM_THREADS=1 .work/unified-cosmology/calibrated-hosts/.venv/bin/python studies/unified_cosmology/code/calibrated_hosts/acquire.py
OPENBLAS_NUM_THREADS=1 .work/unified-cosmology/calibrated-hosts/.venv/bin/python studies/unified_cosmology/code/calibrated_hosts/spectra.py
OPENBLAS_NUM_THREADS=1 .work/unified-cosmology/calibrated-hosts/.venv/bin/python studies/unified_cosmology/code/calibrated_hosts/bands.py
OPENBLAS_NUM_THREADS=1 .work/unified-cosmology/calibrated-hosts/.venv/bin/python studies/unified_cosmology/code/calibrated_hosts/independent_review.py
OPENBLAS_NUM_THREADS=1 .work/unified-cosmology/calibrated-hosts/.venv/bin/python studies/unified_cosmology/code/calibrated_hosts/validate.py
```

`sources.py` records the installed SPARCL client's initialization outcome; this is diagnostic and is not needed by the native retrieval route. Catalogue cones use TAP POST requests to avoid the service's long-URL limit. Failed discovery requests are retained, never classified as scientific nonmatches. `spectra.py` reads target rows from NERSC FITS files with HTTP byte ranges, without downloading the complete HEALPix files. Downloads and regenerated per-host measurements remain under ignored `.work/unified-cosmology/calibrated-hosts/`.

The spectral interface is `spectra/<targetid>.npz`, with `B_`, `R_`, `Z_` arrays `WAVELENGTH`, `FLUX`, `IVAR`, `MASK`, `RESOLUTION`. Adjacent JSON records the complete native row metadata, extension headers, host association, redshift and selected-array hash. The remote URL, ETag and size identify the parent FITS file; its full bytes are **not** claimed hashed. `clean-matches.csv` maps DESI `targetid` to `oz_OzDES_ID`, which joins the original physical host/SN crosswalk.

`band-records.json` contains all 55 signed five-band vectors and their complete 5×5 formal covariance. Use these vectors, not independent Dn4000 and Hδ likelihood factors: their windows share pixels. Ratio-support thresholds control displayed diagnostics only. Forward models should use native resolution and explicit response, dust, aperture and emission assumptions. Keep formal diagonal-IVAR covariance separate from calibration and known-release uncertainty.

Data attribution: DESI Collaboration, public DR1/Iron, [release and paper](https://data.desi.lbl.gov/doc/releases/dr1/), CC BY 4.0; catalogue services provided by [NOIRLab Astro Data Lab](https://datalab.noirlab.edu/data/desi). Publications using these data should include the release's [citation and acknowledgment text](https://data.desi.lbl.gov/doc/acknowledgments/). Source data are unmodified; selected arrays and computed diagnostics are identified separately.
