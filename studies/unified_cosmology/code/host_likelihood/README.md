# Observed host likelihood inputs

The [scientific note](../../notes/host-likelihood.md) distinguishes observations, model-conditioned joint chains, conditional spectral features and failed brightness/model gates. Numerical validation does not certify an age correction.

Run from the repository root using the frozen main environment. First construct the DES–OzDES crosswalk with the neighbouring [survey-selection code](../survey_selection/ozdes.py), which depends on its `acquire.py`, `acquire_followup.py`, `dataprep.py` and `dovekie.py` outputs. Its required interface is `.work/unified-cosmology/survey-selection/followup/dovekie-ozdes-clean.csv` and the matching `results/survey_selection/ozdes-crossmatch.json` provenance record.

```bash
.venv/bin/python studies/unified_cosmology/code/host_likelihood/acquire.py --posterior-archive
.venv/bin/python studies/unified_cosmology/code/host_likelihood/recover.py
.venv/bin/python studies/unified_cosmology/code/host_likelihood/extract.py
.venv/bin/python studies/unified_cosmology/code/host_likelihood/source_audit.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/unified_cosmology/code/host_likelihood/yse_fit.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/unified_cosmology/code/host_likelihood/yse_audit.py
.venv/bin/python studies/unified_cosmology/code/host_likelihood/ozdes_acquire.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/unified_cosmology/code/host_likelihood/ozdes_bands.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/unified_cosmology/code/host_likelihood/ozdes_repeat.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/unified_cosmology/code/host_likelihood/validate.py
```

The main environment includes NumPy, SciPy, pandas, astropy, h5py, sncosmo and iminuit; `pdftotext` and Git are used for source recovery. `yse_fit.py` also needs the existing pinned SALT3 surfaces in `data/des/model`. It downloads native PS1/ZTF passbands through sncosmo and pins the actual arrays. Four workers bound light-curve fitting and spectrum acquisition; no GPU is required. The host posterior archive is about 1.85 GB. Individual spectra and all generated ledgers stay under `.work/unified-cosmology/host-likelihood/`.

The principal interfaces are:

- `host-photometry.csv`: native maggies/errors/masks, filter identifiers, prior foreground correction and error-floor state.
- `joint-host-draws/<name>.npz`: aligned native and derived draws, full joint covariance, observed host fluxes and native filters. Numeric/string arrays only; no author pickle is deserialized.
- `yse-host-crosswalk.csv`: SN/host identity, positional/redshift gates and exact source paths/hashes. `yse-flux-fits.json` separates numerical/sampling success from model adequacy; no object currently passes the scientific brightness-analysis gate.
- `ozdes-band-likelihood.npz`: host IDs, seven signed count-spectrum windows, 7×7 working covariance at three assumed pixel correlations, coverage and ratio-profile grids. The ratios are not calibrated Dn4000.
- `ozdes-Hdelta-records.json` and `ozdes-repeat/*-Hdelta.npz`: signed continuum/contrast, full 2×2 covariance and exact pixel operators.
- `ozdes-repeat-records.json` and `ozdes-repeat/*-exposures.npz`: individual-exposure likelihood vectors, covariance, masks/quality metadata and common-response diagnostics.

Compact provenance and results are in [the result directory](../../results/host_likelihood). Failed acquisition attempts, absent chains, sparse light curves, optimizer warnings, incomplete selection and response failures remain explicit. The design and repeat-test addendum record when outcome inspection occurred. The additional residual screen is exploratory after observing model failures; it is not a preregistered discovery test.
