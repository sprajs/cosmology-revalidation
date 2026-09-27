# Independent host-environment checks

These analyses ask whether independently inferred host properties predict supernova brightness after current observable corrections. They do not turn a host-age correlation into a redshift-dependent cosmological correction.

- `analyse.py` checks the published DustPedia local/global and 1/3-kpc measurements, joins physical SNe to Pantheon+, and compares corrected magnitudes with Cepheid distances. Full released covariance retains repeated light curves and shared Cepheid errors.
- `ztf.py` defines an explicit subset of the public ZTF catalogue and fits joint SALT width/colour/environment relations with grouped held-out prediction. It does not reproduce the unpublished final BayeSN mask.
- `titan.py` joins that sample to the newly located public TITAN host-summary snapshot and tests the incremental predictive value of global stellar age and inferred mean progenitor delay.
- `calibrate.py` tests coverage/power for the conditional Gaussian regression on the observed design. It is a statistical model check, not an injection before survey detection.
- `robustness.py` checks catalogue agreement, a free magnitude intercept, and free redshift-bin intercepts. These were added after inspecting primary fits and are sensitivity tests.
- `validate.py` checks source identities, numerical equivalence and result invariants.
- `coverage.py` audits physical-event overlap with DES/Dovekie and RAISIN; `dovekie.py` measures the small nearby Foundation overlap with its full released covariance. These are descriptive extensions, not preregistered discovery tests. `validate_coverage.py` separately verifies their inputs and numerical checks.

Run from the repository root in the locked environment:

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/environment_validation/analyse.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/environment_validation/ztf.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/environment_validation/titan.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/environment_validation/calibrate.py --draws 300
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/environment_validation/robustness.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/environment_validation/validate.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/environment_validation/coverage.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/environment_validation/dovekie.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/environment_validation/validate_coverage.py
```

`--archive PATH` on the first three programs points to a restored input tree instead of the original local read-only archive. Only the following old-layout paths are required, rather than the whole historical checkout:

```text
sources/updates/2026-09-20-ztf/extracted/kelsey2026-stag1765-supplement/{cigale_local_3kpc,cigale_local_1kpc,cigale_global,photometry_local_3kpc}.txt
sources/updates/2026-09-20-ztf/extracted/ztfsniadr2_lite/tables/{snia_data,globalhost_data,localhost_data}.csv
sources/updates/2026-09-20-ztf/originals/Ginolin25ab_masterlist--903d65d.csv
sources/repos/PantheonPlusSH0ES__DataRelease@7fc6805/Pantheon+_Data/4_DISTANCES_AND_COVAR/README
```

Obtain the environmental tables from the [Kelsey journal supplement](https://doi.org/10.1093/mnras/stag1765) and [ZTF release](https://ztfcosmo.in2p3.fr/download/data). The master list is at [its pinned author commit](https://raw.githubusercontent.com/mginolin/standax/903d65df4f52c7096ff44e9f2dc4acd0bcbd36a8/notebooks/Ginolin25ab_masterlist.csv). The distance README and current main input files come from the [pinned Pantheon+ release](https://github.com/PantheonPlusSH0ES/DataRelease/tree/7fc680548d1ea9fe6ee983a89d8b5635bb5784a8/Pantheon%2B_Data/4_DISTANCES_AND_COVAR). The repository's main acquisition supplies `data/distances/Pantheon+SH0ES.dat` and `Pantheon+SH0ES_STAT+SYS.cov`.

`titan.py --titan PATH` accepts [this exact public CSV](https://raw.githubusercontent.com/SterlingYM/age-of-titans/9b9ed9e5faf8b2f6267c0dce8cea97704d9b95ea/host_props_with_SN_age_good_Mar18.csv), SHA-256 `b827165961c269ec4901061e508778b18e892076567dbf51c4b2981258d41031`. Its default location is the ignored population-transport acquisition folder. This 8,610-row prepaper snapshot differs from the published final 6,983-host sample. Each output records source hashes; a same-named later release is not interchangeable.

Design choices are frozen in `design.json` and `titan-design.json`. The extension was registered after the age-summary catalogue was found and its name overlap counted; the earlier ZTF environment fits were already inspected. None of this is retrospectively presented as a fully blind analysis. Expanded crosswalks, selection failures, fold predictions and simulations go to ignored `.work/environment-validation/`. Only small result summaries are retained in `../../results/environment_validation/`.

The likelihood uses fitted light-curve parameters and SED posterior summaries, not raw photons or the original host likelihood. Errors that are available are propagated; missing joint host posterior covariance, calibration modes and selection models remain explicit. A successful numerical test does not validate these absent physical ingredients. See [the scientific results](../../notes/environment-validation-results.md).

The coverage extension additionally reads the main acquisition's `data/raisin` tables and photometry, plus archive paths `data/des-diffimg/DES-SN5YR_DIFFIMG_HEAD.FITS.gz` and `sources/repos/des-science__DES-SN5YR/4_DISTANCES_COVMAT/{DES-Dovekie_HD.csv,DES-Dovekie_Metadata.csv,STAT+SYS.npz,README.md,DES-Dovekie-SN_Likelihood.py}` from the [DES release](https://github.com/des-science/DES-SN5YR/tree/main/4_DISTANCES_COVMAT). Both new analysis programs accept `--archive` and `--titan` overrides. Input hashes freeze the particular release; covariance is matched to Hubble-diagram order, not metadata order. No downloaded tables are committed.
