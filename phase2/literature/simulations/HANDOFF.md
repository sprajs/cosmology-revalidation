# Forward simulation handoff

The initial P21 pilot is **complete**, with 26,518 generated attempts, 2,678 written light curves and all generated rows in DUMP. See `manifests/PH2_pilot02_P21.json` and `PH2_pilot02_P21-denominator.json`. Historical binary source tree was checked against commit `2fe0f564361a873860661ff61080b4db9c607edf` by the official reproduction agent. The apparently anachronistic source comment is present in the upstream pinned tree; it is not evidence of a local unpinned replacement.

The initial pilot01 failed because the 2023 code builds an invalid zero-width Gaussian map when a 2024 population file specifies only `GENPEAK_SALT2ALPHA`. The official reproduction agent diagnosed this from GDB. Pilot02 removes that single line from a **derived** population file and places the same exact constant alpha=0.15 in the input; this preserves a delta distribution without adding scatter. Original source/PDF/input/log files remain unchanged. See `inputs/constant-alpha-compatibility.json` and `phase2/official/diagnostics/simulation-segfault-gdb.log`.

To run each remaining predeclared pilot from repository root:

```bash
.venv/bin/python scripts/phase2/literature/run_forward.py BS21
.venv/bin/python scripts/phase2/literature/run_forward.py G10
.venv/bin/python scripts/phase2/literature/run_forward.py P21_dmplus010
.venv/bin/python scripts/phase2/literature/run_forward.py P21_dmminus010
.venv/bin/python scripts/phase2/literature/run_forward.py P21_dmz020
.venv/bin/python scripts/phase2/literature/run_forward.py P21_rho000
.venv/bin/python scripts/phase2/literature/run_forward.py P21_rho090
.venv/bin/python scripts/phase2/literature/run_forward.py P21_noisetrue120
```

Each command sets the project-local executable, SNDATA_ROOT and library paths, hashes the executable/input/output/log and refuses to overwrite an existing run. A failed non-nominal input may need a new revision, never an overwritten prior input. Population pilots are predeclared in `inputs/preregistration-pilot02.json`; noise-misspecification pilots in `inputs/preregistration-noise-pilot02.json`. Values are stress brackets, not estimated uncertainty priors.

Do **not** rerun nominal P21 over its existing output. Its raw generated/lightcurve files are in `outputs/PH2_pilot02_P21`. The official reproduction agent has been asked to fit the 2,678 events using the same measured-flux SALT configuration as data/released mocks. Once fitting completes, join measured quality flags back to generated pass rows and estimate total selection conditional on generated variables. Use `audit_forward.py VERSION` to verify generated totals and export a table with unique attempt indices.

Critical denominator rule: the nominal DUMP contains 2,475 duplicated CIDs because early NEPOCH rejects reuse the event counter. The **one-based generated row index** is the attempt ID. HEAD matches exactly the unique DUMP rows with `SIM_EFFMASK=5` and `CUTMASK=4095`; joining every generated row by CID incorrectly aliases failures to selected events. Sentinel epoch/SNR fields for early rejects are not measurements.

Remaining scientific limits:

1. W22 needs SN_age. Nominal HOSTLIB lacks it; the public Wiseman builder references an absent external host-model HDF5 and SN_ages/*.dat. No exact W22 generator or host-cadence crosswalk is established.
2. P21 uses renamed/released 2024 assets; no byte-exact original 2022 mock replay claim. BS21/G10 config choices are stated explicitly. Final nominal counts differ from the old released realization.
3. Constant alpha relocation is exact for that coordinate; the derived input does not add the historical tiny C11 scatter scaffold. This is documented as public-assets nominal construction.
4. Only P21 pilot executed here. Other ready inputs need successful execution and measured-flux refitting. Same seed does not prove paired host/redshift/cadence realizations; compare generated streams before using paired uncertainty estimates.
5. One realization gives a pilot selection response, not repeated-truth frequentist coverage. Full type mixture, final analysis-membership intersection, and uncertainty in empirical calibration maps still need treatment.
6. Noise REDCOV is correlation of **added fudge variance**, not total flux covariance; search and common-template noise are separate. The data-audit agent owns independent real off-signal noise tests and image accessibility.
7. No cosmology fit was run by this workstream.
