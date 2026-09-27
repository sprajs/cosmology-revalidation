# Validation

The [manuscript](../README.md) is supported by a full workflow replay and separate numerical implementations. The independent checks use the same frozen data and standard numerical libraries; they are not independent observations.

## Published result checks

- All **391 input files** match their recorded sizes and hashes.
- The **19 default workflow summaries** exactly matched the previous edition in the September 26 replay.
- All **five supplied alternatives** and **five additional seed runs** completed with their applicable convergence gates.
- The original replay checks passed for their stated scope. The later physical programme found a solver-version defect on a new stellar design and a separate cancellation artefact in one exploratory infrared covariance reconstruction; both are documented below. Physical interpretation and historical reconstruction limits remain explicit.

The [current manifest](manifest.json) verifies the 29 published run records and their outputs. Earlier manifests remain in [provenance/history](../provenance/history/) with their original paths and hashes. The [structure verification](reports/structure-verification.json) records a full fresh campaign after the repository reorganization.

## Independent coverage

| Workflows | Check | Boundary |
|---|---|---|
| `cosmology` | Direct distances, full-covariance likelihood, posterior integration and 200 fixed-truth recovery simulations | Upstream covariance and survey selection not regenerated |
| `bao-shape` | Independently constructed inequalities and constrained projection | Flat geometry and conservative composite-null calibration retained |
| `ages` | Supplement extraction, matching, WLS/GLS, joint latent Gaussian likelihood and influence | No original age PDFs or exact published LINMIX reconstruction |
| `populations`, `dust` | Independent clock/delay integration and mixed-source depth integration | Population checks sample three redshifts; dust example is constructed |
| `des-flux`, `des-predictors` | Fit objectives, held-out densities, gradients, rank diagnostics, paired resampling and extra seed | Frozen accepted epochs and selected sample |
| `calibration` | Design products and independent Gaussian inference | Conditional on modes, priors and local response approximation |
| `raisin`, `timing` | Paired accounting, covariance reconstruction, joins, initializer identity and selection | Missing historical covariance operations remain unresolved |
| `csp-lineage`, `csp-passbands` | Literal row matching and analytic photon integrals | Archive coverage gaps retained; no inferred distance correction |
| HST imaging and dark workflows | FITS extraction, weights, reference variance, read moments, signed slopes and spatial influence | Science WCS not fully independently reimplemented; native CALWF3 products are frozen inputs |
| `signed-baseline`, `sign-selection` | Object-level diagnostics, analytic signed fits and adaptive censored likelihood | Censored integral independently refitted for 12/128 replicates; full synthetic workflow replayed |

Detailed notes: [DES](../docs/validation/des.md), [infrared and selection](../docs/validation/infrared.md), [HST](../docs/validation/hst.md). Core numerical results: [core.json](reports/core.json).

## Commands

From the repository root with the verified input bundle:

```bash
uv sync --frozen --extra predictors
.venv/bin/python research.py verify

# Inspect, then execute a new complete campaign.
.venv/bin/python validation/run.py --name next-campaign --plan-only
.venv/bin/python validation/run.py --name next-campaign
```

The campaign runs all 29 calculations, points the two age-template fits to the newly generated population table, and then runs the four audits. Existing result directories are refused. Its configuration files and `validation.json` are stored under `results/next-campaign/`.

To audit the published reference results without rerunning the workflows:

```bash
.venv/bin/python validation/core_checks.py
.venv/bin/python validation/des_checks.py
.venv/bin/python validation/raisin_checks.py
.venv/bin/python validation/hst_checks.py
.venv/bin/python validation/figures.py
.venv/bin/python validation/record.py
```

Audit programs accept `--results` relative to `results/`; DES and core checks also accept the relevant alternative/second-seed locations. See `--help`. Reports under `validation/reports` are refreshed by these commands, so preserve the dated Git version when comparing campaigns. The figure generator always reads the published reference results.

## What the manifests guarantee

The verifier checks input identities, execution status, posterior gates where present, every recorded output hash, and the current audit-code hashes. Historical runs retain their original code hashes; known changes to path handling and project naming are declared in the layout record. Scientific numerical kernels were not changed by the reorganization.

CSV/JSON records and figures are versioned. Bulk inputs and full sampler arrays remain local and must be restored or regenerated to perform all checks. Runtime timestamps and compressed-array bytes can vary across executions; numerical comparisons and scientific gates matter separately from provenance hashes.

## Additional study sources

The [physical-program report](reports/physical-program.json) checks the new galaxy spectroscopy, local and high-redshift photometry, infrared images, stellar-population bounds and native survey simulations. Its [scientific interpretation](../docs/physical-program-results.md) distinguishes supported observations, conditional injections and failed physical gates. The component instructions restore the new inputs before the aggregate check:

```bash
.venv/bin/python validation/physical_program_checks.py
```

The physical tests include independent spectral integration, exact object joins, passband quadrature, convex/primal optimization, injection coverage, signed-photon interventions, CPU/CUDA agreement and source/output hashes. Passing these tests does not complete production BBC, stellar-model validation, contamination inference or a unified cosmology fit. The [legacy numerical audit](../studies/host_ages/results/galaxy_validation/legacy-nnls-audit.json) checks older BAO and infrared designs separately: the new NNLS residual defect does not occur there, but analytically cancelled infrared mass modes must be removed before fitting covariance weights. Original cosmological summaries are retained unchanged.

The [age-correction execution report](reports/age-correction-execution.json) checks the September 27 studies separately from the 29 earlier reference runs. It covers new source/input/output identities, GLS and latent-Gaussian equivalence, injected-signal recovery, independent flux and distance integration, held-out predictions and population-integration checks. Run `.venv/bin/python validation/age_correction_checks.py` after regenerating the study outputs; `--archive` selects a restored historical input tree. The [scientific report](../docs/age-correction-results.md) explains which observational hypotheses remain unresolved. Passing these numerical checks does not certify unavailable host likelihoods, survey selection or high-redshift transport.

Run `.venv/bin/python validation/studies_checks.py` to verify the restored research source inventory, parse its Python and shell programs, check isolated workspace preparation and input-hash rejection, and repeat the luminosity and HST raw-prefix calculations. [The report](reports/studies.json) records that bounded coverage. It does not certify every native build, classifier or spectral-model investigation. Historical escape-sequence syntax warnings remain visible in the record. The additional source library has its own documented scientific limitations.
