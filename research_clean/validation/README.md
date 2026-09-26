# Revalidation record — 26 September 2026

This directory contains the independent checks supporting the repository's [readable manuscript](../../README.md). The original clean edition is preserved byte-for-byte. New runs, checks and interpretation are additive.

## Scope and outcome

- **391/391 inputs** match their recorded sizes and SHA-256 values.
- **552/552 frozen-edition files** remain unchanged.
- **19/19 default workflow summaries** are exactly equal to the frozen edition.
- **Five supplied alternatives** were executed: joint CPL, flexible expansion bins, R19-priority ages, colour-only DES prediction and the fixed age-template scenario.
- **Five extra seeds** check flat ΛCDM, joint CPL, flexible expansion, the fixed-template fit and DES prediction.
- All **29 run records** report completion; applicable posterior gates pass. Independent audit checks pass. No original workflow code was changed in this pass.

The complete record is [manifest.json](manifest.json). Per-workflow replay and independent-seed comparisons are in [replay-comparison.json](reports/replay-comparison.json). Exact replay is not evidence that the assumed measurement or population model is correct.

## Independent coverage

The audit programs use separate numerical implementations rather than importing the workflow kernels. They still share the frozen data, scientific assumptions and standard numerical libraries; they are not independent observations.

| Workflow(s) | Independent check | Boundary |
|---|---|---|
| `cosmology` | Direct Astropy distances, Cholesky likelihood, intercept profiling, posterior quadrature and 200 fixed-truth Gaussian recovery simulations | Full upstream covariance and selection not regenerated; recovery tested at one truth |
| `bao-shape` | Independently constructed inequalities and constrained projection | Flat geometry and conservative composite-null calibration retained |
| `ages` | Literal supplement extraction, row matching, WLS/GLS, joint Gaussian latent likelihood and leave-one-object influence | No original age PDFs; not an exact published LINMIX reproduction |
| `populations` | Independent ODE clock and adaptive integration for all three delay models at three redshifts | Other redshifts replayed, not separately quadrature-checked |
| `dust` | Independent mixed-source depth integration | Constructed mathematical example; empirical dust population not fitted |
| `des-flux` | All 24 full fits and 36 held-out densities, objective reconstruction, stationarity and numerical sensitivities | Frozen accepted epochs; no historical detector reduction |
| `des-predictors` | NumPy predictive score, compiled likelihood gradients, independent rank diagnostics, paired object/field resampling and second seed | Selected sample; other model families checked algebraically, not all refitted |
| `calibration` | Design products and QR-based Gaussian inference | Conditional on modes, prior and frozen local response approximation |
| `raisin` | Paired correction accounting, resampling, systematic covariance reconstruction and signed-flux inventory | Missing historical covariance operations remain unresolved |
| `timing` | ID joins, truth agreement, initializer identity and selection accounting | Archived author simulations, not newly simulated survey data |
| `csp-lineage`, `csp-passbands` | Literal row multisets and analytic piecewise-polynomial photon integrals | Archive gaps retained; no empirical distance correction inferred |
| `hst-repeat`, `hst-flat` | FITS extraction, coverage/units, constant cancellation and reference variance propagation | One temporal pair per visit; unknown detector covariance remains |
| `hst-geometry` | Full fresh replay; independent dark-aperture quadrature and science weight/coverage checks | Science WCS geometry not independently reimplemented in full |
| `hst-dark-raw`, `hst-dark-calibrated` | Independent time regressions, 32 read-moment matrices, 11,856 slopes and spatial deletion diagnostics | Native CALWF3 outputs are frozen inputs; mixed raw moments do not isolate read noise |
| `signed-baseline`, `sign-selection` | Whole-object baseline diagnostics, analytic signed fits and adaptive censored-likelihood integrals | Censored integral independently refitted for 12/128 evenly spaced replicates; all synthetic workflows replayed |

Detailed source and interpretation notes: [DES](reports/des-notes.md), [infrared and selection](reports/raisin-notes.md), [HST](reports/hst-notes.md). Core calculations: [core.json](reports/core.json). These records distinguish mathematical reproduction from physical identification and historical pipeline reconstruction.

## Commands

From `research_clean`, with the complete input bundle and locked environment:

```bash
uv sync --frozen --extra predictors
.venv/bin/python research.py verify

# Inspect the intended configurations without running or writing outputs.
.venv/bin/python validation/run.py --prefix next-campaign --plan-only

# Full new campaign; select a prefix that does not already exist.
.venv/bin/python validation/run.py --prefix next-campaign
```

`run.py` first replays the 19 workflows, then executes alternatives and second seeds. Both age-template fits point to the newly generated population correction, avoiding a dependency on stale results. Every workflow writes its input/code/configuration/environment/output hashes. Failures remain visible; existing result directories are never overwritten. After the four audits, `record.py` verifies hashes and refreshes the additive manifest.

Reports under `validation/reports` are **working outputs that a new audit overwrites**. Preserve this dated Git version before another campaign; use a separate checkout when retaining both report sets is important. The HST audit compares its independent results to the original frozen edition; the replay comparison separately compares fresh summaries with that edition. HST JSON elapsed time is runtime metadata and will vary.

To rerun just the independent checks against this campaign's existing results:

```bash
.venv/bin/python validation/core_checks.py
.venv/bin/python validation/des_checks.py
.venv/bin/python validation/raisin_checks.py
.venv/bin/python validation/hst_checks.py
.venv/bin/python validation/figures.py
.venv/bin/python validation/record.py
```

The core, DES and infrared scripts accept alternate result prefixes; see `--help`. `figures.py` deliberately uses the dated September 26 campaign, not whichever campaign happens to be newest. Figures show posterior intervals or labelled spatial sensitivity ranges; those quantities are not interchangeable.

## Publication and provenance

JSON and CSV result records, independent numerical code, notes and figures are versioned. New full sampler arrays and bulk acquisition files remain local. Their hashes are retained in each `run.json`; absence from a Git-only clone is not a passed input verification. The [source guide](../SOURCES.md) explains public restoration and derived-input dependencies.

The [original edition manifest](../provenance/edition.json) continues to describe the original edition. The [new manifest](manifest.json) describes this audit and verifies that the original was preserved. The extra [author-release comparison](reports/release-comparison.json) was made against the pinned acquisition outside the clean folder; it records literal column agreement and the identical covariance hash. The four core audit programs themselves require only clean-edition inputs and results.

No scientific workflow defect was found by these checks. This is a bounded conclusion about the inspected calculations, not certification of all upstream observations, published conclusions or historical code paths.
