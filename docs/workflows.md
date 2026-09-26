# Running the research

The manuscript calculations share one locked Python environment and one workflow entry point. Run commands from the repository root. Start with the [manuscript](../README.md), then choose a workflow and its [method](methods/README.md).

## Install and verify inputs

Python 3.12 and `uv` are required:

```bash
uv sync --frozen --python 3.12.14 --extra predictors
.venv/bin/python research.py verify
.venv/bin/python research.py list
```

The input manifest identifies 391 frozen files. Bulk inputs are not included in a Git-only clone. The [data guide](methods/data.md) explains acquisition, derived inputs and verification. Public files can be restored where their recorded URLs remain available:

```bash
.venv/bin/python research.py fetch --dry-run
.venv/bin/python research.py fetch --group bao
```

A missing derived bundle or unavailable historical URL is an acquisition gap. The downloader refuses to replace a frozen input with different bytes.

## Run a calculation

```bash
.venv/bin/python research.py run raisin --name my-analysis/raisin
.venv/bin/python research.py run cosmology --name my-analysis/joint-cpl --config configs/joint-cpl.json
```

`--name` is a relative directory below `results/`. Existing destinations and paths outside `results/` are refused. Configurations are merged with documented defaults; unknown options are rejected. Every run writes a summary, tables and a `run.json` containing configuration, input/code/output hashes, dependency versions and execution status.

| Workflows | Purpose | Method |
|---|---|---|
| `hst-geometry`, `hst-repeat`, `hst-flat` | Aperture design, signed repeat photometry and flat-reference uncertainty | [HST imaging](methods/hst.md) |
| `hst-dark-raw`, `hst-dark-calibrated` | Dark read moments and native slope/error comparisons | [HST dark measurements](methods/hst-dark.md) |
| `csp-lineage`, `csp-passbands` | Photometry lineage and physical filter response | [CSP photometry](methods/csp.md) |
| `des-flux`, `des-predictors`, `calibration` | Calibrated light curves, prediction and shared calibration modes | [DES measurements](methods/des.md) |
| `raisin`, `signed-baseline`, `timing` | Infrared distance accounting, signed baselines and simulated timing | [Infrared measurements](methods/raisin.md) |
| `ages`, `populations` | Published ages, regressions and progenitor-delay distributions | [Ages and populations](methods/populations.md) |
| `dust`, `sign-selection` | Explicit dust geometry and measurement-selection models | [Constructed models](methods/constructed-models.md) |
| `cosmology`, `bao-shape` | Released-distance fits and nonaccelerating BAO shape constraints | [Distances and expansion](methods/distances.md) |

## Run the full calculation set

For the 19 default workflows:

```bash
.venv/bin/python replay.py --output my-analysis/baseline
```

For all 19 defaults, five alternative configurations, five additional seeds and the independent audits:

```bash
.venv/bin/python validation/run.py --name my-validation --plan-only
.venv/bin/python validation/run.py --name my-validation
```

A new campaign creates `results/my-validation/{baseline,alternatives,robustness}` and stores its generated configurations and verification record alongside them. It stops on failure and retains the failed run. Audit reports under `validation/reports` are refreshed; preserve the dated Git version before a new campaign if both report sets are needed.

## Read the outputs

The [published results](../results/README.md) are the reference set behind the manuscript. `summary.json` states each result and its interpretation boundary. `run.json` records execution; a completed process can still have failed scientific convergence gates. A numerical replay does not establish the physical correctness of a population or selection model.

`lib/` contains numerical kernels and shared readers; `workflows/` handles analysis and output contracts; `validation/` contains independent checks. The [provenance guide](../provenance/README.md) distinguishes current paths from historical source identities. The [experimental plan](experimental-plan.md) describes work required beyond the completed workflows, including full selection closure and a unified likelihood.

## Additional research and preparation code

The [study guide](../studies/README.md) covers original acquisition, extraction, native instrumented builds and the wider scientific investigations. Those sources are grouped by topic and prepared in a separate ignored workspace; they are not all covered by the manuscript validation campaign. Scientific findings and failed identification tests remain documented. Downloaded dependencies and generated outputs stay local.
