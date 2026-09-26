# Supernova research — curated executable edition

This folder contains the reusable research extracted from the working investigation as of 26 September 2026. It has its own data bundle, numerical code, locked environment and rewritten methods. It can be moved independently of the parent workspace.

The research asks how supernova distance inferences depend on population assumptions, calibration, photometry, selection and numerical fitting. It distinguishes measured quantities, released fitted products, constructed examples and conditional inference. A repeatable calculation is not automatically an identified physical correction.

## Run

Requires `uv` and Python 3.12. The original workspace is not required.

```bash
uv sync --frozen --python 3.12.14 --extra predictors
uv run --frozen --extra predictors python research.py verify
uv run --frozen --extra predictors python research.py list
uv run --frozen --extra predictors python research.py run raisin --name my-raisin
```

Every run creates a new directory under `results/`. Existing names are refused. `summary.json` contains the result and its interpretation boundary; `run.json` records settings, input and code hashes, package versions, execution status and output hashes. Tables and chains are saved alongside them. The supplied result folders were generated using this edition; see [verification](provenance/validation.json) for what was checked.

## Research workflows

| Commands | Purpose | Method |
|---|---|---|
| `cosmology`, `bao-shape` | Fit released distances; test a flat nonaccelerating BAO shape constraint | [Distances and expansion](methods/distances.md) |
| `ages`, `populations` | Extract published host ages, fit conditional regressions, calculate progenitor-delay distributions | [Ages and populations](methods/populations.md) |
| `dust`, `sign-selection` | Calculate explicit dust-geometry and sign-selection models | [Constructed models](methods/constructed-models.md) |
| `des-flux`, `des-predictors`, `calibration` | Refit calibrated fluxes, compare selected-sample predictors, propagate shared calibration modes | [DES measurements and fitting](methods/des.md) |
| `raisin`, `signed-baseline`, `timing` | Track released distance corrections, signed precursor photometry and simulated timing | [RAISIN](methods/raisin.md) |
| `csp-lineage`, `csp-passbands` | Extract original photometry and calculate physical-filter response | [CSP](methods/csp.md) |
| `hst-geometry`, `hst-repeat`, `hst-flat` | Rebuild aperture design, extract signed repeat photometry and propagate flat-reference variance | [HST pixels](methods/hst.md) |
| `hst-dark-raw`, `hst-dark-calibrated` | Extract dark read moments and compare frozen native dark slopes with quoted variance | [HST dark controls](methods/hst-dark.md) |

Settings are explicit JSON files passed with `--config`. Unknown settings are rejected. Examples are in [configs](configs/); each method describes its defaults and required gates. A convenient complete replay is:

```bash
uv run --frozen --extra predictors python replay.py --prefix replication
```

This runs each workflow once, sequentially, retains failed runs and stops on an execution failure. It is a scientific calculation, including MCMC and aperture construction, rather than a quick installation check.

## What is included

- `lib/`: numerical kernels, separated from execution and file writing.
- `workflows/`: extraction, analysis, fitting and numerical validation commands.
- `data/`: only the frozen inputs used by those commands, including necessary intermediate scientific products.
- `methods/`: purpose, equations, assumptions, expected output and interpretation limits.
- `results/`: one generated result set per workflow, with reproducibility records.
- `provenance/`: input sources/hashes, code origin and the edition's verification record.

The methods and workflow interfaces were rewritten. Tested mathematical kernels were retained where changing the algebra would add risk without improving the research. Historical experiment numbering, abandoned runs, agent handoffs, discussion logs, duplicate scripts and test narratives were excluded.

## Current scientific boundary

The baseline released-distance fit favours acceleration within the specified model. Extra age, dust or calibration terms remain conditional assumptions or sensitivity calculations. RAISIN correction accounting, CSP label lineage and HST repeat statistics identify concrete processing questions, but do not by themselves supply a corrected cosmology.

This edition does not contain every prior investigation. In particular, it does not claim an exact reconstruction of the full native SNANA build/instrumentation experiments, released classifier probabilities, a survey-wide selection-normalized population likelihood, historical HST reductions, or a raw CMB likelihood. Native RAISIN branch ambiguity and unresolved prospective simulation gates remain outside the completed workflows. The dark controls reproduce total repeat moments, without isolating electronic read covariance or calibrating a noise correction. The signed-optical BayeSN adapter has passed bounded implementation checks but has no completed observed posterior or held-out NIR result; it is not included as an inference workflow.

Some DES calculations start from explicitly labelled, locally derived fixed objectives or projected matrices. The [data guide](methods/data.md) explains this reproducibility boundary and how public inputs can be restored. The original research remains the full archival record.
