# Reproduction command inventory (draft, 20 September 2026)

## Status and notation

This is a command and input inventory, not a claim that every result was recomputed in this pass. No fit or long simulation was rerun while preparing it.

- **Help-tested** means the program's `--help` path was executed successfully in the current `.venv`.
- **Source-derived** means the invocation was read from the current Python entry point or the corresponding investigation note. These programs generally expose no argument parser.
- **Manifest-derived** means every option was reconstructed from the current saved `configuration.json`. The fit itself was not rerun.

Run commands from the repository root:

```bash
cd /home/szymon/Documents/ChatGPT/supernova
uv sync --frozen
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
```

`uv sync --frozen` was not executed in this bounded audit because the pinned `.venv` already exists. `pyproject.toml` requires Python `>=3.12,<3.13` and pins every direct numerical dependency; `uv.lock` is the lock source.

## Age-signal branch

`analyse.py --help` was help-tested and exposes exactly the positional actions `prepare`, `diagnostics`, and `linmix`. The remaining commands are source-derived. `linmix`, `covariance_linmix.py`, and `mass_revision.py` can run long MCMC jobs; `fullcov_eiv.py` performs repeated numerical optimisations.

```bash
.venv/bin/python scripts/age_signal/analyse.py prepare
.venv/bin/python scripts/age_signal/analyse.py diagnostics
.venv/bin/python scripts/age_signal/analyse.py linmix
.venv/bin/python scripts/age_signal/mock_audit.py
.venv/bin/python scripts/age_signal/sensitivity.py
.venv/bin/python scripts/age_signal/covariance_linmix.py
.venv/bin/python scripts/age_signal/fullcov_eiv.py
.venv/bin/python scripts/age_signal/validate_algebra.py
.venv/bin/python scripts/age_signal/mass_revision.py
.venv/bin/python scripts/age_signal/summarize.py
.venv/bin/python scripts/age_signal/freeze.py
```

The intended order is preparation, diagnostics, primary/sensitivity fits, summaries, then the evidence freeze. Existing LINMIX result files are reused by some scripts, so a clean numerical replay requires archiving the relevant output directories first rather than silently mixing old and new chains.

## Mapping branch

These commands are source-derived and use no command-line options:

```bash
uv run --frozen python scripts/mapping/controlled.py
uv run --frozen python scripts/mapping/csfh_dtd.py
uv run --frozen python scripts/mapping/hostlib_audit.py
```

`scripts/mapping/common.py` is an imported helper. `scripts/mapping/acquire_sources.py` is a source-acquisition utility with no `__main__` entry point, so it is not listed as a scientific reproduction command.

## Standardization and dust branch

These commands are source-derived and use no command-line options:

```bash
.venv/bin/python scripts/standardization/compare_des.py
.venv/bin/python scripts/standardization/dust_identifiability.py
.venv/bin/python scripts/standardization/pantheon_mass_update.py
```

## Cosmology branch

The help-tested top-level interface is:

```text
run.py [-h] {validate,fit} ...
run.py fit [-h] --name NAME
           [--model {lcdm,wcdm,cpl,kinematic,qbins}]
           [--data {sn,bao,joint}]
           [--amplitude {none,fixed,normal,uniform}]
           [--correction CORRECTION] [--column COLUMN] [--zcolumn ZCOLUMN]
           [--magnitude-table MAGNITUDE_TABLE]
           [--scale SCALE] [--tau TAU] [--steps STEPS] [--burn BURN]
           [--walkers WALKERS] [--seed SEED]
```

The validation and preparation/diagnostic commands below are source-derived. Several scripts execute at module top level and intentionally expose no options.

```bash
.venv/bin/python scripts/cosmology/run.py validate
.venv/bin/python scripts/cosmology/prepare_templates.py
.venv/bin/python scripts/cosmology/reference_chains.py
.venv/bin/python scripts/cosmology/diagnostics.py
.venv/bin/python scripts/cosmology/kinematic_limits.py
.venv/bin/python scripts/cosmology/summarize.py
```

The following 19 fit commands are manifest-derived from the current direct children of `runs/cosmology/*/configuration.json`. They are exact saved configurations, including seeds and sampler lengths, and were not executed in this audit:

```bash
.venv/bin/python scripts/cosmology/run.py fit --name bao-cpl --model cpl --data bao --amplitude none --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 12000 --burn 1500 --walkers 40 --seed 41004
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-bao-cpl --model cpl --data joint --amplitude none --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 10000 --burn 1500 --walkers 40 --seed 41005
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-bao-cpl-c14fixed --model cpl --data joint --amplitude fixed --correction runs/cosmology/templates/c14-cpl63-median.csv --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 12000 --burn 1500 --walkers 40 --seed 41010
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-bao-cpl-c14free --model cpl --data joint --amplitude uniform --correction runs/cosmology/templates/c14-cpl63-median.csv --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 18000 --burn 1500 --walkers 40 --seed 41013
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-bao-cpl-c14lcdm70 --model cpl --data joint --amplitude fixed --correction runs/cosmology/templates/c14-lcdm70-median.csv --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 12000 --burn 1500 --walkers 40 --seed 41018
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-bao-cpl-c14mean --model cpl --data joint --amplitude fixed --correction runs/cosmology/templates/c14-cpl63-mean.csv --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 12000 --burn 1500 --walkers 40 --seed 41016
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-bao-cpl-c14slope --model cpl --data joint --amplitude normal --correction runs/cosmology/templates/c14-cpl63-median.csv --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 16000 --burn 1500 --walkers 40 --seed 41011
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-bao-cpl-c14slope-repeat --model cpl --data joint --amplitude normal --correction runs/cosmology/templates/c14-cpl63-median.csv --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 16000 --burn 1500 --walkers 40 --seed 51221
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-bao-cpl-hostmass --model cpl --data joint --amplitude none --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 12000 --burn 1500 --walkers 40 --seed 41015 --magnitude-table data/derived/standardization/Pantheon_W26_lowz_mass_revision.dat
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-bao-cpl-repeat --model cpl --data joint --amplitude none --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 14000 --burn 1500 --walkers 40 --seed 51222
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-bao-cpl-shortdelay --model cpl --data joint --amplitude fixed --correction runs/cosmology/templates/cut40-cpl63-mean.csv --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 12000 --burn 1500 --walkers 40 --seed 41017
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-cpl --model cpl --data sn --amplitude none --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 10000 --burn 1500 --walkers 40 --seed 41003
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-cpl-c14fixed --model cpl --data sn --amplitude fixed --correction runs/cosmology/templates/c14-cpl63-median.csv --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 12000 --burn 1500 --walkers 40 --seed 41009
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-kinematic --model kinematic --data sn --amplitude none --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 10000 --burn 1500 --walkers 40 --seed 41006
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-lcdm --model lcdm --data sn --amplitude none --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 7000 --burn 1500 --walkers 40 --seed 41001
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-qbins-c14fixed --model qbins --data sn --amplitude fixed --correction runs/cosmology/templates/c14-cpl63-median.csv --column delta_mu --zcolumn z --scale 1.0 --tau 1.5 --steps 12000 --burn 1500 --walkers 40 --seed 41014
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-qbins-tau05 --model qbins --data sn --amplitude none --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 12000 --burn 1500 --walkers 40 --seed 41007
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-qbins-tau15 --model qbins --data sn --amplitude none --column delta_mu --zcolumn z --scale 1.0 --tau 1.5 --steps 12000 --burn 1500 --walkers 40 --seed 41008
.venv/bin/python scripts/cosmology/run.py fit --name pantheon-wcdm --model wcdm --data sn --amplitude none --column delta_mu --zcolumn z --scale 1.0 --tau 0.5 --steps 10000 --burn 1500 --walkers 40 --seed 41002
```

`replay_fits.py --help` was help-tested and exposes `--workers WORKERS` (default 3) and `--archive-first`. The source-derived bulk invocation is:

```bash
.venv/bin/python scripts/cosmology/replay_fits.py --workers 3 --archive-first
```

That command reruns every saved direct-child fit concurrently and writes `replay.log` files. With `--archive-first`, it copies selected current products to `runs/cosmology/preliminary/` only when a destination file does not already exist. It is a large, output-mutating operation and was not run here.

## Directional branch

These commands are source-derived. `audit.py` accepts `all`, `data`, or `synthetic` through `sys.argv` and defaults to `all`; it does not validate any other string. `c2.py` accepts `check` or `run` and defaults to `run`; in the current implementation any value other than `check` enters the fit path.

```bash
.venv/bin/python scripts/directional/audit.py data
.venv/bin/python scripts/directional/input_audit.py
.venv/bin/python scripts/directional/audit.py synthetic
.venv/bin/python scripts/directional/lowz_wide.py
.venv/bin/python scripts/directional/check_optimizer_flags.py
.venv/bin/python scripts/directional/c2.py check
.venv/bin/python scripts/directional/c2.py run
.venv/bin/python scripts/directional/ray_check.py
```

The C2 gradient check is intended to precede the C2 fit. The fit and synthetic commands can be expensive; none was rerun in this audit.

## Independent audit branch

These commands are source-derived and expose no command-line options:

```bash
.venv/bin/python scripts/audit/independent.py
.venv/bin/python scripts/audit/directional_crosscheck.py
.venv/bin/python scripts/audit/c2_crosscheck.py
```

## Manifest shapes and current-input audit

The repository has no single enforced JSON manifest schema. At this audit snapshot there are **69** files matching `runs/**/*manifest*.json`: 68 JSON objects and one JSON-list acquisition ledger. The 68 object manifests occupy ten distinct top-level-key shapes; counting the list shape gives eleven. The observed shapes are:

| Count | Top-level keys |
|---:|---|
| 25 | `code_sha256, configuration, created_utc, git_dirty, git_revision, inputs_sha256, lock_sha256, outputs_sha256, purpose, python` |
| 20 | the preceding keys plus `code_capture` |
| 10 | `code_sha256, config, git_head, inputs, outputs, plan_sha256, platform, purpose, python, time_utc, versions` |
| 4 | `command, inputs, plan, seed, time, versions` |
| 3 | `code_revision, config, extra, input_sha256, name, output_sha256, python, utc, versions` |
| 2 | `code_revision, configuration, environment, experiment, inputs, outputs, plan_sha256, script_sha256, utc` |
| 1 | JSON list acquisition ledger |
| 1 | `code_imports, fit_snapshot_age_cases, inputs_sha256, outputs_sha256, purpose, scope, time_utc` |
| 1 | `code_revision, environment, experiment, inputs, outputs, plan_sha256, script_sha256, seed, utc` |
| 1 | `git_revision, implementations, inputs_sha256, outputs_sha256, purpose, seed, time_utc` |
| 1 | `inputs_sha256, outputs_sha256, purpose, random_seed, time_utc` |

A read-only scan resolved relative paths from the repository root and inspected the fields `inputs`, `inputs_sha256`, and `input_sha256`. All **358 referenced input occurrences exist**; there are **zero missing referenced paths**. This is an input-presence and current-input-hash audit, not a complete manifest validator or output-integrity proof.

There are 36 current-hash mismatches in five manifests:

- `runs/audit/c2-manifest.json`: the recorded `scripts/audit/c2_crosscheck.py` hash predates the current file.
- `runs/audit/manifest.json`: the recorded independent-audit plan hash predates the current plan.
- `runs/cosmology/summary/manifest.json`: 32 saved summary/configuration input hashes predate the current replayed fit products.
- `runs/directional/initial-results/manifest-data.json` and `manifest-synthetic.json`: each records an earlier `scripts/directional/audit.py` hash.

These mismatches were left intact as provenance evidence. A final integrity checker should distinguish immutable historical snapshots from manifests expected to describe current bytes, validate output hashes separately, and report schema-specific missing fields rather than treating every JSON file as one schema.

## Newly acquired primary-source text paths

The exact ZTF v3 text/source paths for methods review are:

```text
/home/szymon/Documents/ChatGPT/supernova/sources/updates/2026-09-20-ztf/extracted/2605.06799v3.txt
/home/szymon/Documents/ChatGPT/supernova/sources/updates/2026-09-20-ztf/extracted/mnras.tex
```

The two September primary papers are available as PDF text and unpacked arXiv source at:

```text
/home/szymon/Documents/ChatGPT/supernova/sources/updates/2026-09-20-ztf/extracted/2609.12083v1.txt
/home/szymon/Documents/ChatGPT/supernova/sources/updates/2026-09-20-ztf/extracted/2609.12083v1-source/
/home/szymon/Documents/ChatGPT/supernova/sources/updates/2026-09-20-ztf/extracted/2609.16972v1.txt
/home/szymon/Documents/ChatGPT/supernova/sources/updates/2026-09-20-ztf/extracted/2609.16972v1-source/
/home/szymon/Documents/ChatGPT/supernova/sources/updates/2026-09-20-ztf/extracted/kelsey2026-stag1765-supplement/
```

Their source-backed methods, sample reuse, claims, data availability, and public-input barriers are recorded in `docs/investigations/ztf-source-update.md`.
