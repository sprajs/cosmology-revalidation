# Reproducing and inspecting this investigation

The supplied local workspace is the evidence bundle. Git preserves lightweight code, plans, results and manifests; original input archives, derived data, large covariance arrays and chains remain local and intentionally ignored. **A fresh Git clone alone cannot run everything.** Preserve/copy `data/`, `papers/`, `sources/` and the ignored numerical/source files under `runs/` with the repository, or reacquire each exact URL/commit and verify the recorded digest. Do not replace a missing author input with the newest similarly named release.

The [current source map](literature-and-data-status.md) and [unsent data requests](author-data-requests.md) distinguish public local inputs from genuinely unavailable exact-author products. The latter remain blocked even with a complete copy of this workspace.

## Environment and non-destructive inspection

Run from `/home/szymon/Documents/ChatGPT/supernova`. The environment is local: Python 3.12.14, `pyproject.toml`, and `uv.lock` pin the numerical stack. No global configuration was changed.

```bash
cd /home/szymon/Documents/ChatGPT/supernova
uv sync --frozen --python 3.12.14
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
uv run --frozen python scripts/verify_record.py
```

The checker reads source/output digests and exact code snapshots, verifies local documentation links and writes `runs/verification/record-check.json`. It distinguishes historical failed/preliminary manifests from current output claims. It does not rerun scientific fits and does not certify physical assumptions. Different branches use documented manifest schemas; path-only output lists are explicitly identified and hashed at final verification rather than falsely described as having original output hashes.

`runs/final-manifest.json` hashes the lightweight final record and records the completed 19-fit replay. Its upstream manifests bind the bulk inputs. `scripts/freeze_record.py` creates that packaging manifest after a successful check; it is a maintainer command, not needed for read-only inspection. The verification output is excluded from packaging hashes to avoid a self-referential cycle.

The calculation requires the local acquisition bundle. Source URLs, commits, hashes and acquisition errors are in `catalog/`, the per-experiment manifests, and dated branch source ledgers. `catalog/integrity.json` is an earlier acquisition check, not the final scientific-output check. The original acquisition state is Git commit `472f13f`; numerical development and final record are on `codex/independent-supernova-investigation`.

## Authoritative outputs and historical attempts

- Final branch reports are under `docs/investigations/`. The [final report](final-report.md) synthesizes them and distinguishes local calculations from source-verified author results.
- Current cosmological posteriors are direct children of `runs/cosmology/<experiment>/`, with `configuration.json`, `summary.json`, `q_history.csv`, local `chains.npz` and `manifest.json`. There are 19 configurations, including two independent-seed repeats. `runs/cosmology/summary/all_fits.csv` is the comparison table.
- Every fit was replayed with code bytes captured **at process import**, stored by SHA-256 in `runs/cosmology/provenance/`. All posterior/scientific summary fields agree with the initial seeded results; differences are runtime and an added zero-row revision diagnostic. `runs/cosmology/replay-comparison.json` records the comparison. Final source hashes, not an earlier dirty Git revision alone, identify executed code.
- `runs/cosmology/preliminary/` retains initial summaries/configurations/manifests, whose old code-at-completion recording was inadequate to prove the executed bytes. They are not authoritative manifests for current output paths. `validation-initial-16/` retains the failed coarse interpolation check; 40 nodes per interval pass the final independent/full-likelihood checks.
- `runs/age_signal/unsuccessful/` preserves unseeded/short pilots. Current age posterior summaries and chains have seeded parallel LINMIX and independent retained-chain diagnostics. Exact historical code/plan bytes are archived only when their digest matches the original manifest. A silent LINMIX-loop issue and ignored serial seed are documented rather than patched in the original upstream source.
- `runs/directional/initial-results/` is historical. Final C1/C2, coordinate and injection records are identified by the directional report/manifests and their independent audit snapshots.
- Empty posterior tails in finite chains are sample outcomes, not exact zero/one probabilities. Profile intervals, posterior intervals and likelihood-ratio improvements are labelled separately.

## Full calculation order

Commands below overwrite their own derived outputs. For a clean replay, work in a **separate full copy** of the bundle, preserving the original. Retain the saved cosmology `configuration.json` files. Archive/delete only that copy's old computed age posterior outputs before running its fits: some age scripts reuse existing result files. Preserve original inputs and the source records. Do not run many MCMC jobs simultaneously with multithreaded BLAS.

The longer [command inventory](reproduction-commands-draft.md) records exact CLI options for every saved fit, which help paths were tested, and a historical pre-freeze presence/hash audit. Its transitional mismatch counts are superseded by the final `record-check.json`. The root actually replayed all 19 cosmology fits after that inventory; the inventory author's statement that they did not rerun them describes only that bounded indexing task.

1. **Standardization inputs and host-mass revision** (produces the revised table consumed by later sensitivities):

   ```bash
   uv run --frozen python scripts/standardization/compare_des.py
   uv run --frozen python scripts/standardization/dust_identifiability.py
   uv run --frozen python scripts/standardization/pantheon_mass_update.py
   ```

2. **Age sample, regressions and tests:**

   ```bash
   uv run --frozen python scripts/age_signal/analyse.py prepare
   uv run --frozen python scripts/age_signal/analyse.py diagnostics
   uv run --frozen python scripts/age_signal/analyse.py linmix
   uv run --frozen python scripts/age_signal/mock_audit.py
   uv run --frozen python scripts/age_signal/sensitivity.py
   uv run --frozen python scripts/age_signal/covariance_linmix.py
   uv run --frozen python scripts/age_signal/fullcov_eiv.py
   uv run --frozen python scripts/age_signal/validate_algebra.py
   uv run --frozen python scripts/age_signal/mass_revision.py
   uv run --frozen python scripts/age_signal/summarize.py
   uv run --frozen python scripts/age_signal/freeze.py
   ```

   LINMIX is imported from its locally preserved immutable upstream source, not a global install. Long regressions can take substantially longer than the scalar baseline fits. Full-covariance latent-age profiles remain simplified sensitivity models, not exact age-posterior/selection inference.

3. **Host/progenitor physics and frozen correction templates:**

   ```bash
   uv run --frozen python scripts/mapping/controlled.py
   uv run --frozen python scripts/mapping/csfh_dtd.py
   uv run --frozen python scripts/mapping/hostlib_audit.py
   uv run --frozen python scripts/cosmology/prepare_templates.py
   ```

4. **Cosmology, diagnostics and figures:**

   ```bash
   uv run --frozen python scripts/cosmology/run.py validate
   uv run --frozen python scripts/cosmology/replay_fits.py --workers 3
   uv run --frozen python scripts/cosmology/reference_chains.py
   uv run --frozen python scripts/cosmology/kinematic_limits.py
   uv run --frozen python scripts/cosmology/diagnostics.py
   uv run --frozen python scripts/cosmology/summarize.py
   ```

   `replay_fits.py` reads the saved direct-child configurations and dispatches each exact model/correction/prior-scale/sampler-length/seed combination. It does not invent missing configurations. For an isolated representative fit:

   ```bash
   uv run --frozen python scripts/cosmology/run.py fit --name reproduced-baseline --model cpl --data joint --steps 10000 --burn 1500 --walkers 40 --seed 41005
   ```

   The reference-chain command processes downloaded DESI chains; it does not rerun a CMB likelihood. The final figures were visually inspected for clipping, captions, uncertainty intervals and explicit differences in vertical scales.

5. **Directional branch:**

   ```bash
   uv run --frozen python scripts/directional/audit.py data
   uv run --frozen python scripts/directional/input_audit.py
   uv run --frozen python scripts/directional/audit.py synthetic
   uv run --frozen python scripts/directional/lowz_wide.py
   uv run --frozen python scripts/directional/check_optimizer_flags.py
   uv run --frozen python scripts/directional/c2.py check
   uv run --frozen python scripts/directional/c2.py run
   uv run --frozen python scripts/directional/ray_check.py
   ```

   These require the recovered, hash-verified C2 covariance. The HEL principal fit is reproduced; a full multiframe/global-optimum/null-distribution analysis is not claimed. The public pickle is not executed; the age curve is reconstructed from safe tabular/notebook evidence.

6. **Independent numerical audit, then integrity check:**

   ```bash
   uv run --frozen python scripts/audit/independent.py
   uv run --frozen python scripts/audit/directional_crosscheck.py
   uv run --frozen python scripts/audit/c2_crosscheck.py
   uv run --frozen python scripts/verify_record.py
   ```

   These audit implementations do not import the investigated scientific likelihoods. A new replay changes timing/manifest hashes and may require refreshing the final inventory; numerical equivalence must be checked independently of wall time. Cross-platform BLAS, floating-point ordering and stochastic diagnostics can introduce small numerical differences even with fixed seeds.

## Interpretation and completion boundary

The [experiment register](experiment-register.md) links each plan to its implementation/output. Reproductions were not blinded to published values. New controlled tests were recorded before their outcomes; amendments and failures remain in the decision log. Several principal Astra Ultra investigators recorded separate initial analyses, with an additional falsification auditor; a faster acquisition agent gathered bounded late sources and command inventories. Agreement between agents is not independent empirical evidence.

Exact author-only age PDF/mapping/selection products, Son's executed full-CMB correction bundle, Hoyt's signed vector and ZTF final posteriors remain unavailable. No runnable command here silently substitutes a different dataset or claims to fill those gaps. The completed deliverable is the bounded, independently checked investigation with explicit barriers, not a claim that every raw-data pipeline has been reconstructed.
