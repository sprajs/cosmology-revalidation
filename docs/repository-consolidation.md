# Repository consolidation — 26 September 2026

`main` is the unified research branch. The three previous local branches formed
one linear history: `master` at `472f13f`,
`codex/independent-supernova-investigation` at `2f19e34`, and
`codex/des-raw-data-adversarial` at `792eb4c`. Each tip is an ancestor of `main`;
no conflicting histories or remote branches existed. Their redundant branch
names were retired after consolidation. All historical commits remain reachable.

The previously uncommitted research, numerical fixes, Markdown, experimental
plan, figures, configurations, source snapshots and result records are now
versioned together. Source snapshots under `runs/` deliberately stay at their
recorded paths because manifests and native refit procedures depend on them.
The [clean executable edition](../research_clean/README.md) is the recommended
entry point for reusable calculations; the dated original research is the full
evidence archive. Consolidation does not promote an exploratory or failed
experiment into a validated cosmology result.

## Storage and cleanup

Large acquired inputs, numerical arrays, native executables, environments and
execution logs remain local under explicit Git exclusions. A Git checkout alone
is still not the complete data bundle. Keep the whole working directory when
transferring it, or restore exact inputs using the acquisition and edition
manifests. The new [local-asset inventory](../catalog/consolidation-2026-09-26/retained-local-assets.json)
records hashes for newly excluded research assets; it supplements rather than
replaces existing acquisition manifests.

Cleanup removed superseded temporary validation copies, the old temporary
Python environment, regenerable bytecode/lint caches, two truncated downloads,
autoconf backup files and stale Git temporary objects. The complete replacement
PDF was hash-verified before removing the truncated bodies. Scientific failure
records, run inputs and native calibration builds remain available. The stable
clean-edition environment is now `research_clean/.venv`, installed from its
unchanged lockfile using `uv sync --frozen --python 3.12.14 --extra predictors`.

The [cleanup ledger](../catalog/consolidation-2026-09-26/cleanup.json) gives exact
paths, reasons and logical byte counts. These are not guaranteed physical disk
savings: package files can share storage. Scope was the supernova workspace and
its specifically identified temporary research artifacts, not unrelated sibling
projects. Git storage was repacked with pruning and reflog expiry disabled.

Acquired third-party repositories remain separate. The one locally modified
SNANA checkout has its exact base commit, header hash and reproducible
[patch](../catalog/consolidation-2026-09-26/vendor-patches/manifest.json) preserved
in the consolidation catalog. Build-tool dependencies and generated compiler
probes are excluded from research source coverage.

## Verification and recovery

- Every original branch tip is retained in the consolidated history. Exact
  original refs are in [branches-before.json](../catalog/consolidation-2026-09-26/branches-before.json).
- All 3,853 previously untracked source/document files checked by the
  [coverage audit](../catalog/consolidation-2026-09-26/source-coverage.json) are
  versioned; the 57 explicit exceptions are third-party dependencies or generated
  compiler probes. Previously ignored acquired archives are outside that count.
- All 552 files in the clean-edition manifest remain byte-identical. Its 19
  workflows completed in the new locked environment; all 19 result summaries
  match exactly. Across the replay products, 78 files match byte for byte; the sole
  remaining difference is the recorded runtime in `hst-geometry/stage_a.json`.
- The original cosmology and slab-attenuation numerical checks and conditional
  sampler convergence-gate check passed. All 304 Python files checked under the
  research script and clean-edition code directories parsed successfully.

The [validation record](../catalog/consolidation-2026-09-26/validation.json)
links these checks to retained replay ledgers. Duplicate replay outputs were
removed after comparison; the frozen reference products remain. This was not a
rerun of the entire archival native-refit, classifier, selection/BBC or inference
programme.

A verified bundle of the original three branches, the pre-consolidation tracked
patch and initial working-tree inventory are stored locally under
`.git/backups/consolidation-20260926/`. The bundle contains the complete original
branch histories. No remote was configured and nothing was pushed or published.
