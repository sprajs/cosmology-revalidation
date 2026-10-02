# Working in Reproducible

Reproducible runs experiments discovered by Prospector using software built in
Irreducible. Own the experiment design, input lineage, fitting/orchestration,
visualization and concise account here. Shared equations, physical models,
numerical kernels and their scientific tests belong in Irreducible. Read
[experiment conventions](docs/experiments.md), [storage](docs/storage.md) and
[development](docs/development.md) before changing the workflow.

## Start with a question

Read the packet and its history before running anything. A Prospector handoff is
a versioned candidate design, not a finished experiment. Retain its source
revision/path/hash and a small immutable JSON snapshot; describe the tested
claim, simplifications and information lost. Inspect the actual Irreducible
schema, compiled discovery and current capabilities. Missing physics or data is
a blocker, never permission to substitute an unrelated model and call it a
reproduction. Keep one current interface and update coupled callers together.

The built-in runner executes one request. Multi-step engine recipes are still a
proposal. More complex fitting may use a small experiment-specific controller
around the compiled interfaces. Do not rebuild the old Python physics engine
here or introduce a generic workflow framework without a concrete consumer.

## Evidence and computation

Distinguish observations, fitted summaries, assumptions and synthetic controls.
Record units, axes, frames, calibration, masks, selection, source versions and
dependence. Shared data ancestry is not independence; unknown covariance remains
unknown. Verify hashes before reuse. Parameters, seeds and numerical policies
must be explicit. Preserve failures and changed-source identities; never edit a
receipt or loosen a tolerance merely to obtain acceptance.

Pin the intended Irreducible revision and rebuild from clean source. The actual
build identity and executable hash belong in the full run record; a source
commit alone is not an executable identity. Keep execution, numerical,
inference and interpretation gates distinct. A CI pass or exit 0 cannot establish
an observational result. Record the checks actually performed and limits.

Use bounded runs, coordinate resources and respect the shared four-job compute
budget when building Irreducible. Use `gpt-6.1-sol` for routine delegated work
when delegation is requested. Inspect scientific comparisons before promoting
findings; qualified claims need the appropriate independent evidence.

## Keep the history small

Normally use README, packet JSON, request JSON and optional plot/controller
source. Append meaningful dated findings and corrections to that packet's
Markdown; do not generate a public receipt file for every run. Git keeps the
edit history. Cite an external durable archive for exact published inputs and
full receipts. Failed and null findings deserve concise entries too.

`data/`, `results/`, `runs/`, `simulations/`, `downloads/`, `notebooks/` and `.work/`
are ignored local storage. Keep plots and notebook sources with their experiment;
write generated figures and executed notebooks locally. Do not force-add
ignored output, third-party archives or downloaded papers. License/publication
permission must be checked before redistributing inputs. Preserve original
inputs; only remove inventoried disposable outputs within authorized scope.

## Contributions and PRs

Follow [CONTRIBUTING.md](CONTRIBUTING.md). Work on a `codex/` branch, stage explicit
coherent paths and inspect the staged diff. Push useful checkpoint commits
promptly when publication is authorized, including before opening a PR or
handover. Open a coherent PR with changes, actual checks and scientific limits.
Separate unrelated work; a coupled migration can be large. Preserve other
people's uncommitted work and coordinate one integration owner.

The owner gives standing authorization to create PRs, review our own changes in
a separate deliberate review pass, and merge ready PRs without asking again.
Inspect the complete final diff, fix actionable findings, run the applicable
checks and require green CI on the latest integrated head and up-to-date base.
Resolve review conversations and verify the exact head before merging. A
separate human reviewer is not required; scientific independence and qualification
still need their own evidence. This workflow applies across Irreducible,
Reproducible and Prospector. Direct pushes to `main`, force pushes and external
messages still require specific authorization; never bypass a failed check.

After an authorized merge, verify head and merge SHA, wait for CI on that exact
`main` commit, fetch/prune and update local `main` fast-forward-only. Prove the
branch is an ancestor of fetched `origin/main` before deleting it with `git branch
-d`. Preserve dirty checkouts, extra commits and branches in active worktrees;
never force-delete. If `main` fails, stop new merges and repair or revert through
the same review/check gates. Report deferred cleanup and remaining blockers.

## Cross-project baseline handoff

For `experiments/lcdm-baseline`, the reviewed experiment-specific native
controller is the only native packet route. Prospector owns the immutable source
map/candidate; Reproducible owns the exact data axes, experiment/reference
comparisons and full attempt records; Irreducible owns all production physics.
Freeze typed native input/output, source/SDK identities, domain and downstream
error allocations before implementation. Coordinate one owner per project and
the four-job total; a native consumer build uses one job. Keep full Planck
massive-neutrino/thermal-drag/CMB blockers distinct from the runnable supplied-drag
massless control. Never substitute a runnable variant for a faithful full-model
claim. Review the candidate hash and consumer revision again after a handoff
changes; do not infer a CLI operation from an experiment name.

Record later engine prerequisites separately from a packet's executed findings.
Native thermal-relic E/H now has a distinct native thermal distance/ruler
consumer with explicit physical-density/species mapping; use its guide and
actual pinned SDK rather than inheriting qualification from another model. Also,
fixed-design estimator variance or a synthetic H0 sampling law does not supply
observational calibration uncertainty. An identified released contrast can be
tested without guessing the remaining host/nuisance axes. Preserve the original
SDK/model/request pins and attempts; changing them requires an explicitly
reviewed new experiment/build identity and its own execution.

For `experiments/released-ladder`, use the reviewed bounded controller and pinned
SDK/source blobs. The unchanged full47 relative Gaussian profile is distinct
from the released MCMC box support and its fixed coordinate. Check the primary
rounded-table host join and retain unresolved axes and selection discrepancies.
The source-fixed44 target uses active46 original coordinates and a separate
closed source-box support check. Its unboxed profile and conditional source-named
constraint-mean shifts are not a normalized posterior, a box optimizer or a
combined systematic uncertainty. Preserve unmatched source rows and selection/
constraint discrepancies; held-out qualification needs resolved event and
cross-covariance lineage. Never turn synthetic constraint-row perturbations into
a physical calibration uncertainty or replace a historical SDK pin without a
reviewed new identity.

Use the dedicated [SN observer/passband controller](experiments/sn-observer-passband/README.md) for its frozen conditional slice. Keep the faithful candidate blocked until event/frame/reduction and calibration-law gaps close; a passing distance or optical integral does not qualify a released likelihood.

The [current campaign](experiments/lcdm-campaign/README.md) is a diagnostic
source/dependency snapshot serving Irreducible's sole active roadmap, not a
second development plan. Its bounded thermal DESI controller has its own
reviewed consumer/build identity and keeps the full-model packet blocked.
Require exact request, model/order, native policy/status and reference identities
before comparison. Retain high-precision reference differences at their declared
precision, including independent refinement axes and tail allocations; binary64
conversion must not erase a failed reference gate. Use a new immutable attempt
after source corrections, preserving the earlier receipt and its limitations.
