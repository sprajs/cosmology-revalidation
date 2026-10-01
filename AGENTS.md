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

Independent review and passing applicable CI on the latest integrated head/base
are required before merging. Authorization to run an experiment does not itself
authorize a merge or an external message. Use standing merge authorization only
when the user has given it for this repository. Direct pushes to `main` and force
pushes require a specific request. No additional reviewer count is imposed.

After an authorized merge, verify head and merge SHA, wait for CI on that exact
`main` commit, fetch/prune and update local `main` fast-forward-only. Prove the
branch is an ancestor of fetched `origin/main` before deleting it with `git branch
-d`. Preserve dirty checkouts, extra commits and branches in active worktrees;
never force-delete. If `main` fails, stop new merges and repair or revert through
the same review/check gates. Report deferred cleanup and remaining blockers.
