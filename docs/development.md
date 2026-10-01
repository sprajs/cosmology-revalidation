# Development and publication

Read [AGENTS.md](../AGENTS.md) and [contributing](../CONTRIBUTING.md). This workflow
follows Irreducible's branch, checkpoint, review and verified-cleanup principles;
its native-engine CI matrix is not required for this experiment workspace.

## Checks

```sh
uv sync --frozen
uv run python scripts/check_repository.py
uv run python -m unittest discover -s tests -v
```

CI runs the packet/storage/link checks and consumer tests on Python 3.11 and 3.13
for PRs targeting `main`, pushes to `main`, and manual runs. It does not download
scientific data or rebuild Irreducible. Runner tests challenge admission,
source/receipt identity, failure retention and path handling using a test engine.
For a consumer change also exercise a real clean, pinned Irreducible binary and
inspect its gates. Use the example command in the README. Plot changes need an
actual render and visual inspection. Documentation changes need truthful examples
and valid local links, not a full scientific campaign.

## PRs

Agents use `codex/` branches. Preserve unrelated changes, stage coherent explicit
paths, review the staged diff, commit useful checkpoints and push them promptly
when authorized. A coupled redesign can be substantial; unrelated work belongs
in another PR. Explain dependencies on Prospector/Irreducible revisions and test
the integrated candidate after a dependency changes. Open a PR and report checks
and limitations accurately. Attach created PRs to the Codex chat.

The owner gives standing authorization for agents to review their own PRs and
merge ready changes without asking again. Do a deliberate review pass over the
complete final diff, fix actionable findings, resolve conversations and require
applicable green checks on the latest integrated head and up-to-date base.
No separate human reviewer is required. Record the review conclusion and actual
checks. Use a merge commit and verify the exact head before merging; perform the
merge directly once ready rather than relying on GitHub's auto-merge setting.
Direct `main` pushes, force pushes and wiki/external messages still require
specific authorization. Scientific qualification remains a separate gate.

## After an authorized merge

Verify the PR state, exact head and merge SHA with GitHub. Wait for post-merge CI
on that exact `main` commit. Inspect status/worktrees, fetch/prune, prove local
`main` is an ancestor of `origin/main`, then update it fast-forward-only. Before
`git branch -d`, prove the feature tip equals the verified PR head and is an
ancestor of fetched `origin/main`. Preserve dirty checkouts, extra commits and
branches used by active worktrees; never use force deletion. Remote branch
deletion is separate and depends on repository settings. Report deferred cleanup.

If `main` fails, stop new merges, preserve logs and fix or revert through a PR
under the same review and fresh-check gates. Do not bypass the failure or weaken
tests merely to make the status green.
