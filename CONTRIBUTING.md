# Contributing

Bring a scientific question we can test and an account someone else can inspect.
Experiment designs, careful data lineage, useful plots, failed tests and clearer
instructions are welcome, including agent-written contributions.

Prospector owns paper discovery and review. Irreducible owns shared physical
models and numerical calculations. Reproducible owns their use in a selected
experiment. Describe the exact paper/version, faithful claim, chosen simplified
test and what the simplification loses. Do not present a synthetic example or a
missing-capability sketch as a reproduced result.

Work on a branch and open a coherent PR. Include the problem, resulting behavior,
checks actually run and unresolved limits. Useful intermediate commits are
welcome; split unrelated changes rather than imposing a line-count limit. Keep
the code, config and explanatory documentation together. Agents use `codex/`
branches and follow [AGENTS.md](AGENTS.md).

Push meaningful checkpoint commits for backup. The owner authorizes agents to
review their own PRs in a separate deliberate pass and merge ready changes
without another confirmation. Fix actionable findings and require green
applicable checks on the latest head and up-to-date base. No separate human
reviewer is required. Direct `main` pushes and force pushes still need a specific
request. CI checks PRs and merged `main`; it does not replace review or scientific
evidence.

Keep bulk inputs, outputs, simulations, rendered figures and executed notebooks
out of Git. Publish short findings in the experiment's Markdown and cite a durable
archive for full scientific evidence. Preserve inputs and failures, and check
third-party terms before copying assets. See [experiment conventions](docs/experiments.md),
[storage](docs/storage.md) and [development](docs/development.md).
