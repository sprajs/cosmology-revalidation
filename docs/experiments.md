# Experiment conventions

## Admit a Prospector candidate

1. Find a reviewed `candidate_design` in Prospector. Preserve its commit, source
   path and byte hash; copy that small JSON into `experiments/<id>/candidate.json`.
2. Write the question, faithful claim, equations/conventions, proposed test,
   simplifications, falsification/check criteria and limitations in README.md.
   Keep measurements, fitted products and assumptions distinct.
3. Inspect the actual Irreducible build/schema. If the required model, data or
   likelihood is missing, record `status: blocked`, `execution: null` and concrete
   blockers. A blocked design is still useful history.
4. For a runnable packet, retain the accepted request unchanged and pin its clean
   engine revision. The validator checks the candidate/request hash and reviewed
   consumer revision. A changed request needs a newly reviewed design. Runtime
   scientific validation remains Irreducible's responsibility.
5. Acquire the identified inputs under `data/`; fill each input's bytes, SHA-256,
   role, source and semantics. Semantics should declare units/axes, frames,
   calibration, selection, covariance/dependence and unresolved gaps. Hash equality
   does not verify those declarations. Input paths in requests are relative to
   the repository root because the engine runs with that working directory.
6. Run a new attempt, inspect all qualification gates, then append a concise
   finding, negative result or failure to the packet's history.

The [schema](../schemas/experiment.schema.json) owns the packet fields. Use
`origin.kind: prospector_candidate` with repository `sprajs/prospector`, revision,
path, snapshot and SHA-256. Development smoke packets identify themselves
explicitly; they are not paper reproductions. Paper models are never replaced
by convenient ΛCDM/CPL requests and described as the original physics.

The built-in runner admits a single complete Irreducible request. It verifies
inputs before and after execution, records discovery/build identity, snapshots
the request/config and preserves stdout, stderr and the engine's immutable store.
It never executes commands from packet JSON, installs models or implements an
engine multi-step recipe language. Large controllers belong in small reviewed
experiment-specific scripts, sharing equation implementations through Irreducible.

## Keep a useful account

A packet normally contains README.md, experiment.json, request.json and optionally
plot.py. Add candidate.json for a literature handoff. Add a small analysis script
or methods section only when needed. Do not copy libraries, papers, full source
snapshots, tables or chains into the packet.

Use one history table or a few dated Markdown entries. Record the question/variant,
meaningful findings, source commit and dirty status, request/input hashes, engine
build ID, actual gates and the interpretation boundary. Distinguish a rerun at
the same design from a new scientific question. Correct a finding visibly and
retain its prior provenance in Git. Record failures and null findings too.

Execution, numerical acceptance, inference adequacy and physical interpretation
are separate. Report an exact, approximate, conditional or blocked reproduction
against a named target rather than a generic “validated” label. No new packet
inherits the historical manuscript's scientific checks.

For a published result, place the complete run store, exact inputs or lawful
acquisition routes, environment and plotting source in a durable research archive.
Record that archive's permanent identifier in the Markdown. The ignored local
store is convenient execution evidence; it is not a preservation service.

## Notebook choice

Plain Python is the default because it runs in CI and has a small, reviewable
diff. Use Jupytext percent-format `.py` files when cells improve exploration; the
example's plot.py includes cell markers. Keep plotting code here because
Irreducible intentionally has no visualization layer. Plot reported values rather
than quietly recomputing the scientific model in plotting code.

Jupyter/Jupytext are optional exploration tools, not runtime dependencies:

```sh
mkdir -p notebooks
uv run --extra plots --with jupytext jupytext --to ipynb experiments/expansion-background/plot.py --output notebooks/background.ipynb
uv run --extra plots --with jupyterlab jupyter lab notebooks/background.ipynb
```

In the notebook, call `plot(run_path)` with the saved run you want to inspect;
the command-line entry cell only runs in script mode. Move useful
edits back into the tracked `.py` source. Record the optional tool versions in an
archive if notebook rendering is part of the published method. Generated notebook
outputs, HTML, figures and widget state stay in ignored storage.
