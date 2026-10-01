# Reproducible

Run a scientific experiment, inspect its assumptions, and leave a useful account
for the next person.

Reproducible is the experiment workspace in a three-project workflow:

| Project | Responsibility |
| --- | --- |
| [Prospector](https://github.com/sprajs/prospector) | Find and review papers; hand off a versioned candidate investigation. |
| [Irreducible](https://github.com/sprajs/irreducible) | Implement and test the physics, numerical calculations and scientific interfaces. |
| **Reproducible** | Design the selected experiment, run the compiled calculations, visualize outputs and keep a compact experimental history. |

New physics belongs in Irreducible. Questions, data lineage, parameter choices,
fitting orchestration and plots for a particular investigation belong here.
A reviewed idea may remain blocked if its model, data or likelihood is missing.

## A small packet for each experiment

An experiment normally has four files under `experiments/<id>/`: a Markdown
account, `experiment.json`, an Irreducible `request.json`, and an optional plot
script. A Prospector-origin experiment also retains the small candidate JSON.
The account explains the claim, simplifications, checks, findings and limitations;
Git records its changes. We do not publish a file for every rerun.

Inputs live in ignored `data/`. Full engine receipts, tables, chains, simulations
and rendered figures live in ignored `results/` or `simulations/`. The public
packet identifies the evidence; a published scientific result should cite a
separate durable archive containing its exact inputs and full run records.

## Run the example

Python 3.11 or newer and [uv](https://docs.astral.sh/uv/) are needed for this
workspace. Build Irreducible separately using its
[build guide](https://github.com/sprajs/irreducible/blob/main/docs/getting-started.md),
at the commit named in the [example packet](experiments/expansion-background/experiment.json).
From this repository root:

```sh
uv sync --frozen
uv run python scripts/packet.py
uv run python scripts/run.py expansion-background --irred ../irreducible/target/debug/irred --name example
uv run --extra plots python experiments/expansion-background/plot.py results/expansion-background/example
```

The [background example](experiments/expansion-background/README.md) evaluates
chosen flat ΛCDM parameters at eight redshifts. It uses no observed data and fits
nothing. The runner verifies the packet, input identities and clean engine source
revision, then calls the current `irred run` interface. It records the actual
build, executable, request and receipts. Each destination must be new; failed
attempts remain local and visible.

This first runner handles one compiled request per packet. Irreducible's proposed
multi-step recipe runner remains proposed. More involved experiments can add a
small, reviewed orchestration script using existing engine interfaces; keep
shared physical equations and numerical kernels in Irreducible.

## Real-data native reference control

[The LCDM baseline audit](experiments/lcdm-baseline/README.md) uses the installed
Irreducible C++ SDK on exact ordered 13-row DESI DR2 Gaussian BAO compression.
Its bounded experiment-specific controller compares distances, ratios and the
normalized density with independent Decimal quadrature/LDLT references and
retains failed attempts. It is a chosen supplied-drag massless-radiation control;
full Planck 2018 massive-neutrino/thermal-drag/CMB reproduction remains blocked.
Numerical success establishes neither a fit nor a posterior result. Follow the
packet's native command and exact source/SDK/data pins.

Subsequent engine work has passed a released-ladder matrix comparison and added
an explicit thermal-relic E/H provider and conditional estimator variance.
The [baseline account](experiments/lcdm-baseline/README.md#engine-prerequisite-progress)
separates that progress from this packet's original massless model and build.
A newer engine does not automatically replace a pinned SDK or change an earlier
run's scientific identity.

[The released ladder experiment](experiments/released-ladder/README.md) now
executes the unchanged 3492-row, 47-column compact SH0ES design through the
installed native library. It verifies all 37 host columns against the pinned
primary Cepheid table and compares the full fit, H0-related contrast variance
and declared synthetic constraint-row responses with independent QR/SVD.
The source MCMC's fixed coordinate, unresolved nuisance identities and
paper/release selection discrepancy remain explicit. Its narrower controller
runs while the full paper-reproduction packet stays blocked; it supplies no
observational posterior or reconstructed systematic H0 uncertainty.

## Plots and notebooks

Keep visualization source next to its experiment and render into the run folder.
Start with a plain Python script. Use a Jupytext percent-format `.py` source when
cells improve exploration; keep executed `.ipynb` copies and their outputs in
ignored `notebooks/`. [Experiment conventions](docs/experiments.md) explain the
handoff and notebook workflow. No notebook service or automatic paper importer
is installed.

## Historical work and data

The old supernova/cosmology code and generated results were removed from the
active tree during the October 2026 reset. The
[historical packet](experiments/legacy-supernova/README.md) retains a short account,
qualification limits and the exact previous Git snapshot. It guides future
experiments; the new example does not reproduce that campaign.

[Storage and restoration](docs/storage.md) describe the retained input identities
and ignored local archives. A fresh clone has source manifests, rather than the
bulk input bundle. Deleting current files has not reduced the remote Git history;
use a shallow clone when you only need the new workspace.

## Contributing

Use a branch and a reviewable PR, including for agent work. See
[CONTRIBUTING.md](CONTRIBUTING.md), [AGENTS.md](AGENTS.md) and the
[development workflow](docs/development.md). Lightweight CI checks packet
integrity, storage rules and runner behavior; it does not qualify a scientific
result. [BSD 3-Clause](LICENSE) covers original repository code and documentation;
third-party inputs retain their own terms.
