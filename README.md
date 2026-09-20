# Supernova populations and cosmic acceleration

Independent local investigation, 20 September 2026. **Scientific analysis is complete; exact-reproduction barriers remain explicitly documented.** This supersedes the original acquisition-only README, preserved in Git commit `472f13f`.

Host age, dust, calibration and selection can materially change supernova cosmology. Standard released distances favour acceleration in the models tested. A Son-like age correction can reverse the present-sign inference, but the public evidence does not uniquely identify that extra correction or establish robust present deceleration.

Start here:

| Record | Purpose |
|---|---|
| [Final scientific report](docs/final-report.md) | Findings, competing evidence, conditional results and limitations |
| [Reproducibility guide](docs/reproducibility.md) | Pinned environment, calculation order, exact configurations and verification |
| [Current literature/data map](docs/literature-and-data-status.md) | Primary revisions through 20 September and exact/approximate/blocked reproductions |
| [Claim-and-evidence ledger](docs/claim-evidence-ledger.md) | Claims, counterevidence and missing discriminators |
| [Experiment register](docs/experiment-register.md) and [decision log](docs/decision-log.md) | Plans, changes, failed attempts and final status |
| [Independent falsification audit](docs/investigations/independent-audit.md) | Different numerical implementations and adversarial scientific checks |
| [First principles](docs/first-principles.md), [causal model](docs/causal-model.md), [correction ledger](docs/correction-ledger.md) | Mathematical assumptions, latent variables and already-applied corrections |
| [Unsent author-data requests](docs/author-data-requests.md) | Precise missing products and what each would decide |

The main inference starts at released corrected distances/covariances and BAO measurements. It is not an independent raw-photometry or CMB-map reconstruction. Measured fluxes, inferred ages, author simulations, corrected distances and reference chains remain labelled separately. New synthetic tests establish mathematical/algorithmic behaviour under their stated generators.

Original `papers/`, `data/`, `sources/` acquisitions were preserved. Derived tables live under `data/derived/`; scientific scripts under `scripts/`; outputs, configurations and manifests under `runs/`. Current source additions are dated and hashed. Large inputs, chains and covariance arrays remain local and are excluded from Git; **a Git-only checkout is not the complete scientific data bundle**. See the reproducibility guide before transferring this record.

The acquisition catalogues under `catalog/` and older `docs/research-map.md`, `data-inventory.md`, `open-gaps.md` and `replication-boundary.md` are historical records. Their counts and “future analysis” language describe the earlier collection, not the completed investigation.

No authors were contacted and no data or results were published or uploaded. Third-party source licenses remain with their authors.
