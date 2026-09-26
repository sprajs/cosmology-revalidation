Prepared, not executed. Root releases each stage separately with JSON fields `protocol_sha256`, `freeze_sha256`, `stages` containing the requested stage.

Protocol SHA: `2b4b1e21f77c50ec8de1d3abcec22118a1aa987d0322c5e4b1cf68572b96b8b6`
Freeze SHA: `3a38d43acfbe63710959459b5d623b5f02da07bee68d49dd9e15b0b7e7b64b9e`
Runner SHA: `a7d34c24c02ae9b9dd7518a6960d746e2b5ddf0fdd32123376907295ee29e456`

First stage only:

```sh
phase2/env-official/bin/python phase2/pt64/run.py generation --release /absolute/path/to/root-release-generation.json
```

After its independent gate review, use the same runner with `joint`, then `nir`, each with its own release. No separate summary command is needed: the NIR stage computes only the frozen estimands after all gates and full null identity pass.

Seed 26092672, fixed 64 attempts, no pooling/reseeding/refill. Two output-identity generation runs represent one set of 64 new draws. Exact original physical inputs are symlinked read-only. Generation parents exist and are writable. Generator registry is private and expected to change; the fitter's distinct registry stays empty. The already validated original-source copies/binaries remain unchanged. Preflight replays existing eight-event NIR/null checks and confirms unchanged function source and threshold logic without fitting.

13 planned fitting jobs: nominal joint12, joint9, actual local first-MINUIT ±2 starts; each true/fitted NIR arm12/9/D−0.2/D+0.2; null12. Every new event must pass. Mean/covariance results are withheld on failure; inspect the saved failure/attempt ledgers before any separately designed follow-up.
