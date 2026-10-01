# Expansion background control

**Runnable example; no observational or paper-reproduction claim.**

This packet exercises Reproducible's execution and plotting path using
Irreducible's `background.evaluate`. It replaces the old repository's broad
computation with a deliberately small starting calculation. Matter density 0.3,
H₀ = 70 km/s/Mpc and eight redshifts are chosen controls. Irreducible owns the
equations, integration and numerical checks; this directory owns their use.

The [packet](experiment.json) states the question, assumptions, limitations and
required engine revision. The [request](request.json) contains the scientific
parameters. There are no data inputs, fitted parameters, uncertainty estimates
or random draws. This is not a supernova likelihood or a CMB calculation.

From the repository root, after building the pinned Irreducible revision:

```sh
uv run python scripts/run.py expansion-background --irred ../irreducible/target/debug/irred --name example
uv run --extra plots python experiments/expansion-background/plot.py results/expansion-background/example
```

The destination must be new. All receipts and rendered figures remain ignored.
At z = 0, the declared physical scale supplies H₀ and zero radial distance;
numerical qualification is reported by the engine. Successful execution proves
this consumer path ran, not that the cosmology describes observations.

## History

| Date | Change | Evidence and scope |
| --- | --- | --- |
| 2026-10-01 | Added initial packet | Native engine smoke example; no historical inference retained as a new result. |

After a substantive run, add a short dated finding here with the source commit,
request hash, engine build identity, gates and limitations. Archive exact full
receipts separately when a result is published; do not commit each rerun.
