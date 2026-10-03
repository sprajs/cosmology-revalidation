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

From the repository root, after building the pinned Irreducible revision and
creating the [plot environment](../../README.md#run-the-example):

```sh
uv run python scripts/run.py expansion-background --irred ../irreducible/target/debug/irred --name example
.work/plots/bin/python experiments/expansion-background/plot.py results/expansion-background/example
```

The destination must be new. All receipts and rendered figures remain ignored.
At z = 0, the declared physical scale supplies H₀ and zero radial distance;
numerical qualification is reported by the engine. Successful execution proves
this consumer path ran, not that the cosmology describes observations.

## History

| Date | Change | Evidence and scope |
| --- | --- | --- |
| 2026-10-01 | Added initial packet | Native engine smoke example; no historical inference retained as a new result. |
| 2026-10-01 | Executed clean committed source | Eight evaluations accepted under the numerical contract; z = 0 returns H₀ = 70 km/s/Mpc and zero radial distance. Inference not applicable; interpretation unqualified. |

After a substantive run, add a short dated finding here with the source commit,
request hash, engine build identity, gates and limitations. Archive exact full
receipts separately when a result is published; do not commit each rerun.

The recorded smoke execution used Reproducible source
[`a8e765515c6343f983c53697c7fea8f79653de9b`](https://github.com/sprajs/reproducible/commit/a8e765515c6343f983c53697c7fea8f79653de9b)
with a clean working tree and Irreducible source
[`ef4cb694addac1ca25c332253cff3fbecd0cbbce`](https://github.com/sprajs/irreducible/commit/ef4cb694addac1ca25c332253cff3fbecd0cbbce).
Its request SHA-256 is `497bbe99b43ebd16b287de2b45b19423058cf28193d24a17716b3c847d5ab897`; the actual engine build ID is
`cef351d8ef3f449e2e7c309dcdb41bd5e30e024323c0b36edfeae28acb709aed`. Full local evidence is in ignored
`results/expansion-background/committed-check/`. This control run establishes
consumer execution only; no observations or historical results were reanalysed.
