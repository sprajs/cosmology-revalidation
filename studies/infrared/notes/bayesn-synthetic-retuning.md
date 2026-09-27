# A separate dense-metric synthetic sampling test

The first synthetic case failed the required zero-divergence diagnostic with both proper distance priors. Its 5 and 14 divergent transitions remain part of the [original result](bayesn-synthetic-validation.md). A separate computational attempt uses a full 47-coordinate mass matrix and a target acceptance probability of 0.99. It retains the identical optical likelihood, proper priors, latent coordinates, spectral operator, integration resolution, four chains, 1,000 warmup iterations and 1,000 retained draws per chain. The maximum tree depth remains 10. No observed light curve is fitted and no observed infrared value is unsealed.

The [design](../code/bayesn-synthetic-retune-design.json) fixes fresh seeds before the new outcomes. The [settings and density validation](../results/bayesn-synthetic-retune-validation.json) constructs the actual NumPyro sampler, verifies its dense 47 × 47 metric layout, and finds bitwise identical old/new joint densities at twelve synthetic linear-kernel test points, including the broad distance-prior tails. Those checks require no SED evaluation or sampling. The original independent physical forward and local gradient checks remain prerequisites; this validation does not establish convergence.

The two proper-distance arms each run on one CPU. The original arm times were 29.86 and 43.58 minutes; the new acceptance target can require more leapfrog steps. Dense matrix algebra at this dimension is small compared with the spectral evaluations, but the runtime improvement or deterioration is unknown. An external process-group limit of 7,200 seconds, with a 30-second termination grace, preserves a bounded attempt. The internal arm timer is checked between complete chains and is not a substitute for that external limit. Complete and partial outcomes are retained in a new directory, including each finished chain's adapted step size and ordered inverse-mass matrix.

The original diagnostic function and all its thresholds apply unchanged. Synthetic infrared prediction is allowed only after both arms pass and their optical chains and diagnostic files are hash-frozen. The joint infrared score then uses the original scorer and its original Monte Carlo precision criterion. The remaining eleven synthetic datasets are not launched automatically.

From the repository root, use the isolated environment described in the [original reproduction instructions](../code/bayesn_heldout/README.md). All output paths below must be fresh:

```sh
"$TASK_PYTHON" studies/infrared/code/bayesn_synthetic_retune_validate.py \
  --original .work/infrared/heldout-synthetic \
  --output .work/infrared/retune-settings-validation.json
"$TASK_PYTHON" studies/infrared/code/bayesn_synthetic_retune.py prepare \
  --original .work/infrared/heldout-synthetic \
  --validation .work/infrared/retune-settings-validation.json \
  --work .work/infrared/heldout-synthetic-dense99
timeout --signal=TERM --kill-after=30s 7200 \
  "$TASK_PYTHON" studies/infrared/code/bayesn_synthetic_retune.py run \
  --work .work/infrared/heldout-synthetic-dense99 \
  > .work/infrared/heldout-dense99-case00.log 2>&1
```

Only after both arms pass, the guarded `prepare-score` action freezes the optical records and copies the original **synthetic** held-out vector into the new attempt directory. It never copies observed photometry. The unchanged scorer can then be run on that new directory. A failed convergence or predictive-precision gate remains a failed result; no divergent draws are removed.
