# Optical convergence passes; infrared score precision remains insufficient

The separate dense-metric attempt completed all eight optical chains without divergences and passed every original sampling threshold. The subsequent six-measurement synthetic J/H prediction failed its predeclared Monte Carlo precision requirement. This establishes a usable optical sampler for this one engineering case, but does not qualify the infrared comparison. No observed light curve was fitted or observed infrared value unsealed. The remaining eleven synthetic datasets have not been run.

The [preserved result](../results/bayesn-synthetic-dense99-case00.json) reports:

| Quantity | External-distance arm | Broad-distance arm |
|---|---:|---:|
| Complete chains × retained draws | 4 × 1,000 | 4 × 1,000 |
| Maximum rank R̂ | 1.00465 | 1.00372 |
| Minimum bulk / tail ESS | 1,070 / 1,190 | 1,350 / 1,292 |
| Minimum energy BFMI | 0.954 | 1.020 |
| Divergences / maximum-depth transitions | 0 / 0 | 0 / 0 |
| Arm wall time | 74.37 min | 93.73 min |
| Joint synthetic J/H log predictive density | −16.0830 | −15.7570 |
| Numerical log-score MCSE | 0.1081 | 0.1109 |
| Predictive likelihood-weight ESS, out of 4,000 draws | 89.8 | 78.8 |

The external-distance minus broad-distance log-score difference is −0.3260, with numerical MCSE 0.1549. That exceeds the fixed 0.05 limit. These errors measure finite-draw integration precision, not astrophysical uncertainty, and the displayed difference is not an accepted scientific preference. A relatively small subset of optical posterior draws carries the joint infrared likelihood. Passing diagnostics for the sampled coordinates therefore did not guarantee precise integration of this held-out quantity.

Both arms share the same 47-coordinate model, signed optical flux likelihood and latent priors except for the proper distance distribution. The external-distance arm combines the stipulated ΛCDM distance with the declared intrinsic gray scatter; the broad arm convolves a proper uniform distance-modulus prior with that same scatter. This simulation cannot independently validate the assumed cosmology, trained spectral model, population selection or an age correction.

The first synthetic case failed the required zero-divergence diagnostic with both proper distance priors. Its 5 and 14 divergent transitions remain part of the [original result](bayesn-synthetic-validation.md). A separate computational attempt uses a full 47-coordinate mass matrix and a target acceptance probability of 0.99. It retains the identical optical likelihood, proper priors, latent coordinates, spectral operator, integration resolution, four chains, 1,000 warmup iterations and 1,000 retained draws per chain. The maximum tree depth remains 10. No observed light curve is fitted and no observed infrared value is unsealed.

The [design](../code/bayesn-synthetic-retune-design.json) fixes fresh seeds before the new outcomes. The [settings and density validation](../results/bayesn-synthetic-retune-validation.json) constructs the actual NumPyro sampler, verifies its dense 47 × 47 metric layout, and finds bitwise identical old/new joint densities at twelve synthetic linear-kernel test points, including the broad distance-prior tails. Those checks require no SED evaluation or sampling. The original independent physical forward and local gradient checks remain prerequisites; this validation does not establish convergence.

The two proper-distance arms each ran on one CPU. Relative to the original 29.86 and 43.58 minutes, the retry took 74.37 and 93.73 minutes. The attempt finished within the external process-group limit of 7,200 seconds and 30-second termination grace. The internal arm timer is checked between complete chains and is not a substitute for that external limit. All outcomes remain in a separate directory, including each finished chain's adapted step size and ordered inverse-mass matrix.

The original diagnostic function and all its thresholds were applied unchanged. The optical chains and diagnostic files were hash-frozen before the synthetic infrared vector was copied into the retry directory. The original scorer sums the likelihood across the complete J/H vector within each shared latent draw before posterior averaging. Its conservative numerical-error estimate uses the larger of between-chain and contiguous-block estimates. The new [saved-record summary](../code/bayesn_synthetic_retune_summary.py) reconstructs the joint densities, flux covariance and score accounting from the prediction arrays, checking their source, input and optical-freeze hashes. It makes no model calls. No divergent draw was removed and no gate relaxed.

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

```sh
"$TASK_PYTHON" studies/infrared/code/bayesn_synthetic_retune.py prepare-score \
  --work .work/infrared/heldout-synthetic-dense99
"$TASK_PYTHON" studies/infrared/code/bayesn_heldout/score.py \
  --work .work/infrared/heldout-synthetic-dense99 --case 0
"$TASK_PYTHON" studies/infrared/code/bayesn_synthetic_retune_summary.py \
  --work .work/infrared/heldout-synthetic-dense99 --external-outcome completed \
  --full-output .work/infrared/heldout-synthetic-dense99/case00-final-geometry.json \
  --output .work/infrared/dense99-case00-summary.json
```

Further computation would need to improve precision of the joint predictive integral while preserving the fixed optical posterior and original error criterion. The present result does not justify beginning an observed-data comparison or claiming predictive coverage from a single synthetic case.
