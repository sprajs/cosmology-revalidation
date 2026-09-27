Changing the CMB calibration nuisances or the linear supernova luminosity coefficient does not require a new cosmological background or CMB spectrum. The native Cobaya dependency graph confirms this for the declared CPL linear target. All six fast parameters caused zero theory recalculations, and cached versus uncached likelihood and prior vectors agreed exactly. The smooth01 target passed the corresponding five-nuisance checks. The smooth03 target has the same dependency structure; it was not separately timed.

The new sampler therefore places H0, ωb, ωc, logA, ns, τ, w and wa in the slow block, with A_planck, P_act, Tcal, Ecal, A_fg and the sampled linear coefficient ε in the fast block. Fixed parameters are omitted. Cobaya's triangular covariance transformation lets slow proposals move correlated fast parameters, while fast proposals leave every slow coordinate unchanged. The complete covariance is retained. No prior transformation, altered likelihood, dragging approximation or convergence relaxation is introduced.

One proposal per coordinate per cycle is the default: eight slow-block and six fast-block proposals for the linear target. At the cubic-reference point, changing the whole slow block cost a median 0.04663 CPU seconds; changing the whole fast block cost 0.01639 seconds. Using these measured costs with four independently seeded, 20,000-attempt correlated-Gaussian controls gave:

| Fast-block oversampling | Minimum slow-coordinate ESS per modeled CPU, relative to one block | Minimum fast-coordinate ESS per modeled CPU, relative to one block |
| --- | ---: | ---: |
| 1 | 1.43× median; 1.11–1.91× across seeds | 1.39× median |
| 4 | 0.80× median; 0.36–1.00× across seeds | 3.83× median |

These comparisons use whitened coordinates and the minimum ESS within each block, with the first quarter discarded solely for the performance diagnostic. They are **a known-Gaussian calculation with modeled cosmology costs**, not a measurement of cosmological sampling efficiency. Whole-block timings are used; the faster timings of individual nuisance changes would understate the cost of changing several nuisances together. Gaussian controls retain all rejected states and their source arrays. Factor 4 spends more effort on fast directions than these costs and the slow-direction objective justify.

Six real CPL-linear pilots, two seeds and 400 attempts per layout, separately checked execution. One block required 351/373 theory recalculations, unit-factor blocks 229/232, and fourfold fast oversampling 101/102. Their slow-coordinate effective sample sizes after the short diagnostic discard were only 1.38–5.54. Four expensive native-CAMB fallbacks also made their timing ratios uneven. They cannot establish stationarity, mixing, posterior estimates or an empirical cosmological efficiency gain. The final invocation completed within its four-spectrum cap. An earlier invocation saved its first pilot and then failed during relative-path serialization; its source, result and draws are preserved under `.work/unified-cosmology/inference/blocked-benchmark-01`. The path handling was fixed before replay, making eight native calls across both invocations. No capped evaluation was converted into a rejected point.

`modern_blocked_sample.py` is a separate entrypoint for fresh runs. It defaults to fast factor 1, proposal scale 1.6, R−1<0.005 and the existing interval-convergence threshold 0.05. Automatic speed measurement and dragging remain off; all repeated states are retained without oversampling thinning. Proposal learning, finite-support initialization, MPI metadata protection and the existing convergence and native-correction requirements remain. Each manifest binds the sampler, imported initialization helper and metadata helper, as well as the existing complete target identity. Resume requires identical target, helper bytes, settings and MPI size. Unknown parameter dependencies fail until audited.

The optional `--gpu` flag selects the separately identified `modern_gpu` factory and uses a distinct run directory. It is off by default. GPU likelihood closure and the chosen surrogate's holdout validation are separate launch requirements. The CPU cubic-reference benchmark does not establish an optimum for a quartic/GPU implementation. Both remain numerical proposal chains, requiring independent convergence and native-CAMB correction before scientific use.

The configuration can be inspected without starting a run:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/modern_blocked_sample.py \
  --model cpl --evolution linear --surrogate PATH_TO_PINNED_MODEL --describe
```

The configuration validator requires no likelihood or spectrum calls and checks 16 model/evolution/factor combinations, the native proposal's zero slow-coordinate change during fast steps, and the eight planned linear/smooth01 MPI starting points:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/blocked_sampler_validate.py
```

The preserved performance experiment is repeatable in new ignored directories:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/blocked_benchmark.py \
  --surrogate .work/unified-cosmology/inference/surrogate/cubic-0509.npz \
  --work .work/unified-cosmology/inference/blocked-benchmark-replay \
  --output .work/unified-cosmology/inference/blocked-benchmark-replay.json

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/blocked_benchmark_control.py \
  --surrogate .work/unified-cosmology/inference/surrogate/cubic-0509.npz \
  --work .work/unified-cosmology/inference/blocked-control-replay \
  --pilot-results .work/unified-cosmology/inference/blocked-benchmark-replay.json \
  --output .work/unified-cosmology/inference/blocked-control-replay.json
```

The first report's individual-parameter cost model is superseded for interpretation by the second report's whole-block timing model. The second report reuses the first report's hashed Gaussian draws, so this correction does not select or rerun stochastic outcomes. Source review used the installed Cobaya model dependency analysis, `check_blocking`, native `BlockedProposer` and its Metropolis implementation; their identities are preserved in the result records. No live sampler or scientific target source was edited.
