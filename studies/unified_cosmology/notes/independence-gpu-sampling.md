# Frozen independence proposals for the GPU luminosity targets

A separate driver prepares new independence-Metropolis chains for the existing linear and smooth01 luminosity models. It reuses the validated CPU sampler, normalized Gaussian/Student proposal, rejection-weight accounting and convergence tests without editing their sources. This is a sampling change, not a new cosmological likelihood, luminosity prior or observed result.

The driver reconstructs each source with `modern_gpu.configuration` and requires exact equality of its complete GPU target identity. The full parameter/prior, likelihood and theory configuration equals the existing CPU configuration after replacing only the spectral matrix multiplication class. New manifests retain `gpu=True` and `fast_lensing=True`, so the existing exact-correction and reporting consumers select the correct backend.

Preparation copies complete chain rows from each of four source chains, expands their integer frequency weights, and trains on the middle 40% after discarding the first 50%. The final 10% is withheld. These source chains are provisional and are used only to learn a proposal. Neither their estimates nor the new proposal means are qualified measurements. Each proposal is frozen, has support over all real sampled coordinates, and is never adapted during the new run.

Production uses the unchanged explicit native CAMB initialization before evaluating its first sampled point. This addresses the recorded initial cold/warm numerical-state discrepancy. It consumes no proposal random draws and changes no physical prior. Configuration equivalence alone is not a performance or posterior-convergence claim. The separately recorded efficiency pilot below evaluates the unchanged GPU target; it does not qualify a posterior.

The configuration validator forbids CAMB spectra, backgrounds, model construction and likelihood evaluation. It checks both actual GPU source identities, unchanged priors, exact agreement with the CPU driver's sampler-options dictionary, independent initialization streams and rejection of incompatible parents. The original sampler's synthetic detailed-balance, rejection-weight and convergence validation remains applicable because that sampler is imported unchanged. Preparation alone does not start a new sampling campaign.

All four chains must pass the original multivariate and confidence-bound convergence conditions and the additional unchanged rank-Rhat/ESS checks. Exact native posterior correction and its overlap/stability gates remain mandatory. A maximum-sample or interrupted run is not reclassified as converged. The driver refuses resumption and existing output folders.

From the repository root, preparation and inspection are:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1
INDEPENDENCE_PYTHON=.work/unified-cosmology/external-probes/.modern-venv/bin/python
INDEPENDENCE_DRIVER=studies/unified_cosmology/code/inference/independence_gpu_run.py

"$INDEPENDENCE_PYTHON" "$INDEPENDENCE_DRIVER" freeze SOURCE_GPU_CHAIN_FOLDER NEW_PROPOSAL_FOLDER
"$INDEPENDENCE_PYTHON" "$INDEPENDENCE_DRIVER" describe NEW_PROPOSAL_FOLDER
```

The prepared proposal folders are `.work/unified-cosmology/inference/independence-proposal-gpu-linear-v1` and `independence-proposal-gpu-smooth01-v1`. Their snapshot identities, proposal hashes and configuration checks are recorded in `results/inference/independence-gpu-validation.json`. The linear proposal has 14 coordinates and the smooth01 proposal 13. Different default sampling/initialization seeds are fixed for each model in `independence-gpu-design.json`.

After allocation and any requested density/throughput pilot, a fresh four-rank run can use:

```bash
.work/unified-cosmology/external-probes/.modern-venv/bin/mpiexec -n 4 \
"$INDEPENDENCE_PYTHON" "$INDEPENDENCE_DRIVER" sample NEW_PROPOSAL_FOLDER NEW_CHAIN_FOLDER
```

The `sample` action is the only action that performs native initialization or launches chains. No existing chain, scientific target, covariance, or likelihood input is modified by this driver.


## Linear efficiency pilot

The fixed 400-request linear-model pilot used proposal/acceptance seeds 273280/273281 and one CPU thread plus one GPU context. It produced 364 finite candidates, 36 prior rejections and no other invalid target evaluations. Its acceptance rate was 179/399 = **44.86%** after initialization. The 400 occupied draws contained 180 distinct states; the longest hold was 19 draws. The raw independent-proposal importance ESS was 123.22, with maximum normalized candidate weight 4.84%. These short-run quantities diagnose proposal efficiency, not convergence, mode coverage or a cosmological result.

Evaluations took 225.42 seconds. Seven native fallbacks accounted for 192.92 seconds; the median requested evaluation took 0.0935 seconds. Fixed native initialization (30.32 seconds), model setup (36.72 seconds), input hashing and reporting are separate costs. Concurrent activity on the same host was not controlled, so this does not establish a CPU-versus-GPU speedup. The subsequent full four-chain run still needs every convergence and native-correction gate.

The unchanged independent saved-pilot validator reconstructs all 400 proposal densities, prior terms, Metropolis decisions and occupied-state holds. The maximum independent log-proposal discrepancy is 1.53×10⁻¹² and the prior-plus-likelihood-component closure error is 8.53×10⁻¹⁴. Full records remain in the ignored pilot folder; compact results are `results/inference/independence-gpu-linear-pilot.json` and `independence-gpu-linear-pilot-validation.json`. The source is `code/inference/independence_gpu_pilot.py`.

```bash
"$INDEPENDENCE_PYTHON" studies/unified_cosmology/code/inference/independence_gpu_pilot.py \
.work/unified-cosmology/inference/independence-proposal-gpu-linear-v1 \
NEW_LINEAR_PILOT_FOLDER --seed 273280 --acceptance-seed 273281

"$INDEPENDENCE_PYTHON" studies/unified_cosmology/code/inference/independence_pilot_validate.py \
.work/unified-cosmology/inference/independence-proposal-gpu-linear-v1 \
NEW_LINEAR_PILOT_FOLDER
```

Use explicit production seeds, distinct from pilot streams. The allocated linear production uses 273282/273283 (sampling/initialization); the frozen driver's default seeds are not used for that run.


## Smooth-model efficiency pilot

The separate fixed 400-request smooth01 pilot used proposal/acceptance seeds 273284/273285, one CPU thread and one GPU context. It produced 367 finite candidates, 33 prior rejections and no other invalid target evaluations. The first request was outside the prior; the second initialized the chain. Subsequent acceptance was 203/398 = **51.01%**. The 399 occupied draws contained 204 distinct states, with longest hold 11 draws. The raw independent-proposal importance ESS was 186.25 and the maximum normalized candidate weight 1.97%. These are proposal diagnostics, not a qualified posterior or evidence of complete mode coverage.

Evaluation took 351.83 seconds, of which 323.23 seconds (91.9%) was spent in the 11 requests that required native fallback. The median request took 0.0681 seconds. Fixed native initialization (38.24 seconds), model setup (47.68 seconds), hashing and reporting are separate costs. Host contention was not controlled. The efficiency evidence supports considering a separate longer run, but it cannot replace the four-chain convergence, native-correction and weighted-stability requirements.

The unchanged independent validator checked every one of the 400 requests, including the initial prior rejection, and reproduced all acceptance decisions and occupied-state holds. Maximum independent log-proposal disagreement was 1.76×10⁻¹²; prior-plus-likelihood-component closure was 1.14×10⁻¹³. Compact results are `results/inference/independence-gpu-smooth01-pilot.json` and `independence-gpu-smooth01-pilot-validation.json`. Full candidate records and frozen request arrays remain in the ignored pilot folder. No production chain is started by the pilot.

```bash
"$INDEPENDENCE_PYTHON" studies/unified_cosmology/code/inference/independence_gpu_smooth_pilot.py \
.work/unified-cosmology/inference/independence-proposal-gpu-smooth01-v1 \
NEW_SMOOTH_PILOT_FOLDER --seed 273284 --acceptance-seed 273285

"$INDEPENDENCE_PYTHON" studies/unified_cosmology/code/inference/independence_pilot_validate.py \
.work/unified-cosmology/inference/independence-proposal-gpu-smooth01-v1 \
NEW_SMOOTH_PILOT_FOLDER
```
