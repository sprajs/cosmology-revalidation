# Proposal cost and the interpolation boundary

A fixed-seed diagnostic supports reducing the proposal scale while preserving the spectral approximation, its boundary and the eventual exact cosmological target. It does **not** measure posterior convergence or effective samples per second.

The test uses the native Cobaya random-direction/radial-mixture proposal and a snapshot of its 13-parameter covariance. It cycles through six accepted parameter vectors from three early chains. The same 192 proposal displacements are rescaled in each comparison; these are paired sensitivity calculations, not independent Markov chains.

| Proposal scale | Prior rejection | Fast calculation | Exact fallback required |
|---|---:|---:|---:|
| 2.4 | 31 | 156 | 5 |
| 1.6 | 17 | 174 | 1 |
| 1.2 | 10 | 182 | 0 |

Every fallback is caused by exceeding the existing maximum absolute whitened coordinate of 4. None of the tested polynomial spectra is nonfinite or nonpositive in the guarded TT, EE or lensing channels. The largest excursion is 4.727. At scale 2.4 the boundary crossings are in the acoustic-angle direction (three), cold-dark-matter density (one) and w₀ (one).

Nine ordinary uncached fast posterior evaluations have a warm median of **0.0308 seconds** in the curated rerun, compared with **0.0441 seconds** in the original diagnostic. Both timing records are retained; host contention changes wall time. No full native CAMB spectra were calculated: the diagnostic explicitly disables `camb.get_results` in its own process and checks the surrogate's exact-call count remains zero.

If an exact fallback costs 30 seconds, the observed fractions imply about 0.81, 0.19 and 0.031 seconds per attempted proposal at the respective scales. This is a conditional cost scenario, conservatively charging prior rejections a complete fast evaluation. It is not a measured chain average. Five tail events are too few to estimate the fallback rate precisely, and zero events at scale 1.2 does not establish zero probability. The calculation also does not tell us whether the expensive proposals would be accepted.

Scale **1.6** is therefore a reasonable throughput trial without changing the target. Smaller steps may increase autocorrelation, so the choice must ultimately be judged by convergence and effective samples per wall time. Enlarging the interpolation boundary would instead require additional approximation validation; these proposal counts do not justify it.

## Repeat the diagnostic

The [compact result](../results/inference/proposal-throughput.json) includes the six numerical centres, complete proposal covariance, fixed seed, model/source identities, classification counts and the historical result/source hashes. Generated chain snapshots and all 576 proposal records remain in ignored storage. To replay the published numerical design without the original chain files, prepare the recorded modern environment and `cubic-0509.npz`, then run from the repository root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/proposal_throughput.py \
  --frozen-design studies/unified_cosmology/results/inference/proposal-throughput.json \
  --model-file .work/unified-cosmology/inference/surrogate/cubic-0509.npz \
  --output .work/unified-cosmology/inference/proposal-throughput-repeat \
  --summary .work/unified-cosmology/inference/proposal-throughput-repeat/summary.json
```

The output directory must be new. This frozen-design route was executed and reproduced all classification counts; its timing varied as expected. To diagnose a different chain snapshot, replace `--frozen-design` with `--chain-folder PATH`; the script records source and snapshot hashes before evaluating the proposals. It never modifies the input chains, covariance, model or active likelihood code.
