This diagnostic checks numerical likelihood sensitivity at **32 fixed points from an already qualified native-CAMB corrected measurement**. It compares the stored declared-accuracy likelihood with a new native calculation at identical physical and nuisance parameters, doubling only `AccuracyBoost`, `lAccuracyBoost` and `lSampleBoost` from 1 to 2. It retains the full CMB, BAO and supernova likelihood, including the BAO response to a changed drag sound horizon. The physical model, priors, data, lensing accuracy, lens margin and requested multipole ranges remain fixed. No spectrum emulator or acoustic-coordinate reinversion is used.

The parent must pass `measurement_summary.summarize_run`, including its independently recomputed chain, correction-weight and source-identity checks. Each of four chains contributes eight chronological strata. Within its existing native-correction points, positions are fixed at `floor((j+0.5)*N/8)` for j=0,…,7, after sorting by the original expanded-chain position. Selection uses neither the accuracy differences nor extreme importance weights. Original weights, their mass in the full parent and the full chain weight masses are retained as context only. Repeated physical points are retained and their distinct count reported.

For each point, the code records

\[
\Delta_i=\sum_k[\log L_{k,\,\mathrm{accuracy}\,2}(\theta_i)
                  -\log L_{k,\,\mathrm{declared}}(\theta_i)],
\]

with every component difference preserved. It checks that the inferred stored prior sum equals the new prior sum and that the likelihood and posterior differences account for one another. The unweighted mean of these 32 differences supplies an arbitrary common offset for centering. **It is not an estimate of a posterior expectation.** A common offset cancels from normalized weights, whereas variation across parameters can change inference.

The numerical screen, declared before evaluating qualified posterior points, flags centered total RMS above 0.05 log-likelihood units or a maximum centered absolute difference above 0.10. It also flags a component's centered maximum above 0.20 to expose cancellations. A total maximum of 0.10 would limit pairwise relative-likelihood changes to exp(0.20) among the tested points only. No threshold is placed on a purely common offset. These are diagnostic thresholds, not calibrated confidence guarantees or replacements for the existing posterior gates. Earlier fixed-reference/holdout accuracy comparisons were known when the screen was designed.

A background-only reconstruction at declared accuracy checks the stored rdrag, matter density and q/j values. It then compares H(z), the acoustic angle and recombination quantities with the higher-accuracy background. This adds no native spectrum calls. Higher-accuracy spectra, complete likelihood components, finalized numerical settings and separate CAMB maximum multipole/provider-array lengths are retained. Each native evaluation has an isolated log, preserving warnings without conflating warning counts with counts of invalid cosmologies.

Two accuracy levels and 32 deterministic points cannot establish asymptotic numerical convergence, cover every tail or rule out localized errors. The code performs **no 32-point posterior reweighting**, significance calculation or measurement qualification. A failed point remains failed; it is never replaced. A previous incomplete log prevents a silent retry, and the report distinguishes known native evaluations from incomplete attempts whose call count is unknown. Cached records require both prior manifest hashes and canonical payload checks.

The accompanying validation uses only synthetic arithmetic, parameter configurations and cache/qualification failures. It checks deterministic selection, unchanged CPU/GPU native configurations, signed density accounting, common-offset invariance, cancelling-component flags, independent CAMB parameter copies, and rejection of altered caches or unqualified parents. It performs zero native background, spectrum or likelihood evaluations:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/native_posterior_precision_validate.py
```

Once an actual parent qualifies, prepare and inspect its immutable selection without executing native calculations:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/native_posterior_precision.py \
  --chain-folder QUALIFIED_CHAIN_FOLDER \
  --correction-summary QUALIFIED_CORRECTION_SUMMARY \
  --cache .work/unified-cosmology/inference/native-posterior-precision/TARGET_NAME \
  --output studies/unified_cosmology/results/inference/native-posterior-precision-TARGET_NAME.json \
  --plan-only
```

For execution, repeat the command without `--plan-only`, optionally adding `--workers 4`. The default is one worker; the maximum is four, with one OpenMP/BLAS thread each and at most one native point evaluation per selected index. Existing completed records are checked and reused. Source, data, environment or parent changes require a new cache. A separate result name is required for each physical/luminosity target.
