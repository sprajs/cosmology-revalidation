# A separately gated native-accuracy refinement screen

This implementation can test whether genuine accuracy-1 to accuracy-2 likelihood changes shrink at accuracy 3. **No physical refinement evaluation has been performed.** The [synthetic validation](../results/inference/native_precision_refinement_validation.json) checks the numerical accounting, stage gates, exact point selection, source bindings and immutable caches without constructing a cosmological model or calling CAMB.

The [design](../code/inference/native_precision_refinement_design.json) and [executor](../code/inference/native_precision_refinement.py) preserve the original 32 chronological strata, including repeated physical points. They first invoke the pure [precision-review consumer](../code/inference/native_precision_review_consumer.py), which freshly qualifies the original posterior and verifies the original screen, thermal review, native records, source/data/version identities, spectra and logs. A genuinely surviving total or component density-variation flag is required. A repaired passing screen, incomplete points or an unqualified parent cannot start this programme.

Execution has two separate stages:

1. **Replay:** evaluate the original first selected point at accuracy 2 in a fresh process. Compare every likelihood component, the prior sum, total log density, all stored derived quantities and all five native spectra with the original accuracy-2 record. Record the elapsed time. Component, prior-sum and total-density replay tolerances are 10⁻⁸; scaled derived differences are limited to 10⁻⁸ and each spectrum's maximum difference, divided by that spectrum's maximum absolute reference amplitude, to 10⁻¹⁰. This is a reproducibility gate, not a change to the original numerical-spread thresholds.
2. **Refinement:** only after the sealed replay and complete hash ledger pass, evaluate the same 32 points at accuracy 3, with one fresh process per point. Only `AccuracyBoost`, `lAccuracyBoost` and `lSampleBoost` change from 2 to 3. Physical and nuisance coordinates, likelihoods, data/covariances, priors, neutrinos, halofit, lensing conventions, multipole requirements and calibration remain fixed. There is no spectral emulator or acoustic-angle reinversion.

Each worker also reconstructs an accuracy-2 auxiliary background from a copy of the finalized native parameters, using the existing thermal adapter. It must match the original accuracy-2 H, rdrag, matter density and q/j values under the unchanged background tolerances. This protects the reference-background path without modifying either full native likelihood.

The accuracy-3 minus accuracy-2 log-likelihood differences use the original centered RMS limit 0.05, total maximum centered difference 0.1 and component maximum centered difference 0.2. Centering is over all 32 points, not separately by chain. Component tests therefore retain evidence of cancellations in the total. Prior and density-accounting failures remain fatal. The original unrepaired and thermal-reviewed accuracy-1 to accuracy-2 summaries are retained beside the new result.

Plans, per-point payloads, spectra, process logs, one-use attempt tickets and process claims are sealed and hashed. A previous log or claimed slot without a completed record is an incomplete attempt with an unknown call count; it is never silently retried. Failed replay blocks refinement. There are at most 33 requested native `model.logposterior` invocations in this campaign, with no explicit warmup call. This is an invocation count, not instrumentation of CAMB's internal solver. Completed records retain the separate auxiliary-background count; incomplete attempts remain visible. Existing source and parent bindings are rechecked before and after each worker and stage.

The [validator](../code/inference/native_precision_refinement_validate.py) tests LCDM and CPL configuration dictionaries, replay mutations, nonfinite backgrounds, constant offsets, cancelling components, and total density variation. Its 33-record synthetic cache exercises exact original IDs, no-call replay, failed/missing replay refusal, one-use claims, interrupted-attempt handling, parent/source deletion and numeric tampering. Parent qualification and native outputs in those cache fixtures are explicitly mocked; actual physical calculations are forbidden. It is implementation validation, not validation of accuracy 3.

After root review and a genuine qualifying parent screen exist, the staged commands are:

```sh
TASK_PYTHON=.work/unified-cosmology/external-probes/.modern-venv/bin/python
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1
"$TASK_PYTHON" studies/unified_cosmology/code/inference/native_precision_refinement.py prepare \
  --screen /absolute/path/to/original-32-screen.json \
  --thermal-review /absolute/path/to/thermal-review.json \
  --work .work/unified-cosmology/inference/native-refinement-new
"$TASK_PYTHON" studies/unified_cosmology/code/inference/native_precision_refinement.py replay \
  --work .work/unified-cosmology/inference/native-refinement-new
```

Review `replay-summary.json`, including its timing and every gate, before allocating the next stage. A numerical replay success alone is not a runtime forecast for accuracy 3. Refinement defaults to one worker and permits at most four:

```sh
"$TASK_PYTHON" studies/unified_cosmology/code/inference/native_precision_refinement.py refine \
  --work .work/unified-cosmology/inference/native-refinement-new --workers 1
```

A `report --stage replay` or `report --stage refine` command inspects the same cache without starting workers. Validation itself is reproducible using `native_precision_refinement_validate.py --output NEW_VALIDATION_PATH`; the installed execution guard requires the source-bound public validation record.

Passing would mean only that no large 2→3 variation was detected at these 32 points. It would neither repair accuracy 1 nor prove asymptotic convergence, represent a weighted posterior or justify cosmological shifts from 32-point reweighting. Any later accuracy-2 target correction must have its own identity, all 2,000 native evaluations and every original importance-weight/chain/stability gate. The current implementation makes no such correction and changes no active target.
