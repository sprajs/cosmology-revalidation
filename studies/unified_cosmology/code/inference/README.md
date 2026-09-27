# Shared supernova, BAO and CMB inference

Scientific definitions and interpretation are in [joint inference](../../notes/joint-inference.md). The [study overview](../../README.md) distinguishes qualified posterior measurements from fitted points, failed convergence checks and reconstructed published results.

## Inputs and environments

Run from the repository root. Prepare the [survey](../survey_selection/README.md) and [external-probe](../external_probes/README.md) inputs first. Their downloads, pinned third-party implementations and environments live in ignored `.work/`. For Pantheon+, run `data_variants.py` with the primary external-probe Python environment. It uses corrected apparent magnitudes, a free absolute offset and zHD > 0.01; it does not introduce a Cepheid absolute calibration.

The normalized supernova interface contains `zHD`, `zHEL`, `MU`, `covariance` and `precision`. The numerical `MU` column can be a corrected apparent magnitude or a released distance modulus: the unknown constant is integrated out. The complete release covariance determines the likelihood. Classification weights and plotting uncertainties are not applied a second time. Dovekie, Pantheon+ and DES3YR are alternatives, not independent factors.

Use `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1` for concurrent workers unless deliberately benchmarking a different allocation. CAMB is compiled numerical code; fixed linear-response compression accelerates lensing without changing the calculation. Extra threads in every concurrent process can reduce throughput.

## Numerical checks and late-time fits

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.venv/bin/python \
  studies/unified_cosmology/code/inference/validate.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.venv/bin/python \
  studies/unified_cosmology/code/inference/late_geometry.py --validate

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.venv/bin/python \
  studies/unified_cosmology/code/inference/late_geometry.py \
  --model lcdm --steps 20000
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.venv/bin/python \
  studies/unified_cosmology/code/inference/late_diagnostics.py \
  .work/unified-cosmology/inference/late-geometry/lcdm-none-dovekie \
  --output studies/unified_cosmology/results/inference/late-lcdm-none-dovekie.json
```

Models are `lcdm`, `wcdm` or `cpl`. `--evolution linear` adds a free magnitude drift; `smooth01` and `smooth03` analytically integrate four shape coefficients with declared Gaussian scales. They test assumptions about luminosity evolution; they are not measured age corrections. `--no-sn` gives a separate BAO-only target. Existing chains require `--resume`, which extends the recorded chains and appends run history.

Four interacting-walker ensembles are independently seeded. Diagnostics compare the same fixed walker index across independent ensembles and report the worst convergence result; they do not pretend that walkers within one ensemble are independent chains. The broad CPL prior has a genuine low-matter-density secondary ridge. `late_tail.py`, `late_nested.py` and `late_nested_refine.py` provide separate profile and nested-sampling checks. Their designs, failed initial comparison and isolated dependency lock are retained.

## Direct CMB calculation

```bash
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 CLIPY_NOJAX=1 \
  .work/unified-cosmology/external-probes/.venv/bin/mpiexec -n 4 \
  .work/unified-cosmology/external-probes/.venv/bin/python \
  studies/unified_cosmology/code/inference/run.py --model cpl
```

This is the separately declared full Planck 2018 target. `modern_run.py` supplies the contemporary cropped Planck + ACT + SPT configuration. Both share one cosmology with supernovae and BAO. The sound horizon is calculated by CAMB, not freed independently. `ExpansionDiagnostics` measures q and jerk from that provider's expansion history.

The contemporary target is not labelled an exact author reproduction. Its calibration and dependence assumptions are explicit. CAMB's analytic CPL interface restricts w₀+wₐ≤0; the late-time distance calculation has no such restriction. A likelihood reference point is not a posterior measurement.

## Accelerated contemporary calculation

First run the frozen training acquisition in the [external-probe guide](../external_probes/README.md). The 512 attempted training points yield 509 finite spectra; three fail the declared physical prior. The separate 96 attempted holdout points never enter the fit.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/surrogate_fit.py \
  .work/unified-cosmology/inference/surrogate/cubic-0509.npz
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/surrogate_check.py \
  .work/unified-cosmology/inference/surrogate/cubic-0509.npz \
  --output studies/unified_cosmology/results/inference/surrogate-holdout.json
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CLIPY_NOJAX=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/mpiexec -n 4 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/modern_sample.py --fast-lensing --proposal-scale 1.6 \
  --surrogate .work/unified-cosmology/inference/surrogate/cubic-0509.npz
```

The training fit refuses to overwrite an existing model. The spectral model is only a numerical proposal; CAMB supplies exact background distances, and points outside its interpolation envelope use native CAMB. `modern_fast.py` separately applies an exact matrix reassociation of the lensing response, independently checked at physical and synthetic spectra. This operation is not part of the spectral approximation.

The sampler prints its output directory, including a prefix of the actual model-file hash. Set `chain_dir` to that printed directory before the following commands. The published-chain mean and covariance initialize sampling only; they are not prior information. Source, actual likelihood bytes, covariance, priors and numerical model are hashed in each rank's manifest.

Use separate invocations with `--evolution linear` and `--evolution smooth01` for independent luminosity-sensitivity targets. Their output directories and manifests are distinct. A preliminary importance comparison of the baseline samples failed the required effective-sample-size gates; those provisional weighted intervals are not measurements. The [sensitivity note](../../notes/luminosity-sensitivity.md) specifies the integration, overlap checks and remaining broader `smooth03` alternative.

Per-rank JSON manifests also provide the parameter configuration for diagnostics. Cobaya 3.6.2 was observed to race when different MPI ranks wrote its common YAML metadata; one baseline metadata file contains a trailing fragment. Its chain files and immutable rank manifests remain separate. `mpi_metadata.py` restricts the three common metadata output streams to rank zero while preserving all MPI communication and per-chain output. The guard passes a four-rank, 12-round write check. It changes no sampled density. New runs use the guard; the original baseline metadata is retained as observed and is not used to infer its parameter list.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/diagnostics.py "$chain_dir" \
  --output "$chain_dir/diagnostics.json"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CLIPY_NOJAX=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/exact_correction.py "$chain_dir" \
  --diagnostics "$chain_dir/diagnostics.json" --points 2000 --workers 4 \
  --output studies/unified_cosmology/results/inference/modern-cpl-exact-correction.json
```

Exact correction refuses a failed diagnostic, changed chain cohort, changed target or fewer than 2,000 selected points. Each selected row is re-evaluated with both native CAMB and the frozen proposal, including a check against the stored sampling density. All failures and all **untrimmed** importance weights are preserved. Pareto smoothing diagnoses overlap but does not alter reported weights. Required gates are independent-chain Rhat ≤ 1.01, bulk/tail ESS ≥ 400, importance ESS ≥ 400 and Pareto k < 0.7. Each chain must retain importance ESS ≥ 50 and 5–50% of the pooled weight. Weighted chain means must lie within 0.3 pooled posterior standard deviations of the pooled mean, and estimated mean Monte Carlo errors must remain below 0.1 standard deviations for 10, 20 and 40 batches per chain. Each batch contains at least ten selected points. These explicit stability thresholds were fixed before any cosmological density correction. Constant weights have no Pareto tail to fit; a zero-weight chain fails. Finite diagnostics cannot prove the absence of a remote mode.

The exact calculation is substantially slower than proposal sampling. Cached point records allow resuming only with matching source, model, chains and selection identities. Each record seals its numerical payload, and the final summary records every native file hash; resuming cannot silently accept a changed payload or reseal a previously summarized file. If a larger or different native-point selection is needed, use a distinct `--name` such as `exact-correction-4000`; preserve the original attempt. A failed proposal or failed convergence check stays a failed numerical attempt; it is not reported as a cosmological uncertainty interval.

`measurement_summary.py --run CHAIN_FOLDER CORRECTION_SUMMARY --output REPORT.json` consumes a qualified correction. Repeat `--run` for separately sampled targets. It rechecks the parent diagnostic, all source/dependency identities, selected chain bytes and every native record, then recomputes the weight and stability gates. Its tables retain the weighted parameter covariance, finite scalar-prior and coupled CAMB-support diagnostics, and conditional acceleration/jerk fractions. Zero sampled tail events do not become a claim of certainty or a Gaussian significance.

## Interpreting qualified samples

The following consumers require the qualified native correction; they cannot turn an unfinished chain into a measurement.

```bash
python studies/unified_cosmology/code/inference/quantile_precision.py \
  --chain-folder "$chain_dir" --correction-summary "$correction_summary" \
  --output "$quantile_report"
python studies/unified_cosmology/code/inference/expansion_history.py \
  --chain-folder "$chain_dir" --correction-summary "$correction_summary" \
  --cache "$background_cache" --output "$expansion_report"
python studies/unified_cosmology/code/inference/luminosity_history.py \
  --chain-folder "$chain_dir" --correction-summary "$correction_summary" \
  --output "$luminosity_report"
```

Use the same modern environment and thread settings as the sampled target. [Expansion history](../../notes/expansion-history.md) reports past H(z), q(z), jerk and cosmic age from that target's CAMB background. Its bands are pointwise, and its age is conditional on the cosmological model. [Luminosity history](../../notes/luminosity-sensitivity.md) retains the full conditional uncertainty of analytically integrated brightness coefficients. A zero-width baseline brightness curve is an imposed assumption, not a measurement of no evolution. [Quantile precision](../../notes/quantile-precision.md) diagnoses uncertainty in estimated interval endpoints; its block-resampling spread is not an extra astrophysical error or a calibrated numerical confidence bound.

The [posterior precision audit](../../notes/native-posterior-precision.md) compares native CAMB accuracy settings at 32 predefined points from a qualified posterior. It is a bounded numerical sensitivity screen, not a second posterior or a proof of accuracy throughout the parameter space. The [conditional supernova quadratic check](../../notes/sn-predictive-check.md) verifies normalized likelihood accounting and compares the aggregate residual quadratic with its Gaussian replication distribution. Its posterior-averaged tail fractions are not calibrated frequentist p-values. Both checks retain failures and require the qualified parent before accessing observations.

The [joint-lensing comparison](../../notes/joint-lensing-bridge.md) separately replaces the current two lensing factors with the released ACT–Planck–SPT MUSE likelihood, including cross-experiment lensing covariance. Its baseline and extended ACT ranges share a new native spectral cache but retain separate weights and qualification results. This changes the SPT estimator as well as the covariance. [Probe omissions](../../notes/probe-omission.md) instead remove one declared component group using stored native likelihoods, with no further spectrum calculation. Both comparisons require adequate overlap; neither supplies missing cross-probe covariance or measures an age correction.

## Broader numerical support for luminosity alternatives

The initial linear and smooth-luminosity pilots visited regions outside the original interpolation envelope often enough that exact fallbacks dominated their runtime. They remain unqualified pilots. A new [numerical design](broad-spectral-design.json) keeps the original centre and physical priors, but expands the two conditional dark-energy directions of the training covariance by a factor of three in standard deviation. Its 768 training attempts and 128 independent holdout attempts use distinct fixed random seeds. Invalid requests are recorded without replacement. The polynomial degree and fitting method were chosen before those new outcomes.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/broad_spectral_training.py --workers 8
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/broad_spectral_fit.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/broad_spectral_check.py
```

The broader cubic model lives at `.work/unified-cosmology/inference/broad-spectral/surrogate-cubic-v1.npz`. Of 896 requested points, 689 training and 116 holdout evaluations were finite. Its untouched holdout includes likelihood errors as large as 140 in log likelihood, so it is retained as a numerical diagnostic rather than used to report cosmology. The largest errors occur at low native likelihood near the early-time CPL support boundary, but this does not establish adequate posterior overlap. The original baseline model stays unchanged. Any replacement numerical model requires a newly declared design, fresh independent validation, independent-chain convergence and native density-correction gates. Training and holdout spectra are not cosmological observations.

## Quartic numerical model

The fixed quartic model uses all 509 original and 689 broader finite training spectra, with 495 polynomial features. Its independent validation contains 128 new requests split equally between the original and broader coordinate distributions; 117 are finite. Earlier holdouts informed this repair and are explicitly not its independent validation. Run `quartic_spectral.py acquire --workers 6`, then `quartic_spectral.py fit`, `quartic_spectral.py check` and `quartic_spectral_audit.py` with the modern environment. Acquisition and fitting can run concurrently because the new holdout never enters the fit. Existing model outputs are not overwritten.

The mixed-width holdout has median log-likelihood error +0.014, with central 68% range [−0.758, +0.511] and much larger errors in some tails. The resulting model remains a proposal requiring native correction. Fresh independent sensitivity chains use:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CLIPY_NOJAX=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/mpiexec -n 4 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/modern_blocked_sample.py \
  --evolution linear --seed 272812 --gpu \
  --surrogate .work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz
```

The separately sampled `smooth01` target uses seed 272813. The same convergence, native density correction and weighted-stability gates apply. The GPU option requires the dependencies and validation below; omitting it retains CPU matrix evaluation with a different recorded numerical identity.

## Optional GPU calculation

`modern_gpu.py` dispatches only the spectral proposal's dense matrix multiplication to the GPU in float64. It inherits the original polynomial, exact background, interpolation-envelope checks, native fallbacks and unit conversions. It changes neither the CMB physics nor the required native-CAMB correction. The GPU identity records the wrapper source, CuPy/CUDA package versions, runtime, driver and device.

Install the pinned optional [GPU dependencies](gpu-requirements-lock.txt) with `uv pip install --no-deps --python .work/unified-cosmology/external-probes/.modern-venv/bin/python -r studies/unified_cosmology/code/inference/gpu-requirements-lock.txt`. The recorded installation added packages while leaving all 54 existing package versions unchanged. CuPy documents these bundled CUDA components in its [installation guide](https://docs.cupy.dev/en/v14.2.0/install.html).

`gpu_validate.py` compares the complete CPU/GPU likelihood at eight fixed interior points for each of three brightness models. Maximum log-posterior difference is 6.83×10⁻¹³. Its component timing includes feature transfer and result retrieval; it does not measure full-chain speed or validate the polynomial against native spectra. The independently validated [blocked sampler](../../notes/blocked-sampling.md) accepts `--gpu` for fresh runs with separate identities and directories. Existing CPU chains are not converted in place.
