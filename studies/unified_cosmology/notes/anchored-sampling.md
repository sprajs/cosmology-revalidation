# Dedicated inference with the released absolute calibration

The separate anchored target uses the complete Pantheon+SH0ES covariance for 77 calibrator and 1,580 noncalibrator light-curve rows. It replaces the Dovekie distance factor. It adds neither a separate Hubble-constant prior nor another Cepheid covariance. The [released calibration audit](calibration-interface.md) describes the retained covariance and its unresolved internal decomposition. This target preserves the CMB and BAO factors, physical priors and cross-probe independence approximations of the main measurement.

Fresh sampling provides a route if the separately tested [anchored importance bridge](anchored-bridge.md) lacks overlap. No anchored chains or cosmological likelihood evaluations have been run by this scaffold's validation. Either an explicitly audited anchored-bridge weighted mixture or a Dovekie-trained Gaussian/Student-t mixture can supply a normalized proposal: its training distribution is not treated as an anchored posterior. The unchanged Metropolis–Hastings sampler evaluates the actual anchored density, including the reverse/forward proposal-density ratio, at every trial.

The numerical proposal uses the frozen accuracy-one spectrum surrogate and background. The declared native target can use CAMB accuracy one or two. Accuracy two changes only `AccuracyBoost`, `lAccuracyBoost` and `lSampleBoost`; lensing accuracy remains four and all physical settings, likelihood factors and priors are retained. Its exact/proposal correction evaluates native accuracy two directly, without requiring a second accuracy-one correction cohort. Proposal-density reconstruction must still match the recorded chain, so a changed numerical library state cannot silently redefine the sampled density.

All four chains must terminate through the original sampler criteria: both consecutive multivariate R−1 statistics below 0.005, bound statistic below 0.05 and the sampler's independent rank gate. The checkpoint, convergence ledger, four terminal records and candidate ledgers are checked and hashed. The unchanged independent chain diagnostics then require rank Rhat ≤ 1.01 and bulk/tail ESS ≥ 400 on the common tail after discarding 30% of the saved post-burn chain. A run stopped at a finite sample budget cannot pass merely because marginal diagnostics look acceptable.

At least 2,000 stratified points undergo native correction. All original raw-weight, Pareto-tail, independent-chain and batch-stability gates remain in force. Rejected or nonfinite evaluations remain in the record; weights are never clipped. The dedicated measurement consumer reconstructs the selection and all weighted summaries from sealed native records, rechecks live source/data/version identities, and reports no measurement after any failed gate. Original generic consumers reject the new target and record schemas.

The synthetic validator checks the two model families at both native accuracy settings, unchanged priors and proposals, independent starting streams, the actual sampler settings, nonconstant integer chain-weight expansion, sealed cache replay, original sampler-stop gates, altered-record rejection and a complete 2,000-point synthetic measurement path. These tests establish the plumbing and arithmetic, not sampling efficiency, posterior convergence, native accuracy or coverage of unseen modes. An explicit initialization call is recorded for each correction worker; a fully cached replay launches no workers and makes no initialization call.

The proper execution order is:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1
PY=.work/unified-cosmology/external-probes/.modern-venv/bin/python
CODE=studies/unified_cosmology/code/inference
MPIEXEC=.work/unified-cosmology/external-probes/.modern-venv/bin/mpiexec

# Inspect and freeze the separate target before any launch.
$PY "$CODE/anchored_sampling.py" describe --model lcdm \
  --surrogate .work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz \
  --native-accuracy 2 --output .work/unified-cosmology/anchored-target-lcdm-a2.json

# This bridge report must already exist with a qualified parent and complete
# finite source/target closure. Failed child overlap is allowed for training only.
$PY "$CODE/anchored_proposal.py" \
  --bridge .work/unified-cosmology/inference/anchored-bridge-lcdm.json \
  --reference-proposal .work/unified-cosmology/inference/independence-proposal-lcdm-v1 \
  --output .work/unified-cosmology/inference/anchored-proposal-lcdm-v1

# Choose an available CPU explicitly; this starts no production chain.
$PY "$CODE/anchored_pilot.py" --model lcdm --native-accuracy 2 \
  --surrogate .work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz \
  --proposal-folder .work/unified-cosmology/inference/anchored-proposal-lcdm-v1 \
  --output .work/unified-cosmology/inference/anchored-pilot-lcdm-v1 --cpu 4

# Review the bounded efficiency diagnostic before scheduling four chains.
# Every output directory must be fresh.
"$MPIEXEC" -n 4 "$PY" "$CODE/anchored_sampling.py" sample --model lcdm \
  --surrogate .work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz \
  --native-accuracy 2 \
  --proposal-folder .work/unified-cosmology/inference/anchored-proposal-lcdm-v1 \
  --output .work/unified-cosmology/inference/anchored-lcdm-a2

$PY "$CODE/diagnostics.py" .work/unified-cosmology/inference/anchored-lcdm-a2 \
  --output .work/unified-cosmology/inference/anchored-lcdm-a2/diagnostics.json
$PY "$CODE/anchored_correction.py" .work/unified-cosmology/inference/anchored-lcdm-a2 \
  --diagnostics .work/unified-cosmology/inference/anchored-lcdm-a2/diagnostics.json \
  --points 2000 --workers 4 \
  --output .work/unified-cosmology/inference/anchored-lcdm-a2/native-summary.json
$PY "$CODE/anchored_measurement.py" .work/unified-cosmology/inference/anchored-lcdm-a2 \
  --correction-summary .work/unified-cosmology/inference/anchored-lcdm-a2/native-summary.json \
  --output .work/unified-cosmology/inference/anchored-lcdm-a2/measurement.json
```

These commands are a separate future execution path. They do not replace the active main measurement or authorize interpreting an unconverged anchored result. The released selection and absolute-calibration assumptions remain conditional, and the improper flat supernova-magnitude measure does not define normalized model evidence.

A separate `anchored_pilot.py` can screen a genuine bridge-trained proposal before committing four chains to it. Its 400 independent requests and acceptance uniforms are frozen before model construction. It keeps every rejected hold, evaluates the unchanged anchored numerical proposal, and reports acceptance, holding concentration and descriptive ESS without cosmological summaries. The pilot stops after four successful native fallbacks or 600 seconds checked after a completed call; an external timeout is needed to cap an in-flight calculation. Its one native initialization call is recorded separately. A stopped pilot remains a partial diagnostic and cannot satisfy any measurement gate.

The native-accuracy label in this pilot identifies the eventual correction target. The pilot itself tests proposal efficiency at accuracy one; it does not establish the adequacy or posterior overlap of accuracy two. Its synthetic controls independently check the normalized Gaussian/t5 density, detailed balance, a known Gaussian target and explicit prior-rejected holding times. No real anchored pilot has been evaluated by those controls.
