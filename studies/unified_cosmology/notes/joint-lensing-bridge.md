# Sensitivity to the joint ACT–Planck–SPT lensing release

This is a prepared, separately identified sensitivity calculation. **No posterior result exists yet.** It evaluates the [public joint lensing alternative](joint-lensing-availability.md) without modifying or replacing the current cosmological measurement.

Both alternatives are fixed before any posterior evaluations:

- **Baseline:** 35 lensing bandpowers, retaining the current ACT window support L=40–763, replacing Pan's SPT 2018 temperature reconstruction with SPT 2019–2020 polarization MUSE and using the released joint covariance.
- **Extended:** 38 bandpowers, additionally extending ACT to L=1300, as in the joint paper.

They use identical newly computed native spectra and retain separate weights, identities and qualification results. Neither is selected according to its outcome. The global Hartlap precision factor also changes when the new data-vector dimension changes; the replacement is not merely appending independent SPT data.

At each previously selected, qualified source point, the unchanged source model is evaluated once with native CAMB. Its full likelihood components and total prior must reproduce the saved native result within 10⁻⁷ in log density. Raw potential and lensed CMB spectra are then archived with parameter, source, data, numerical-environment and payload identities. There is no spectral emulator in this step. The recorded `native_evaluations` counts full-source `model.logposterior` invocations; CAMB may reuse an identical cosmological state, so this is an upper bound on actual new solver evaluations, not an instrumented solver count. The two new lensing likelihoods use exact algebraic compression of the released ACT/Planck response operators plus the SPT windows; this is checked against the unmodified author's likelihood.

The untrimmed posterior weight is the original exact/proposal weight multiplied by the ratio of the new joint likelihood to the old ACT–Planck and Pan-SPT likelihoods. Primary CMB, BAO, supernovae, physical theory and priors remain identical. The Pan-only foreground amplitude `A_fg` retains its normalized uniform [0,2] prior as an auxiliary variable; it is included in stability diagnostics and integrated out of physical reporting. This avoids an implicit change of parameter-space measure.

Fixed Gaussian normalization constants are recorded using the actual precision, including Hartlap scaling. They cancel when posterior weights are normalized. Since these alternatives use different observed data vectors and conventions, these offsets are **not** evidence ratios or significances.

A source posterior is accepted only through the qualified measurement consumer, which recomputes its original correction and chain checks. Both new variants must independently pass the unchanged raw-weight, Pareto-tail, independent-chain and contiguous-batch stability gates. Failed evaluations and inadequate overlap remain explicit and withhold posterior results. These finite checks do not establish coverage of unseen modes. If reweighting fails, separate sampling is needed; trimming weights or relaxing gates is not part of this design.

The new cross-covariance treats dependence between lensing experiments. It still supplies no primary-CMB/lensing, lensing/BAO or lensing/supernova cross block. It uses the released Gaussian MUSE compression, not the original joint delensed-EE likelihood. These conditional assumptions carry into both sensitivity measurements.

## Preparation and execution

The [design](../code/inference/joint-lensing-bridge-design.json), [implementation](../code/inference/joint_lensing_bridge.py) and [validation](../results/inference/joint-lensing-bridge-validation.json) are separate from active samplers and frozen targets. Validation uses synthetic densities, qualification/cache failures and existing stored spectra; it makes no new native CMB or GPU calls.

With a qualified source correction available, first create its immutable plan:

```bash
.work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/joint_lensing_bridge.py \
  --chain-folder /absolute/path/to/qualified/source/chains \
  --correction-summary /absolute/path/to/qualified/correction-summary.json \
  --cache "$PWD/.work/unified-cosmology/joint-lensing-bridge/source-name"
```

The default stops after preparation. Add `--execute-native --workers 1 --output /absolute/path/to/new-result.json` only when resources and the qualified parent are available; up to four workers are supported. `--summarize-only` reuses a complete, verified native spectral cohort. Existing partial attempts without sealed results are blocked from automatic retries, and completed or partial numerical caches are verified before reuse. Generated spectra remain ignored; authored validation and compact results are retained.
