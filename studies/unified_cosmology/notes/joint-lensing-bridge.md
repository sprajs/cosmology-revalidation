# Sensitivity to the joint ACT–Planck–SPT lensing release

This is a prepared, separately identified sensitivity calculation. **No posterior result exists yet.** It evaluates the [public joint lensing alternative](joint-lensing-availability.md) without modifying or replacing the current cosmological measurement.

Both alternatives are fixed before any posterior evaluations:

- **Baseline:** 35 lensing bandpowers, retaining the current ACT window support L=40–763, replacing Pan's SPT 2018 temperature reconstruction with SPT 2019–2020 polarization MUSE and using the released joint covariance.
- **Extended:** 38 bandpowers, additionally extending ACT to L=1300, as in the joint paper.

They use identical native spectra and retain separate weights, identities and qualification results. Neither is selected according to its outcome. The global Hartlap precision factor also changes when the new data-vector dimension changes; the replacement is not merely appending independent SPT data.

At each previously selected, qualified source point, a verified spectral sidecar from the original native correction is reused when available. The optional [spectral correction companion](../code/inference/spectral_correction.py) retains that calculation's raw potential and lensed CMB spectra, tied to its complete native record by content hashes. It leaves the original density-correction code, selection, weights and gates unchanged. The exact and proposal models have separate provider states; a native comparison tests that the captured spectra still belong to the exact model after the proposal evaluation. Missing or failed captures never authorize using an unrelated provider state.

Without a valid sidecar, the unchanged source model is evaluated once more with native CAMB. Its full likelihood components and total prior must reproduce the saved native result within 10⁻⁷ in log density. There is no spectral emulator in either path. The recorded `native_evaluations` counts additional full-source `model.logposterior` invocations: zero for verified reuse, one for reevaluation. CAMB may cache an identical cosmological state, so this is an upper bound on new solver calls. The two new lensing likelihoods use exact algebraic compression of the released ACT/Planck response operators plus the SPT windows; this is checked against the unmodified author's likelihood.

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
