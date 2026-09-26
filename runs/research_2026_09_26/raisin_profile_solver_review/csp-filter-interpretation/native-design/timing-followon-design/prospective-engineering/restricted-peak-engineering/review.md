# Restricted peak estimator: source and execution handoff

The private build and static checks pass. No native fit or new photon generation has run with this estimator. The original eight-draw unrestricted branch remains failed: valid final coordinates do not undo its unsupported intermediate physical means. The new estimator is explicitly distinct and must pass disabled-binary identity and active noiseless recovery before the same eight photons can be evaluated.

The binary is `phase2/pte/restricted-peak-engineering/build/bin/snlc_fit.exe`, SHA256 `4ba9c40660fd1d1df42dd300b64f56c0aaa42e1d6e5c891c6d4b56c60737cf14`. Its single build took 15.445649625 seconds of the authorized 120 seconds. The prior support binary remains byte-identical at SHA256 `d3ed187ec43968894ec651a6f94a8a9e1258058e2155dc9f523e03cf9fd56ade`.

## Domain and actual source paths

For fixed heliocentric redshift, each exposure requires

`phase_min <= (MJD_i - peak_absolute)/(1+zHEL) <= phase_max`.

Intersecting these requirements over **all117 native R8 exposure times**, including the optical rows when fitting only NIR, gives

`max_i[MJD_i-(1+zHEL)*phase_max] <= peak_absolute <= min_i[MJD_i-(1+zHEL)*phase_min]`.

The phase interval is the intersection of the loaded SNooPy table `[-20,70]` and KCOR table `[-20,85]`. At the fixed native redshift `0.453000009059906`, the absolute peak interval is `[57671.342999365806, 57715.067000181196]` MJD. It uses no brightness, fitted residual, fitted peak, or true peak. The initializers and the original final ±4-day engineering recovery requirement remain separate. Tabulation limits do not establish empirical model accuracy or a physical population timing prior.

The changes are limited to two existing source files and one new private header. `restricted-peak.patch` records the exact delta from the prior support build.

- `snlc_fit.car:8424` calls the domain helper after first-iteration user and simulation initial-value overrides. It precedes the subsequent covariance initialization and initial peak-grid adjustment.
- `snlc_fit.car:8129` calls it before the later-iteration early return. The helper asserts exact full metadata, redshift, MJDOFF, phase tables and native bounds, then retains the same bounds.
- `snlc_fit.car:8451` passes all `SNLC8_MJD` values, the actual fixed fitted-redshift initializer, native MJDOFF and loaded table limits. The native `MNPARM` call at unchanged `snana.car:36586` consumes `INIBND`.
- `snlc_fit.car:4443` checks the actual FCN peak after native prior-only/out-of-bound early returns and before physical means. This also covers covariance/sigma calls that bypass native priors.
- `snlc_fit.car:5390` checks every USRFUN phase, shape, redshift and finite nuisance value before nearest-filter selection, SNooPy, extinction or KCOR. Auxiliary rest-frame calls are checked at their actual phase/redshift. This refusal remains active even if the optional output-only `PROSP_FIT_SUPPORT` flag is absent.

`PROSP_HARD_PEAK_DOMAIN` absent or exactly `0` returns before native writes, checks, or new output. Exactly `1` enables the new estimator. Other values abort. The executor explicitly removes inherited values before setting the intended mode. The hook never clips a coordinate, deletes an epoch, adds a likelihood penalty, changes random draws, or substitutes an objective. Unsupported physical calls terminate with an explicit record and status86.

Bounds are represented in native peak coordinates by subtracting MJDOFF. The helper verifies endpoint phases in both the exact R8 FCN operation order and the separate R4 mask arithmetic. Only inward representable rounding is permitted if numerically necessary; no adjustment is required for this fixed cadence. Raw and safe bounds are both logged. The native covariance iteration, state-dependent weights and recentered five-day peak prior are otherwise unchanged. This is not a Gaussian maximum-likelihood replacement.

## Validation completed without fits

Thirteen standalone C tests pass: absent/zero mode leaves even deliberately invalid inputs and bounds untouched; all current endpoints pass; a nonzero MJDOFF gives the same absolute interval; tighter KCOR limits narrow the domain; incorrect row count, floating redshift, changed metadata, out-of-domain FCN/mean, invalid shape, nonfinite nuisance and invalid environment setting terminate as intended. The first standalone compilation caught a diagnostic format count and test-harness indentation, preserved under `first-static-format-failure/`, before any native compilation. A provenance helper initially lacked the private compiler PATH; that draft is also preserved there.

Independent endpoint arithmetic yields raw endpoint phase extrema `[10.0922233605,70]` and `[-20,39.9077766395]` in R8. The corresponding R4 mask extrema are `[10.0922231674,69.9999923706]` and `[-19.9999980927,39.9077758789]`.

A separate synthetic log fixture tests the execution checker: complete schemas/endpoints pass; missing bound records, inconsistent mean-call counts, a hard-abort record and an invalid final peak are refused. This fixture is explicitly synthetic and is not a native fit result.

AstraScience independently reviewed patch/base identity and source placement. Its review is under `runs/research_2026_09_26/astra_design/raisin_timing_assets/prospective_noise_review/physical_peak_domain/patch-review.json`. Static inspection cannot replace the pending disabled-binary replay. The broad prebuild directory inventory included the generated `ypatchy.log`; only that build log changed during compilation. Its prior bytes are retained as `ypatchy.before-build.log`, and all actual frozen source files remain unchanged.

## Frozen execution stages, not executed

`execution-protocol.json`, `execution-freeze.json` and `run_restricted.py` define three independently released stages:

1. `disabled`: absent and explicit0 mode on the saved noiseless joint input. Every CSP entry/measurement/model/error/phase/weight/objective token, FITRES science row and support record must equal the successful old `fits-readme/noiseless_joint` run exactly. All HEAD columns and the same PHOT file are checked.
2. `noiseless`: enabled joint12 and NIR12. Original thresholds remain: distance0.001mag, peak0.01day, maximum standardized residual0.02, exact masks, positive definite covariance, fixed-C amplitude stationarity, iteration stability and all-call support.
3. `noisy`: the exact existing eight photons, with joint12/9 and ±2-day starts; NIR truth/fitted-peak arms each12/9 and ±0.2mag post-initialization starts; then a full null adapter. Every original gate remains, including all-draw membership, 117/6 exact rows, distance/peak/C stability and final ±4day recovery. Only HEAD PEAKMJD changes for the peak intervention. The full null model/W/objective/FITRES comparison and bounds checks precede creation of the paired-effect table.

The new outputs are root-level `P/fits-restricted/` and `P/derived-restricted/`, preserving the established relative NML path depth. A **separate120-second cumulative fitter budget**, at most40seconds per process and one worker, includes all17 planned new processes. Earlier failed-estimator time is preserved in its original ledger and is not erased or recycled. Generation is absent from this executor. There is no retry, refill, best-Q choice, source support relaxation or automatic64-draw continuation.

`PROSP_DOMAIN_SUMMARY` reports search extrema and minimum endpoint distances over all evaluated FCN trials. Final solution margins are instead calculated from the final CSP_OBJECTIVE peak against the fixed logged bounds. A trial touching an endpoint does not make an interior optimum a boundary fit.

Invocation, only after the root supplies the stage-specific frozen release:

```
phase2/env-official/bin/python phase2/pte/restricted-peak-engineering/run_restricted.py disabled --release <root-release.json>
```

Replace `disabled` by `noiseless` or `noisy` only for the corresponding separately released stage. A successful engineering comparison would establish a conditional native timing mechanism for this cadence/template/noise/weight estimator, not a historical or survey timing correction and not corrected cosmology.
