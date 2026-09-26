# Start-only MINUIT hook — prepared, not built or fitted

The exact patch is `start-only.patch`, SHA256 `6d1b6ce87d632224c411e6d8cb80b167cae93485e3e72db0977d362285d5cd05`. It is an additive proposal on the already frozen restricted-domain source. The existing binary, source tree, failed original8 branch and partial restricted runs remain unchanged. No new native fitter was built or executed here.

The failed `joint_minus` case applied a −2day **pre-grid initializer**, but native adjustment selected the same actual MINUIT entry for all eight objects. It remains useful as an initialization-grid sensitivity and does not satisfy the original displaced-entry gate. See `../start-diagnostic/diagnosis.md`.

## Implementation

In `MNFIT_DRIVER`, copy each source `INIVAL(IPAR)` to a local scalar. Pass that scalar to `MNPARM`. For the parameter named `PKMJD` only, a C helper optionally adds exactly−2 or+2days in **covariance iteration1 only**. The driver’s first parameter name must be `ITER`. The same nominal NML initializer is used in all branches: the shift is relative to each object's own native post-grid entry. It is not the old NML−2 branch and is not anchored by a fitted-outcome choice of offset.

The environment variable `PROSP_MINUIT_PEAK_SHIFT` is absent/zero by default. In that path the helper changes no arguments and writes no records; `MNPARM` receives the exact copied native value. The executor must remove inherited values for every default, noiseless and NIR job. Nonzero values are restricted to±2, require `PROSP_HARD_PEAK_DOMAIN=1`, a free peak, finite inputs and an entry inside the existing fixed bounds before MINUIT sees the coordinate. Invalid cases terminate with status87. There is no clipping, new penalty, epoch change or mask change.

**Common INIVAL is never written.** It is a const pointer in the C helper; only the separate local optimizer scalar and private readback flag are writable. Thus the first prior center, grid history and already prepared initialization state stay native. The existing CSP_ENTRY record continues to mean the native initializer/prior state; it is explicitly no longer a record of the modified MINUIT start. `PROSP_MNPARM_PEAK` names both quantities.

After successful parameter definition, enabled mode calls native `MNPOUT` into separate local name/value/error/bound/index buffers. Source `minuit.F:9641–9677` retrieves U/ALIM/BLIM and writes only those output arguments. It performs no objective evaluation or minimization. The new readback record therefore exposes MINUIT's actual stored external coordinate, not just the proposed argument. It checks the peak name/index and unchanged bounds. Stored minus common source must equal the declared first-iteration±2 to the existing<0.004day precision gate; later iterations must equal the unshifted source exactly. Readback remains within the fixed domain. The existing all-mean support guard remains authoritative for every subsequent physical evaluation.

The same objective **function and prior center** are preserved at the first iteration. This does not require identical first-iteration W at different trial coordinates: native diagonal model variance depends on those coordinates. Compare prehook/native state or objective components at common coordinates and verify identical first CSP_ENTRY/prior state. Later covariances and recentered priors may follow different fit histories. This test diagnoses iterative-estimator path stability; it is not a proof of one globally fixed Gaussian objective or a unique optimum.

## Standalone evidence

Eleven tests pass in `test-result.json`. The harness links the exact existing native `minuit.o` (SHA256 `9bbc88ab0568b229a5e4bad888f37fb58383c5fa2755beae2ef9a7336dc2cea7`) and invokes only MNINIT, MNPARM and MNPOUT. It links no photometry, light-curve objective, minimization, covariance or RNG. Its test-only INTRAC stub reports a noninteractive terminal.

Absent, zero,−2 and+2 cases each inspect48 actual MINUIT stored values: four parameters over12 iterations. Only the first-iteration peak changes, by exactly the declared±2; all input arrays/prior centers, other parameters and later entries remain exact. Refusal tests cover shifting a fixed peak, an out-of-domain start, an undeclared offset, NaN, missing hard-domain mode, deliberately missed MNPOUT shift and deliberately altered later-iteration readback. The latter corruptions occur only in the test fixture, never in the patch. Ordinary disabled scientific replay remains mandatory after a separately authorized private build.

## Additive execution proposal, not authorized here

1. Build a separate private tree only after root approval; retain the restricted source and binary. Freeze source/build identity first.
2. With offset absent and then0, replay the already completed nominal joint12 and joint9 scientific states. Require exact full saved model/error/weight/objective/CSP_ENTRY/FITRES/support identity. These are identity controls supporting reuse of the original nominal results, not new scientific baselines or draw selection.
3. Use nominal NML bytes for new first-iteration±2 local-entry jobs. Verify MNPARM and MNPOUT records, unchanged first prior center/native CSP_ENTRY/masks and identical state at common coordinates before interpreting any differences. Keep all original numerical, support and final±4day engineering gates; no lowest-Q winner.
4. Only if these gates pass resume the previously declared NIR true/fitted-peak/null sequence using the exact existing photons and nominal joint12 peak estimate. The new environment variable must be absent in every NIR/default job. Preserve the full null checker and same-PHOT adapter checks.

Carry the existing11.420924458seconds into the same restricted-estimator120-second native budget, with40seconds maximum per process. No reseed, new photons, deleted CID, relaxed tolerance, automatic64-draw continuation or observed-population correction follows. The initializer-minus run and every earlier failure remain in the ledger.
