# Prospective native timing recovery: engineering record

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

This experiment tests a specific mechanism: changing the peak time used by an
infrared distance fit from the injected truth to a peak estimated from the same
optical/infrared photons. It does not estimate the survey's population bias or
validate its existing cosmological corrections. The earlier archived-precision
experiment remains stopped under its original failed gate; see the
[asset-recovery report](raisin-timing-asset-recovery.md).

## Fixed scientific scope

The frozen protocol uses one DES16E1dcx-derived observing cadence, 117 signed
measurements in grizJH, eight fixed noisy draws and one noiseless draw. The
cadence itself was inherited from a selected archive, so this is a conditional
mechanism experiment. It is not a fresh simulation of epoch or object selection.
The truth has stretch 1, host AV 0, RV 1.5180000066757202, redshift
0.453000009059906 and peak MJD 57707.80078125. The distance modulus is
42.00489057043259. The two fitting arms must retain identical infrared photons.
No draw may be replaced because of an unfavorable fitted result.

The generator and fitter derive from SNANA v11_04d, commit
`10ec91297e4482d593cb5d3d055b10d4aa915071`. Two generator bounds repairs described
below are explicit deviations from that historical source. The fitter uses
twelve covariance iterations, checked against nine and displaced starting
points. It is a defined stable conditional estimator, not a reproduction of the
historical three-iteration estimator. Shape and host extinction remain fixed;
the joint arm fits distance and peak. Native peak-prior recentering, quoted-error
dependence on the realized flux, support boundaries and all iteration states
are retained in the audit.

The artifact root is
`runs/research_2026_09_26/raisin_profile_solver_review/csp-filter-interpretation/native-design/timing-followon-design/prospective-engineering/`.
The scientific protocol is `protocol-v5.json`, SHA-256
`376a8abb276f8ec917e8c9313a527e2bd9f44c4370a38bf6a50fd26496ddfb80`.
Additive freezes preserve each failed attempt and subsequent administrative
repair; no failed scientific outcome is deleted or relabelled as a successful
replicate.

## Generation result

Generation passed independent verification. There are exactly eight noisy
attempts, eight written events and 117 retained measurements per event. Their
negative-flux counts are 11, 6, 6, 5, 4, 6, 4 and 7: all 49 negative observations
remain in the 936-row data. The noiseless event contains 117 observations.
Output-only instrumentation leaves every native HEAD and PHOT table column
exactly unchanged. All 1,053 physical rows reproduce the native mean, noise
shift, error and serialization equations. Six non-generated peak placeholders
per event are distinguished from physical measurements.

The reported-to-true variance ratio spans approximately 0.975864–1.023867 in
these draws. This is the known realization-dependent quoted-error rule, not a
measured physical noise miscalibration. The native error-map applicability audit
found no matching configured map for this cadence; it did not remap fields or
silently introduce a different noise model. See the
[noise-source review](raisin-prospective-noise-review.md).

Independent evidence: `v5-generation-review/`, `generation-gate.json`, and
`generation-v5/{original,ledger,noiseless}/`. FITS observation MJD is double
precision; flux, quoted error and selected header fields have their native
single-precision formats. The earlier text-export rounding model must not be
imposed on this new FITS input path.

## Noiseless recovery passed; first noisy branch stopped

The two completed joint fits have exactly identical numerical states and
FITRES scientific rows with and without support instrumentation. The independent
root parser checked 24 iteration callbacks, exact 117-row input membership,
covariance positivity, objective arithmetic, fixed nuisance values and phase
support. Both recover distance with error +0.0000603123 mag and peak with error
+0.00160196 day. The largest model residual divided by its quoted measurement
error is 0.01369394, within the predeclared 0.02 limit.

The infrared-only noiseless fit also passes: distance error +0.0000387753 mag,
exact fixed peak and maximum standardized model residual 0.01114779. The complete
independent audit covers 36 callbacks. The source-matched all-call support
instrumentation reports no invalid mean evaluations in these noiseless fits.
Evidence: `fits-readme/`, `fits-resume/noiseless_NIR/` and
`root-noiseless-complete-review/`.

The eight noisy joint fits then failed the original all-call support gate.
Two events encountered unsupported template phases during optimizer trials:
CID4 has 304 invalid calls and evaluated phases spanning approximately
−83.819 to +107.915 rest-frame days; CID7 has 147 invalid calls and reaches
+107.878 days. The mean model's audited common support is −20 to +70 days.
No subsequent peak-intervention arm was executed and no paired timing effect
is certified from this attempt. The broad native peak bounds allowed invalid
intermediate evaluations even though a final fitted point might be in range.
Convergence at that final point cannot erase the unsupported calculations.

This is a failure of the specified computational estimator under these noisy
inputs. It is not evidence that a supernova's real luminosity differs from the
model. The separate physically restricted estimator below preserves all eight
events and this original failed branch. It cannot inherit a successful status
from an unsupported run.

## Restricted estimator: recovery and actual start tests pass

The private source patch intersects the loaded model and KCOR phase support
over all117 observation times at the fixed redshift. Its absolute peak bounds
are `[57671.342999365806,57715.067000181196]` MJD. These bounds come from the
cadence and tables, not the true peak or measured brightness. They define a new
constrained numerical estimator; changing MINUIT's bounded coordinate mapping
can also change interior optimization trajectories. Every physical mean call
has an independent support guard. The original ±4day final recovery requirement
remains a separate gate.

The compiled binary passes exact disabled identity in both absent and zero
modes: each comparison matches all2,832 CSP state records, FITRES science rows,
support summaries and native input tables. With the bounds active, noiseless
joint/NIR distance errors remain +0.0000603123/+0.0000387753mag. An independent
root audit checks all24 callbacks against the physical FITS measurements.

Three noisy joint cases completed before a start-intervention gate stopped the
executor. The nominal12 and9iteration cases have exactly identical final
distances and peaks for all eight draws. The minus-initializer case agrees to
4.58e−9mag and1.90e−7day. Independent checks of all264 callbacks pass membership,
covariance, objective, final stationarity and all-call support. These successful
state checks do **not** establish a displaced-start test: every minus case
actually enters the optimizer at the same peak as its nominal counterpart.

Source inspection explains the difference. The specified −2day initializer is
present before `FITINI_ADJUST`, but that native routine searches a five-point
peak grid and overwrites the initial value before MINUIT receives it. The
gate correctly refuses to count this as a −2day optimizer-entry displacement.
Further, the native peak prior uses `INIVAL` as its center, so directly shifting
that shared variable would change both a start and a prior. A separate,
explicitly logged local MINUIT-value intervention now moves only the optimizer
start while preserving the prior center. Four disabled controls reproduce all
saved scientific states exactly. Actual MINUIT readback confirms first-iteration
shifts of minus/plus two days for every draw, with later starts unmodified.
The two completed start tests agree with the original nominal branch to
4.87e-9mag and1.91e-7day. No failed case is discarded.

The infrared stage then stopped before running a fit because its null HEAD
adapter changed unrelated string padding. An independent byte comparison finds
exactly712 spaces rewritten as NULs; every numerical column is unchanged.
A narrower adapter now preserves the complete file outside the 32 peak-time
bytes and passes exact-file identity for its null. The failure and partially
written adapters remain preserved; the completed infrared result is below.

Evidence is in `restricted-peak-engineering/`, `fits-restricted/`,
`root-restricted-noiseless-review/` and `root-restricted-joint-review/`. The
original scientific thresholds and each failed execution remain recorded.

## Completed same-photon infrared engineering comparison

The byte-preserving adapter passes exact-file null identity and changes only the intended peak-time fields. All nine infrared jobs pass the unchanged gates. Two independent reviews check all 816 callbacks against the six generated J/H measurements, assigned fixed peaks, positive covariance, objective arithmetic and numerical controls. The largest distance difference across numerical controls is 1.65e-9 mag. All 98,316 recorded infrared mean calls remain inside model support. The complete null model, weights, objective, FITRES and HEAD bytes are exact.

| Draw | Fitted minus true peak (days) | Paired infrared distance change (mag) |
|---|---:|---:|
| 1 | -0.25000000 | -0.003783499 |
| 2 | +0.04296875 | +0.000634137 |
| 3 | +0.23437500 | +0.003531197 |
| 4 | -0.16406250 | -0.002478561 |
| 5 | -0.06250000 | -0.000941816 |
| 6 | +0.28906250 | +0.004366774 |
| 7 | -0.05078125 | -0.000760085 |
| 8 | +0.14453125 | +0.002261087 |

The mean paired change is **+0.000353654 mag**, sample SD **0.002877872 mag** and ordinary Monte Carlo standard error **0.001017481 mag** across these eight generated draws. Fitted-peak errors have mean +0.02295 day and SD 0.19008 day. This engineering sample does not resolve a nonzero mean timing response. It is not a bound on the observed population: only one selected cadence, one fixed template/dust state and the stated native noise model were generated. The peak is estimated from the same optical and infrared photons, so its correlation with the infrared distance noise is retained.

Evidence: `restricted-peak-engineering/start-only-hook/engineering-result.json`, `independent-nir-review/`, `full-null-identity.json` and `root-nir-engineering-review/`. The latter was produced by the independent root script `scripts/research_2026_09_26/verify_prospective_nir_states.py`. Cumulative fitter activity for this restricted-estimator sequence is 35.23 seconds; source builds and earlier failed branches retain separate ledgers. A separately frozen 64-draw conditional study is being prepared; these eight engineering outcomes will not be pooled into it.

## Preserved engineering failures

These failures limit reproducibility until repaired; they are not evidence for
extinction, population evolution or new cosmology.

1. An initial LIBID lock of 1 did not match the only cadence, LIBID 11. It ran to
   the native retry limit without generating a physical light curve. The
   corrected identifier was frozen before the next attempt.
2. The historical generator overflowed two arrays while computing cut summaries,
   even with cut application disabled. One array lacked the extra slot required
   by an inclusive loop; two temporary time/sort arrays used one-based indices
   with insufficient allocation. The first event's flux/noise had been generated
   before this crash, but no event was written. The exact local repairs alter
   storage bounds, not the mean, noise draws or cuts. Replaying that seed is not
   an additional independent draw.
3. Missing output parent directories caused a later pre-generation failure.
   Existing failure files were preserved before the parents were created.
4. The fitter's shared simulation registry counted multiple paths with the same
   version name despite a private data path. A private reference-data root with
   an empty simulation registry resolved lookup without changing observations.
5. The legacy README checker copied an overlong provenance-path token into a
   small buffer. A private metadata adapter inserts one blank line after
   `DOCUMENTATION:`. All tokens and parsed YAML values are identical; numerical
   data and fitter binaries are unchanged.
6. The first successful joint fits exposed an audit-code error:
   `CUTFLAG_SNANA=3` means both native pass bits are set, whereas the checker
   required zero. `ERRFLAG_FIT=0` and the native cut summary confirm success.
   Independently, a checker inherited a single-precision MJD expectation from
   the old text path; new FITS times retain double precision. Both interpretation
   errors must be corrected transparently against source and existing outputs,
   without rerunning or selecting different fits.

The actual displaced starts and same-photon infrared intervention now pass
the unchanged numerical gates. The next step is a separately frozen conditional
simulation with new noise draws. A larger conditional
study requires its own frozen design and must retain all original failures.
Native activity and the separate build budgets remain recorded.
