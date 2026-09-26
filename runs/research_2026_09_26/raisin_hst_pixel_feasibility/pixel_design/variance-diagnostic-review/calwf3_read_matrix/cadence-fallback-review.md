# Cadence fallback: what the available darks could identify

The exact-match avenue is **not yet exhausted**. The first40 date-ranked long-exposure headers at each visit give0 full-frame SPARS50/NSAMP≥8 matches. The existing ±30-day CAOM lists contain79 and73 such long-exposure roots, leaving39 and33 nearer headers unchecked. The ±180-day expansion is CAOM metadata only and cannot establish their read patterns. Root has authorized a separate metadata-only continuation for those remaining nearby roots; the original search/caps/results remain intact. This review neither changes Sol's search nor authorizes an exposure download.

Independent aggregation of the saved80 headers finds8 full-frame STEP50/NSAMP16 ramps near the search visit and6 near the template visit. Each visit also has20 SPARS25 subarray ramps; the other full-frame ramps use STEP100/200/400 or SPARS200. The metadata do not prove that any alternate sequence contains the exact science read times. RAW and IMA products are available for the inspected candidates. Representative full-frame STEP50 RAW and IMA sizes are34,027,200 and168,258,240 bytes respectively; matching a product URI is not a reason to acquire the larger product now. Full ledgers and source hashes are in `fallback-metadata-result.json`.

## Ranked options

| Rank | Candidate | What can be measured | Additional transport assumption |
|---|---|---|---|
| 1 | Repeated full-frame SPARS50 with the same eight-read sequence | Matched detector temporal/spatial covariance and its actual slope/aperture contraction | Detector state, persistence and processing must still be comparable; dark photons differ from science sky photons. |
| 2 | Longer full-frame SPARS50 whose original RAW first8 reads exactly match | A prefix measurement using those original reads; no temporal covariance interpolation is needed | Verify exact TIME/reset/read-history identity and process the prefix without later-read rejection or calibration leakage. A whole-ramp FLT is not the prefix product. |
| 3 | Same sequence from a more distant epoch | The same sampling operator on a different detector state | Explicit temporal stability across calibration cycles, amplifier behavior and persistence. Keep epoch groups separate. |
| 4 | Full-frame STEP50, including a verified late seven-read near50-second segment | Its own covariance; a falsification test of cadence-independent white-read or stationary-lag assumptions | Transfer to SPARS50 requires invariance to prior reads/reset age and exact lag matching or a declared interpolation model. This is not an exact prefix. |
| 5 | Full-frame other STEP/SPARS sequences | Their own lag/scale dependence, common electronics and spatial structure | Extrapolation/interpolation of unmeasured time lags and cadence history; not an identified SPARS50 variance. |
| 6 | SPARS25 subarrays | Region/layout-specific noise and a check for gross failure of a universal detector-noise claim | Timing, readout layout, detector coverage and reference-pixel corrections all change; no automatic full-frame aperture calibration. |

The [official WFC3 timing documentation](https://hst-docs.stsci.edu/wfc3ihb/chapter-7-ir-imaging-with-wfc3/7-7-ir-exposure-and-readout) distinguishes the early rapid/logarithmic reads in STEP from SPARS and notes that subarray timing differs even for the same sequence name. Its rounded table is a design guide; actual per-read TIME/reset metadata must establish identity. In particular, a STEP50 tail with approximately50-second spacings has a different preceding read history from the first seven nonzero SPARS50 reads. Later-read DQ decisions also cannot be silently inherited into a supposedly shorter exposure.

## Why a different cadence cannot determine the target covariance

For target slope coefficients h, the quantity needed is hᵀK_target h, with spatial covariance included for an aperture. A true matched prefix supplies a principal submatrix of the **same early measurement process**. A different sequence instead measures K_other. Without a relation between the two, positive-semidefinite covariance matrices can agree on every alternate-sequence measurement and differ in the unobserved target slope direction. More pixels or repeated alternate ramps do not remove that structural gap.

A stationary covariance K(t_i,t_j)=k(t_i−t_j), independent of read history and detector state, could allow a STEP50 tail to constrain the lag values relevant to SPARS50. That is a new physical assumption to test. Equal spacing alone does not establish it. Comparing separate predeclared windows within repeated STEP50 ramps can challenge time-since-reset stationarity; comparing cadences can reject a universal independent-read model. Agreement is weaker evidence and does not prove transport to a cadence never measured. No resampling or interpolation should be labelled newly observed independent reads.

For constant-rate Poisson accumulation, an initial accumulated random count is a common offset and cancels from a fixed-intercept slope. That exact algebra does **not** make electronic read noise, persistence, nonlinear response or adaptive weighting invariant to the earlier ramp history. The same distinction applies to common-zero-read subtraction.

## Stop rules and remaining scope

Require at least four comparable repeats for the planned bounded repeated-ramp control; retain fewer as metadata/engineering feasibility only. Stop exact-match inference if actual read times/reset conventions are unavailable, prefix history cannot be verified, or the source/product replay cannot retain the declared reads and calibration units. Missing optional metadata remain unknown, not permission to infer a match from EXPTIME. Do not search indefinitely or expand resource caps silently. Current root authorization is metadata-only.

If exact sampling remains unavailable after the frozen continuation, the useful next proposal is a separately named **cross-cadence falsification control**, preferably the already observed full-frame STEP50 groups. Freeze its native-sequence or tail-window estimand before pixels. Retain signed RAW pairs, fixed detector masks, master-dark training dependence, temporal/spatial blocks and actual source-operator checks from `dark-control-design.md`. It cannot supply a SPARS50 error correction. No dark result alone identifies the cause of the science repeat-variance deficit or justifies changing SN distances or cosmology.

This note used saved metadata and one official timing page only. No duplicate MAST query, native fit, image/pixel read or outcome-selected cohort change occurred.
