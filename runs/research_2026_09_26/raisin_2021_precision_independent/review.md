# Independent review of 2021 rounded-input timing precision fits

This review parses the native `CSP_ROW`, `CSP_OBJECTIVE`, `CSP_WROW`/`CSP_WDIAG`, FITRES, source DAT, and baseline instrumentation outputs independently of Astra's runner/parser. It executes **no native fit**. The [protocol](protocol.json) and [initial executable hash](execution-freeze.json) were saved before this verifier ran; two source-convention corrections and their first outputs remain in [`v1-peak-comparison/`](v1-peak-comparison/) and [`v2-native-cap/`](v2-native-cap/). They corrected the verifier, not the scientific gates: native peak is exactly the DAT header after R4 storage ([amendment](peak-storage-amendment.json)), and the data-Q cap multiplier is the original archived-2021 FITRES Q rather than the new baseline native Q ([amendment](cap-anchor-amendment.json)). The final verifier's hash is fixed in [v3-execution-freeze.json](v3-execution-freeze.json).

**Outcome:** all **248 endpoint cases** pass; **6 of 40 corner-stage cases** fail the frozen data-Q cap. Those six are all among the 32 directed corners. The eight candidate SIMLIB-time cases pass. No fit or case was removed and no tolerance changed. The independently failed cases are:

| New CID | Original CID | ΔD (mag) | Δdata Q | Frozen cap |
|---|---:|---:|---:|---:|
| 200007 | 2 | +0.0000430 | +0.023190 | 0.013789 |
| 200008 | 2 | −0.0000336 | −0.024359 | 0.013789 |
| 200011 | 3 | −0.0000351 | +0.012243 | 0.010000 |
| 200012 | 3 | +0.0000197 | −0.011886 | 0.010000 |
| 200031 | 8 | −0.0000434 | +0.017598 | 0.010000 |
| 200032 | 8 | +0.0000202 | −0.016450 | 0.010000 |

The [case ledger](case-results.json) retains **all 288** endpoint/corner cases and each gate. Every native run has an `ERRFLAG_FIT=0` FITRES row, source DAT hash and declared-only input changes, exact accepted physical `(MJD, band, data flux, error)` multiset after R4 flux/error storage, fixed shape=1/AV=0/RV=1.518/peak, and `NDOF=n−1`. The eight endpoint clone controls have exact non-CID FITRES and full three-iteration native-state equality to the unperturbed absent instrumentation. All final W and inverse C matrices are positive definite: minimum eigenvalues across endpoint/corner cases are **0.062593** for W and **0.303388** for C. Recomputed `(data−model)^T W(data−model)` closes native objective Q within **5.33×10⁻¹⁵**; maximum final-iteration D step is **2.074×10⁻⁵ mag** and positive-amplitude optimum gap in frozen final C is **6.52×10⁻¹⁰ mag**, both far below 0.001. Maximum absolute ΔD is **6.24×10⁻⁵ mag** for endpoints and **1.31×10⁻⁴ mag** for corners, so all distance caps pass. Endpoint maximum |Δdata Q| is 0.007446; corner maximum is 0.024359.

The source [precision ledger](../astra_design/raisin_timing_assets/rounding_2021/precision-ledger.json) contains 32 MJD tokens. An independent libc `snprintf("%.3f")` enumeration finds **one** compatible float32 state for each; every corresponding inverse open midpoint cell matches the frozen bounds exactly and is **0.00390625 observer day** wide ([timestamp ledger](timestamp-cells.json)). No extra 0.0005-day allowance was applied after inverting the float32 print state.

The post-independent [cross-check](author-result-crosscheck.json) agrees with Astra's saved outputs for all 288 cases: ΔD and archived Q caps exactly, and raw Q/ΔQ within 5.4×10⁻¹⁵. It also reproduces their endpoint PASS and corner FAIL. This validates the implementation of the *fixed rounded-input estimator and gates* in these saved local runs. Six directed corners failing data Q mean the precision perturbation stage does **not** meet its full frozen engineering gate; no paired timing or population-likelihood conclusion follows.
