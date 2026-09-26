# Independent saved global-profile review

The completed DES16C1cim fixed-covariance numerical pilot passes independent
arithmetic, native-stream, membership, support and declared resolution checks.
The fitted distance is nevertheless ambiguous across nearly tied shape
branches. No native fit was run by this review.

## What was verified

`verify_profile.py` imports no Astra fitting/export implementation. Its only
local helper is this reviewer's earlier raw-log parser. In 4–5 seconds on one
BLAS thread, it checks all 36 frozen output hashes, then reconstructs:

- All 86,813 saved native mean vectors and all four covariance/row metrics.
  Maximum disagreement is 9.095e−13 in the quadratic, 7.106e−15 mag in DLMAG,
  and 1.422e−14 in amplitude.
- All 93,159 actual native stream replies against exact coordinates and saved
  arrays, including the initial/final unchanged-state replay and 3,172 native
  amplitude checks. Maximum pointwise amplitude-scaling discrepancy is
  4.633e−9. The stream returns exactly the initial data, rows, C and precision.
- Both covariances against their raw native precision exports. Primary C is
  the B reference; alternate C is evaluated on all B rows at the pre-existing
  A/start1 coordinates. The 65 A rows are exactly the positive subset of the
  74 B rows, with identical native-ingested flux/error coordinates. Each A
  metric inverts the marginal C[A,A]. No precision principal block or error
  replacement is substituted.
- Every 261-shape/33-AV coarse and 521-shape/65-AV fine grid coordinate, all
  CSV branches and quadratic values, every sampled fine-profile local branch's
  continuous-shape follow-up, and all 1,548 native algorithm-boundary probes.
  The native float32 step is 0.00461538415402174. Stored knots are checked
  against the released FITS file; algorithm boundaries and ±1e−8 offsets match
  exactly. All means and tested profiled amplitudes remain in the declared
  native mean and DLMAG support.

The largest coarse-to-final changes are 1.012e−6 in Q and 1.666e−7 mag in
DLMAG, passing both the executed 0.001 thresholds and the earlier review's
stricter 0.0001 suggestions. Every reported minimum also equals the best of
*all* saved probe vectors to numerical precision; no better cached candidate
was discarded by the result summarizer. This establishes the declared finite
search gate, not a proof of a global minimum over the continuum.

## Result and its limitation

| Frozen covariance anchor | A positive-row minimum DLMAG | B signed-row minimum DLMAG | B−A |
|---|---:|---:|---:|
| B reference | 42.753754402 | 42.588805188 | −0.164949213 mag |
| Pre-existing A reference | 42.756846530 | 42.591207757 | −0.165638773 mag |

The anchor changes the paired minimum difference by −0.000689559 mag. This is
only the specified covariance-state sensitivity; it does not test covariance
feedback, the free-peak fit, historical MW94, or a different noise model.

Under the primary B metric there are ten sampled fine-profile local branches,
all within ΔQ=0.551 of the minimum. In particular:

| B branch | Stretch | Empirical AV | DLMAG | Q−Qmin |
|---|---:|---:|---:|---:|
| Lowest saved Q | 1.032307647 | +0.051309317 | 42.588805188 | 0 |
| High-stretch competitor | 1.189230708 | −0.054760828 | 42.767926739 | 0.115738988 |

The two distances differ by **0.179121551 mag**, much more than the numerical
resolution. The high-stretch signed-data branch is only +0.014172338 mag from
the positive-row global minimum. Thus the −0.16495 mag minimum-to-minimum
change largely reflects selection of a different nearly tied branch.

At the two B branch optima, exact marginal/conditional block algebra gives

`Q_B = Q_positive,marginal + Q_negative|positive`.

Moving from high stretch to the lowest-Q branch worsens the positive-row
term by 1.987576318 and improves the nine-negative-row conditional term by
2.103315306, leaving a net improvement of only 0.115738988. The block sum
closes within 2.85e−14. Separately, the parent's four-point decomposition is
independently reproduced in `parent-decomposition-check.json`; the fixed
A-shape/AV amplitude response is +0.012454 mag, followed by a −0.177403 mag
shape/AV branch response under B.

The report also saves descriptive fixed-Q level sets at ΔQ=1,4,9. For each
saved shape/AV coordinate, it includes the **exact analytic amplitude** range
`a_hat ± sqrt((Qmin+ΔQ−Qprofile)/(hᵀWh))`, then transforms to DLMAG. At ΔQ=1,
the primary B distance envelope is **42.471390–42.865636**, versus
42.638328–42.861942 for A. These are finite-saved-coordinate envelopes, not
posterior or calibrated confidence intervals. Unsampled nuisance coordinates
could enlarge them. At ΔQ=9 the declared shape boundaries are reached: the
minimum shape-edge costs are 4.249695 for A and 8.119063 for B. That broad
level set is also domain-limited; the ΔQ=1 and 4 sets do not reach those
shape edges. AV edges are more than 130 above the minimum.

## What this permits

The same frozen numerical diagnostic can be extended to the predeclared
cohort, preserving per-object failures, competing branches and domain limits.
Duplicate-row mapping must preserve multiplicity and verify covariance
exchangeability before treating identical keys as interchangeable. The pilot
does not support an identified distance correction or a population bias
estimate. A and B have different row sets, so their minimum Q values must not
be compared as a likelihood ratio. The present objective has fixed state C,
fixed peak, empirical AV (including negative values), and no sign-selection
normalization. It does not establish a physical dust inference or cosmology
correction. A generative recovery experiment must use a declared censoring or
selection likelihood, and retain rather than collapse these distance branches.

## Artifacts

- `result.json`: complete independent checks and quantitative results.
- `all-local-modes.json`: every fine-profile local branch for all four metrics.
- `level-sets.json`: saved-coordinate distance envelopes and boundary flags.
- `block-decomposition.json`: marginal/conditional decomposition of every B branch.
- `parent-decomposition-check.json`: separate check of the parent's four points.
- `verify_profile.py`: reproducible checker; no native execution.
- `manifest.json`: final hashes of review artifacts.

Input protocol SHA256:
`9e128dfc60be4579f0c1825d49b06a9e2102867f2ab2793ccd744b77286e73e6`.
Final pre-score amendment SHA256:
`a77086e65e357388df0bad735d62434ebc906bcbe685c9aa0a1426063653283b`.
Native profile array SHA256:
`a7fd4cb87e0fb825e0135d99ca0f3cc86b5369f654742e376b8f619f444d3f70`.
