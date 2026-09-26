# Public training-input overlap diagnostic

Freeze before examining residuals stratified by public training-list membership.
The nominal DES release documentation identifies K21 as its training sample,
with recalibration and observer-U exclusions. A newly acquired archive from the
pinned public SALTShaker commit supplies config-driven input lists and headers.
This is not an execution-linked roster of accepted DES5YR training observations.

First identify DES header IDs from the K21 input lists, verify exact same-ID sky
positions against the released DES HEAD, and join the existing 43 discovery and
1,020 validation IDs. Retain ambiguous or mismatching joins explicitly; no
outcome-based correction of membership. Write membership and hash it before
computing any stratified residual statistic. The fixed original 43-derived
observer vector remains unchanged, even if discovery objects have training
overlap. Do not reselect or retrain the discovery model.

Report original projected fixed-vector I, M and G, amplitude M/I and field
contributions separately for validation objects present and absent from the
public K21 DES input list. Every original validation object must appear exactly
once. These descriptive subsets are not independent-training tests: common
calibration and global trained surfaces remain, non-DES physical aliases have
not all been ruled out, and the precise accepted training execution is unknown.
Do not report a new global p-value or choose a subset as the primary result.

This tests whether the previously transferred residual is numerically confined
to objects reused in the identifiable public DES training input. It does not
certify a fully independent sample, prove model overfitting, or infer cosmology.
