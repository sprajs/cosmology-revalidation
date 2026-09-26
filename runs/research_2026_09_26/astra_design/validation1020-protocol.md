# Fixed-template transfer to previously unscored high-Ia DES objects

Frozen after the 64-object exploratory signal, the 43-object high-Ia secondary,
and its matched modern SNANA-mean/covariance check. This is an adaptive
validation design, not a pristine survey preregistration. Select every released
DES object with PROB_SNNV19>.999 that is outside all 64 discovery IDs: 1,020
objects, fixed from metadata before inspecting their residual scores. The
published mask is retained. No object is excluded on its residuals or fit score.
CID order, published parameter initialization, fields, epoch counts, coefficient
vector and covariance are saved and hashed before executing the analysis.

Primary model: exact exported current pinned SNANA accepted-epoch mean and
frozen covariance, with OPT_MWCOLORLAW=99 (exact F99), same SALT3, passbands,
mask and fit options as exact43. Use independent native local shape/colour/time
Jacobian at exported SNANA coordinates, exact amplitude derivative, and exact
SNANA-mean observer templates. Project all four fitted SALT directions. Retain
this tangent approximation explicitly. One-to-one match every epoch by band,
time and observed flux, retaining duplicate multiplicity; fail rather than
silently drop mismatches. Retain frozen noise and selection limitations.

The only primary prediction is the signed observer-band vector inferred from
43 discovery objects with exact SNANA mean/C and the previously primary ridge
scale 0.02 mag. In griz zero-sum gauge it is approximately
[-.00522845988, .00783720569, .01926460352, -.02187334933] mag of extra dimming.
Freeze its full-precision three-coordinate mean and shared posterior covariance
before validation. Do not estimate, rescale, reverse, select ridge scales or
choose a different template family from the new objects. This template is a
calibration-shaped empirical prediction, not a physically identified zero-point
correction.

For each held-out object let r be its nuisance-projected whitened residual,
T its three observer templates, v=T*c_discovery, a=r dot v, I=v dot v.
Primary predictive log-score gain relative to zero template is G=sum(a-I/2).
The fixed-conditional-Gaussian signed matched filter is Z=sum(a)/sqrt(sum(I)),
with one-sided tail; G is not converted by sqrt(2G). Report absolute projected
chi-square, epoch count and coverage, per-object and per-field (a,I,G), and the
information-weighted transfer amplitude a/I only descriptively, never used to
retune the primary prediction. Enumerate all 2^10 whole-field residual sign
flips and use the one-sided tail of sum(sign_field*a_field), keeping the fixed
prediction/information. All field contributions are retained. Conditional
independent-field sign symmetry can fail under common SALT training,
calibration, published clipping or selection, so it is not a universal p-value.

Secondary: integrate the single frozen three-dimensional discovery posterior
once over the joint 1,020-object validation likelihood, preserving its shared
cross-object uncertainty. Report that joint predictive gain and fixed
conditional-Gaussian calibration; do not sum marginal object predictive gains
that discard shared covariance. Do not select between this result and the
fixed-vector primary. For transparency retain all discovery and earlier
measurement-only failures, existing historical-F99 stress, and limits on
physical noise validation. No corrected cosmology, BBC or SALT retraining is
licensed by this residual transfer alone.

## Duplicate-epoch secondary, frozen before validation residual inspection

The existing audit found exact duplicate groups in 68 DES objects. Tag exact
(CID, band, MJD, observed flux, quoted error) duplicates in the accepted export,
report group/epoch/object counts and their full-sample contributions. As a
secondary sensitivity keep the first occurrence in the stable epoch ordering,
take the corresponding marginal covariance submatrix, and recompute the SALT
nuisance projection and the same fixed-template scores. This does not assert
that every duplicate is erroneous, change the primary published mask, or
authorize deleting data. No tuning or exclusion based on residuals is allowed.

## Frozen rest/phase comparison, added before validation residual inspection

While the 1,020-object native export is running, and before validation residual
scores are read, freeze a second prediction from the exact43 discovery design:
the three previously defined rest/phase coefficients, same fixed 0.02 ridge.
Use native host RV=3.1, host RV=2.0 derivatives and colour*tanh(phase/20), the
same definitions as discovery. Save this mean and its full 3D posterior
covariance separately without altering the frozen observer coefficient file.
Compute its fixed gain and single shared posterior predictive gain on the same
validation objects. Report the paired fixed observer-minus-rest gain, its
conditional Gaussian scale and field-sign distribution using the actual
correlation of both predictions. Retain both outcomes and no family winner
selection or retuning. A predictive preference only distinguishes these two
specific regularized empirical forecasts, not all rest-frame versus observer
systematics. Near-collinearity and shared-model errors remain material.
