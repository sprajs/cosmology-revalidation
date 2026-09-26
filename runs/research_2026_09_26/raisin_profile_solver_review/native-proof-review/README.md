# Independent native proof review

The saved one-object modern-native proof passes independent reconstruction.
No native fit or new observed-data profile was run by this reviewer.

`verify_native_proof.py` imports no Astra fitting or export code. It reparses
the raw native logs, joins source photometry by unique band/MJD, and checks
the saved arrays and hashes. All 8 protocol inputs and 44 outputs match;
the original current binary and source also match the original refit
manifest. All 65 target FITRES fields and all 893 target LCPLOT rows match
the output-only audit run exactly. The LCPLOT flags account for 74 accepted
data rows, 419 rejected data rows and 400 plotted model rows. Exact exposure
identity uses high-precision PHASE2 MJD and raw source MJD; the rounded
LCPLOT MJD is not used as an exact timestamp.

Every probe has the same 74 unique accepted exposures. Their flux/error
values equal float32 ingestion of the author-raw values; 65 are positive
and 9 negative. All common A/B flux/error strings match exactly. Fixed,
sorted and reversed native mean vectors agree exactly after mapping.

The actual positive/negative DLMAG offsets are +.150001525879 and
-.149997711182 mag after float32 seeding. Maximum **pointwise** multiplicative
mean error is 2.059e-9 relative, or 3.058e-8 quoted-error units. All other
seed coordinates close exactly. This verifies the tested distance interval,
not every AV/shape/distance combination in a future broad profile.

The full native quadratic plus prior closes within 2.274e-13 for all six
runs; the peak/other prior contribution here is zero. Independent Cholesky
inversion differs from exported C by at most 9.095e-13 flux-units squared;
all covariance matrices are positive definite. The distance probes change
C by .28351%/.37373% in Frobenius norm, confirming that the **reference** C
must remain fixed in the proposed diagnostic.

An independent Cholesky-whitened one-column least-squares solve gives
amplitude 1.0002038776176378 and DLMAG 42.57600394601274, with quadratic
93.98758585949474. It matches the executable fixed-C profile kernel within
8.89e-16 in amplitude, 7.11e-15 mag and 2.85e-14 in the quadratic. These
are arithmetic checks at fixed shape/AV/peak, not a new distance correction
or a globally minimized light-curve result.

The initial checker expected one final native export block; each log instead
contains two. `checker-amendment.json` records that failure. The corrected
checker requires both complete blocks to be exactly identical for flux,
inverse covariance, parameters/prior and setup before reading one. It does
not silently overwrite disagreeing records or loosen numerical tolerances.

## Required A/B covariance treatment

For B's 74-row anchor covariance, A is the verified 65-row positive subset
with identical shared raw errors. Use `C_A=C_B[A,A]` and then invert it.
Equivalently,

`W_A = W_AA - W_AN solve(W_NN,W_NA)`.

These independent constructions agree within 2.78e-17. The principal block
`W_B[A,A]` is a conditional precision and differs from the needed marginal
precision by 1.55e-7 in relative Frobenius norm here. Its small numerical
effect for this particular matrix does not make the operation valid.

No adjustment of common-row errors is required for the verified author A/B
pair. For a different error-convention arm, explicitly replace the changed
measurement variance on the diagonal while retaining the declared common
model/MW covariance components, verify units and positive definiteness, and
label the change. Do not substitute that arm's separately refitted native C:
it can include a different parameter-dependent anchor and changes the
experiment. A genuine flux-unit transformation must transform y, model and
the entire covariance together.

The covariance restriction is the marginal metric of the underlying frozen
Gaussian model used to replay the native fitting rule. It is **not** the
conditional noise covariance/mean after positive-flux selection. A/B minima
cannot be multiplied as independent evidence, compared as equal-dimension
likelihoods, or interpreted as a selection-corrected dust posterior.

The broad A/B protocol still requires a common fixed peak and native mean
oracle, all supported shape cells and nuisance branches, complete boundary
checks, amplitude-domain/active-prior handling, and saved lower-envelope
convergence. Native all-fixed mean probes can change prior-control flags:
`FCNCHI2_PRIOR` conditionally overwrites its AV guard when `INISTP_AV!=0`.
Therefore the mean oracle's incidental reported prior must not replace
the declared reference prior, especially near AV>5. The present proof at
AV about .054 is unaffected. Protocol acceptance is recorded separately
once the native-engine agent freezes that protocol.

## Final pre-score protocol disposition

The native-engine agent froze amendment `a77086e65e357388df0bad735d62434ebc906bcbe685c9aa0a1426063653283b` before global profile scoring. It addresses the amplitude/support, exact float32 boundary, covariance-asymmetry and completion-versus-convergence guards. The relative-output-prefix change follows a preserved initialization failure, with no observed/model/objective change. This amended design is accepted for the declared restricted pilot. No profile outcome file has been read in this review; acceptance of the design is not certification of its eventual numerical result. `final-review.json` identifies the exact accepted runner and limitations.
