# Independent ten-object fixed-C cohort verification

**All ten saved numerical profiles pass independent verification.** This
certifies the restricted fixed-state calculation and its recorded gates;
it does not identify a physical distance correction. No new native fits
were launched by this review.

## Provenance and arithmetic

The final cohort manifest has **560 files and14 executed source files**;
every hash verifies. Its SHA256 is
`bfa57a6caa3aae1a9af3f05685988d0978c65591ba2bc1a800804800831d3688`.
The active-attempt index selects the immutable original C1 profile, corrected
C3/v2 profile, and remaining v3 profiles. The original wrong-object C3
initialization is retained and is not treated as a scientific result.

Each provisional per-attempt manifest hashed `runner.log` before the engine
printed its final result. For all ten, removing **exactly** the final
`json.dumps(result, indent=2)+'\n'` restores the recorded prefix hash.
The independent per-object ledgers preserve that explanation and both
hashes; the final manifest certifies the complete closed logs. No unaccounted
array, source or measurement-file mismatch remains.

The checker independently verifies:

- Original source-row membership with multiplicity and native float32
  flux/error ingestion; exact reference/fixed-reference/stream-final
  band, MJD, flux and error identity; exact fixed parameters and independent
  pre-query native-prefix mean/precision identity. The alternate C anchor
  retains the same B rows. Current/audit nominal reference fits close at the
  declared tolerances.
- Each full covariance against the raw native precision export and its
  objective quadratic. A uses `C_B[A,A]` or the corresponding alternate
  marginal covariance before inversion. Identical duplicate keys retain
  multiplicity; all required covariance-swap tests return exact zero error.
  Covariances from different native parameter states are not required equal.
  Swap invariance establishes numerical exchangeability under the supplied
  covariance, not physical independence of duplicated measurements; a future
  noise generator still needs their exposure-identity/covariance ledger.
- **670,988 saved native means × four metrics**, using independently
  reconstructed amplitude, DLMAG and direct residual quadratic. Maximum
  absolute discrepancies are **5.457e−12 in Q, 7.106e−15 mag in DLMAG,
  1.599e−14 in amplitude**.
- All coarse/fine shape and AV coordinates, explicit float32 algorithm
  boundaries and neighboring probes, branch CSV arithmetic, finest minima,
  amplitude/support/normal-equation gates, edge flags and cached minima.
  The largest coarse/fine changes are **0.000456797 in Q** and
  **0.000701446 mag**, within the frozen0.001 thresholds.
- Counts of **712,924 native replies**, with every one of the **20,958
  competitive actual-amplitude reply pairs** independently compared with
  its cached reference mean and analytic scale. Maximum relative discrepancy
  is **4.633e−9**, below2e−8. The full raw-vector pilot audit remains linked;
  the remaining cached mean vectors are independently scored here without
  redundantly parsing every raw full-vector reply. All baseline replays are
  exact.
- All ten branch/level-set outputs and cohort summary fields. Every saved
  ΔQ1/4/9 amplitude interval is reconstructed from the quadratic, with the
  distance and native-mean domain intersections and each boundary flag.

The initial checker reused the pilot's two-callback assumption. C3's fitted
reference has one complete final callback, while its fixed probes have two
identical copies. The revised parser requires one or two complete blocks and
verifies exact duplicate equality when present. The initial code and explicit
checker amendment are retained; no numerical tolerance was loosened.

## Conditional profile shifts

These are differences between minima of the B signed-row and A positive-row
objectives, under the common B-anchor covariance rule. They are not estimated
population biases or corrections.

| Object | Added negative rows | B−A minimum DLMAG, mmag |
|---|---:|---:|
| DES16C1cim | 9 | −164.949 |
| DES16C3cmy | 1 | −2.569 |
| DES16E1dcx | 1 | −0.015 |
| DES16E2clk | 0 | 0 exactly |
| DES16E2cqq | 0 | 0 exactly |
| DES16S1agd | 5 | −49.376 |
| DES16S1bno | 6 | −19.045 |
| DES16S2afz | 2 | −1.431 |
| DES16X3cry | 4 | −4.753 |
| DES16X3zd | 5 | −3.156 |

Both zero-change controls reproduce **identical entire A/B metric arrays**
within each anchor, not merely nearly equal selected minima. The nearest
fine-grid B-branch comparator is a different numerical quantity: it can be
slightly displaced from the continuously refined A minimum even when the
objectives are identical. DES16E2cqq's −0.366787 mmag matched-branch offset is
such a grid-resolution artifact, not a processing effect. Its actual global
minimum difference is exactly zero.

Branch ambiguity remains central. DES16C1cim has ten sampled B local branches;
the high-stretch competitor has almost the positive-arm distance at only
ΔQ≈0.116 above the B minimum. The frozen nearest-shape comparison gives
+14.172 mmag, versus −164.949 mmag for the global-minimum comparison.
DES16S1agd's nearest-shape comparison is −30.856 mmag and costs only
ΔQ=0.0822 above its B minimum. These are nearly tied empirical solutions,
not unique dust/distance determinations.

Every ΔQ1 level set remains inside the declared shape box, but DES16E1dcx
already reaches a shape edge at ΔQ4. At ΔQ9, shape edges are reached for
DES16C1cim, DES16E1dcx, DES16S1agd, DES16S1bno and DES16S2afz. The reported
broader envelopes are therefore restricted by that domain; no extrapolated
support or confidence calibration is supplied.

All forty metrics pass the executed minimum/edge gates, yet this is still
finite-grid numerical validation rather than a continuum-global proof.
The covariance is held at a native state, the peak is fixed, AV remains an
empirical parameter that may be negative, and the underlying Gaussian
quadratic lacks sign-selection normalization. Minimum Q values across A/B
row sets must not be interpreted as a likelihood ratio. The next scientific
gate is the separately frozen generative/likelihood recovery experiment,
followed by native-template recovery that preserves competing branches and
measurement/selection uncertainties.

## Reproduction

`verify_cohort.py` checks active object artifacts, raw measurement membership,
covariances, saved vectors and geometry. `final_certify.py` verifies the
closed final manifest, aggregate fields, native amplitude replies and zero
controls. Neither imports the owner's fit implementation or launches a
native executable. `final-certification.json` is the concise machine
certificate; `DES*.json` files hold per-object evidence. `protocol.json`,
`checker-amendment.json` and `manifest.json` preserve the independent review's
design and versions.
