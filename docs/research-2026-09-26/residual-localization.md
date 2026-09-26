# Where the frozen residual pattern is measured

The pattern is **not confined to UV-exposed passbands or to one side of peak
brightness**. Removing the UV-exposed passbands lowers its descriptive
amplitude and especially its fixed-amplitude gain, but leaves most of its
information and a positive matched product. This is a post-discovery mechanism
diagnostic at fixed nominal covariance and nuisance tangent, not a new
significance test or a physical explanation.

The [protocol](residual-localization-protocol.md) was saved before partitioned
scores. Design-only masks/information preceded any observed-flux reads. A
passband passes an optical threshold when at most 1% of its positive
`lambda*T(lambda)` throughput lies below that rest wavelength, using measured
zHEL and a maximum 5-Angstrom observer integration step. This tests blue tails
in addition to the native fit's existing 3500–8000 Angstrom representative
wavelength cut. It is not an SED-weighted photon fraction.

Every subset uses its marginal covariance, then a fresh four-coordinate
nuisance projection. The exact43 observer vector is unchanged. Let M be its
matched residual product, I its information and G=M-I/2 its fixed-amplitude
log-score gain. M/I is a descriptive amplitude relative to that original vector.

| Retained epochs | Epoch count | Information I | M/I | Fixed gain G |
|---|---:|---:|---:|---:|
| All | 39,606 | 330.336 | .58244 | 27.233 |
| At most 1% throughput below 3000 A | 38,467 | 318.870 | .57729 | 24.644 |
| At most 1% throughput below 3500 A, primary | 33,499 | 265.710 | .51523 | 4.048 |
| At most 1% throughput below 4000 A | 27,646 | 165.990 | .49163 | -1.389 |
| Before or at peak | 10,995 | 86.236 | .65575 | 13.431 |
| After peak | 28,611 | 227.309 | .55401 | 12.277 |

All 1,020 objects remain represented with rank four after the primary 3500 cut;
1,017 have surviving rows/rank four at 4000. Pre-peak has 1,004 objects with
rows and 988 rank-four fits; post-peak has 1,020 and 1,018 respectively. The
rank-deficient objects' information is at floating-point roundoff and is kept
in the ledger. A negative fixed gain near amplitude .5 means the original
unit-amplitude template overpredicts; it does not mean the direction vanished.

UV-exposed 3000/3500 complements on their own carry effectively zero information
about this band-offset vector: the surviving band's offset can be absorbed into
the object's amplitude. The 4000 complement has only I=.00214. Dividing by such
tiny information would be misleading. Separate subset scores do not add to the
full result because their fitted nuisance directions differ.

There is nevertheless an exact nested information decomposition under the
conditional linear model. Embed the subset's efficient score weights in the
full epoch vector. Their covariance with full-data weights equals I_subset;
the difference is covariance-orthogonal to the retained score. For the primary
3500 threshold, the removed information is **64.627**, its matched product is
**55.499**, and its descriptive amplitude is **.85876**, compared with .51523
for the retained optical information. Its fixed gain is 23.186. This component
includes cross-wavelength constraints supplied by the removed epochs; it is
**not** the UV-only light-curve score. No physical-null probability is assigned.

An independent least-squares projection reproduces all 9,180 object/partition
score rows to 2.27e-13. The nested covariance identity closes to 1.46e-15, and
4,114 recorded input/output hashes verify, using preserved source/protocol
snapshots for the documented pre-outcome reporting amendment. The original
all-epoch score closes exactly. Saved artifacts include
[design information](../../runs/research_2026_09_26/residual_localization/design/result.json),
[scores](../../runs/research_2026_09_26/residual_localization/score/result.json),
[independent verification](../../runs/research_2026_09_26/residual_localization/independent-verification.json),
and the [script](../../scripts/research_2026_09_26/residual_localization.py).

These results do not refit nonlinear light curves after cutting bands, propagate
shared calibration/training uncertainty, regenerate classification/selection,
or identify dust versus intrinsic spectral evolution. They motivate retaining
broad optical spectral alternatives as well as UV-specific ones.
