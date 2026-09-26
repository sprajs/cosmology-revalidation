# Signed pre-explosion RAISIN ancestor measurements

The recovered signed DES16 ancestors make a direct baseline-noise check possible
for this RAISIN photometry cohort. Across the ten provenance-selected objects,
3,476 measurements precede the released peak date by more than 180 observer days;
1,736 are negative. Their aggregate scatter after fitting a separate constant
baseline to each object and band is close to the quoted diagonal error scale,
but this does not validate independent Gaussian errors or every baseline.
No photometry or error estimates were changed.

This uses the actual author DIFFIMG ancestors identified in the
[sign audit](raisin-flux-sign-audit.md), not the different DES SMP reduction.
It does not close the separate [SMP provenance gap](offseason-noise-protocol.md).
The cohort was fixed by source lineage before this calculation; it is not an
unselected survey population. All ten are also in each released 79-object
nominal optical, optical+NIR and NIR cosmology table; the
[membership ledger](../../runs/research_2026_09_26/astra_design/raisin_signed_refit/released-cosmology-membership.json)
pins those table hashes. The sign result and native-fitter instability
were already known. The new [protocol](../../runs/research_2026_09_26/raisin_signed_baseline/protocol.json)
was frozen before these baseline outcomes, SHA256
`fbaac7dca06161c53b6c180af4281ddbe1bcf98e91b81efdb73a9d1e92085d3d`.

## Quantities and uncertainty boundary

Keep all finite signed griz measurements with positive quoted errors and
`MJD < PEAKMJD - 180`. The declared sensitivity changes 180 to 365 days.
There is no SNR clipping or PHOTFLAG cut. Full source-line and flag records
remain in the outputs. The cuts are in observer time, not rest-frame phase.

For each object and band, write `x_i=f_i/sigma_i`, `v_i=1/sigma_i` and fit
`b=(v'x)/(v'v)`. The residual is `r=x-v*b`; report `Q=r'r` divided by `n-1`.
The tabulated baseline uncertainty `1/sqrt(v'v)` is conditional on the diagonal
quoted-error model. Its ratio to the fitted baseline is not a calibrated
survey-level significance.

The source can have a constant or changing offset from host subtraction or
reference-template contamination even before the SN explosion. A nonzero
baseline is therefore not automatically a photometric zero-point error,
astrophysical luminosity, or a correction to subtract from the SN.

For pairs in the same object and band, also report
`(f_i-f_j)^2/(sigma_i^2+sigma_j^2)` in fixed lag bins. This cancels a constant
baseline. Its reference expectation is one only under the stated independent,
correct-variance error model. Pair products of fitted residuals account for
the artificial correlation caused by fitting that baseline: for distinct
indices, `E[r_i*r_j]=-v_i*v_j/(v'v)` under the diagonal Gaussian reference.
The reported product excess subtracts that expectation. Pair counts are not
independent sample sizes.

## Descriptive results

| Pre-peak cut | Measurements | Negative | Object-band groups | Pooled Q/dof | Median group Q/dof | Group range |
|---|---:|---:|---:|---:|---:|---:|
| More than 180 observer days | 3,476 | 1,736 | 40 | 1.05553 | .99830 | .57150–2.77054 |
| More than 365 observer days | 2,856 | 1,428 | 40 | 1.02868 | 1.00730 | .61121–1.77290 |

The weighted raw baseline estimates have median −.0994 FLUXCAL and range
−2.0206 to +5.7711 for the primary cut; the sensitivity gives median −.0996
and range −1.4263 to +5.6659. These have differing measurement uncertainties
and flux scales. They must not be interpreted as a single magnitude shift.

The highest primary within-band dispersion is DES16X3zd/z: 68 epochs,
`Q/dof=2.77054`. All groups, including this one, remain in the result; no
outcome-based cut was added. A pooled statistic near one does not erase such
heterogeneity or certify the tails.

For the primary cut, same-integer-MJD pairs separated by less than half a day
have mean normalized squared difference .84838 and mean corrected residual
product +.23489. This uses 233 pairs from only five object-band groups. The
365-day sensitivity gives .84580 and +.22626 from 232 pairs in four groups.
The short-lag comparison therefore suggests shared measurement structure under
the diagonal reference, but does not estimate a universal correlation coefficient.
Distinct-night pairs have mean squared differences .99534, .98799, .99988 and
1.09853 in the fixed .5–7, 7–30, 30–180 and ≥180-day bins. Their corrected
products are +.07863, +.07121, +.06748 and −.03447. Cadence, season, error regime
and repeated records can affect these aggregates.

The signed baseline supports using an approximately correct aggregate noise
scale as an initial synthetic reference. It does not validate near-peak source
photon noise, flux-dependent quoted errors, template covariance, or the fitted
SN-population likelihood. In particular, a positive-only subset would remove
about half of these near-baseline measurements; a simulation using only the
surviving times would not inherit their positive conditional noise mean.

The [group table](../../runs/research_2026_09_26/raisin_signed_baseline/baseline-groups.csv),
[pair table](../../runs/research_2026_09_26/raisin_signed_baseline/pair-summaries.csv),
[flag strata](../../runs/research_2026_09_26/raisin_signed_baseline/flag-strata.csv)
and [source-row ledger](../../runs/research_2026_09_26/raisin_signed_baseline/baseline-source-rows.csv)
preserve the full denominators. The
[script](../../scripts/research_2026_09_26/raisin_signed_baseline.py), executed
source snapshot and [result](../../runs/research_2026_09_26/raisin_signed_baseline/result.json)
retain input/output hashes. The
[independent review](../../runs/research_2026_09_26/raisin_signed_baseline/independent-review/review.md)
passes: a separate raw-file parser and scalar calculation reproduce all 6,332
mask-expanded source rows, 80 group calculations, 331 lag/night pair bins and
242 flag strata. This verifies the computation within the stated measurement
assumptions; it does not promote the diagonal model to a validated noise law.
