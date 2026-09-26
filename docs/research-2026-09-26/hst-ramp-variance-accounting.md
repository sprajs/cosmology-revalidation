# HST repeat-noise result: accounting for the quoted variance

The [current masked-aperture experiment](raisin-hst-signed-pixel-noise.md) found
repeat second moments smaller than the sum of the quoted diagonal pixel
variances. This follow-up examines possible explanations after that result; it
is exploratory and does not change the original mask, sample or error arrays.
It is not yet an identified calibration correction.

## Source and execution boundary

The public images declare CALWF3 `3.7.3 (Jan-07-2026)`. The pinned source commit
`6a1147d7bdccb7e2a7b73af276f4597c6420fbb8` declares exactly that version. Its Git
blob identities and downloaded bytes are preserved under
`runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_variance_source/`.
This matches the declared source version, not the archive executable's exact
build. The current upstream main is already version3.8.0 and is retained only
as a separately identified source comparison.

CALWF3 first initializes read/Poisson errors, applies nonlinearity and dark
corrections, and then fits the read ramp. For ordinary accepted multi-read
intervals, `linfit` constructs a new slope error from its read-noise expression
and estimated source/dark count increments. The input ramp errors also affect
SNR-dependent weighting and rejection. The final flat-field division then adds
flat-reference uncertainty explicitly. Consequently, the handbook's general
error-propagation description cannot simply be treated as an exact decomposition
of every final FLT ERR pixel. Dark-reference uncertainty in intermediate ramp
arrays is not automatically the same additive term in the final slope error.

The source links are the pinned
[ramp fitter](https://github.com/spacetelescope/hstcal/blob/6a1147d7bdccb7e2a7b73af276f4597c6420fbb8/pkg/wfc3/lib/wf3ir/cridcalc.c),
[step ordering](https://github.com/spacetelescope/hstcal/blob/6a1147d7bdccb7e2a7b73af276f4597c6420fbb8/pkg/wfc3/lib/wf3ir/doir.c),
and [image arithmetic](https://github.com/spacetelescope/hstcal/blob/6a1147d7bdccb7e2a7b73af276f4597c6420fbb8/pkg/wfc3/lib/wf3ir/math.c).

## The supplied diagonal flat component is small

The exact header-named PFLT and two visit-specific DFLT reference files were
acquired from CRDS. Their 1024-square planes have LTV=(5,5), while the trimmed
1014-square FLTs have LTV=(0,0); unit LTM gives the explicit five-pixel reference
offset. No alternative flat or empirical flat amplitude was fitted.

For each unchanged masked operator, let `w` be its full signed aperture/annulus
weight and `y` the calibrated SCI value. The pipeline's diagonal flat term is
`Σ(w y)²[(ERR_P/P)²+(ERR_D/D)²]`. Under independent reference pixels, two
exposures using the same flat also have covariance from overlapping **detector**
pixels, with coefficient `(w2 y2)(w4 y4)` multiplying that same relative
reference variance. Coincident sky positions are not automatically coincident
detector pixels.

| Visit | Flat marginal contribution / quoted pair variance | Fraction removed by same-pixel shared-flat covariance |
|---|---:|---:|
| Search | 0.0000525810 | 0.00000340231 |
| Template | 0.0000514067 | 0.00000274097 |

These particular diagonal terms are far too small to account for the roughly
0.28 fractional repeat-variance deficit. This is not a bound on unknown
cross-pixel flat covariance or other calibration systematics. The supplied
reference ERR planes do not specify those covariances. Full arithmetic,
per-exposure and per-position records are in `flat_component/`. Independent
arithmetic checks reproduce all1,360 exposure operators and340 paired
covariances, with maximum numerical gap2.84e−14; the coordinate mapping and
source formula were checked separately.

The same support audit records SAMP/TIME without changing eligibility. More
than99.5% of valid operator-pixel occurrences have SAMP8; the others are retained
and include shorter accepted ramps. Native source includes the zeroth read in
the reported SAMP count even though the ordinary slope fit omits it. Thus a
seven-point nonzero-read calculation is relevant to the dominant path, but it
is not an exact reconstruction of every selected pixel's read history.

## A conditional mismatch to the independent-read model

The pinned fitter labels21electrons as a correlated-double-sample (CDS) noise
and uses it in the slope fit's read-variance normalization. STScI defines CDS
noise as the noise of a **difference of two reads**. Under the additional
assumption that the two reads have independent equal noise, each read has
variance21²/2; subtracting a common zeroth read adds a common covariance term
that cancels from a slope fitted with an intercept. The CDS definition alone
does not establish that independence, especially across different read
separations. [STScI detector documentation](https://hst-docs.stsci.edu/wfc3ihb/chapter-5-wfc3-detector-characteristics-and-performance/5-7-ir-detector-characteristics-and-performance)

An [independent covariance challenge](../../runs/research_2026_09_26/raisin_hst_pixel_feasibility/pixel_design/variance-diagnostic-review/calwf3_read_matrix/review.md)
makes this ambiguity explicit. Different illustrative temporal covariance
matrices can share the same adjacent CDS variance but produce reported/actual
slope-read-variance ratios2.000,1.154 or0.648. None was fitted to the detector.
Even the direction of an error correction therefore cannot be deduced from the
CDS convention alone.

For seven equally spaced reads and each fixed native weighting-power branch,
the [synthetic source audit](../../runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_variance_source/linfit_algebra/review.md)
finds the reported read-variance term to be2.00–2.20 times the exact weighted
slope variance under that independent-read model. Its endpoint-count Poisson
term is0.933–1.000 times the exact cumulative-Poisson term. These ratios compare
specified algebraic models. They neither measure the real detector covariance
nor justify dividing the observed ERR array by a fitted factor.

Real temporal and spatial read correlations, reference-pixel subtraction,
SNR-selected weight branches, cosmic-ray/read rejection and uncertain ramp
means remain relevant. The decisive next input is repeated matched-cadence
ramps whose covariance can be measured without fitting away the slope-noise
direction. Independently detrending every ramp would erase precisely the
quantity needed to test the slope-error model.

This concern has relevant prior instrument research.
[WFC3 ISR2009-23](https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs/_documents/2009/WFC3-2009-23.pdf)
discusses cadence-dependent noise and detector-history effects. Its effective
noise calculation uses formal uncertainties from individual ramp fits. That
description does not identify the full repeated-ramp slope covariance needed
here; treating it as validation of the current error plane would additionally
require estimator, processing-version and covariance closure. This is an
inference about validation scope, not evidence that the historical measurement
was wrong.

The [metadata feasibility review](../../runs/research_2026_09_26/raisin_hst_pixel_feasibility/pixel_design/variance-diagnostic-review/calwf3_read_matrix/cadence-fallback-review.md)
keeps matched sampling separate from approximate substitutes. The first bounded
search inspected80 nearby long-exposure headers without finding full-frame
SPARS50 ramps; it had not exhausted the nearby candidate list. The separately
recorded [additive metadata continuation](raisin-dark-ramp-metadata.md) found
six search and eight template full-frame SPARS50 ramps with16 reads, compared
with8 science reads. A longer SPARS50 RAW ramp
could supply a valid early prefix if its actual read times, reset convention
and preceding history match. A STEP50 tail or a SPARS25 subarray instead
requires an additional transport assumption about detector covariance.

In particular, the desired variance is `hᵀ K_target h`. Measurements of a
different cadence determine its own covariance, not `K_target`, unless a
relation between the two is established. The
[official readout tables](https://hst-docs.stsci.edu/wfc3ihb/chapter-7-ir-imaging-with-wfc3/7-7-ir-exposure-and-readout)
and [FITS sample metadata documentation](https://hst-docs.stsci.edu/wfc3dhb/chapter-9-wfc3-data-analysis/9-5-specific-tools-for-the-analysis-of-wfc3)
distinguish proposal NSAMP (excluding the initial read) from FITS NSAMP
(including it). Rounded exposure time alone is therefore a
search aid, not evidence of an exact ramp match. No mismatched dark exposure
has been treated as a measurement of the supernova cadence's read noise.

The historical RAISIN photometry and simulation-error bridge remain separate.
In particular, its traced NIR simulation helper uses an empirical image
background estimate; the presence of a discrepancy against current FLT ERR
would not prove that the same term was used in the published supernova errors.
No distance shift, selection correction or cosmological refit is inferred here.
