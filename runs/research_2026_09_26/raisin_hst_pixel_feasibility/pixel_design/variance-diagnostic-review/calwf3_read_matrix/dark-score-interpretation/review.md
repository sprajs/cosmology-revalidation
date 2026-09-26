# What the completed RAW moments identify, and the next calibrated control

**The RAW experiment measures total digitized repeat-ramp variability. It does
not resolve the suspected factor-of-two read-noise normalization. The next useful
experiment is a source-faithful native first-eight-read prefix replay, starting
with the prespecified NORMAL search pair.** A fitted noise scale, lower-tail
selection or robust-mixture interpretation of these totals would skip the
important detector-processing step.

This interpretation is explicitly post-outcome. `compare_saved.py` opens only
saved sufficient arrays/tables and the already frozen six operators; it opens no
RAW image, changes no mask or read weight, and fits no parameter. Its input
hashes are in `protocol.json`. An initial path-resolution failure happened before
loading the arrays; that source is retained separately. The source-reference
contrasts are arithmetic diagnostics, not newly preregistered tests or a
calibration estimate.

## What the actual totals show

The native `cridcalc.c` `linfit` sets `rdns=21/gain` and its caller passes
`wf3->mean_gain`. For the declared seven nominal read times and mean gain2.35,
the equal-weight *reported read-only* term is .00114077 DN²/s². The hypothetical
model of independent reads with single-read noise `21/sqrt(2)` electrons gives
.000570384 DN²/s². Their factor2 relation assumes the CDS interpretation and
independent reads; it is not a measurement of the real detector covariance.

The 28 primary pair/quadrant totals instead range .0242580–.520840 DN²/s²,
**21.3–456.6 times the reported read-only term**. The separate NORMAL pair gives:

| Quadrant | Total DN²/s² | Total / reported read-only reference |
|---|---:|---:|
| B | .0231900 | 20.33 |
| C | .0564980 | 49.53 |
| A | .0459473 | 40.28 |
| D | .0207433 | 18.18 |

The amplifier-gain electron columns in the original experiment are valid
coordinate transformations, but they are not the same operation as the native
mean-gain argument. This comparison keeps that distinction explicit. It omits
Poisson, cosmic-ray, unstable dark, persistence, bias/reference and other
components; it must not be labelled observed/full-ERR. A large total does not
refute a smaller read term being overestimated.

The saved mean-vector decomposition shows that one quadrant-wide mean slope
accounts for at most .6611% of a primary total, and .114–.307% for the NORMAL
pair. That rules out *that particular constant spatial mean* as the main total
component. It does not rule out spatially varying repeat-mean changes, variable
hot pixels or localized ramp steps. The fixed aperture contractions show strong
heterogeneity: for example, the search_1 mean half-squared contrast is21.13 while
the median is.0353 DN²/s², with one value5046.66. These are descriptive summaries
of the fixed positions, not justification to remove the large values or treat
a median as read-noise variance.

There are zero selected RAW DQ or digitizer-endpoint contributions in the saved
ledger. **The cosmic-ray count is unknown, not zero:** unprocessed RAW DQ does not
supply the native ramp's later CR detection. The externally frozen bad-pixel
union does not remove every time-variable hot pixel or transient event. A CR
step has a nonzero slope projection and can dominate precisely the second
moment that was deliberately retained.

## What Gamma and the six powers constrain

For repeat read vectors Y_A,Y_B, the saved quantity is

`Gamma = E[(Y_B-Y_A)(Y_B-Y_A)^T]/2`.

For potentially different means/covariances its expectation contains
`(K_A+K_B-K_AB-K_BA)/2 + (m_B-m_A)(m_B-m_A)^T/2`.
Only independent, equal-mean, equal-covariance repeats make it K. Spatial
averaging adds another mixture of detector states. Subtracting the *global*
mean-vector outer product leaves spatially varying means in that mixture.

For every declared fixed h, `h^T Gamma h` is identified for the measured masked
population. Thus the full matrix preserves the slope direction that per-ramp
linear detrending would destroy, and it supports the prespecified alternative
weights without rereading pixels. It does not uniquely decompose electronic,
Poisson, CR, dark-state and reference components. A CDS measurement also fixes
only a difference projection of the read covariance; temporally correlated
models can share the same CDS and have very different slope noise, as the prior
synthetic AR(1) counterexamples demonstrate.

The actual power10/power0 ratios across saved pairs/quadrants span .634–1.718;
the independent-white-read model gives1.532. These fixed projections demonstrate
operator dependence, not a reason to choose whichever power produces the most
convenient variance. Native CALWF3 chooses its power from a measured SNR and may
split/reject reads. Its estimator is therefore not any one globally fixed h,
and a mixture of our six marginal contractions cannot recover that adaptive
estimator without the joint read/selection state.

The aperture/annulus contrasts retain spatial cross-terms that diagonal pixel
variances omit. They are measurements of particular spatial operators on these
ramps; they do not give all pixel-pixel covariances, independent temporal repeat
counts or an independent sample size of247. Quoted finite-field uncertainty
must retain pairs and spatial dependence. No Gaussian significance or ERR
rescaling follows from millions of pixels or overlapping NORMAL sensitivity.

## Additional source detail for the next comparison

`cridcalc.c` explicitly calls `EstimateDarkandGlow` even for DARKCORR=OMIT; its
Oct2010 source note explains the unconditional call. The routine uses a fixed
.036electron/s dark and zero amp glow. `linfit` adds the measured dy and that
fixed ddg to the reported Poisson term. Under DARKCORR=OMIT, dy can already
contain dark signal. Thus the source-faithful dark-mode ERR comparison must
retain separate read, measured-dy and fixed-dark contributions. It is not a pure
electronic-noise comparison. At the current nominal span the fixed-dark term is
about1.9% of the equal-weight reported read term, insufficient to explain the
large RAW totals but relevant to a percent-level follow-up. This is a statement
about the executed formula, not proof that the physical dark equals its fixed
value or a proposed source repair.

The fixed-weight Poisson contraction is also distinct from the source endpoint
approximation: `h^T min(t_i,t_j)h = .00357140 s^-1` for power0, versus
`1/(t_last-t_first)=.00333331 s^-1`. No observed source/dark rate was estimated
here to turn those coefficients into a claimed variance budget.

## Decisive next experiment: native prefix, first NORMAL pair

The design can be established **without using full16 later reads or asserting
unmeasured actual read times**. It is a conditional nominal-time characterization
of the source estimator.

1. Make new eight-read RAW copies from chronological first8 groups, including
   zero plus7 nonzero reads. Preserve each SCI/ERR/DQ/SAMP/TIME representation and
   pixel bytes. Renumber the reversed EXTVER sequence to8…1 and change only the
   source-required NSAMP/NEXTEND/exposure bookkeeping. Keep each recorded
   SAMPNUM, SAMPTIME, DELTATIM, TIME and relevant ROUTTIME. Derive the nominal
   prefix exposure end from its own retained metadata, not the discarded
   read16 time. Record every header edit. Keep the original EXPFLAG/TDF evidence.
   Do not use a full16 IMA, FLT, CR mask, slope or fit result.
2. Gate the constructor by applying it to an existing eight-read science RAW:
   it must preserve the raw observation arrays and reproduce the previously
   verified native SCI/ERR/SAMP/TIME exactly. Archive DQ4096 mismatch stays a
   separate unresolved workflow issue; it is not a numerical SCI/ERR failure
   and must not be silently imported into the new fixed reference mask.
3. Use the exact already-built CALWF3 3.7.3 and its frozen references, initially
   on idbx41onq/idbx43p7q. Retain the **dark headers' original switches**:
   BLEV/ZOFF/NLIN/DQI/CR/UNIT as specified, DARK/FLAT omitted and ZSIG as recorded.
   Source-match every required reference and the gain/unit output. This is a
   native dark-prefix estimator, not the F160W science recipe. No master dark is
   subtracted in this primary, avoiding that particular calibration-self-use
   channel. Reference masks and nonlinear calibration can still be shared.
4. Keep the existing externally frozen union mask and spatial operators. Report
   every finite/failed/saturated/no-fit count over that same denominator. Preserve
   all native CR decisions and output DQ as outcomes/strata, not a reason to
   manufacture a low-variance population. No new residual clipping or exclusion
   based on the resulting repeat contrast. Retain the unprocessed primary.
5. For the native slope outputs report `sum(d^2)/sum(V_A+V_B)` and
   `sum[d^2-(V_A+V_B)]`, signed repeat means, and the already declared spatial
   contrasts by pair/quadrant/block. Also retain `sum[d^2/(V_A+V_B)]` as a
   descriptive, non-Gaussian-pivot statistic. The ERR estimates and rejection
   decisions use the same ramp, so their endogeneity and conditional selection
   remain explicit. Do not compare a native adaptive slope numerator to the old
   fixed-h denominator.
6. An output-only state hook, with disabled/active SCI+ERR exact-identity tests,
   should record the native intervals/read counts, selected power and the read,
   measured-dy and fixed-dark variance terms at the frozen spatial support.
   Aggregate counts over the whole fixed mask; record detailed states at the
   pre-existing aperture supports if a full state dump is too large. It must not
   alter rejection, weights or floats. This attributes changes to declared
   processing operations without identifying their physical correctness.

The initial NORMAL pair is an engineering and conditional-repeat control, not
an outcome-selected replacement for all14. After unit/state/mask gates, extend
the same operator to the frozen3+4 pairs and preserve EXPFLAG strata. The11 TDF
warnings remain unresolved; the units are DN²/nominal-s² unless further clock
information is supplied. Their printed prefix times match, but that does not
certify actual physical clocks. NORMAL-only still has one search pair and no
template pair, so it cannot produce replicated epoch-specific noise calibration.

Propose no new downloads: all inputs exist. A separately frozen initial stage
can cap native execution at120seconds for the constructor control and NORMAL
pair, with1GB new outputs, one worker and partial results on failure. Estimate
full-cohort IMA/FLT disk cost before extension rather than delete provenance or
relax caps after outcomes. This proposal is not an execution release.

## What this resolves and what remains next

If native processing removes most RAW excess while its reported variance does
not follow the repeated native slopes, the discrepancy is localized to a
specified estimator and detector state. That is much more informative than the
current total-RAW/white-read comparison. It still mixes photon/dark noise,
correlated electronic noise, rejection conditioning and reference uncertainties;
a numerical ratio alone cannot identify the factor-two hypothesis.

A subsequent *separately labelled* F160W science-switch/reference control would
measure transport across ZSIG, DARK and FLAT processing. It must not be called
historical RAISIN reconstruction or use a self-trained dark reference as
independent evidence. Reference-pixel contrasts can help isolate electronic
common modes, but reference pixels are a different detector population and do
not by themselves identify active-pixel slope covariance. A read-level CR/noise
mixture would need a frozen physical step model and recovery checks; fitting
one post hoc to recover the expected read variance is lower priority than the
now-feasible native prefix test.

The machine-readable proposal is `prefix-control-design.json`. Source `unitcorr.c` and `doir.c` explicitly give COUNTS/S for the dark-default FLATCORR=OMIT output, hence DN per nominal second; only the science flat branch multiplies by mean gain and uses ELECTRONS/S. The executor must assert those actual units. Root independently re-extracted all8 RAW pairs after this review: its `dark_ramp_execution/root-raw-review/result.json` reproduces every matrix/mean exactly and all saved quadrant/block/aperture contractions within its floating tolerances. That validates extraction arithmetic, not physical attribution.
