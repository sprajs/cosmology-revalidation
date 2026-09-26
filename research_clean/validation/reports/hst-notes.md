# HST detector and repeat-photometry revalidation

The five HST workflows have internally reproducible numerical results, and a separate FITS-level implementation recovers the main numbers without importing their workflow or library functions. The evidence supports a discrepancy between observed repeat moments and the supplied diagonal variance model on these fixed samples. It does **not** yet identify an incorrect detector model, a correction to historical supernova photometry, or a cosmological shift.

Run the independent audit with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 research_clean/.venv/bin/python research_clean/validation/hst_checks.py
```

The audit records input hashes, its own source hash and numerical tolerances in [hst.json](hst.json). Its 2,363 checks passed in approximately 11 seconds. Comparison files are the original frozen result tables; the audit starts again from the FITS pixels, frozen sky-geometry arrays and reference tables. It independently integrates the dark circle/pixel overlaps by numerical quadrature and constructs temporal slope weights by a weighted design-matrix pseudoinverse. It does not rebuild the full science WCS operator; that remains the separate `hst-geometry` replay. Native CALWF3 execution is also outside this clean-edition replay: the two native dark FLTs are hashed starting inputs.

## What is actually being measured?

Science exposures are eight current WFC3/IR F160W FLTs, four in each of two visits. Exposures 1 and 3 define the source masks; exposures 2 and 4 supply the signed repeat contrast. The 170 retained sites share a fixed 6 arcsecond sky lattice, 0.4 arcsecond aperture, and 1.2–2.0 arcsecond background annulus. All sites pass the masked coverage rule in all eight exposures. Minimum measured aperture and annulus coverage are 0.900756 and 0.751479. The strict all-pixels-valid branch has no support. No amplitude or positivity cut is imposed, although finite-value, quoted-error and native DQ eligibility do use held-out product quality. The result therefore describes this selected support.

Weights cancel a spatially constant image to at most 1.24e-14. FLT SCI units are electrons/s, and the relative PHOTFLAM factor places the measurements on the first exposure's scale. For each held-out pair, the numerator is the sum of squared signed differences; the denominator is the sum of the two propagated **diagonal** ERR variances. No full detector covariance has been measured.

| Quantity | Search visit | Template visit |
| --- | ---: | ---: |
| Squared repeat / quoted variance | 0.7082532853 | 0.7287245596 |
| Weighted-intercept centred ratio | 0.7106205175 | 0.7375782633 |
| Supplied flat marginal fraction of pair variance | 0.0000525810 | 0.0000514067 |
| Fraction removed by same-detector-pixel shared-flat covariance | 0.0000034023 | 0.0000027410 |
| Ratio after removing only this modelled shared-flat term | 0.7082556950 | 0.7287265571 |

The independent extraction agrees with the frozen science and flat summaries at numerical precision. Reference-to-FLT mapping is verified as a five-pixel offset. Shared-flat covariance joins **detector indices**, not coincident sky positions. It assumes independent reference pixels and no cross-reference covariance. Those measured reference-error terms are far too small to account for the roughly 27–29% variance-ratio deficit. The calculation leaves unknown spatial and temporal covariance unresolved.

These within-visit FLT contrasts do not reconstruct the eventual supernova template subtraction. When two photometric epochs share a noisy template, its contribution is correlated and can cancel in their difference; adding two marginal errors without the covariance would overstate that difference variance. No shared-template covariance amplitude has been measured here.

There is one temporal repeat pair per visit, not 170 independent repeated exposures. The existing tile bootstrap is descriptive and relies on an approximation about independence between tiles. Its nominal interval is not a demonstrated coverage guarantee for a detector-wide calibration parameter. In particular, the template centred-ratio interval already reaches about 1.001.

## An additional influence diagnostic

We added fixed spatial-group deletion as a statistical stress test. The group boundaries come from the existing geometry; no outcome selects which group to remove. These are post-hoc sensitivity ranges, **not confidence intervals**, and no deletion replaces the full primary statistic.

| Full-sample measurement | Full ratio | Range after leaving one fixed spatial group out |
| --- | ---: | ---: |
| Science search, 16 sky tiles | 0.708253 | 0.681285–0.744120 |
| Science template, 16 sky tiles | 0.728725 | 0.597862–0.760653 |
| Native dark pixels, 256 detector blocks | 4.201565 | 0.720173–4.216400 |
| Native dark apertures, 16 broad detector regions | 0.547440 | 0.523321–0.563747 |

The template's most influential tile contains 25.61% of the total squared signal, compared with 11.13% for search. This makes the template result more spatially concentrated, even though neither science sample is driven by a single tile.

## Why the dark pixel and aperture results differ

The reference-only dark mask has 993,750 pixels and 247 aperture sites. The primary full-pixel ratio is 4.2015651045, whereas the independent quadrature aperture ratio is 0.5474400637. These remain distinct spatial measurements: the aperture subtracts a local background and samples only specific detector locations.

One raw-coordinate pixel `(506,945)` contributes 82.8641% of the full squared pixel difference and has zero weight in every fixed aperture. The ten largest pixels contribute 85.3588%. Four pixels with native DQ bit 32 in either image contribute 84.1749%. As a labelled counterfactual, removing those four pixels would give 0.664909; this is **not** a revised estimate because the flags are processing outcomes. It explains why a single spatial average is an inadequate summary of this image pair. The block-deletion diagnostic independently makes the same point without redefining pixel eligibility.

The [STScI calibration handbook](https://hst-docs.stsci.edu/wfc3dhb/chapter-3-wfc3-data-calibration/3-3-ir-data-calibration-steps) identifies DQ 32 as unstable and explains that it can mark ramps with at least four identified jumps. It also distinguishes count-rate products with flat/gain processing (electrons/s) from those without it (counts/s). The dark products here are in **DN per nominal second**. A DQ label is useful processing evidence, but does not by itself establish a physical origin for the outlier.

## Raw reads and temporal covariance

All fourteen unsigned RAW files were reopened, promoted to float64 before subtraction and checked against their recorded read times. Seven disjoint primary pairs and one overlapping NORMAL-pair sensitivity calculation were kept separate. Eleven exposures have `EXPFLAG=INDETERMINATE`; nominal read coordinates do not establish correct physical clocks.

For all 32 pair/quadrant combinations, an independent centred-covariance plus mean decomposition recovers the original seven-by-seven read-moment matrices. Direct slope and matrix contractions agree to 6.38e-16 in the reported squared-slope units. All 11,856 signed aperture slope values (eight pairs, 247 sites, six fixed time-weight choices) match within 1.55e-13 DN per nominal second using independently integrated geometry.

A further diagnostic compares `h' C h` with the result obtained after setting temporal off-diagonal entries of the **centred spatial moment matrix** to zero. The full/diagonal ratio ranges from 0.6683 to 2.3137 across the 32 pair/quadrant cases. Thus temporal cross terms can change even the direction of a diagonal approximation's error. These are spatial moments of paired detector realizations containing shot noise, events and detector-state changes; they are not isolated electronic read covariance or transportable ERR multipliers. The small coherent-mean fraction does not isolate the remaining processes.

## Comparison with published calibration and remaining experiments

The current science FITS headers all identify CALWF3 3.7.3 and `NLINFILE=a2412448i_lin.fits`. [Huynh et al. (2026)](https://arxiv.org/abs/2602.12110) identify that file as the February 2026 update to the pixel-based nonlinearity calibration. This is evidence that the supplied products use a recent calibration, not evidence that they reproduce the historical RAISIN reduction. The paper's reported bright-source performance does not validate the present faint-sky repeat variance.

The pinned local CALWF3 [`cridcalc.c`](../../../runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/source/pkg/wfc3/lib/wf3ir/cridcalc.c) uses an SNR-dependent weighted fit, a `21/gain` read-noise parameter labelled as correlated-double-sampling noise, and separate source/dark terms. Its initial [`noiscalc.c`](../../../runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_native_replay/source/pkg/wfc3/lib/wf3ir/noiscalc.c) uses amplifier-specific reference read-noise values. These are different operations. The six fixed raw time-weight operators are not equivalent to the native adaptive weighting, segment rejection and calibration chain. A factor-of-two argument obtained by interpreting a double-sampling noise value as independent single-read noise is a **conditional source-formula hypothesis**; the mixed RAW moments above do not yet supply the physical covariance needed to prove or refute it.

The next decisive work is to obtain several independent, timing-verified dark and science repeats; preserve native sample rejection and reference versions; estimate temporal and spatial covariance on a training set; and predict held-out aperture contrasts with the same operators. Compare reference versions in controlled native reruns to establish which changes are historical-calibration differences. Include injected faint sources and shared-template covariance before transferring any result into supernova flux or bias-correction likelihoods. Repeated visits are essential: more apertures within the same two exposures cannot substitute for independent detector realizations.

No arithmetic discrepancy in these five workflows was found. Remaining gaps concern identification, sampling and transport to cosmology, rather than a demonstrated mistake by either the clean implementation or published calibration.
