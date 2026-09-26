# HST dark read moments and calibrated variance

## Purpose and inputs

These two completed calculations investigate the detector side of the signed-pixel variance question. They preserve the original fixed cohort, time weights and reference-only masks. Neither estimates a correction to supernova photometry.

`hst-dark-raw` begins with fourteen original MAST WFC3/IR RAW FITS files. It uses seven disjoint primary pairs (three search, four template), plus the prespecified overlapping NORMAL search pair `idbx41onq` → `idbx43p7q`. The latter is a sensitivity comparison, not another independent pair. Eleven exposures have `EXPFLAG=INDETERMINATE`; their physical read timing remains uncertain. Nominal SPARS50 times are coordinates for the calculation, not verified physical clocks.

`hst-dark-calibrated` begins with the two frozen FLTs produced by locally built, pinned CALWF3 3.7.3 from the first eight reads of that NORMAL pair. Those local derivatives are bundled and hashed. This command repeats extraction and analysis, **not** native calibration. The native build, reference environment and original prefix-construction/execution hashes are identified in [native provenance](../../provenance/hst-dark-native.json). Restoring these derivatives requires the bundle or the original native procedure; downloading an archive FLT is not an equivalent substitute.

The native dark processing kept DQI/ZOFF/BLEV/NLIN/UNIT/CR corrections and omitted ZSIG/DARK/FLAT/PHOT, retaining original nominal times. Its output is DN per nominal second, unlike flattened science FLTs expressed in electrons per second. Native rejection, segment selection and variance terms have not been separately identified by these comparisons.

## Fixed support

Both commands rebuild the mask from the union of the two original BPIXTAB references, before reading outcome pixels. The CALWF3 coordinate convention maps one-based reference positions to zero-based RAW coordinates by subtracting one and adding the five-pixel reference border. The active rectangle is `[5:1019, 5:1019]`; it contains 993,750 eligible pixels. Native output DQ flags and residual magnitudes never change this mask.

The 256 candidate centers are `(32+64i, 32+64j)` in zero-based RAW coordinates. Exact circle/pixel overlap uses radii 3.1189083821 pixels for the aperture and 9.3567251462–15.5945419103 pixels for the annulus. At least 90% aperture and 75% annulus coverage are required, retaining 247 positions. For reference-good mask `m`, aperture area `a` and annulus area `b`, weights are

`w = m*a − sum(m*a)/sum(m*b) * m*b`.

Area closure and `sum(w)=0` are checked. No full-aperture renormalization, output-flag rejection, or residual clipping is used. Five broad regions and a 16×16 grid of detector blocks retain spatial variation.

## Raw extraction and numerical quantities

FITS physical unsigned integers are promoted to float64 **before subtraction**. EXTVER 15 down to 9 gives the first seven nonzero chronological reads. The zero read is omitted. For later-minus-earlier read vector `d` at each pixel, retain

`mu = mean_pixels(d)` and `Gamma = mean_pixels(d d^T)/2`.

At fixed nominal times `t`, use `u_i = abs((i−3)/3)^p`, weighted mean `tbar`, and `h_i = u_i(t_i−tbar) / sum[u(t−tbar)^2]`. Powers are 0, 0.4, 1, 3, 6 and 10; power zero is primary. The slope operator must satisfy `sum(h)=0` and `h^T t=1`.

Every region checks the identity `mean[(h^T d)^2]/2 = h^T Gamma h` to absolute 1e−12 plus relative 1e−10 tolerance. Both reference-masked and geometry-only quadrant moments are retained. The saved NPZ includes quadrant/block mean vectors and matrices. CSVs contain every weighting power, signed aperture contrast, and input DQ/digital-endpoint contribution without removing those pixels. A secondary gain conversion uses the physical amplifier entries in the pinned CCDTAB, not the nominal 2.5 header gain.

The primary masked score is decomposed as `0.5*(h^T mu)^2` plus its centered spatial remainder. The coherent fraction ranges from about 0.00030% to 0.661% over all 32 pair/quadrant rows. Even the centered remainder mixes unflagged events, dark-current changes, shot noise, digitization and detector-state variation. It is not isolated electronic read noise.

## Calibrated comparison and expected outcomes

On the same fixed support, compute `d = SCI_later − SCI_earlier`, `V = ERR_earlier² + ERR_later²`, and `R = sum(d²)/sum(V)`. This diagonal reference assumes independent exposure errors; it does not assert a complete covariance model. All fixed-mask values must be finite, quoted errors nonnegative, and pair variance positive.

The expected full-pixel ratio is **4.201565**. Quadrants B/C/A/D give **0.669390 / 0.616820 / 14.723955 / 0.799958**. For aperture contrasts, compare `sum[(w^T d)²]` with `sum[sum(w² V)]`; the expected ratio across 247 sites is **0.547440**. These are different spatial quantities and must not be interchanged.

The descriptive influence table retains the ten largest squared pixel differences without altering any result. RAW pixel `(x=506, y=945)` contributes **82.8641%** of the full squared sum, has later native DQ bit 32, and has exactly zero weight in the fixed apertures. Four distinct pixels with bit 32 in either exposure contribute **8,305.450** of **9,866.900** total squared difference. These facts explain the strong spatial dependence; they do not identify a detector mechanism or justify dropping pixels to select a preferred noise scale.

## Repeat

```bash
uv run --frozen python research.py run hst-dark-raw --name my-dark-raw
uv run --frozen python research.py run hst-dark-calibrated --name my-dark-calibrated
```

Settings are frozen; both workflows reject additional configuration keys. Inputs are checked against the manifest before analysis. Raw read matrices, quadrant/block statistics and all signed aperture contractions were compared with the original executed outputs. The calibrated result was checked against the separate independent extraction and circle quadrature review. The edition's [validation record](../../provenance/history/curation-validation.json) records the actual replay and numerical tolerances.

The fourteen full RAWs account for most of the added bundle size. Keeping the original bytes permits fresh extraction and hash-checked public restoration. Native IMA histories, build trees, trial runs and review narratives are excluded. Unresolved physical covariance, telemetry and historical reduction prevent any error rescaling or cosmological propagation from these controls.
