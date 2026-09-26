# Proposed current-product HST aperture-noise experiment

2026-09-26. Design only, before any SCI/ERR values or aperture outcomes were read by this reviewer. Root owns implementation/release. The fixed inputs are the acquired eight F160W FLTs (four first-visit, four template-visit) and two DRZ products for DES16E1dcx. This tests an explicitly defined extraction operator on public2026 products. It cannot reproduce private RAISIN photometry, establish a detector count-rate correction, or by itself validate near-peak SN uncertainties.

**Primary recommendation:** signed, same-sky aperture differences between independently exposed native FLTs within each visit, with candidate locations/source masks determined from other exposures. This distinguishes erroneous variance scale and spatial correlation from a fitted SN mean. Drizzled-image noise and common-template differences are secondary operator diagnostics, with their dependence retained.

## Metadata and resource boundary

[Independent structural validation](/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/raisin_hst_pixel_feasibility/pixel_acquisition/independent_validation/review.md) finds1014×1014 FLT SCI/ERR/DQ/SAMP/TIME arrays. SCI/ERR use electrons/s; FLT WCS is TAN-SIP. The DRZ products have SCI/WHT/CTX and final TAN WCS, with four listed FLT contributors, square kernel, pixfrac1 and0.12825arcsec output pixels. WHT says UNITLESS; final drizzle masks named by the headers are absent. No read of these arrays occurred in this design.

Use one worker, initially ≤120s of extraction computation and ≤100MB of new arrays, after a synthetic small-array/unit gate. No redrizzling, image fitting, new survey fields or external catalogs in the minimal primary. Pin packages/scripts/masks/operators before outcomes. Stop rather than widen masks or choose a more favorable aperture if support or mapping fails. A changed implementation with scientific choices requires a dated pre-score amendment.

## Fixed cohort, geometry and masking

1. Order the four exposures of each visit by EXPSTART, breaking ties by filename. Exposures1&3 form the **design pair**; exposures2&4 form the **primary measurement pair**. This split interleaves time/dither changes. The exposure-order reversal is one declared, dependent sensitivity; never select the better split or treat it as a fresh sample.
2. Construct a6-arcsec square lattice in a tangent plane centered at the published coordinate(8.9685°,−43.35812°), aligned east/north, origin at that coordinate. Keep positions whose full2-arcsec annulus lies inside all eight FLT footprints with an additional2-pixel detector-boundary margin. Exclude a fixed5-arcsec radius around the SN in every branch. Use the same sky coordinates throughout; no recentring on test-image peaks.
3. The primary physical aperture has radius0.4arcsec; sky annulus1.2–2.0arcsec. Use native-pixel fractional overlaps and the complete SIP/WCS, including LTV/subarray conventions. Require a verified official PAM/orientation or a separately certified equivalent area Jacobian; Sol is auditing that geometry without pixel values. Do not silently replace it with a constant. Aperture/annulus overlap integration must pass doubled-subpixel-resolution synthetic flux/variance closure to0.1%; otherwise refine the numeric quadrature only, retaining both outputs.
4. Design-image source masks use only exposures1&3 from each visit. In each native image subtract a64×64-pixel tile median background (tiles at detector origin; partial edge tiles retained). Flag a pixel when the absolute3×3 box-summed residual exceeds5 times the square root of the summed ERR². Do not fit an empirical noise renormalization. Map the union of flagged pixels from all four design exposures to sky and dilate by2.4arcsec; reject lattice centers falling within that union. This is a deterministic conservative source/artifact mask, not a calibrated detection significance. Broad diffuse emission may remain; record design-image background strata, do not relabel all accepted positions as physically empty sky.
5. Reject a measurement at a position only for invalid geometry, nonfinite SCI/ERR, nonpositive ERR, or DQ≠0 in any native pixel with nonzero aperture/annulus coefficient. Apply the same DQ rule to design support. This is an intentionally conservative pipeline-good-pixel estimand, not validation of rejected epochs. Preserve every exclusion reason and DQ-bit count. Because native DQ itself is data-dependent, the test remains conditional on pipeline acceptance; no claim of unconditional detector noise follows.
6. No cut uses the signed aperture flux, test residual, measured test SNR or extreme test statistic. Retain negative flux and all tails. No sigma clipping or manually selected empty patches. Background/source-brightness strata derive only from the design pair (fixed terciles after eligibility), and are descriptive.

First release the frozen candidate/mask/support ledger **before computing measurement-pair aperture sums**. Target at least100 accepted positions spanning at least12 of16 equal tangent-plane footprint tiles, with at least3 positions in each occupied tile used for resampling. If less, report the shortfall;≥30 positions permits only a limited descriptive ledger, and<30 ends the noise-scale comparison. Do not loosen the predeclared selection. Footprint partition edges depend only on WCS bounds.

## Native photometric operator and quoted variance

For exposure `e`, let `P_ep` be the verified PAM, `a_ep` the aperture fractional overlap and `b_ep` the annulus overlap. Estimate background as the unweighted-by-ERR fractional-area mean of the **original flat-fielded SCI** in the annulus:

`B_e = Σ b_ep SCI_ep / Σ b_ep`.

Then define the signed aperture sum

`f_e = κ_e Σ a_ep P_ep (SCI_ep−B_e)`.

Here `κ_e=PHOTFLAM_e/PHOTFLAM_ref`, using one fixed first-visit reference exposure. This puts all epochs on the same present-day flux-density convention; it is not a CRNL correction. Build the whole linear weight vector explicitly: aperture coefficients `κ aP`, annulus coefficients `−κ(ΣaP)b/Σb`. It sums to zero before multiplying SCI and therefore removes a constant native-SCI sky exactly. The same annulus pixels must not be separately treated as independent fitted background noise.

The diagonal pipeline comparator is

`V_e,diag = Σ_p w_ep² ERR_ep²`.

PAM and κ act on SCI and ERR consistently. The test evaluates this diagonal covariance approximation; it does not assume native pixels are independent in reality. Shared dark/flat/astrometric uncertainties and correlated read/ramp noise can violate it. FLT ERR is already the native ramp-fit uncertainty; do not add source/sky Poisson or a textbook read-noise term again. Conversely DRZ IVM background weights need not include source Poisson. The [WFC3 calibration documentation](https://hst-docs.stsci.edu/wfc3dhb/chapter-3-wfc3-data-calibration/3-3-ir-data-calibration-steps) assigns final FLT ERR to the ramp fit. The [file-structure documentation](https://hst-docs.stsci.edu/wfc3dhb/chapter-2-wfc3-data-structure/2-2-wfc3-file-structure) describes SCI/ERR/DQ roles. Report the inherited possible dependence of ERR on the same measured signal; `d/sqrt(V)` is not an exactly Gaussian pivot merely because the pipeline outputs V.

Synthetic closure must verify constant-SCI cancellation, a known signed injected aperture signal, PAM/κ unit scaling, and direct versus assembled-weight variance including the annulus. A fixed-flux injection checks extraction arithmetic only; it cannot validate detector CRNL, photon noise, persistence or empirical covariance.

## Primary estimands: two visits, not hundreds of independent observations

At each retained position `j`, separately for each visit,

`d_j=f_2j−f_4j`, `V_j=V_2j,diag+V_4j,diag`.

Record all `d_j,V_j,z_j=d_j/sqrt(V_j)`. Primary outputs are signed standardized mean, uncentered `Σz²/N`, and a centered variance ratio `Σ(d−a)²/V /(N−1)` with `a=Σd/V /Σ1/V`. The fitted intercept distinguishes a coherent pair offset from fluctuation scale; it is not silently subtracted from the raw ledger. Compare paired aperture contrast covariance versus sky separation in fixed6–12,12–24,>24arcsec bins, and give detector-position/background strata. Report full distributions/tails without selecting thresholds after scores.

Give a16-spatial-tile delete-one-tile range and block-bootstrap descriptive interval (seed fixed in protocol;2000 resamples). These quantify sensitivity to the sampled field, not physical-null p-values: two visits cannot calibrate general temporal covariance, and shared detector calibration errors need not obey independent tile sampling. Do not call `Σz²` a chi-square test with N independent degrees of freedom. First/template results remain separately visible, alongside a count-weighted descriptive pooled value only.

Declare radius0.2 and0.8arcsec sensitivities on the same frozen locations/masks and annulus, plus the design/measurement split reversal. They are dependent diagnostic summaries, not a family of discoveries. No uncertainty rescaling, error-floor fit, rejected-position recovery or final cosmology correction follows from a variance mismatch.

## Shared-template control with the same measured photons

On the common primary mask retain the vector `x=(S2,S4,T2,T4)` and diagonal comparator `D=diag(VS2,VS4,VT2,VT4)`. Every contrast is a fixed matrix `H x` with covariance **`H D Hᵀ`**, retaining the off-diagonal entries. A common reference `Tbar=(T2+T4)/2` gives differences `(S2−Tbar,S4−Tbar)` with predicted off-diagonal `(VT2+VT4)/4`; separate references `(S2−T2,S4−T4)` lack that term under independent-exposure noise. The same-template difference cancels T exactly and returns the primary search-pair contrast.

Report a secondary cross-product diagnostic

`(S2−Tbar)(S4−Tbar)−(S2−T2)(S4−T4)`

whose independent-noise expectation is `(VT2+VT4)/4` when within-visit means match. Save its algebraic decomposition `.5(S2−S4)(T4−T2)+.25(T2−T4)²`. It is a contrast of the same four measurements, not independent covariance validation. Within-visit PSF/astrometric/flat drift or cross-exposure correlations change the expectation; between-visit static-scene mismatches are reduced by this contrast but must not be assumed absent in raw science-minus-template fluxes. This test cannot infer the full cross-SN template covariance or establish a survey-wide correction.

## Drizzle and source/host limitations

[Drizzle's noise documentation](https://hst-docs.stsci.edu/drizzpac/chapter-3-description-of-the-drizzle-algorithm/3-3-weight-maps-and-correlated-noise) explains why output pixels are correlated and why pixel RMS scaled by aperture area misses covariance. [AstroDrizzle's weighting description](https://hst-docs.stsci.edu/drizzpac/chapter-5-drizzlepac-software-package/5-2-astrodrizzle-the-new-drizzle-workhorse) distinguishes exposure, ERR and inverse-background-variance weighting. Therefore no `1/WHT` variance claim is permitted until the actual weight type and contributing masks are linked. If that gate remains open, DRZ may supply only a separately labelled signed-aperture spatial diagnostic, with no variance calibration or exact FLT→DRZ closure assertion.

Public DRZ uses all four exposures, including mask-design data, so its mask is not independent of its own noise. A DRZ comparison is secondary and selection-conditioned. Reconstructing disjoint split drizzles requires a later source/config/mask freeze, not an undocumented interpolation in this primary.

Blank-aperture tests probe the chosen background regime. They cannot establish SN photon-noise validity, host subtraction residuals, aperture corrections or centroid uncertainty at the actual SN. A separately declared known-SN-position diagnostic can retain four signed per-exposure measurements without SNR cuts, but one object/four exposures has little covariance information; PSF/breathing, residual alignment, real short-time light-curve change and host gradients compete with noise. Do not use blank-field success to calibrate those terms. WCS/PAM, common calibration references, persistent detector patterns and background subtraction remain explicit nuisance sources; no parameter is retuned until a supported mechanism is identified.

The minimum successful deliverable is a frozen mask/geometry ledger, complete signed aperture/variance tables, primary pair contrasts with spatial dependence, and the common-template covariance check. It is a present-day extraction-noise measurement, not historical reproduction, a CRNL test or a bias-corrected distance likelihood.
