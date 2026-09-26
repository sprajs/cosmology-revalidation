# Pre-outcome amendment: declared masked-pixel secondary

2026-09-26. Root requested this amendment before SCI/DQ support counts, aperture sums or noise scores. It is additive to `protocol.json` (SHA-2566768b7db03776d575d695efe1c4bd2da53256e2a10b788ae9cce323bc7a91806); that file and its README remain unchanged. No array values or outcomes were inspected to choose thresholds.

The original **strict full-coverage primary**, source-mask rules, fixed grid, exposure split, support target and stopping rule remain in force. An inadequate strict sample is reported as an inadequate primary. The separately declared branch below may run if its own support passes; it never substitutes for or retroactively relaxes the primary.

## Secondary operator

Use only radius0.4arcsec, annulus1.2–2.0arcsec, original design exposures1&3 and measurement exposures2&4. The same candidate sky grid, fixed SN exclusion and design-image source mask are shared by both branches. Do not add masked variants of the radius or split-reversal sensitivities.

For each exposure and candidate define a fixed validity indicator on native pixels:

`v_p = (DQ_p==0) AND finite(SCI_p) AND finite(ERR_p) AND (ERR_p>0)`.

PAM and geometry must also be finite and positive on the footprint; a failed geometry/calibration check rejects the position rather than imputing a value. DQ/ERR/finite status can depend on the same observed data: this is a pipeline-accepted-pixel estimand, not an unconditional noise test. Do not cut on signed flux, measured SNR or aperture residual.

With original fractional overlaps `a_p,b_p` and the verified official PAM `P_p`, compute

`A_full=Σa_p P_p`, `A_good=Σv_p a_p P_p`,

`B_full=Σb_p`, `B_good=Σv_p b_p`.

Require **A_good/A_full ≥0.90** and **B_good/B_full ≥0.75**, inclusive, separately in every one of the eight fixed exposures used by the common comparison. Aperture coverage is geometrical PAM-weighted area; annulus coverage uses fractional native-pixel area, matching its original flat-fielded-SCI mean. Save all eight coverage fractions and rejection reasons. Nonpositive denominators fail. There is no empirical renormalization or threshold tuning.

The retained-pixel background and complete signed weight vector are

`Bhat_e=Σv_p b_p SCI_p / B_good`,

`w_ep=κ_e[v_p a_p P_p−(A_good/B_good)v_p b_p]`,

`f_e=Σw_ep SCI_ep`, `V_e,diag=Σw_ep² ERR_ep²`.

Evaluate sums only on valid entries, so invalid NaN/Inf values do not leak through a numerical `0*NaN`. The coefficient of the annular subtraction uses **A_good**, not A_full; using the original coefficient would fail constant-SCI cancellation. Aperture and annular weights must be assembled before squaring. Carry the same κ/PHOTFLAM and PAM convention as the primary. Do **not** divide by aperture coverage, impute missing pixels, apply a point-source aperture-loss correction or copy the full-coverage variance.

This operator integrates the retained pixels, and removes a constant native-SCI sky exactly. It is **not an estimate of total point-source flux**. Different masks in two dithered exposures produce different responses to unresolved sources, host gradients, PSF wings and residual flat structure. Therefore repeat-contrast excess can contain real static-scene mismatch as well as noise covariance. The independent design mask reduces that risk without proving it absent. Report the mask/coverage distribution and spatial/background stratification alongside the residual statistics.

## Support, dependence and execution gates

Keep the original geometry footprint and source-mask construction independent of the branch. Pixel-quality holes must not be silently dilated into the common astrophysical source mask merely to reproduce the strict full-coverage rejection; they belong in each branch's explicitly recorded validity/coverage rule. Any necessary deterministic treatment of invalid design pixels during source detection must be specified in the frozen executor before support counts, shared by both branches and independent of test aperture outcomes.

The secondary retains the same100-position,12-of16-spatial-tile target;≥30 positions permits the same explicitly limited descriptive treatment, and<30 stops its score. If the strict primary stops, record that stop before showing the secondary. Neither branch's outcome chooses the other's support or parameters.

Use the same signed pair estimands, centering ledger and spatial dependence treatment as the original design. Common-template comparisons retain their full `H D Hᵀ` covariance. Within-visit equal mean response is now an additional approximation when masks differ; do not interpret a cross-product departure uniquely as shared stochastic template noise. No extra source/sky Poisson variance is added to FLT ERR.

Before any observed aperture sums, the implementation must verify on synthetic values and the actual **geometry-only** coefficient vectors: constant-SCI cancellation, equivalence of direct background subtraction and dot product, masked weights and variance, and κ/PAM unit scaling. The original doubled-resolution0.1% quadrature gate applies. Numerical quadrature details and any deterministic resolution refinements must be frozen in the executor and preserved; no scientific threshold changes are implied. Root releases geometry/support first and measurement sums separately. This amendment does not authorize execution.

## Geometry prerequisite now available

Sol's official-map audit in `pixel_acquisition/pam_validation/` acquired `ir_wfc3_map.fits` and the exact header IDCTAB `w3m18525i_idc.fits`. All eight FLTs are full-frame1014² with zero LTV offsets and identity LTM. Official PAM orientation/support therefore maps directly onto the native pixels. Its normalized area agrees with the native TAN-SIP WCS Jacobian to a maximum relative discrepancy7.82×10⁻⁴ on the predeclared11×11 grid; halving the derivative step changes the Jacobian by at most4.56×10⁻⁹. Use the official PAM, retaining that audit and its source hashes. This certifies geometry for the proposed operator, not photometry/noise outcomes.
