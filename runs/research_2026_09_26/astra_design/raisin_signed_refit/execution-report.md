# Modern controlled RAISIN signed-epoch refits: numerical gate

The source-defined paired inputs and native reading gates pass. The numerical stability gate fails. No fitted-distance correction, uncertainty interval, or sign-selection bias is inferred from these runs.

The frozen design is `docs/research-2026-09-26/raisin-signed-refit-design.md`; its original bytes are retained in `input-manifest.json`. This execution note corrects interpretations found by subsequent source inspection without rewriting that preregistration. The authoritative numerical audit is `qualified-stability/`, produced by `correct_objective_audit.py`. Earlier `result.json`, `stability-result.json` and CSV files are preserved and hashed, but their first-fit summaries and fields called “best_start_by_objective”, “comparable_objective”, and “conditionally_preferred” must not be treated as scientific rankings.

## Inputs and native closure

The ten DES16 objects are selected by exact author-product/raw-ancestor provenance. They all occur in each released nominal 79-object optical, optical+NIR, and NIR cosmology table; this membership was checked without dropping or reselecting objects (`released-cosmology-membership.csv`). The four arms are released R, author-positive A, author-signed B, and released-plus-author-negative H. A and B use the same raw error convention, so their input difference is exactly measured nonpositive epochs. H is an explicitly labelled convention bridge. All released R bytes, unchanged nonoptical rows, duplicate multiplicities and raw error/flag provenance are preserved.

The source audit subsequently found that all 3,133 unique positive author-to-product matches obey `product_error = round(raw_error * 1.086 * 0.4 * ln(10), 3)` exactly. This numerical identity is consistent with a magnitude-error round trip; the executed conversion script is not established. The frozen arms were not silently rescaled after this discovery. The error factor is 1.0002429644, and a future mapped-H branch would require an explicit amendment.

The initial modern SNANA run required only a DOCUMENTATION wrapper around the original byte-identical vpec table. Its failed initialization is retained. The first complete object took 5.58 seconds including initialization, with about 0.46 seconds reported for fitting. The 80 primary fits completed in about 77 seconds. Released input and byte-identical copy give exactly equal fitted values, reported χ² and LCPLOT bytes. Parent's independent parser verified all frozen row values, header fields, hashes and multiplicities (`raisin_refit_input_review/input-check.json`). Native accepted data were matched to the original rows with explicit print-precision tolerances and ambiguous-date multiplicities retained.

These are **modern explicit MW-law 99** comparisons. The archived author FITRES declares SNANA v11_04k; parent independently acquired that source and found its original default law 94 with the release's law-99 line commented. Historical version/law reproduction is a separate experiment. It must not be silently conflated with these results.

## Numerical readiness fails despite nominal success flags

All 240 fits—80 primary cases and two deterministic shape starts for every case—report ERRFLAG=0 and positive reported parameter covariance. Nevertheless:

- 28 of 80 cases have a DLMAG range above 0.001 mag across starts; 20 exceed 0.01 mag.
- The largest ranges are 0.224064 mag in DLMAG, 0.132682 in AV, 0.193845 in stretch, and 1.52734 days in peak date.
- The accepted-row mask is equal across starts for 76 cases. The four changing-mask cases are free-peak DES16C3cmy in each data arm. No objective ranking is made across those masks.
- Time-sorted and reverse-time permutations of the first object's B rows preserve the accepted rows but change DLMAG by +0.000118 and +0.000076 mag, respectively. The input information is unchanged; small optimizer path dependence remains.

The 20 declared local fixed-shape profiles all complete with ERRFLAG=0 and the same accepted mask as their original free-B fit. Across the local offsets, DLMAG spans are 0.012192 mag for DES16C1cim, 0.008862 for DES16S1agd and 0.006187 for DES16E1dcx. These are numerical mapping diagnostics, not likelihood intervals. The original resource prose said 21 points; one E1dcx endpoint lay outside the tabulated grid and was omitted in the frozen job list before execution, giving 20.

Independent source-extracted grid integration confirms piecewise stretch-slope changes, including changes up to 1.95 mag per stretch in individual rest-band/phase cells, while limiting value jumps are below 5.1e-7 mag. Such knots undermine local Hessian interpretations but do not alone establish the cause of the 0.224-mag multistart spread. See `raisin_profile_solver_review/grid-result.json`.

## What the native objective actually means

Pinned modern `snlc_fit.F90` states at lines 3941–3945 that `FCN_FITCHI2` excludes `CHI2INI`, although minimization includes it. `qualified-stability/fit-ledger.csv` records the reported data χ² and the separately printed prior term. The latter is rounded native log output (range 0–0.0012 in these runs); their sum is only a display-precision minimand diagnostic.

The 5-day peak term is **recentered on the previous iteration's fitted peak**. `FITINI_PARVAL` updates INIVAL for iterations after the first (around lines 7919–7950), and `FCNCHI2_PRIOR` uses peak minus INIVAL (around lines 4842–4848). The released header sets the initial peak; it is not an immutable Gaussian prior center throughout the native algorithm. The frozen design's phrase “inherited fitted timing constraint” needs this qualification. Adding `(final_peak − header_peak)^2 / 25` would reconstruct the wrong minimand.

Likewise, OPT_COVAR_FLUX=0 does not freeze all weights. `FITINI_COV` builds diagonal model errors and the MW rank-one covariance from iteration-reference predictions, then adds quoted measurement errors and inverts the matrix (around lines 10660–10730). Iterations 2/3 use that stored matrix inside FCNSNLC, while iteration 1 has its separate weighting behavior. Thus matching accepted rows is necessary but insufficient to establish a shared likelihood across starts or fixed-shape runs. OPT_CHI2_SIGMA is zero here; these fits do not provide a normalized parameter-dependent Gaussian likelihood by default.

Consequently, no start is selected as a global likelihood maximum, even if its printed data χ² is smaller. The parameter spread is directly observed and remains a valid failure of numerical readiness. The pairwise first-start distance changes remain archived as diagnostics, not promoted as corrected distances.

## Next controlled gate

Parent authorized a single-object proof with DES16C1cim, author-signed B, fixed released header peak. The native exporter first reproduces the original modern fixed-peak result, then saves exact accepted rows, mean, inverse covariance, covariance, coordinates, current prior center and covariance-reference state. Fixed-coordinate mean probes test DLMAG multiplicativity and row permutation. The failed generic SNCID_LIST_FILE engineering attempt is retained: its reader stores SALT x0/x1/c names only, so it does not supply SNooPy DLMAG/STRETCH/AV. The replacement uses the native SNooPy INIVAL keys with explicit seed precision rather than pretending these models share parameter conventions.

If these mean/export gates close, define a **new conditional objective** with one declared reference covariance and mask. For fixed shape, AV and peak, prove the native relation `F(D)=a F(Dref)`, `a=10^[-0.4(D−Dref)]`; then the positive interior amplitude profile is `a=(hᵀC⁻¹y)/(hᵀC⁻¹h)`. It must not be used if the optimum amplitude is nonpositive or lies beyond the verified distance/prior support. A flat-distance marginal likelihood is not obtained by integrating with a flat-amplitude measure.

The global fallback must scan every tabulated shape cell, retain endpoints and distinct AV branches, and refine until objective and distance-profile changes meet declared tolerances. Fixed peak comes first. A subsequent free-peak sensitivity must freeze epochs for the entire admitted peak range and use one explicit stationary peak prior or none, with that change labelled. Retain competing distance modes and profile envelopes rather than reporting local Hessian errors as global uncertainty. A normalized parameter-dependent covariance likelihood would be a separate declared model, not a hidden repair of the native estimator. Neither approach establishes physical noise covariance or a population bias correction.

The matched-noise sign-selection simulation is therefore still stopped at its fitter-readiness gate. Its toy mechanism remains valid; the proposed native simulations cannot yet deliver reliable bias estimates through this optimizer.

## One-object native operator gate now passes

`native-profile-gate-v2` reproduces every checked modern reference FITRES value exactly for fixed-peak DES16C1cim B. Its 74 accepted epochs include nine negative fluxes and 65 positive fluxes. All exported flux/error values equal their raw author's exact native float32 representation. The initial assertion of equality to unrounded decimal input failed and is retained in `native-profile-raw-check-failure.json`; the largest expected storage differences are 3.60e-6 flux units and 1.13e-6 error units. No observations or uncertainties were changed to obtain closure.

Native DLMAG shifts near ±0.15 mag reproduce multiplicative flux scaling to 2.06e-9 relative, using the actual NML float32 offsets. Fixed-coordinate means are exactly invariant to sorted/reversed input rows. Native C changes by 0.284% and 0.374% in Frobenius norm under these two distance changes; these variant matrices are diagnostics and are not substituted into the frozen-C objective. Independent Cholesky and direct scalar amplitude profiles agree in objective to 3.76e-11, with exact native total objective reconstructed to 7.11e-14. See `native-profile-independent-check.json` and the independent collaborator review.

An isolated opt-in stream hook calls the original full observer-frame USRFUN after the native reference has loaded. It does not implement new K-corrections or interpolate an approximate mean surface. Its source patch, build log and binary hashes are in `fixed-c-profile/`; the existing binaries remain unchanged. All seven separate-process native vectors (nominal, two distances, two shapes and two AV values) match the stream exactly. One hundred design vectors repeated in randomized order, then a nominal replay, also match exactly; returned reference rows and C/W remain unchanged. One hundred vectors take 0.0384 seconds after initialization. This removes the five-second-per-process overhead that would otherwise prevent a useful global search.

The next A/B fixed-C experiment is frozen in `fixed-c-profile/global-profile/protocol.json`. Its additional covariance anchor is fixed before outcomes from the original A/start1 coordinates, exported over all B epochs. The common B-anchor A covariance is the marginal underlying-Gaussian submatrix; it is not the distribution of errors conditioned on positive selection. Even a numerically stable A/B distance difference remains a conditional processing sensitivity rather than a bias-corrected likelihood.
