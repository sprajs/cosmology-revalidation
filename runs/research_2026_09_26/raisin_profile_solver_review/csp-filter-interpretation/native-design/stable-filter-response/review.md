# Stable native CSP filter-processing response

The frozen five-job comparison passes all declared numerical gates for all 42 selected CSP objects. Changing the 188 securely identified WIRC-J observation labels from released `J` (RC1) to released `j` (RC2) shifts the unweighted mean fitted distance by **+0.000483432007 mag** across all 42, and **+0.000634504509 mag** across the 32 objects selected as affected from the input lineage. The individual responses range from −0.003213875645 mag (2007jg) to +0.005489917673 mag (2004gu). All ten metadata-defined unchanged controls are exactly unchanged. Holding the high-redshift sample fixed, the corresponding unweighted high-minus-low distance contrast changes by **−0.000483432007 mag**.

This is a conditional response of the current released inputs and fixed historical empirical model to a source-motivated filter-operator change. It is not a physical WIRC calibration correction, a regenerated bias correction, an uncertainty bound on the full calibration problem, or a corrected cosmology. RC2 is the paper-motivated approximation to WIRC-J; the actual WIRC throughput and the empirical model's inherited training conventions remain separate questions.

## Frozen design and execution

The protocol was frozen before the full-cohort changed-filter outputs. It retained all 42 objects, all quoted fluxes/errors/header peaks, all cuts, the same SNooPy_B18 surfaces and KCOR file, the historical v11_04k native objective and covariance update, fixed stretch 1 and empirical AV 0. Four ambiguous 2006kf J/Jdw rows remain unchanged. The 32-object affected roster is defined by the prepared raw-row ledger, never by the fitted response.

The five jobs were changed-filter 9 iterations, changed-filter 12 iterations, changed-filter 12 with actual post-initialization first-iteration distance starts shifted by −0.2 and +0.2 mag, and a complete unchanged nominal 12-iteration copy. The primary estimator is the default 12-iteration result; no fit was chosen by chi-square or effect size. The 3-iteration response is a separately labelled historical stopping-recipe sensitivity.

All jobs returned zero process status and all 210 object/job FITRES rows have ERRFLAG 0 and the prescribed fixed parameters. The changed 9-iteration states exactly equal the first nine states of changed 12. The unchanged 12-iteration copy has exact native callback states, science rows and LCPLOT bytes. The actual ±0.2 first-entry shifts meet the frozen 1e−12 coordinate tolerance.

The native jobs used 15.349783268 seconds, bringing the shared cumulative native budget to 43.947051949 of 600 seconds; no invocation failed in this stage. The additional 120-second and per-job 60-second caps were respected. Native originals and previous failed gates remain untouched.

## Numerical and row-identity gates

All 2,394 saved callbacks satisfy the direct residual quadratic check. The maximum absolute difference between reconstructed data quadratic plus printed prior/sigma terms and the native objective is 1.54614e−11. The model-dependent covariance-map relative Frobenius discrepancy is at most 8.00424e−8, below the frozen 5e−6 cap. Every accepted native mean and covariance passes the positive-mean/positive-definiteness checks.

For all changed objects, D9 equals D12 exactly at the exported double precision. Across the 12-iteration jobs, the last two distance increments and covariance changes are exactly zero. The maximum final distance difference across the three tested starts is 2.44134e−8 mag; the largest fixed-C positive-amplitude optimum gap is 1.84195e−8 mag. These are far below the unchanged 0.001-mag gate. Tested-start agreement establishes the bounded numerical check, not arbitrary global uniqueness or a physical likelihood theorem.

Every changed callback contains the same physical measurements as the corresponding nominal callback after the exact declared J→j row map; duplicate multiplicities are preserved. There are 1,812 accepted final measurements, including 170 relabelled observations in 31 objects. The metadata-affected object 2008gp has only out-of-final-window relabelled rows and remains in the 32-object estimand. The ten controls match their corresponding nominal start at every callback, FITRES science row and LCPLOT row.

All iteration-2-and-later accepted phases lie within the empirical grid −20 to +70 days; shape is within the source grid; full nonzero observed-filter throughput mapped to rest wavelength satisfies each native wavelength screen. The same nine objects have unsupported first-iteration late phases: 2007S, 2008ar, 2006et, 2007as, 2008bc, 2008hv, 2009aa, 2006kf and 2009al. Their source-defined tail extrapolation is retained. Convergence from the tested starts does not validate that initial empirical extrapolation.

## Primary and secondary response

| Quantity | Stable 12-iteration primary | Original 3-iteration secondary |
|---|---:|---:|
| Mean changed minus nominal, 42 objects | +0.483432 mmag | +0.482095 mmag |
| Mean changed minus nominal, metadata-affected 32 | +0.634505 mmag | +0.632749 mmag |
| High-minus-low contrast change, high-z held fixed | −0.483432 mmag | −0.482095 mmag |
| Maximum absolute unchanged-control response | 0 | 0 |

The secondary minus primary mean is −0.001337 mmag across all 42. The largest per-object difference between the two processing responses is 0.072372 mmag (2008ia). The original nominal 3-iteration convergence failures are not erased by this small difference between *paired responses*. The separate nominal continuation changed 2008ia's distance by 1.231436 mmag. Both arms must be numerically checked before their cancellation is interpreted.

Native chi-square values remain descriptive fit diagnostics. The original two-object baseline already had data quadratics 105.490 for 14 degrees of freedom (2004ef) and 64.888 for 38 (unchanged 2005hc). Exact reproduction and numerical convergence do not make this empirical model a validated description of these measurements. Different filter arms also change model-dependent covariance, so comparing their raw fit quadratics is not a likelihood-ratio test or evidence for a physical calibration winner. No alternative Gaussian maximum-likelihood estimator was fitted; the earlier distinction between the native fixed-point estimating equation and a normalized parameter-dependent Gaussian objective remains.

The current-input 2007A comparison still fails the earlier archived chi-square reproduction gate because the released and archived peak headers differ. No header was repaired to force historical agreement. No object or high-chi-square row was removed. No bias-correction, selection, training or systematic-covariance product was regenerated, and no cosmology was fitted.

## Independent review and provenance

Root independently checked the raw paired photometry lines and native logs without importing this agent's mapping or arithmetic: 1,890 changed physical callback comparisons and 450 exact control comparisons pass, with the same response means/range and objective closure. See `runs/research_2026_09_26/csp_native_filter_response/root-state-review/stable-response/result.json`.

Primary machine artifacts in this directory are `result.json`, `summary.json`, `responses.csv`, `state-histories.json`, `control-gates.json`, `accepted-row-accounting.json`, `execution-ledger.json`, and the five complete native log/FITRES/LCPLOT directories. Earlier failed 3-iteration gates and successful nominal continuation remain at their original paths.

Frozen identifiers:

- Protocol: `0ebd121b8b842aebd5d9d7c7969b14d1d85fd6e376ca81db2141b72c1994eb91`.
- Active runner v2: `c81bd4f5c0527fd65e84296a797a8a90414a345ef7394910cbc2261aa1f2b7fd`; its only amendment caches existing shifted nominal log parsing.
- Runner cache amendment: `7d3c43cef3451080450f889cb49d0faf419cb0b793dc420a0096cd7660834fa4`.
- Pre-output checker: `d944da49d2dd7971eebbd650934dcdfbd24572f5a65ad6b383d98f4fbaf6d4ef`.
- Checker freeze: `877bf4a6d01a8b4b631506f2f8a6d8a841e5c9bc668deca56e3fc13aba365a98`.
- Root execution release: `212ff20cc8b2135f0c51a09e5774f59830bbd06e75e301a2df80f564fcfb3d50`.
- Unchanged instrumented historical binary: `716d5c235468329e878bd56e678cfa55b2a316f333c36ec411e9a320c1770938`.

A final manifest records all owned files, external review linkage and the data symlink target. No source or data was overwritten to obtain these results.
