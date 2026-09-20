# Independent directional cosmography investigation

Initial findings recorded 2026-09-20 before reading other investigators' scientific conclusions. Plan: `docs/experiments/directional-plan.md`, parent COS-04. The original author snapshot remains unchanged. This analysis starts from public Pantheon+ corrected-magnitude or SALT2-summary tables; it does not refit photometry.

## Claims, method, and literature baseline

Sah, Rameez and Sarkar, [arXiv:2606.09650v1](https://arxiv.org/abs/2606.09650), [MNRAS DOI 10.1093/mnras/stag844](https://doi.org/10.1093/mnras/stag844), argue that a progenitor-age magnitude correction makes the monopole of their directional cosmographic deceleration coefficient positive, while leaving its local dipole almost unchanged. Their main Table 1 is **C2**, which fits SALT2 magnitude, stretch and colour jointly with redshift- and survey-independent Gaussian population means and widths. C1 uses the supplied already-standardized magnitudes and total covariance. These are two likelihoods of extensively shared observations, not independent discoveries.

The cosmographic model is

\[
D_L={c\over H_0}z\left[1+{1-q\over2}z-{1-q-3q^2+J\over6}z^2\right]{1+z_{HEL}\over1+z},
\quad q=q_m+q_d\cos\theta\,e^{-z/S},
\]

where the fitted cubic parameter is the jerk-curvature combination J=j0-Omega_k in FLRW conventions. The fixed direction is equatorial RA=168 deg, DEC=-7 deg; the LG version uses RA=162.95389715 deg, DEC=-25.96734154 deg. The main cut is .00937<zHEL<.8. The redshift-dependent substitution into the isotropic cubic series is a phenomenological model, rather than a unique covariant reconstruction of a directional spacetime expansion history. In particular the q-squared term produces higher angular powers, and a coefficient inferred over a finite redshift range is not automatically a model-independent measurement of a derivative at z=0.

The authors already report a zHEL<.5 check with no qualitative change, and their predecessor says the LambdaCDM cubic distance approximation is within about 2% through .8. Those are counterevidence to treating any nonzero Taylor error as a refutation. Their use of a fixed 0.030 mag/Gyr correction amplitude and an externally supplied redshift-age evolution is an astrophysical assumption, not measured by this directional likelihood. The safe spline reconstruction below reproduces Sah's supplied CSV curve, **not an independently obtained exact Son correction array**.

Live arXiv metadata on 20 September still lists only v1; repository HEAD is unchanged at `8ba821e2fce205b8a760bfbef9d5637b1df19c09` and GitHub releases is empty. The full main paper, including C1 Appendix A, and predecessor [2411.10838v2](https://arxiv.org/abs/2411.10838) methods, tables and covariance appendix were read. A source-acquisition colleague resolved the main version-of-record PDF from Oxford ORA after our initial HTTP403: `sources/updates/2026-09-20-ztf/originals/searches/sah2026-stag844-version-of-record.pdf`, SHA256 `1db9c966d8707c3aac9e045c01275218287f832587dd472cd609494654580490`. Its complete six pages were read, confirming the same equations, main Table 1, C1 Table A1 and z<.5 sensitivity statement. Source snapshots and retrieval outcomes are under `runs/directional/sources/`.

## Input and implementation audit

The C1 `Z_mbcorr.csv` columns exactly equal the corresponding original Pantheon+ table columns in original order: 1701 light-curve rows, 1543 distinct CID values (a string-count diagnostic, not an independently resolved physical-SN count). Repeated measurements are preserved. `statsys_mbcorr.npy` exactly equals the released Pantheon+ STAT+SYS covariance in value, so it is a covariance, not an inverse, and no additional tabulated statistical errors should be added. It is positive definite; its tiny maximum asymmetry is 3e-8 from source rounding. Our symmetric realization uses the lower triangle, exactly matching the authors' lower-triangle Cholesky likelihood. Row sorting applies the identical permutation to both covariance axes. The author's C1 model adds a separately fitted sigma_M squared to this already-total covariance; we reproduce that choice, without interpreting it as a fresh measurement-error term. The paper's method prose says zHEL<=.8, but the public code and Figure 2 cut use strict <.8: this excludes the row exactly at .8 and gives 1563 C1 rows after the lower cut, versus Ray's inclusive 1564.

The README-linked [Zenodo15629813](https://zenodo.org/records/15629813) is publicly resolvable through its API even when its browser page is blocked. The acquired `cov_final.npy` is 203939336 bytes, MD5 `09996e2b37009aa7d0c7de13b79c90c7` (matches published checksum), SHA256 `ac4d002f3c6bf7dd2defe79e286ca6294fbfb9761c8fcc9d0a9796f81100e0eb`. It is a symmetric positive-definite 5049x5049 matrix, paired with 1683 distinct original row indices. The indexed rows are sorted by zHEL; 1653 survive zHEL<.8 and 1547 survive the full main-paper cut. The 18 omitted rows, including one above .8, are preserved in `c2-input-audit.json`: 15 distinct CID strings overall, 14 below .8. This differs slightly from the predecessor prose's 17 omitted rows/15 unique SNe below .8, while the released row counts and covariance dimension agree. The original C2 summary columns exactly match the full Pantheon table. The supplied matrix includes additional systematics, so its x1/c diagonals exceed the SALT2-only error squares; comparing them as if they should be identical would be incorrect. This validates the released product and its documented ordering; it does not reconstruct every Lane covariance component from raw calibration inputs.

The likelihood normalization is N ln(2pi)+ln|C|+r^T C^-1 r for C1 and 3N ln(2pi)+ln|C|+r^T C^-1 r for C2. C2 includes the full repeated 3x3 population covariance block with alpha/beta-induced off-diagonals. Its mean profiling uses an invertible change from (M,X,C) to (M-alpha X+beta C,X,C), at fixed alpha and beta; determinant terms remain unchanged.

Public script variants contain hard-coded private `/Storage/animesh/` paths and pickle loads; one current C1 variant (`analysis_c1.py`, lower-case a) contains pasted runtime-warning output and fails syntax parsing, with an additional undeclared ray call. This is a public reproducibility defect. It does **not** demonstrate an error in the code actually used for the published fits: our independent implementation recovers the predecessor C1 table. We never execute those scripts or deserialize pickle data. The supplied `Progenitor_age.ipynb` constructs `CubicSpline` directly from sorted `median_deltaage.csv`; we repeat precisely that safe construction.

## Reproduction results

C1 fixed-direction, reversed-selection-bias, no-age fits reproduce the rounded predecessor Table 5 entries. We use H0=70 as a distance scale absorbed by a free magnitude intercept. All dipole starts converge to the same reported minima to numerical precision for the data fits.

| Frame | q_m | q_d | S | J | -2 ln L | Improvement in -2 ln L over q_d=0 |
|---|---:|---:|---:|---:|---:|---:|
| HEL | -0.36933 | -6.276 | .024499 | .1224 | -1507.3733 | 30.5104 |
| CMB | -0.39242 | 20.543 | .011239 | .2635 | -1518.2713 | 28.2515 |
| LG | -0.41151 | -39.569 | .012357 | .2883 | -1501.1963 | 119.7230 |
| HD | -0.48889 | 11.399 | .010440 | .7102 | -1520.7310 | 6.3430 |

The same safe age curve moves q_m to -0.01388, -0.03425, -0.05479 and -0.12538 respectively. Its shift is about +.355 to +.364, while q_d barely changes. This reproduces the substance of the main paper's **own C1 Appendix A**, which explicitly reports near-zero monopoles (-.01,-.03,-.05,-.12), rather than its C2 positive values. Treating the C1/C2 distinction as a numerical disagreement would be a category error. Some last published digits are not matched: for example C1 HEL age-corrected S=.024875 versus .024 printed, HD q=-.12538 versus -.12, and CMB likelihood improvement 27.822 versus 27.9. We therefore classify the full C1 table as a **faithful released-input, close numerical reproduction**, with residual rounding/optimizer discrepancies rather than claiming every printed digit agrees. Keeping Pantheon selection corrections instead of reversing them gives yet another explicitly labelled C1 sensitivity, saved in `c1-fits.json`.

We report likelihood improvements directly. Scale S is unidentified when q_d=0, and several fits hit the imposed lower scale limit; therefore the usual regular Wilks chi-square calibration is not established. We do not relabel those improvements as a verified number of Gaussian sigma. The paper labels its tabulated statistic Delta ln L, but the public code and reproduced values correspond to differences in the -2 ln L objective.

C2 **full HEL reproduction is complete** after safe input recovery. The analytic likelihood gradient, including covariance derivatives and profiled means, agrees with central finite differences at maximum scaled error 7.7e-9. That validation was recorded before its fit results were inspected. The complete free-population likelihood gives:

| C2 HEL | q_m | q_d | S | J | -2 ln L dipole | -2 ln L isotropic | Improvement |
|---|---:|---:|---:|---:|---:|---:|---:|
| No age correction | .0095173 | -31.7734 | .00938 | -.646482 | -184.912667 | -148.022189 | 36.890478 |
| Supplied age correction | .3539906 | -32.1337 | .00938 | -.377093 | -184.735624 | -147.794707 | 36.940917 |

These reproduce all four headline HEL entries of published Table 1 at its printed precision: (.01,-31.8,.0094,36.9) and (.35,-32.1,.0094,36.9). The monopole shift is +.344473. Both dispersed optimizer starts agree to <=1.1e-10 in the objective and <=3.1e-7 in q_m. All nuisance parameters are included in `c2-fits.json`; alpha=.1583, beta=3.1752, sigma_M=.15672, sigma_X=.96299 and sigma_C=.05522 for the uncorrected dipole. The scale is at the imposed lower bound in both fits, so the unconstrained S-gradient need not vanish. This is an **exact-equation, released-input numerical reproduction of the HEL point estimates and likelihood ratio**, not an independent raw-data reconstruction, validated posterior probability of deceleration, or verification of the quoted 5.7-sigma calibration. Other C2 frames, tomographic bins and contour surfaces were not refitted in this bounded branch.

An independent Astra audit evaluated the four C2 final likelihoods using a different whole-observable transformation `(m+alpha*x-beta*c,x,c)` with the correspondingly transformed **full** measurement covariance and diagonal intrinsic variances. The maximum -2 ln L difference is 9.1e-12, and profiled means agree to 4.3e-14 (`runs/audit/c2-crosscheck.json`). Its separate direct-Cholesky C1 check agrees with the eigensystem implementation to 1.1e-11, and an independent ODE luminosity distance plus Nelder-Mead fit recovers the LambdaCDM synthetic q bias within 9.3e-9 (`runs/audit/directional-crosscheck.json`). This checks implementation consistency; investigator agreement is not independent observational evidence.

## Controlled Taylor truncation tests

At z=.8 the cubic distance modulus error at the **true** cosmographic coefficients is -0.02947 mag for flat Omega_m=.3 LambdaCDM, -0.04964 mag for Einstein-de Sitter, and -0.04846 mag for flat constant-q=0 coasting. These correspond to -1.35%, -2.26% and -2.21% distance errors. LambdaCDM quadrature agrees with independent Astropy distances to <1e-12 relative. Einstein-de Sitter and coasting use exact closed forms. Thus a universal claim of <2% accuracy across arbitrary models is too strong, but the LambdaCDM-specific statement is satisfied.

Crucially, fitting a free intercept, q and J absorbs much of this error. Exact noiseless models injected at actual Pantheon+ zHD/zHEL and weighted with the full C1 covariance (.01<zHEL<.8; 1552 rows) give:

| Injected exact model | True q0 | Cubic-fit q0 | Bias | Forecast sigma(q0) | Bias/sigma |
|---|---:|---:|---:|---:|---:|
| LambdaCDM, Omega_m=.3 | -.5500 | -.52363 | +.02637 | .08397 | .314 |
| Einstein-de Sitter | .5000 | .52409 | +.02409 | .07195 | .335 |
| Flat constant-q coasting | .0000 | .02748 | +.02748 | .07802 | .352 |

For zHEL<.4 the same biases fall to +.00620, +.00968 and +.01078 respectively, each <.07 forecast sigma. They approach zero under tighter cuts. At .8, 100 paired correlated-noise realizations per model give mean q values -.515, .531, .035, with ensemble SDs .077,.065,.071. The cubic self-injection control recovers q=-.55, J=1 and q_d=0 with chi2 <3e-23. An exact isotropic LambdaCDM injection with dipole freedom yields q_d=.0084 and improvement only .0023 in chi2: the series error plus this actual angular selection does not create the observed large dipole in this controlled model.

These are conditional synthetic checks, not a demonstration that LambdaCDM or the adopted covariance/population model is true. In particular they use C1 weighting and do not yet quantify a joint C2 nuisance-population misspecification. Nevertheless they are strong counterevidence to the hypothesis that cubic truncation by itself explains the approximately +.35 shift caused by the supplied age correction. Truncation slightly biases q towards deceleration in these examples but does not turn the accelerating LambdaCDM injection into a decelerating result.

The initial low-z noisy ensemble hit the arbitrary J bound in 62-66% of z<.1 draws and 2-4% of z<.2 draws. We preserve it, flag it as an identifiability/bounds diagnostic, and rerun widened bounds before using low-z ensemble means. The primary low-z statement above concerns noiseless recovery, not those bounded means.

The completed widened-bound rerun has no q or J boundary hits. Its q standard deviations are about 1.7 at z<.1 and .47-.50 at z<.2, rather than the much smaller values artificially imposed by the original J bound. Low-redshift truncation improves while the derivative estimate becomes very weak. Six total solver calls across original/wide runs return non-success termination despite finite nearby optima (all calls and flags retained). Independent Nelder-Mead reruns of those points all terminate successfully, improve the objective by at most 8e-12, and shift q by at most 5.2e-7. Thus those flags do not affect the numerical findings. The low-z ensemble is a diagnostic of lost precision, not a calibrated posterior.

## Subsequent coordinate dispute and revised Ray analysis

The [Sah response 2608.02484v1](https://arxiv.org/abs/2608.02484) (3 August) objects specifically to Ray et al.'s use of Galactic CMB coordinates (264,48 deg) as though they were equatorial RA/DEC. Its claim does not identify this error in Sah's own code or in our reproduction. Independently applying the wrong coordinates to .00937<zHEL<=.8 (1564 rows) gives exactly the reported 724/840 split. Astropy transforms the direction to ICRS (167.78661,-7.14539 deg); both that and the rounded (167.8,-7.1) give 538/1026. Sah's response table labels (167.8,-7.1) but prints 539/1025. An exploratory precision check finds that using the public Sah main-code direction (168,-7), or (167.94,-6.94), gives 539/1025. A slightly different reference direction is therefore a plausible explanation, not a verified identification of the response's actual executed pipeline. Its printed direction/count remain inconsistent by one net row; this does not challenge the demonstrated large Galactic/equatorial coordinate error.

Crucially, the current [Ray2607.20570v2](https://arxiv.org/abs/2607.20570v2), submitted 3 August 2026 17:55:09 UTC, **withdraws the earlier numerical deceleration conclusion**. It acknowledges the coordinate mistake and says that it cannot reproduce its original full-sample baseline either, which cannot be explained by a sky-coordinate change. It reports q=-.490 without and -.267 with the supplied age correction, and the corrected hemisphere counts 538/1026, agreeing with our count. Its reported hemispheric q values are -.527/-.464. Therefore the v1 deceleration numbers must not be presented as the current result or independent support for Sah2026.

We read the revised method and discussion in full. Ray v2 fixes J=1 and uses an equal .15 mag diagonal uncertainty while fitting q, alpha, beta and an intercept (its stated H0 and M are exactly degenerate). This differs substantially from C2's full covariance and fitted J/population model. Our direct implementation of that stated simplified estimator, safely applying the same Sah spline, gives q=-.46311 then -.24911 for zHD, -.43304 then -.21142 for zCMB, and -.43603 then -.21484 for zHEL. Thus the qualitative accelerating result and approximate .21-.22 age shift reproduce, but **the exact printed q values do not**. The printed method does not state an unambiguous cosmological-redshift field or enough additional pipeline details to resolve this remaining discrepancy. This is an approximate reproduction, not a claimed numerical correction to Ray v2.

The hemispheric means are a coarse, differently weighted statistic from the local exponential dipole; their agreement alone does not rule out a scale-dependent low-z dipole. Neither Ray's simplified error treatment nor multiple decompositions of the same observations make this an independent measurement with an established covariance-calibrated significance.

## Strongest counterevidence and limits

The strongest evidence against a broad numerical dismissal is the reproduced C1 table, correctly recovered public C2 covariance, exact self-injection recovery, and modest q bias under exact-model injections. The strongest limits on the paper's cosmological interpretation are its assumed luminosity-age correction, constant C2 population model, omission/reversal of selection terms, peculiar-velocity/frame dependence, nonregular dipole significance calibration and the phenomenological nature of its directional cubic coefficient. A dipole in an observer-frame Hubble diagram does not by itself exclude a cosmological constant: peculiar motions and an accelerating background can coexist. This branch has not fitted a physical velocity-field model or raw-survey selection likelihood, so it cannot adjudicate that alternative observationally.

The C1 age shift is a conditional response to an imposed redshift correction. It is not fresh evidence for the physical progenitor-age law and cannot be combined as independent evidence with other Pantheon+ reanalyses. A rigorous directional significance measurement would need a null simulation ensemble incorporating the fitted scales, sky mask, duplicate observations, covariance, selected populations and velocity treatment. A discriminating astrophysical interpretation additionally needs a validated luminosity-evolution model or measurements that break its distance-redshift degeneracy.

Reproduction commands:

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/audit.py data
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/input_audit.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/audit.py synthetic
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/lowz_wide.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/c2.py check
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/c2.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/directional/ray_check.py
```

All numerical outputs and input/code hashes are retained under `runs/directional/`. Final per-experiment manifests include the actual executed code snapshots, input and output hashes, configurations and seeds. Initial C1/synthetic outputs are retained in `initial-results/`; after adding an optional widened-bound helper argument, the original data/synthetic runs were repeated and their outputs remained byte-identical. The original C2-imported helper is reconstructed by removing only that later unused optional argument/branch and verified against the recorded initial code hash; its exact snapshot is preserved separately. No author contact, external message or publication occurred.
