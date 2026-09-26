# Independent calibrated-flux likelihood: 12-SN adversarial pilot

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The independent engine closely reproduces the calibrated mean fluxes and fitted brightness, stretch and colour. It also exposes a non-equivalent model-uncertainty prescription and a serious local peak-time Hessian artifact. The held-out observations do **not** establish a preferred physical colour correction. This is an observation-level check of released calibrated fluxes, not a detector-pixel reconstruction or a selection-complete population analysis.

The pilot was fixed before this work: 12 evenly spaced redshift ranks from the original DES sample, 1,394 supplied raw calibrated-flux rows, and 528 published accepted rows. Every published accepted row matches the exact runtime row by original photometry row ID. `phase2/independent_flux/PLAN.json` records the analysis before independent outcomes; `AMENDMENTS.md` records the evidence-led corrections. Phase 1 files were not edited.

[Independent flux diagnostics](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/independent_flux/independent_flux_diagnostics.png)

## Independent construction and baseline gate

`scripts/phase2/independent_flux/engine.py` contains explicit photon-weighted passband integration, the SALT polynomial colour law and linear extrapolation, observer-frame F99 Galactic extinction, redshift/time dilation, and calibration-primary normalization. One path uses independent SciPy tensor splines for the SED; a separate path uses sncosmo's native SALT3 interpolation and integration. The optimizer, objective, conditional predictions, derivative checks and Hessian calculations are implemented independently. No SNANA fitter or forward-model library is called by these scripts.

Both paths use the actual original DES5YR M0/M1 surfaces, DECam transmission curves, primary SED and primary magnitudes. They retain `MAG_OFFSET=0.27`, `mB=-2.5 log10(x0)+10.635`, and `FLUXCAL` zero point 27.5. The effective exported Galactic E(B−V) is used once, without another factor of 0.86. Fitted redshift is heliocentric. Phase is `(MJD−t0)/(1+zHEL)`. The masks already impose the original phase and rest-band-centre cuts; no cosmological distance relation enters the light-curve fits.

At all 528 exact accepted epochs, native sncosmo minus SNANA has median fractional flux difference −0.0000250, 95th absolute percentile 0.000488, and maximum 0.001221. The largest difference is 0.0957 observational standard deviations. The fully independent SciPy interpolation/integration path has corresponding absolute 95th percentile 0.000549 and maximum 0.001796, or at most 0.1338 observational standard deviations. These are quantified approximate reproductions; they do not pass a strict all-point 0.001 relative tolerance or establish bitwise equivalence.

The explicit integral using the sncosmo interpolator agrees with sncosmo's own integral to about 6×10⁻⁶ at the 95th percentile of the plotted model grid. Refining independent integration from 2 Å to 1 Å changes flux by 1.3×10⁻⁶ at the 95th percentile. Interpolation conventions dominate the remaining difference. One non-accepted grid point has negative sncosmo model flux while SNANA prints zero; it is preserved separately and not divided by zero in fractional summaries.

The frozen-covariance refits converge from two starts to objectives agreeing within 5.6×10⁻¹¹. Maximum differences from SNANA are 0.000987 mag in mB, 0.02546 in x1, 0.000502 in c, and 0.0856 observer days in t0. The first three differences are below 0.030 of their reference marginal standard errors. The independent SciPy SED path gives similar results. All 24 full-data fits are interior with positive observed Hessians. The time prior is the actual preceding-iteration `INIVAL` centre, not the original search-peak estimate. The independent prior is a smooth Gaussian; SNANA uses its interpolated table.

## Covariance discrepancy, with the representation made explicit

An initial portable-bundle mistake was caught by comparing the runtime log with the supplied files: SNANA loads `salt3_lc_variance_0/1` and `salt3_lc_covariance_01`, while the initial bundle supplied the distinct `salt3_lc_model_*` products. Both sets are now preserved and identified. Actual SNANA fits always used the correct full-release directory. Initial confounded comparisons and exact source snapshots are retained under `phase2/independent_flux/initial_wrong_error_maps/`; they are not evidence for the final numerical claims.

After using the verified runtime maps and previous-iteration error-stretch parameter, an independently transcribed SNANA sigma formula matches exported model magnitude errors: median relative difference 0.0000213, 95th absolute percentile 0.000478, maximum 0.001213. Reconstructing the **full model flux covariance** with the exported preceding fit and the nonlinear magnitude-to-flux conversion matches the actual frozen matrix to relative Frobenius errors 0.0000245–0.000192. Galactic and photometric terms are checked separately. Thus there is no unexplained calibration or units compensation in this comparison.

Native sncosmo with the same runtime maps has model magnitude errors lower by a median 34.7%, with maximum absolute discrepancy 52.8%. This is primarily a different transport of the uncertainty, not a difference in mean flux. The relevant derivation is:

Let `M=M0+x1 M1`, `C=10^(−0.4 c CL)`, and let `A_i` contain the calibrated passband normalization. The observer flux is

`F_i = A_i x0/(1+z) ∫ λobs T_i(λobs) MW(λobs) C(λrest) M(p,λrest) dλobs`.

For a genuine rest-basis covariance `K_ij(λ,λ′)`, linear propagation gives

`Cov(F_i,F_j) = A_i A_j x0²/(1+z)² ∬ [λ T_i MW C] K_ij [λ′ T_j MW′ C′] dλ dλ′`.

Time dilation changes the phases at which both the mean and covariance are evaluated; it supplies no additional independent noise factor. Multiplicative flux transformations act on both the mean and its uncertainty. The scalar effective-wavelength approximation to the basis variance is `V=V00+2 x1 V01+x1² V11`. With normalized photon weights and ignoring Galactic extinction just for this comparison, native sncosmo uses approximately `V/⟨M⟩²` for relative variance. Executable SNANA instead uses `V/[⟨CM⟩/(1+z)]²`; its denominator retains colour and the redshift factor despite the nearby comment describing a no-colour quantity. The separate same-band colour-dispersion term is added in both engines.

The SALT3 formalism describes variance in the underlying flux surfaces and includes the squared colour multiplier in its band-variance expression. [Kenworthy et al. 2021, equations 3–4](https://arxiv.org/pdf/2104.07795), locally preserved as the explicitly identified v1 paper, supports this rest-basis interpretation. The [sncosmo SALT3 implementation](https://sncosmo.readthedocs.io/en/stable/api/sncosmo.SALT3Source.html) is independently inspected in the pinned installed package. In a constant-spectrum test with directly simulated rest-basis perturbations, propagated fractional variance remains invariant under a common colour/redshift rescaling; the SNANA-denominator expression instead changes by `(1+z)²/C²`. This demonstrates non-equivalent prescriptions under that stated interpretation. It does not establish the authors' intent, fully reconstruct the DES model training, or prove that the native approximation describes real supernova residuals better. Neither approximation supplies the full wavelength covariance in the double integral.

The source-level difference looks large in model-only sigma, but measurement noise dilutes it. Across objects, the median native/official **total** per-point sigma ratio is 0.948–1.008. Model variance contributes a median 1.3%–26.9% of total variance in eleven objects and 95.0% in the lowest-redshift object. The same-band colour-dispersion term contributes object-median 0.4%–13.5% of model diagonal variance, so it is not the dominant diagonal term here; it remains correlated and is retained. The Galactic term contributes at most 0.102% of median total variance in these objects.

Holding the native covariance fixed changes fitted mB by at most 0.00493 mag, x1 by 0.0691, and c by 0.00512. Corresponding standard-error ratios relative to the independent fit with official frozen covariance range from 0.861–0.985 for brightness, 0.913–1.006 for stretch and 0.890–0.995 for colour. Error changes exceed mean-parameter changes in this pilot.

## A positive Hessian is not sufficient validation

For mB,x1,c, the independent-to-SNANA generalized covariance eigenvalues span 0.9461–1.1997; absolute correlation differences reach 0.0521. Some joint variance directions therefore differ by 20%, although marginal standard errors differ by only about 4%. A downstream three-dimensional summary likelihood should retain this covariance sensitivity, not infer exact joint agreement from the diagonal alone.

The full four-dimensional matrix has a stronger failure. For CID1896213, SNANA's local Hessian gives σ(t0)=0.123924 days, while the independent smooth Hessian gives 0.564431 days. Profiling the other three parameters gives Δχ²=0.9975 and 1.0028 at displacements of ±0.564431 days. SNANA's separately reported MINOS error is 0.564311 days and the published scalar error is about 0.5644. The profile and MINOS agree; the local Hessian time error does not. Three other pilot objects have smaller versions of this discrepancy. The evidence is consistent with local curvature sensitivity to the tabulated interpolation, but this work does not isolate every internal Minuit derivative step. No covariance was silently repaired.

Halving the independent Hessian differencing steps changes its matrix by at most 3.6×10⁻⁶ in relative norm on the native interpolation path. Using independent SciPy SED interpolation retains the time-width discrepancy. Full-matrix positive definiteness alone cannot certify that all fitted uncertainties are trustworthy.

## Held-out mean and noise tests

The fixed epoch split hashes CID and integer observing night; all bands of a night stay together. It holds out 114 observations across twelve objects. The secondary leave-one-band-out experiment has 35 folds across eleven objects with at least three bands. It does not treat overlapping folds as independent experiments. Each score uses the same held-out rows, full train/test covariance conditioning and local observed-Hessian parameter propagation. Plug-in scores are also saved. All fits needed for these paired comparisons succeeded without a boundary or nonpositive Hessian.

The alternatives were declared before their outcomes: a fixed Rv=3.1 F99 rest-frame colour curve replacing SALT's effective colour law, with its value at the B reference wavelength subtracted to remove the amplitude gauge; and empirical `c(p)=c+d tanh(p/20)` with d~Normal(0,0.10), bounded to ±0.30. F99's fitted signed coefficient is an effective colour-shape coordinate. It is **not** a measured nonnegative dust column, a dust-population decomposition, or evidence for a physical Rv.

| Paired mean comparison | Epoch Δ log predictive density | Band-held-out Δ log predictive density |
|---|---:|---:|
| F99 − SALT | +0.312 ± 0.632 | +2.104 ± 2.995 |
| Phase colour − SALT | −1.304 ± 0.689 | −0.785 ± 3.989 |

The ± values are object-cluster standard errors of the summed differences, not posterior odds or significance certificates. Shared calibration uncertainties are not represented by treating objects as independent clusters. F99 improves epoch prediction in only 5/12 objects; CID1896213 supplies +0.591 while CID1330817 supplies −0.133. The flexible phase law has its strongest band-held-out gain at CID1480038 (+3.328), despite its negative aggregate epoch result. The heterogeneity is retained.

Replacing official frozen noise with native frozen noise gives Δ log predictive density +0.704±0.692 for the same SALT mean. Conditional held-out χ² is 117.33 for 114 points with official noise and 123.19 with native noise. Both remain plausible at this resolution; the pilot cannot empirically choose the transport convention. F99 minus SALT is +0.328 under native noise, again small. A normalized residual distribution alone is not a proof of independent Gaussian observations, because conditional correlations and parameter uncertainty are present.

These scores condition on the published SALT-selected mask, full-data reference parameters used to construct frozen noise, and the preceding full-data time-prior centre. Consequently they contain nuisance/mask information from the full data. They are controlled conditional sensitivity tests, **not** a clean prospective validation or a selected-population likelihood. No released distance correction, LCDM fit, or acceleration criterion is used to choose a colour law.

## Recomputed covariance and synthetic counterevidence

Recomputing parameter-dependent noise is a different objective from the original frozen iteration. All 48 separately registered recomputed fits converge, but omitting versus including the Gaussian log determinant matters:

| Covariance convention | Maximum absolute ΔmB, recomputed χ² only | Maximum absolute ΔmB, recomputed normalized Gaussian |
|---|---:|---:|
| Native sncosmo | 0.03452 mag | 0.00372 mag |
| Independent SNANA formula and nonlinear flux conversion | 0.05467 mag | 0.00881 mag |

These shifts are relative to the reference fitted parameters. They do not mean that the original frozen SNANA objective is secretly the recomputed χ²-only arm. A determinant constant within a frozen fit cannot change that fit; it can change a fit whose covariance varies with trial parameters. Reporting these arms separately prevents a silent likelihood replacement.

Synthetic draws use the actual twelve cadences with the exported frozen covariance as the stipulated generating law. Noiseless fits recover the generating parameters to 3.6×10⁻¹² in the stated scaled coordinates. An independent Cholesky generator with 10,000 linearized draws per cadence reproduces analytic sampling standard deviations to 1.5%. The distinction between inverse posterior Hessian and fixed-prior repeated-sampling covariance is calculated explicitly.

Ten nonlinear noisy realizations per cadence give 120/120 interior converged fits with positive Hessians under each noise arm. Under the generating frozen covariance, pull standard deviations for ln(x0),x1,c,t0 are 1.071, 1.008, 1.099, 0.953; nominal 95% coverage is 93.3%, 95.8%, 92.5%, 97.5%. Using the smaller native covariance on those same draws broadens the pull dispersions to 1.177, 1.071, 1.179, 1.036 and coverage becomes 90.8%, 95.0%, 91.7%, 95.0%. These are conditional recovery diagnostics with limited repetitions. Generating with official covariance cannot establish that official covariance is physically correct.

The strongest counterevidence to treating the real-data phase-colour loss as exclusion is the deliberate synthetic injection: even when d=0.06 is present, the same twelve-object held-out test scores the flexible model **below** SALT by 3.83±2.47. Sparse training data, uncertain shape/colour/time and the extra predictive uncertainty limit power. This pilot cannot reliably exclude that amplitude of phase-colour behavior. More events, stronger independent phase coverage and genuinely training-only noise/mask construction are needed before a physical interpretation.

## Reproducibility and remaining boundary

Run from the repository root with the isolated environment already pinned by the official reconstruction:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/phase2/independent_flux/fixed_check.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/phase2/independent_flux/exact_check.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/phase2/independent_flux/fit.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/phase2/independent_flux/noise.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/phase2/independent_flux/validate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/phase2/independent_flux/plot_results.py
```

The portable fixture must include the exact objective, prior-setup and covariance-component exports prepared by the official reconstruction. Code/package/input/output hashes accompany each run; loaded source snapshots are retained under `phase2/independent_flux/provenance/`. The downloadable diagnostic figure has PNG and PDF forms. Main machine-readable evidence is `fixed_summary_exact.json`, `independent_refits.json`, `heldout_predictions.json`, `fit_summary.json`, `noise_results.json`, `noise_summary.json`, `validation_results.json` and `validation_summary.json`.

This work validates a bounded independent flux likelihood, identifies concrete covariance limits, and actually tests competing mean/noise prescriptions on observations. It does not retrain SALT, remove shared calibration/training assumptions, reconstruct selection/classification, infer a unique dust or age mechanism, or break luminosity-evolution versus distance degeneracy. Scaling to the full sample must retain these failures and the uninformative counterexamples rather than promote this pilot to a physical correction ranking.
