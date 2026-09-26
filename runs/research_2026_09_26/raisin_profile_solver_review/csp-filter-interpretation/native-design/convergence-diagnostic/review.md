# CSP native covariance iteration: numerical convergence, not a new distance likelihood

The nominal **42-object** continuation passes the frozen numerical gates without changing the old .001-mag thresholds, dropping objects, altering photometry or selecting a different covariance model. Every default state in iterations 1–3 exactly reproduces the preserved three-iteration run; every default state in 1–9 is identical between the nine- and twelve-iteration executions. All 168 fitted object/job combinations have ERRFLAG0 and retain fixed stretch1, AV0 and the released header peak.

This does not repair the separate **2007A input-lineage discrepancy**: the released header peak differs from the archived raw fit header, and its original archived-Q gate remains failed. Nine objects also retain first-iteration epochs outside the tabulated template grid. The result is a converged numerical recipe on the current released inputs, not an exact reproduction of every historical input or empirical validation of the model.

## What the historical code does

The source chain is preserved in [source-proof.json](source-proof.json). `SNANA_DRIVER` loops until `ITER == NFIT_ITERATION`; nominal NML requests three iterations, and the compiled `MXITER` is12. No test of consecutive distances or flux-covariance matrices controls this stopping condition. `FITPAR_ANA` can repeat a fit when the MINUIT parameter-error covariance status is poor; that is different from convergence of the observation-space flux covariance.

At each iteration `FITPAR_PREP` calls `FITINI_COV` using the previous native state before minimizing. With nominal `OPT_COVAR_FLUX=0`, model errors still enter the diagonal and Milky Way uncertainty still contributes correlated flux covariance. The native matrix is frozen during that minimization. Its diagonal model and correlated MW terms depend on previously fitted flux, hence on distance. `OPT_COVAR_FLUX=2` would update the matrix during objective evaluation, but the corresponding full model-covariance route is not provided for SNooPy; toggling it is not a supported convergence repair.

For these fixed shape, AV, peak and final epoch masks, source formulas imply, apart from saved single-precision rounding,

`f(D) = exp[−K(D−Dref)] fref`,

`C(D) = E + exp[−2K(D−Dref)] M`,

where `K=ln(10)/2.5`, E is the distance-independent native data/fudge covariance and M contains the model-error and MW covariance at the reference amplitude. The source squares the stored data/fudge errors in single precision; the algebra check preserves that convention. The independently reconstructed C2→C3 transition and all continuation transitions close with maximum relative Frobenius error **8.0043e−8**, beneath the separately frozen 5e−6 representation check. No covariance is projected onto a positive-semidefinite approximation or otherwise repaired.

The original full-cohort maximum relative changes were **0.1990230 for C** and **0.2028272 for W=C⁻¹**, both at2009ab. These are distinct metrics; the larger number previously labelled covariance was the weight-matrix change.

## Frozen experiment and result

[protocol.json](protocol.json), SHA256 `85766863c80e41a4de7f057a7b43bf256cc6451c27cea71ee2d865e8185e9061`, specifies four nominal-only jobs: nine iterations, twelve iterations, and twelve iterations from actual first-iteration post-initialization distance shifts of −.2 and +.2 mag. [pre-execution-amendment.json](pre-execution-amendment.json), SHA256 `af8587c90bbf145f27eca440c618f9e04ce172482b1029cd79f6ac0b2bbbdd64`, explicitly separates initial empirical support from the final objective and freezes the runner. Root reviewed66 hashes and released only these four jobs in [nominal-continuation-release.json](../../../../csp_native_filter_response/root-state-review/nominal-continuation-release.json). The saved release SHA is `7e286a9ecf125a3bd1f7ae6d8898bd5efa84231fdb24e87947f40490a2e00797`.

| Check | All42 result |
|---|---:|
| D12−D9 | Exactly0 at exported precision |
| Last two distance increments, all three starts | Exactly0 at exported precision |
| Maximum final distance range across starts | 2.529e−8 mag |
| Maximum last-two covariance update, whitened operator norm | 1.405e−7 |
| Maximum reconstructed native objective discrepancy | 1.547e−11 |
| Mean D12−D3 | +0.0000413886 mag |
| Maximum absolute D12−D3 | +0.001231436 mag, 2008ia |

The next largest absolute changes from three to twelve iterations are 2006kf +.000251817 mag, 2009ab +.000131280, 2004ey +.000097646 and 2007ca −.000090812. These are changes of numerical recipe on the same data, not inferred physical corrections. Reaching equal distances at the saved precision does not mean exact real-arithmetic convergence; the remaining covariance changes are at approximately single-precision storage scale. The [case-level results](result.json), [compact summary](summary.json) and [complete state histories](state-histories.json) retain every object and each tested start.

The four native jobs took **12.0312 seconds**, increasing the carried total to **28.5973/600 seconds**. There were no new execution failures. The unmodified historical binary and the previously validated instrumented binary remain unchanged. The exclusive native worker was released before postprocessing.

## Initial extrapolation and final support

Source `genmag_snoopy.c` clamps the bracketing index at the grid edge and then uses the unrestricted interpolation ratio from the final two packed **magnitude and magnitude-error** nodes; its out-of-range ratio abort is commented out. This defines a linear extrapolation in those tabulated quantities. It is not a claim that the template was trained or validated at those phases.

The nine first-iteration cases are 2007S,2008ar,2006et,2007as,2008bc,2008hv,2009aa,2006kf and2009al, reaching **+93.9151 rest-days** against a source grid ending at+70. They remain in every original and continuation initialization. No rows were removed or clamped. All accepted epochs from iteration2 onward lie inside the source grid, and actual nonzero passband support remains inside the native wavelength screen. Final means and covariance matrices remain positive.

The three starts and9/12 comparison show independence from this **tested set of amplitude initializations** at the quoted precision. They do not prove global basin uniqueness or make the extrapolated initial templates empirical evidence. Because only distance varies after fixing peak/shape/AV and selecting the final mask, the unsupported initial rows feed later iterations through the earlier amplitude/covariance state, not as retained final measurements. Initial empirical support remains explicitly false for these nine cases even though numerical convergence and final-objective support pass.

## Fixed point versus a normalized Gaussian objective

Let `r=y−f(D)` and `W=C(D)⁻¹`. The converged native scheme solves the mean estimating equation

`f_Dᵀ W r = 0`.

For a separately posited Gaussian measurement model with the same parameter-dependent C, the normalized objective would instead be

`L(D)=rᵀWr + log(det C)`

(up to constants and any specified prior). Its derivative is

`L_D=−2 f_DᵀWr − rᵀW C_D W r + tr(W C_D)`.

Here `f_D=−Kf` and `C_D=−2K(C−E)`. The [case-level decomposition](result.json) reports all three terms separately, reconstructing C at the final D rather than confusing the iteration's previous-state matrix with a self-consistent C. The largest absolute mean-score term is .0001545 per mag, while the sum including covariance and determinant terms ranges from4.04 to5318.96 per mag. This demonstrates an **equation/objective distinction**, not a measured distance bias or a recommendation to maximize the alternative objective. No alternative maximum-likelihood distance was fitted.

The nominal historical `DOCHI2_SIGMA` is false. If enabled, its source term is `2 sum log(error_i/current_previous_error_i)`, using the diagonal data/model/fudge error convention. It omits the determinant of the full correlated MW covariance, and at option0 it does not make the covariance in the quadratic vary continuously within minimization. Turning it on would change the estimating procedure; it is not equivalent to adding the full Gaussian normalizer and is not part of this diagnostic.

The difference between these equations does **not** alone establish bias of the native estimator. Under a correct conditional mean, the mean estimating function has zero expectation at the generating parameter; that statement is weaker than finite-sample unbiasedness of its root. A correctly specified normalized Gaussian model would likewise have cancellation, in expectation, between its covariance quadratic and determinant scores. Neither expectation is established by an observed residual or by these empirical error tables. Selection, model misspecification, trained error interpretation and finite-sample recovery need their own generative checks.

## Meaning for the filter response

A subsequent filter processing comparison should compare matched **converged nominal versus converged relabelled** branches, retaining the original three-iteration outputs as the historical numerical recipe. A three-iteration response can be reported as a deterministic response of that recipe, but it must not be silently labelled a converged correction after its gate failed. Comparing a twelve-iteration changed branch against a three-iteration nominal branch would mix numerical and filter changes.

All42 objects, including2007A and the nine initial extrapolation cases, remain in the accounting. The current experiment does not erase the2007A historical-header mismatch, improve model goodness of fit, establish physical WIRC/RC2 equivalence, regenerate bias simulations or provide a corrected cosmology covariance. No full42 changed-filter fit or alternative Gaussian-MLE fit was executed here.

Root independently reproduced all2016 callback quadratics, exact prefixes, covariance norms, masks, starts and final-phase gates: [independent continuation review](../../../../csp_native_filter_response/root-state-review/continuation/result.json). Maximum independent quadratic discrepancy is1.54614e−11.
