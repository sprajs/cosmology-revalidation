# Numerical proposal for the anchored likelihood

The calibrated CPL replacement could not be measured by reweighting the existing unanchored posterior: its raw effective sample size was only 3.74 out of 2,000 evaluated points, and its posterior is withheld in the [bridge record](../results/inference/anchored-bridge-cpl.json). This motivates a separately sampled anchored target. It does not justify interpreting the few high-weight points as a posterior or replacing the complete calibrated likelihood with an independent Hubble-constant constraint.

The first direct numerical search stopped at its declared allocation of eight successful native spectral fallbacks. It recorded 1,145 target evaluations: 1,131 finite values, 14 ordinary nonfinite/domain rejections, and no recorded CAMB or evaluation failures. Two Powell starts completed their 450-evaluation allocations; the third was interrupted. Neither gradient/Hessian stencil ran, and no proposal was frozen. The eight fallback calls took 244.3 seconds combined; other calls had a median time of 0.173 seconds. These are computational measurements, not cosmological results. The stopped manifest, ledger, initialization and result are hash-bound in the [new design](../code/inference/anchored-refined-proposal-design.json).

The separate [refinement](../code/inference/anchored_refined_proposal.py) starts from the deterministic best finite point of that saved search. It freshly verifies the same target and all input/source identities, then re-evaluates that one point and requires every component density and prior to agree within an absolute tolerance of 10⁻⁷. It runs one Powell search with at most 450 evaluations, followed by the original two central gradient/Hessian stencils. Its independent allocation is at most 1,200 target calls, 16 successful native fallbacks and 1,800 seconds after model initialization. One explicit native numerical initialization is recorded separately. Every caught CAMB failure stops construction, including failures that the frozen theory exposes as a logged failure followed by a nonfinite density.

This is a new one-start construction, not a resumed or completed version of the original multistart experiment. The original covariance supplies only affine numerical coordinates. The fixed search box does not restrict the scientific prior or the later proposal's support. The Hessian agreement, stationarity checks, limited curvature flooring and broad-covariance fallback are unchanged. The final normalized proposal is a 90% Gaussian and 10% Student-t mixture with full support; the true anchored target must still be evaluated at each sampled point. No local covariance is a measured cosmological uncertainty, and one refined location cannot establish mode completeness.

The [refinement controls](../results/inference/anchored-refined-proposal-validation.json) check the actual stopped-ledger selection without computing cosmology, a separate known Gaussian target, density replay failures, source/target tampering, resource stops and caught failures. The [bounded pilot controls](../results/inference/anchored-refined-pilot-validation.json) independently replay 400 Gaussian-test MH decisions using SciPy mixture densities, including all rejections and repeated occupied states. These controls certify numerical contracts only. The physical optimizer and any later pilot remain separate recorded computations.

The actual refinement completed 901 target calls in 217.9 seconds, with one successful native fallback and one separately recorded native initialization. The starting-density replay was bitwise identical. Powell converged after 170 evaluations. The largest checked gradient was 0.0800, but the two Hessians disagreed by 0.630 in relative Frobenius norm, exceeding the declared 0.20 threshold. The rule therefore froze the broad reference-covariance fallback, with covariance multiplied by four. This is an explicit proposal-construction limitation, not a fitted uncertainty.

The user then requested closeout while the [pilot](../code/inference/anchored_refined_pilot.py) was unfinished. Its owned process was interrupted and exited with code 130. The final preserved ledger contains **43 of 400 requested candidates**: 24 finite target evaluations and 19 nonfinite/domain rejections. It records three completed native fallbacks and 21 surrogate evaluations; a further native calculation was interrupted, so the in-flight invocation count is unknown. No completed pilot efficiency result, production chain, native correction or anchored posterior was produced. The [closeout record](../results/inference/anchored-refinement-closeout.json) binds the saved files and retains both the original resource stop and the completed refinement.

The declared pilot requests exactly 400 candidates under separately declared seeds and caps. It reports acceptance, holding, raw importance-weight concentration and runtime. It produces no posterior. Four independently initialized chains, unchanged convergence requirements, at least 2,000 native correction evaluations and all original overlap/stability gates are still required before a measurement. Released calibration covariance, shared-probe assumptions and native numerical precision remain distinct questions.

From the repository root, the two inexpensive validations are:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1
TASK_PYTHON=.work/unified-cosmology/external-probes/.modern-venv/bin/python
$TASK_PYTHON studies/unified_cosmology/code/inference/anchored_refined_proposal_validate.py \
  --output .work/unified-cosmology/anchored-local-validation/refinement-recheck.json
$TASK_PYTHON studies/unified_cosmology/code/inference/anchored_refined_pilot_validate.py \
  --output .work/unified-cosmology/anchored-local-validation/refined-pilot-recheck.json
```

The refinement requires the exact pinned original stopped-ledger inputs; it refuses substituted histories. Its actual execution uses:

```bash
$TASK_PYTHON studies/unified_cosmology/code/inference/anchored_refined_proposal.py \
  --model cpl --native-accuracy 2 \
  --surrogate .work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz \
  --reference-proposal .work/unified-cosmology/inference/independence-proposal-cpl-v1 \
  --parent-folder .work/unified-cosmology/inference/anchored-direct-cpl-v1 \
  --output .work/unified-cosmology/inference/anchored-refined-cpl-v1 --cpu 14
```

The output directory must be fresh. The original stopped inputs and all new evaluations are retained in ignored work storage; no downloadable spectra or generated arrays are added to the repository.
