# Predicting the absolute calibration from other observations

This check asks whether the supernova absolute calibration agrees with the expansion model inferred from the CMB, BAO and noncalibrating supernovae. It withholds the **77 calibrator measurement rows** and uses the **1,580 noncalibrator rows** of Pantheon+SH0ES in place of Dovekie. The released cross-covariance remains in the prediction for the withheld rows. No extra H₀ or Cepheid likelihood is multiplied into an overlapping supernova sample.

This is a conditional check of the released Gaussian likelihood. Its upstream photometry, selection and covariance were already constructed using shared calibration information. Excluding likelihood rows does not recreate an independent blind reduction of raw observations. The [unresolved shared-host covariance construction](calibration-covariance-estimands.md) remains a limitation.

## The correlated prediction

Write N for noncalibrators, C for calibrators, and S for the full released covariance. Residuals exclude the common unknown absolute magnitude M. The fitted noncalibrator factor uses the **marginal covariance Sₙₙ**, not a subblock of the full precision matrix. Define

$$
A=S_{CN}S_{NN}^{-1},\quad b=\boldsymbol1_C-A\boldsymbol1_N,
\quad a_N=\boldsymbol1_N^TS_{NN}^{-1}\boldsymbol1_N,
\quad \widehat M_N=\frac{\boldsymbol1_N^TS_{NN}^{-1}r_N}{a_N}.
$$

After integrating M under the same flat magnitude measure, the conditional calibrator residual and covariance are

$$
t=r_C-Ar_N-b\widehat M_N,\qquad
V=S_{CC}-AS_{NC}+\frac{bb^T}{a_N}.
$$

The final term propagates the uncertainty in M inferred from N. The exact normalized decomposition is

$$
\log L_{\rm full,flat\ M}=\log L_{N,\rm flat\ M}
+\log\mathcal N(t;0,V).
$$

The first factor has 1,579 residual dimensions after the magnitude integral; the proper conditional factor has 77. Subtracting another intercept from the latter would change the prediction. The undefined common normalization of the original improper M measure cancels in this conditional density; no model-evidence ratio is asserted.

The predeclared scalar test concerns a common calibration shift. Its fixed vector is

$$
v=\frac{V^{-1}\boldsymbol1_C}{\boldsymbol1_C^TV^{-1}\boldsymbol1_C},\qquad
\sigma=(\boldsymbol1_C^TV^{-1}\boldsymbol1_C)^{-1/2}.
$$

For a cosmological point θ, the observed contrast relative to its prediction is vᵀt(θ), in magnitudes. Adding a constant k to the noncalibrator distance moduli moves t by $k\boldsymbol1_C$ and this contrast by k. The vector b instead describes propagation of latent-M uncertainty; using it as the contrast direction would test a different mode. The vector v depends only on the fixed covariance and row selection, and is not chosen to maximize an observed discrepancy.

If wᵢ are qualified weights from **CMB + BAO + N only**, the predictive CDF at the observed contrast is

$$
F=\sum_i w_i\,\Phi\!\left(\frac{v^Tt(\theta_i)}{\sigma}\right).
$$

The reported two-sided equal-tail area is 2 min(F, 1−F). Averaging the individual two-sided areas would answer a different question and can be badly misleading for a mixture. The full 77-dimensional conditional predictive density is also integrated. These are Bayesian predictive diagnostics, not calibrated frequentist significances, and the one-dimensional contrast cannot detect every residual pattern.

## Numerical requirements and repetition

The [design](../code/inference/calibration-holdout-design.json) precedes observational evaluation. A currently qualified Dovekie parent is mandatory. One CAMB background per selected point evaluates both releases in their original row orders. The old supernova likelihood must reproduce its stored native value before its factor is replaced. All original importance-weight and chain-stability requirements remain in force, with at least 2,000 selected points. Withheld calibrator values and their conditional likelihood never enter the new posterior weights.

Each predictive integral has a separate precision assessment. Its normalized contributions must have effective sample size at least 100, retain 5–50% mass in each chain, and have relative Monte Carlo error at most 0.2 both between the four independent chains and for 10, 20 and 40 contiguous batches per chain. Splitting each chain cannot manufacture precision when the chains persistently disagree. Logarithmic CDF and survival calculations preserve very small tails without reporting underflow as zero probability. A failed predictive check withholds that number even if the noncalibrator posterior qualifies. A failed replacement posterior withholds all posterior and predictive claims.

After restoring the inputs and qualifying a compatible parent, use a fresh cache and report:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/calibration_holdout_bridge.py \
  --chain-folder QUALIFIED_PARENT \
  --correction-summary QUALIFIED_NATIVE_CORRECTION.json \
  --cache .work/unified-cosmology/calibration-holdout/NEW_RUN \
  --output .work/unified-cosmology/calibration-holdout/NEW_REPORT.json
```

The separate [combined calibrated-sample fit](anchored-sampling.md) includes the calibrators in its posterior. Its weights cannot be substituted into this withheld-calibration prediction. Native CMB accuracy and absent cross-probe covariance remain separate qualifications.
