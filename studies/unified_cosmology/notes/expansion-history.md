This calculation turns a **qualified native-CAMB corrected posterior** into the past expansion history. It reports H(z), the deceleration parameter q(z), and the jerk j(z) at 25 fixed redshifts between 0 and 2.33. The sign q<0 denotes positive acceleration of the scale factor. The jerk describes its third time derivative; j<0 means decreasing scale-factor acceleration, and is not the same quantity as the time derivative of q.

The identities used are

\[
q=(1+z)H'/H-1,\qquad
j=1-2(1+z)H'/H+(1+z)^2[(H'/H)^2+H''/H].
\]

Five-point derivatives use central stencils except at z=0, where a forward stencil avoids future redshifts. Repeating with half the step checks numerical sensitivity. The calculation checks the stored native q and j values at z=0, 0.5 and 1, along with the native matter density and drag sound horizon. Radiation and massive-neutrino evolution remain in the shared CAMB background; the calculation does not replace them with a matter-only approximation. [CAMB’s background definitions](https://camb.readthedocs.io/en/stable/results.html) specify H in km/s/Mpc and the component densities needed for the independent Friedmann closure.

Four saved full-native controls exposed a numerical-path discrepancy before any qualified history was produced: an ordinary background calculation differed in the drag sound horizon by up to 1.23 × 10⁻⁶ relative, exceeding the unchanged 10⁻⁷ gate. CAMB's [full nonlinear-lensing calculation](https://github.com/cmbant/CAMB/blob/1.6.6/fortran/camb.f90#L104-L106) temporarily enables `WantTransfer`, then restores the original flag. That changes the [thermal starting time](https://github.com/cmbant/CAMB/blob/1.6.6/fortran/cmbmain.f90#L801-L818), its sampling grid and the finite-tolerance drag-redshift root bracket. Copying the returned parameters or performing one native initialization did not remove the discrepancy. The background adapter now mirrors that temporary flag under the same nonlinear-lensing condition. All four stored native sound horizons then agree bitwise; no spectrum or matter-transfer calculation, tolerance relaxation or physical-prior change is involved. The [fixed controls and source hashes](../code/inference/expansion-history-native-controls.json) preserve the previous failures and the exact CAMB source version.

The reported median and 68%/95% equal-tail intervals are **pointwise conditional posterior intervals**, not simultaneous confidence bands. Acceleration fractions use the original untrimmed exact/proposal weights, with separate chain fractions and batch Monte Carlo errors. A fraction of one is not certainty. “Any acceleration” refers only to at least one of the listed past redshifts. No future extrapolation, Gaussian significance conversion or model probability is provided. Parent convergence, exact-weight stability, input identities and additional stability of these derived histories must all pass.

Cosmic age is the integral of da/[aH(a)] since the hot big bang under the assumed flat GR CMB model and its fixed early-universe and neutrino physics. It is **not an independently measured stellar-population age**. A synthetic test found that CAMB 1.6.6’s default age integral was lower than converged integration by 0.00028147 Gyr (about 0.281 Myr). Its [default Romberg tolerance](https://github.com/cmbant/CAMB/blob/1.6.6/fortran/results.f90#L675-L684) is 10⁻⁴ at AccuracyBoost=1. The initial failed precision check is retained in the validation record. This summary therefore integrates the same unchanged H(a) with adaptive quadrature, checks a separate 256-node Gauss–Legendre rule, and records the default-CAMB age difference. Strict native Romberg integration also agrees in the synthetic validation. No likelihood or physical model is changed.

Fifteen synthetic power-law, radiation/matter/CPL and modern-physics cases test the equations, derivative precision, Friedmann closure and age integrals. Four independently saved full-native controls additionally test the thermal adapter, and seven flag-condition checks preserve unaffected parameter settings. Known weighted quantiles and sign fractions are checked separately. These tests establish implementation accuracy; they supply no measured cosmological history or age. The validation is rerunnable with:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/expansion_history_validate.py
```

Once a correction summary passes every registered gate, run:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/expansion_history.py \
  --chain-folder PATH_TO_QUALIFIED_CHAIN_FOLDER \
  --correction-summary PATH_TO_QUALIFIED_CORRECTION_SUMMARY \
  --cache .work/unified-cosmology/inference/expansion-history/TARGET_NAME \
  --output studies/unified_cosmology/results/inference/expansion-history-TARGET_NAME.json
```

The same command reuses checked background records. A changed input, source, environment or target requires a new cache; named correction folders are read from the supplied summary. The consumer has no provisional or failed-parent bypass. Larger point records and their hash ledger remain under `.work`; the result retains their lineage and all qualification failures.
