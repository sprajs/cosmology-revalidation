# Replacing the supernova sample with absolute calibration

The [anchored target](anchored-target.md) can be evaluated as a separate conditional CMB+BAO+Pantheon+SH0ES alternative using already qualified Dovekie samples. This is a replacement of the entire SN likelihood. It adds neither an H₀ Gaussian nor a second Cepheid or geometric-anchor factor. The release already contains the calibrator information and its covariance.

The new [bridge](../code/inference/anchored_bridge.py) is implemented and [synthetically validated](../results/inference/anchored-bridge-validation.json). **No observational bridge has been evaluated and no anchored joint cosmological measurement follows from these tests.** Its [design](../code/inference/anchored-bridge-design.json) was fixed before such evaluation.

For each selected point, the untrimmed log weight is

$$
\log W_{\mathrm{anchored}}=
\log W_{\mathrm{native/proposal}}
+\log L_{\mathrm{Pantheon+SH0ES}}
-\log L_{\mathrm{Dovekie,native}}.
$$

Both SN factors include their complete covariance and flat-magnitude marginalization terms. Their normalizing constants are retained; the improper flat-M measure does not define a model-evidence ratio. The actual ratio uses the **stored native Dovekie density**. A separately reconstructed Dovekie density must agree to an absolute log-likelihood tolerance of 10⁻⁵ at every point. It is a consistency check, not a replacement for the stored source density.

The helper admits only the existing CPU `modern_fast` LCDM/CPL parents with zero luminosity evolution and official Planck calibration. It first reruns `measurement_summary.summarize_run`, checking parent convergence, native correction, point records, current source, assets, versions and identities. Provisional parents are rejected before any background calculation or cache creation. It then verifies that the old and new configurations are identical after removing only the SN block. CMB, BAO, physical-domain restrictions, normalized priors and nuisance coordinates remain unchanged.

One background calculation per point supplies angular-diameter distances in Mpc at the union of old and new redshifts. The original order and repeated measurements are restored separately. The anchored factor uses 77 released calibrator distances and 1,580 cosmological distances within the 1,657-row covariance. The background uses the previously validated native thermal-initialization adapter; native spectra and transfer-function calls are forbidden. The source and target retain the same H₀ convention, km s⁻¹ Mpc⁻¹. Cosmological derived values remain those of the native parent, while the SN chi-square used in weight diagnostics is replaced by the anchored value.

Every original overlap criterion is retained: at least 2,000 native points, raw weight ESS of 400, four chains with per-chain ESS of 50 and weight fractions between 0.05 and 0.5, Pareto k below 0.7, and the declared chain-mean and batch-error checks. All prior coordinates participate in stability checks. A failed density check or overlap criterion withholds posterior summaries, covariance, sign fractions and paired parameter shifts. Diagnostic weighted summaries are explicitly not measurements. Cache files preserve each point, the parent record hash, target identity and source/data dependencies; altered inputs or payloads are rejected.

Validation recovers a known Gaussian-factor replacement from 20,000 synthetic proposal points using an independent conjugate calculation. It also checks null replacement, additive-density constants, prior/configuration mutations, deliberately concentrated weights, insufficient points, a 400-point fake-background consumer with immutable replay, and redshift/heliocentric-distance conventions using a synthetic provider. No CAMB background, CMB spectrum or observational bridge point is evaluated by these tests. Actual source-density closure and overlap remain necessary after the parent qualifies.

Reproduce the synthetic checks from the repository root:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/anchored_bridge_validate.py --output studies/unified_cosmology/results/inference/anchored-bridge-validation.json
```

The evaluation interface takes `--chain-folder`, `--correction-summary`, a fresh `--cache` under `.work`, and `--output`. A cache may be replayed only with identical identities; previously recorded failures are preserved. The original correction/report entrypoints must not be used to qualify the new target as though it were another Dovekie chain.

Finite importance diagnostics cannot establish coverage of unseen modes. If the original H₀ distribution gives poor support for the anchored alternative, the appropriate next step is separate anchored sampling. The thresholds and weights must not be relaxed to obtain a result. Overlapping SN measurements do not require an independence assumption for this replacement, because the two SN compilations are never multiplied together. The retained CMB/BAO factorization, released ladder and survey assumptions, and unknown physical luminosity evolution remain limitations. A paired mean change is not the significance of a disagreement between independent experiments.
