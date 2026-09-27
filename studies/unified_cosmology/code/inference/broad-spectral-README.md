# Broader numerical spectrum proposal

This is a computational extension for the existing native modern CMB/BAO target. It changes the support of the polynomial spectrum approximation, not the scientific parameter priors, likelihood inputs or cosmological assumptions. The motivation was frequent expensive native fallbacks in exploratory luminosity-evolution proposal runs. Those earlier models and pilot outputs remain separate.

The [frozen design](broad-spectral-design.json) keeps the original acoustic-coordinate centre and changes the numerical Cholesky factor to

$$L_{\rm broad}=L_{\rm original}\,\operatorname{diag}(1,1,1,1,1,1,3,3).$$

Only the last two conditional directions, corresponding to w0 and wa, are widened. The covariance among the first six coordinates and their cross-covariance with the last two are unchanged. This is not a uniform threefold increase in the marginal physical-parameter standard deviations.

There are 768 new training attempts and 128 independent holdout attempts, with seeds 2727701 and 2727702. Every requested point is retained. Native prior/theory failures, including CAMB's high-redshift w0+wa restriction, are recorded without replacement draws. A cubic polynomial with 165 coefficients per spectral element was chosen before new spectra or holdout errors; at least 330 finite training rows are required. It uses all finite training rows and no holdout rows. The unchanged native likelihood uses the already validated exact reassociation of fixed lensing response matrices.

After the [modern external inputs](../external_probes/README.md) are installed, run from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/broad_spectral_training.py --freeze
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/broad_spectral_training.py --workers 8
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/broad_spectral_fit.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/broad_spectral_check.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/broad_spectral_audit.py
```

All new arrays and per-attempt manifests are under ignored `.work/unified-cosmology/inference/broad-spectral/`. Native acquisition can resume matching completed attempts after verifying their source and data identities. A partial `--limit` run is explicitly labelled partial and cannot enter the final fitter. The fitter refuses to overwrite an active model; the first holdout record is likewise preserved. Re-evaluation or a new model choice requires an explicitly named version and an honest record of which validation outcomes have already been inspected.

The `.npz` interface is unchanged: `centre`, `coordinate_cholesky`, `exponents`, `coefficients`, `output_scale`, `length`. It can be consumed by the existing `SpectralSurrogate` without modifying that class or any active target. Outside its numerical envelope, the provider still calls native CAMB. The new fit helper never imports or opens holdout spectra. The audit independently reconstructs all requested points, verifies unchanged target-source identities, checks the separation of train and holdout, and compares the SVD coefficients against a rank-revealing QR solution on predetermined spectral outputs.

Holdout likelihood errors are measured at fixed reference nuisance values. Even good holdout results do not establish posterior accuracy, convergence or importance-weight overlap. Every scientific use still requires the declared independent-chain convergence diagnostics and untrimmed native exact-posterior correction. No cosmological result follows from the training size, computational speed or numerical validation status alone.
