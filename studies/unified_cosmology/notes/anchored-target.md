# An absolute-calibration alternative to the Dovekie target

The public Pantheon+SH0ES release permits a separate joint CMB+BAO+SN target with absolute distance calibration. The new adapter uses **1,657 measurements: 77 calibrator rows and 1,580 other rows**, retaining repeated light curves and the full covariance. It replaces the Dovekie likelihood. It does not add an H₀ Gaussian, another SN compilation, geometric-anchor data, or our separately reconstructed 37-host Cepheid factor.

The [calibration-interface audit](../results/distance_ladder/calibration-interface.json) supplies the observational model. Calibrator rows use their released `CEPH_DIST`; other rows use the shared CAMB background as

$$
\mu=5\log_{10}\!\left[(1+z_{\rm HD})(1+z_{\rm HEL})D_A(z_{\rm HD})/{\rm Mpc}\right]+25.
$$

One common absolute magnitude M is integrated with the declared flat measure. The likelihood retains the full covariance determinant, the intercept-information term and the **N−1** Gaussian normalization. The flat-M measure has no proper evidence normalization. Unlike an unanchored distance-only sample, the mixture of calibrator and Hubble-flow rows leaves absolute H₀ information after this integration. The underlying covariance is retained with the calibration audit's explicit symmetric averaging of tiny printed antisymmetries; raw bytes and the small literal-official projection discrepancy remain documented. This is a conditional released-data target, not a new reduction of Cepheids or SN photometry.

`anchored_adapter.configuration` builds the existing `modern_fast` configuration, then replaces only `released_sn`. The LCDM and CPL alternatives preserve every CMB/BAO likelihood, sampled prior, fixed nuisance parameter, neutrino setting, flatness assumption and physical implementation-domain restriction. Epsilon remains fixed at zero. Additional luminosity-drift models are not included in this preparation. An optional existing spectral surrogate changes only the numerical proposal and would require its own validated native correction for this new target.

The Cobaya wrapper requests distances at unique noncalibrator redshifts, restores their full original order, and delegates the Gaussian calculation to the independently audited `ReleasedCalibration` class. It exposes the same `sn_chi2` derived diagnostic as the existing SN interface. Calibrator rows never receive cosmological distances. Nonzero epsilon is rejected before requesting a background.

The [wrapper validation](../results/inference/anchored-validation.json) compares native and numerical-proposal configurations for both LCDM and CPL. After removing the replaced SN block, each configuration must equal the original target exactly. It checks provider ordering and likelihood normalization against an independent full-matrix calculation using eight synthetic distance curves, verifies current external likelihood assets and package versions, and rejects changed source/data/calibration records and unsupported settings. Native CAMB, background evaluation, model construction and surrogate initialization are forbidden during these checks. These are interface checks, not cosmological posterior measurements.

## Why the original postprocessors cannot yet be reused

| Component | Reuse status |
|---|---|
| Modern CMB/BAO configuration, fast lensing and numerical spectrum class | Reused unchanged by the new adapter; identities include the new SN source and complete release lineage. |
| Chain diagnostic arithmetic | Potentially reusable once a separately validated anchored sampler produces the same chain format. No anchored chains exist yet. |
| Original sampler drivers and `sample_path` | Not compatible: their dispatch accepts the original unanchored samples and factories. |
| Original native correction and `measurement_summary` | Not compatible: backend selection reconstructs the original factories and identities. They must not consume the new preparation manifest. |
| Luminosity histories, SN predictive checks, omissions and related descendants | Require an explicitly anchored parent qualifier and review of calibrator handling; changing a sample label is insufficient. |

The new driver therefore has only a **prepare** action. It writes immutable configuration and target-identity files into a fresh ignored directory, requiring a current passing wrapper validation. It offers no sampling or likelihood-evaluation command. A new correction and qualification route must be reviewed before this alternative becomes an operational measurement. Nothing modifies the four active Dovekie-based targets or their frozen source files.

From the repository root:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/anchored_validate.py --surrogate .work/unified-cosmology/inference/quartic-spectral/surrogate-quartic-v1.npz --output studies/unified_cosmology/results/inference/anchored-validation.json
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/anchored_run.py prepare --model lcdm --output .work/unified-cosmology/anchored-preparation/lcdm-native
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/anchored_run.py prepare --model cpl --output .work/unified-cosmology/anchored-preparation/cpl-native
```

The source acquisition and covariance checks must first pass through the separate calibration interface. Existing preparation directories are never overwritten. The [design](../code/inference/anchored-design.json) fixes the limited scope; no additional independence or covariance-decomposition claim follows from preparing this alternative.
