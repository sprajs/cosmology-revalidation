# Cosmology from the calibrated Pantheon+ supernova sample

The released Pantheon+SH0ES distances give **H₀ = 73.550 ± 1.017 km/s/Mpc** and **Ωₘ = 0.33245 ± 0.01806** in flat matter-plus-Λ cosmology. These are posterior means and standard deviations. The equal-tail 95% intervals are **[71.576, 75.563] km/s/Mpc** and **[0.29770, 0.36847]**. The implied present deceleration parameter is **q₀ = −0.50133 ± 0.02708**, with 95% interval **[−0.55345, −0.44730]**. This model therefore implies present acceleration. It fixes the expansion-history shape and does not test arbitrary dark-energy evolution.

This calculation uses **1,657 released measurements**, including **77 calibrator light-curve rows**. Those rows refer to 42 physical supernovae in 37 hosts; they are not 77 independent calibrating galaxies. The full released covariance, including calibrator-to-other-supernova correlations, enters once. The [interface audit](calibration-interface.md) documents the selection, host associations, tiny printed covariance asymmetries, and unresolved details of the embedded Cepheid covariance. No separate H₀ Gaussian or extra Cepheid factor is added.

## Relation to published results

[Brout et al., Table 3](https://arxiv.org/html/2202.04077v2#S4.T3) report **73.6 ± 1.1 km/s/Mpc** and **0.334 ± 0.018** for Pantheon+SH0ES in flat ΛCDM. Our central values agree closely. Our H₀ standard deviation is about 8% smaller than the printed 1.1, so rounding alone is not an adequate explanation.

We recovered the official 160,000-row flat-ΛCDM chain and its embedded configuration. Its full unfiltered H₀ standard deviation is **1.068**. It contains a clear initialization transient: the first walkers start near M = −19, far from the later distribution near −19.244, with much lower posterior density. Fixed removal of the first 10%, 20%, 30% or 50% gives H₀ standard deviations **1.014, 1.016, 1.019 and 1.014**, respectively, and means **73.527–73.538**. All these widths are consistent with our deterministic **1.017**. The original chain is preserved intact; these are descriptive sensitivity checks, not a selection of the prefix that best matches our result.

This comparison does **not** establish which burn-in or summary procedure produced the paper's table. The embedded historical configuration also uses finite bounds on M, different broad parameter bounds, and a CAMB background. It declares neutrino settings whose actual mapping through the historical interface has not been fully recovered. Our bounded late-time calculation omits radiation. Consequently this is not an exact replay of the historical posterior density. The available evidence supports numerical consistency, not a claim that the paper's uncertainty is erroneous. [Released-chain analysis](../results/distance_ladder/calibration-chain-review.json).

A separate [background sensitivity check](calibration-background.md) tests three explicit neutrino/radiation interpretations at nine points around the solution. The largest distance-modulus change is **0.0003245 mag**, and the changes in log likelihood range from **−0.01085 to +0.01370**. Local curvature estimates change the H₀ uncertainty by only about **0.0004–0.0006%**. This supplies no indication that the late-time background approximation explains the quoted uncertainty difference near the solution; it is not a newly integrated posterior or a bound on remote tails.

The result differs from the reconstructed **73.043 ± 1.007 km/s/Mpc** [distance-ladder matrix](distance-ladder.md) because that target uses a different selected Hubble-flow set and a fixed low-redshift expansion prescription. These are overlapping observations under different analyses, not two independent H₀ measurements to average.

## Model and calculation

For noncalibrators,

$$D_L=\frac{c(1+z_{\rm HEL})}{H_0}\int_0^{z_{\rm HD}}\frac{dz}{\sqrt{\Omega_m(1+z)^3+1-\Omega_m}}.$$

Calibrators instead use their released Cepheid distance moduli. One common absolute magnitude M is integrated with a flat measure. The priors are uniform in **H₀ ∈ [50, 90] km/s/Mpc** and **Ωₘ ∈ [0.01, 0.99]**. The flat-M measure does not define an absolute model evidence. No CMB, BAO, stellar-age constraint or age-dependent brightness correction enters.

At fixed Ωₘ, let η = 5 log₁₀(H₀/70). After eliminating M, the residual quadratic has the form c − 2bη + aη². Integrating in a flat-H₀ measure requires the Jacobian

$$\frac{dH_0}{d\eta}=\frac{\ln 10}{5}\,70\exp\!\left(\frac{\ln 10}{5}\eta\right).$$

It shifts the conditional Gaussian mean in η by (ln 10/5)/a. The resulting finite-bound H₀ distribution is a truncated lognormal. Its exact moments and CDF are then integrated over Ωₘ using deterministic quadrature; there is no sampling uncertainty.

The production calculation passes direct likelihood comparisons, 64/128-node distance refinement, 256/512-node Ωₘ refinement, adaptive normalization checks and direct H₀ integration. Expanding the priors to H₀ ∈ [30,110] and Ωₘ ∈ [0.001,0.999] changes the reported H₀ summaries by less than **2 × 10⁻¹⁰ km/s/Mpc**. This is numerical boundary insensitivity within the model, not robustness to different physics or calibration.

A separate implementation solves the full two-parameter (M,η) generalized least-squares system by LU decomposition, uses different distance and Ωₘ quadratures, and checks normalized likelihoods. H₀ means, standard deviations and quantiles agree within **5.1 × 10⁻¹⁰ km/s/Mpc**; Ωₘ summaries agree within **2.9 × 10⁻¹⁴**. These small differences test arithmetic, not observational accuracy. [Measurement](../results/distance_ladder/calibration-lcdm.json); [independent review](../results/distance_ladder/calibration-lcdm-review.json).

The subsequent [published-table and construction search](calibration-construction-search.md) confirms that the four near-doubled sibling covariance entries are not explained by the printed Cepheid-only errors. An [independent covariance audit](calibration-covariance-estimands.md) rules out omitted nuisance marginalization and common-$M$ projection as explanations, while retaining the possibility of an unseparated host-shared SN term. This does not identify a justified covariance correction; the measurement here remains conditional on the released matrix.

## Reproduction

First acquire and verify the calibration inputs as described in the [interface instructions](calibration-interface.md#reproduction-and-reuse). Use a fresh output path for the calculation; existing results are not overwritten.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/distance_ladder/calibration_lcdm.py --output .work/unified-cosmology/calibration-lcdm-reproduction.json
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/distance_ladder/calibration_lcdm_review.py --producer .work/unified-cosmology/calibration-lcdm-reproduction.json --output .work/unified-cosmology/calibration-lcdm-reproduction-review.json
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/distance_ladder/calibration_chain_review.py --acquire --quadrature .work/unified-cosmology/calibration-lcdm-reproduction.json --output .work/unified-cosmology/calibration-chain-reproduction-review.json
```

Input files and the original author chain are downloaded and hash checked, rather than committed. The compact numerical records retain input and source hashes. This is an independent calculation of shared released observations, conditional on their corrections and covariance, rather than a new independent data set.
