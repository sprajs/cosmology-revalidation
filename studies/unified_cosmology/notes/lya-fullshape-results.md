# Replacing the overlapping Lyα distance measurement

For flat ΛCDM with the existing CMB, Dovekie supernova and other DESI distance data, the new Lyα information changes the conditional mean of $H_0$ from **68.164 to 68.103 km s⁻¹ Mpc⁻¹**. Its conditional standard deviation is **0.247 km s⁻¹ Mpc⁻¹**. This is a small shift within the existing model. The result uses an explicit Gaussian approximation to the rounded published distance summary; it is not a reconstruction of the forest observations or the authors' complete likelihood. The separate native CMB numerical-precision screen is not established by this calculation.

The [DR2 Lyα full-shape paper, Eq. 26](https://arxiv.org/html/2607.27410v3) reports, at $z=2.33$,

$$
(D_M/r_d,D_H/r_d)=(39.32,8.600),\qquad
(\sigma_M,\sigma_H)=(0.33,0.066),\qquad \rho=0.225.
$$

Its broadband Alcock–Paczynski information and BAO information have already been combined. We therefore replace the two overlapping $z=2.33$ rows in the existing 13-row BAO likelihood. We retain the other 11 rows and verify that their released covariance with this pair is exactly zero. Neither the previous Lyα factor nor a separate broadband-AP factor is added again. No growth-rate, $H_0$ or BBN constraint is added.

| Quantity | Conditional mean | Conditional standard deviation | Mean change from the parent |
|---|---:|---:|---:|
| $H_0$ / km s⁻¹ Mpc⁻¹ | 68.10251 | 0.24702 | −0.06130 |
| $\Omega_m$ | 0.304579 | 0.003360 | +0.000852 |
| $q(0)$ | −0.542976 | 0.005040 | +0.001278 |
| $q(0.5)$ | −0.105017 | 0.005726 | +0.001455 |
| $q(1)$ | +0.167280 | 0.004108 | +0.001045 |

These are posterior summaries conditional on the unchanged flat-ΛCDM model, priors, released supernova treatment and CMB likelihoods. They retain present acceleration and earlier deceleration within that model. They do not test all dark-energy or luminosity-evolution alternatives. A zero observed tail count is not a calibrated exclusion probability.

All 2,000 stratified native-corrected parent points were retained. The new untrimmed weights are the parent's native/proposal weights plus the new-minus-old full BAO log likelihood. Their effective sample size is **1,606.5**, Pareto $k=0.206$, and each independent chain contributes 24.7–25.2% of the weight. All unchanged parent and child overlap/stability gates pass. The old BAO density reconstructed from the new background calls agrees with the saved native density within **$2.60\times10^{-12}$** in log likelihood; density accounting closes within **$5.60\times10^{-14}$**. There were 2,000 background calculations and no additional CMB-spectrum calculations.

The preregistered 32 corners of the five printed-number rounding intervals also pass separately. Across these finite scenarios, the $H_0$ means span **68.09997–68.10494**, and the $\Omega_m$ means span **0.304545–0.304614**. These corners have no assigned probabilities. They are not a bound on all interior cases, unreported non-Gaussian tails or astrophysical systematics.

## Evolving dark energy

The same replacement in the constant-brightness CPL model also preserves the accelerating inference:

| Quantity | Conditional mean ± standard deviation | Equal-tail 95% interval |
|---|---:|---:|
| $H_0$ / km s⁻¹ Mpc⁻¹ | 67.408 ± 0.561 | [66.368, 68.462] |
| $\Omega_m$ | 0.31346 ± 0.00546 | [0.30302, 0.32397] |
| $w_0$ | −0.8204 ± 0.0553 | [−0.9260, −0.7155] |
| $w_a$ | −0.6640 ± 0.2083 | [−1.0904, −0.2778] |
| $q_0$ | −0.3450 ± 0.0613 | [−0.4612, −0.2301] |
| $j_0$ | −0.1281 ± 0.3138 | [−0.7222, +0.5170] |

Relative to the original Lyα compression, the mean $H_0$ moves by −0.0526 km s⁻¹ Mpc⁻¹, $w_a$ by +0.0582 and $q_0$ by −0.0140. These are paired changes between correlated analyses, not independent-measurement discrepancies. The present-acceleration interval remains negative. Both signs of jerk remain allowed; the sampled fraction with negative jerk decreases from 75.6% to **66.0%**, with block Monte Carlo errors of **0.014–0.015**. This sensitivity does not establish the sign of the change in acceleration.

All 33 printed-number scenarios separately pass the unchanged importance and stability gates. The nominal raw effective sample size is **1615.8**, Pareto $k=0.0641$, and the four chains contribute **24.0–25.5%** of the weight. The old BAO density closes to **$1.94\times10^{-12}$** and total density accounting to **$5.42\times10^{-14}$**. Across the rounding scenarios, the $H_0$ means span **67.40607–67.40936**, $q_0$ spans **−0.34543 to −0.34453**, and $w_a$ spans **−0.66599 to −0.66191**. These are finite rounding sensitivities, not an uncertainty distribution. The calculation uses all 2,000 original slots, 2,000 background calculations and zero new CMB spectra. [CPL result with all scenarios and provenance](../results/inference/lya-fullshape-cpl.json).

Both model comparisons remain conditional on the original numerical CMB target. The ΛCDM native-accuracy assessment has detected genuine likelihood variation requiring further work; the CPL assessment is separate. The likelihood replacement and its successful weight diagnostics do not resolve that numerical issue.

## Compression and validation limits

The [validation paper](https://arxiv.org/html/2607.27411v2) tests metal absorption, high-column-density systems, UV response, continuum distortion and small-scale modeling. Those assumptions and nuisance priors remain embedded in the published compression; this replacement does not independently refit them. The growth-rate estimate was excluded by the authors because of mock bias. A Gaussian likelihood for the original correlation-function measurements does not by itself prove Gaussian tails for the compressed distances. The [source inventory](../results/external_probes/lya-fullshape-sources.json) did not locate an author machine likelihood in the checked public interfaces; it does not establish absence everywhere. Our CMB configuration also differs from the paper's, so no reproduction of its reported cosmological significance is claimed.

The [complete result](../results/inference/lya-fullshape-lcdm.json) contains all 33 separate target identities, parameter covariances, gate diagnostics, source hashes and cache seals. Independent [kernel controls](../results/inference/lya-fullshape-validation.json) cover dense quadratics, row ordering and the installed BAO interface. The [consumer controls](../results/inference/lya-fullshape-bridge-validation.json) cover Gaussian reweighting, preserved non-BAO factors, failed parents, altered inputs and sealed-cache replay.

Reproduction requires the qualified, hash-identical parent and acquired likelihood assets. From the repository root:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1
TASK_PYTHON=.work/unified-cosmology/external-probes/.modern-venv/bin/python
$TASK_PYTHON studies/unified_cosmology/code/inference/lya_fullshape_validate.py \
  --output .work/unified-cosmology/lya-fullshape/kernel-recheck.json
$TASK_PYTHON studies/unified_cosmology/code/inference/lya_fullshape_bridge_validate.py \
  --output .work/unified-cosmology/lya-fullshape/consumer-recheck.json
$TASK_PYTHON studies/unified_cosmology/code/inference/lya_fullshape_bridge.py \
  --chain-folder .work/unified-cosmology/inference/modern-lcdm-none-dovekie-official_planck-independence-seed273240-7316f91342ee \
  --correction-summary studies/unified_cosmology/results/inference/modern-lcdm-independence-exact-correction.json \
  --cache .work/unified-cosmology/lya-fullshape/lcdm-v1 \
  --output .work/unified-cosmology/lya-fullshape/lcdm-recheck.json
```

Existing cache records are verified before reuse. The [frozen design](../code/inference/lya-fullshape-design.json) and [acquisition code](../code/external_probes/lya_fullshape_sources.py) specify the approximation and source recovery; generated cache records remain outside Git.
