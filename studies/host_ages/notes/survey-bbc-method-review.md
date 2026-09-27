# Native BBC correction targets and normalization

The released BBC option `4336` can be used for the larger fixed R_V=3.1 simulation, but its literal correction target **subtracts the known simulated grey age shift**. A residual age slope after that correction does not establish that a correction trained to absorb an unknown age effect would fail. These are different questions and require separately labeled comparisons. This review checks the actual downloaded source, not only the option name or the original BBC paper.

Source hashes, equations, numerical checks on the existing native outputs, and the completed larger-run review are in [scaled-bbc-method-review.json](../results/galaxy_validation/scaled-bbc-method-review.json). The final review supports the conditional arithmetic and the distinction between age slopes and redshift-mean shifts; it does not establish a correction for observed supernovae. The primary method references are [Kessler & Scolnic (2017)](https://arxiv.org/abs/1610.04677) and [Kessler, Vincenzi & Armstrong (2023)](https://arxiv.org/abs/2306.05819).

## What the native map estimates

`4336 = 16 + 32 + 64 + 128 + 4096`: it requests the multidimensional map, variance scaling, survey/field subsets, a direct distance-modulus correction, and an additive variance floor. The source additionally enables the MAD option. With the dust flag active, `set_DUST_FLAG_biasCor` overwrites training `SIM_ALPHA` and `SIM_BETA` with the input `p1` and `p2`. The released values used here are 0.15 and 2.87. The dust population's intrinsic colour-luminosity parameter is therefore not treated as the fitted effective colour coefficient.

Writing the native amplitude convention as \(m=-2.5\log_{10}x_0\), the direct map averages

\[
 b_\mu = m-M_{0,\mathrm{default}}+p_1 x_1-p_2 c
         -\mathrm{SIM\_gammaDM}-\mu_{\mathrm{true},z}.
\]

The injected grey term is \(\mathrm{SIM\_gammaDM}=-0.030(a-3)\) mag. It is absent from the stored `SIM_mB` normalization but present in the photons, as verified in the [independent photon audit](../results/galaxy_validation/survey-physics-review.json). The old age-training FITRES satisfies that identity within \(1.50\times10^{-8}\) mag. More than five distinct gamma values explicitly force one gamma-map bin. A continuous age population is supported, but this creates no resolved age dimension.

For the alternative *unknown brightness mismatch* question, a derived training table may set only the BBC-facing `SIM_gammaDM` to zero while retaining its original value in `TRUE_AGE_MAGSHIFT`. The photons, fitted quantities, true cosmological distance, and evaluation tables must stay identical. This retains the injected term in the correction target. It is an explicit change of statistical target, not an exact reproduction of the released pipeline. It also changes the native gamma-centering constant and potentially its variance calibration; it must not be described as a mean-only intervention.

Even when evaluation \(\alpha,\beta,\gamma\) are fitted, a dust map whose training coefficients were overwritten to one \(p_1,p_2\) cell is not continuously rebuilt at those fitted values. Changing the initial `p1,p2` and rerunning BBC would rebuild the map. To isolate nuisance absorption, use the same nominal map and reconstruct a second distance with the nominal fitted nuisance coefficients held fixed; compare that with the fully refitted native output on common retained occurrences.

## Distance and host-step signs

For the present single-slope, constant-host-step configuration,

\[
 g(M)=\gamma\left[\tfrac12-\operatorname{logistic}((\log M-10)/0.001)\right]
       -G_{\mathrm{offset}},
\qquad
 \mu_{\mathrm{out}}=m+\alpha x_1-\beta c-g(M)-b_\mu-\overline{M_0}.
\]

Thus positive fitted \(\gamma\) adds \(+\gamma/2\) to a massive host's standardized distance. Full-map setup sets \(G_{\mathrm{offset}}\) to the midpoint of the training gamma range. It is zero for the nominal table; the existing small age-training table would imply −0.0706725 mag. The earlier one-dimensional BBC branch does not initialize that gamma grid. The common arbitrary magnitude zero point must be removed once, using the same common-support weighting when comparing arms.

Native `MURES` is centered by the fitted redshift-bin intercept, whereas output `MU` is centered by the global intercept. Despite the FITRES header wording, the relevant identity for `uave=1` is

\[
 \mathrm{MU}-\mathrm{MUMODEL}=\mathrm{MURES}+\mathrm{M0DIF}.
\]

Using `MURES` alone to measure redshift drift would remove the effect being tested. Reconstructing `MU` from the original `x0`, rounded fitted coefficients, host step, and correction agrees within \(6.22\times10^{-5}\) mag in all three existing positive-dust 1D outputs. The correct residual identity agrees within \(1.01\times10^{-4}\) mag, consistent with their printed precision. The misleading header identity misses by as much as 0.449 mag.

BBC's map uses truth evaluated at observed redshift:

\[
 \mu_{\mathrm{true},z}=\mathrm{SIM\_DLMAG}
 +5\log_{10}\!\left[\frac{d_L(z_{HD},z_{HD};\theta_{BBC})}
 {d_L(z_{\mathrm{sim}},z_{\mathrm{sim}};\theta_{BBC})}\right].
\]

`SIM_LENSDMU` is separate photon scatter; it is not an extra term in this conversion. Report the primary true-distance residual `MU-SIM_DLMAG` and an observed-redshift-truth sensitivity. Do not silently remove lensing or fit a separate intercept in every redshift bin.

## Variance and support

The scatter-calibration function uses the initial host-step prescription rather than subtracting each object's true age gamma. An age term excluded from the mean-map target can consequently remain in the scatter-calibration residual. A positive additive MUCOV floor replaces the scale path; otherwise the scale multiplies the variance excluding the peculiar-velocity contribution. Native BBC intentionally suppresses adding `biasCorErr_mu` separately. Adding all these pieces again would double-count uncertainty.

The native `.COV` file contains the fitted M0-bin Hesse block, conditional on the realized correction table. It is not transformed by subtracting `M0avg`; using it for distance shape therefore requires a free global normalization. Unfitted bins receive a large placeholder variance. It does not supply correction-training Monte Carlo uncertainty, the correlation of paired arms, or survey-systematic covariance. Retraining resamples or independent training shards are needed for those contributions. Existing native 1D outputs contain 169 nominal and 162 age objects, with respectively four and five uninformative redshift bins. Their full 20-dimensional covariance matrices are finite and positive definite; that does not make the excluded bins informative.

For the larger run, retain a unique short merged CID plus the immutable `(shard, original CID, SIM_LIBID)` occurrence key and `SIM_NGEN_LIBID`. Repeated cadence or host IDs identify dependencies, not duplicate supernova occurrences. Verify identical paired latent coordinates, independent training/evaluation seeds, preselection generated-attempt denominators, classifier and quality cuts, native SNR>60 calibration support, all field/cell exclusions, and common-support fractions by age, redshift, mass, and colour. Present both population-level selected-sample differences and common-occurrence differences. A support-changing comparison is not a pure correction comparison.

The physical population remains a mock-age hypothesis with perfect host-mass measurements, fixed R_V=3.1, a reconstructed classifier, and only Type Ia events. Successful native arithmetic or larger cell counts cannot establish the empirical age effect, classifier contamination control, or a cosmological correction.

## Findings while enlarging the simulation

An independent photon check on the first 15,000-attempt paired training shard verifies all nine checked latent coordinates. Across 2,164 common accepted objects and 147,104 noiseless epochs, the imposed signed magnitude shift agrees within \(3.71\times10^{-6}\) mag. Nominal and age-injected selection retain 2,253 and 2,184 objects, with 109 attempt-level selection differences. The imposed fixed-R_V screen remains nonnegative on the checked 1000–30000 Å grid. This verifies the intervention, not its empirical population plausibility.

At four training shards, BBC executes but retains only 19 nominal and 15–16 age objects. Every retained host has \(\log M>10\). The fitted host-step column is consequently constant and exactly degenerate with the free M0-bin intercepts: the standardized design loses one rank. Native gamma errors reach 23 mag despite a successful exit and positive-definite output covariance. Gamma is not measured here. This particular null direction is a global magnitude offset, so globally centered shape can remain identifiable; a null direction changing relative redshift-bin offsets would instead obstruct shape inference. The declared 30-object age-slope gate is not met in this stage.

The native correction grid also depends on each evaluation field's observed-redshift range. Identical training objects and `p1,p2` do not ensure identical correction geometry across separately selected arms. A fixed-map claim requires a common pre-BBC cohort with identical per-field observed-redshift coordinates and a check of actual grid boundaries and training-cell counts. In the original selected-sample comparison, the difference between refitted and externally frozen nuisance distances on the *same age-arm output* still isolates nuisance refitting, because that output's correction cancels exactly.

## Final independent review

The final calculation contains 240,000 generated attempts per arm. Independent token comparisons verify that all 13,724 age-training rows differ between literal and retained targets only in `SIM_gammaDM`; every other token, including `TRUE_AGE_MAGSHIFT`, is preserved. The common-input alternatives have identical evaluation occurrence IDs and observed redshifts, with matching native map geometry. The high-mass restriction independently removes exactly 20 of 161 pre-BBC common candidates. Its host-step factor is exactly −0.5 at machine precision, so fixing gamma chooses an otherwise arbitrary global zero point within this selected stratum. It does not measure a host-mass step.

Independent pivoted-QR ranks agree with the native design audits. The final unrestricted gamma fits are formally full rank but reach the gamma boundary and lack two-sided mass support; their slopes remain withheld. The fixed-coefficient control has 84 common objects, and the high-mass alpha/beta refit has 82. Their inverse-variance weighted effective sample sizes are 64.17 and 63.17, respectively. Direct native distance reconstruction agrees within (6.15\times10^{-5}) mag, consistent with printed coefficients. Weighted within-bin covariance calculations reproduce every reportable point slope within (2\times10^{-15}) mag/Gyr.

In the high-mass stratum, retained-target training leaves a within-bin slope of −0.02791 mag/Gyr while changing the centered 0.7–0.9 redshift-bin mean from +0.05741 mag under frozen nominal calibration to −0.00027 mag. This is not contradictory: the slope regression removes between-bin information. The high-minus-low redshift shift changes from +0.08276 to −0.00674 mag. The residual slope fraction therefore cannot be interpreted as the surviving fraction of a cosmological bias.

All 600 final bootstrap replicate records were independently hash-checked and their percentile summaries and covariance matrices reconstructed. Eight SHA-selected supported replicates per variant and resampling mode were additionally reconstructed from raw native FITRES rows, without importing the producer's estimator. Separate raw-token checks cover 167,322 cloned rows and confirm that cloning changes only occurrence CID. The all-case age slopes use 139/200 joint and 65/100 training-only replicates in each variant. The exploratory high-minus-low contrast has substantially less support: 51/200 joint replicates for the fixed control and 39/200 for the high-mass refit. Full joint covariance across four contrasts and four bins uses only 34 and 25 complete replicates. These are distributions conditional on native support, not validated confidence coverage. In the high-mass case the retained-minus-frozen high–low change is −0.08949 mag, with supported-replicate percentile range [−0.16775, −0.04574] mag; this post-pattern contrast is not a cosmological significance measurement.

The isolated scan-storage patch preserves the original numerical estimator and scientific gates. Independent byte comparisons confirm unchanged science rows and covariance files for four successful regression cases. Of six pilot failures, only one was a storage overflow: the repair exposes its original −0.3 scatter lower-limit outcome after 121 steps, rather than finding an interior solution. Five genuine zero-MAD failures persist. They and the final bootstrap failures must remain visible when interpreting the results. No material arithmetic or target-identity defect was found in the completed review; the scientific and sparse-support limits remain.
