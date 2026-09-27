# The released Pantheon+SH0ES calibration interface

The public Pantheon+SH0ES release supplies an executable, positive-definite joint Gaussian likelihood for calibrated supernova distances. It can therefore be analysed as an **alternative to Dovekie**, with its own released calibration assumptions. Our checks do not certify the original covariance construction, and they do not justify adding our reconstructed Cepheid factor to that likelihood. No cosmological parameters were fitted in this interface audit.

The data and official CosmoSIS implementation are pinned to [DataRelease commit c447f0f](https://github.com/PantheonPlusSH0ES/DataRelease/tree/c447f0fea703fcd0fff57de5000947b5ca81286b). The numerical record is [calibration-interface.json](../results/distance_ladder/calibration-interface.json); the reusable implementation is [calibration_interface.py](../code/distance_ladder/calibration_interface.py).

## What is actually calibrated

The official selection is `zHD > 0.01 OR IS_CALIBRATOR`. It retains **1,657 light-curve measurements: 77 calibrator measurements and 1,580 other measurements**. Ten calibrator rows also pass the redshift cut, but remain calibrators. Repeated measurements are not independent supernovae and are not averaged away.

R22 Table 6 identifies the 42 calibrating supernovae and their 37 hosts. Matching those published names to the release links all 77 rows to our independently identified Cepheid hosts. Explicit aliases are `2005df_ANU → 2005df`, `2008fv_comb → 2008fv`, `N0105 → N105A`, and `N0976 → N976A`. The two `2005df` names also have identical released sky coordinates and peak dates differing by 0.1 day. The host map uses names, not fitted brightnesses or distances. All calibrator host-coordinate fields are missing, so these are documentary host associations rather than newly measured host positions. This does **not** establish an exact measurement-row mapping to the older 77-row ladder matrix.

The released `CEPH_DIST` values closely follow the supernova-free host solution. Relative to our reconstruction they are lower by 0.000598–0.001208 mag, with RMS 0.000934 mag; after a common offset the RMS is 0.000134 mag. In contrast, differences from the full supernova-informed ladder host means have RMS 0.06042 mag and maximum absolute difference 0.17688 mag. R22 Table 6, column (b), explicitly describes distances excluding every supernova. These facts support a Cepheid-only interpretation of the released central distances, but they do not demonstrate byte-for-byte identity with our reconstructed factor. The small offset's exact origin remains unresolved. [R22, Table 6](https://arxiv.org/abs/2112.04510)

## The likelihood and the magnitude convention

Use the observed **`m_b_corr`**, not the already calibrated `MU_SH0ES` column. The theoretical distance modulus is

$$
t_i=\begin{cases}
\mu_i^{\rm Cepheid}, & i\text{ is a calibrator},\\
5\log_{10}\!\left[(1+z_{{\rm HD},i})(1+z_{{\rm HEL},i})D_A(z_{{\rm HD},i})/{\rm Mpc}\right]+25,&\text{otherwise}.
\end{cases}
$$

The predicted apparent magnitude is $t_i+M$. Calibrator redshifts do not determine their predicted distances. The full selected `STAT+SYS` covariance already contains Cepheid information. The official formulation adds the supernova and Cepheid covariance contributions before inversion; it does not add an independent Gaussian constraint on $H_0$. [Brout et al., equations 14–15](https://arxiv.org/html/2202.04077v2#S2.SS3), [official CosmoSIS implementation](https://github.com/PantheonPlusSH0ES/DataRelease/blob/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/5_COSMOLOGY/cosmosis_likelihoods/Pantheon%2BSH0ES_cosmosis_likelihood.py)

For $r=m-t$, $a=\mathbf1^TC^{-1}\mathbf1$, and $\widehat M=\mathbf1^TC^{-1}r/a$, the implementation integrates the single flat magnitude parameter:

$$
\log L=-\frac12\left[(r-\widehat M\mathbf1)^TC^{-1}(r-\widehat M\mathbf1)
+\log|C|+\log a+(N-1)\log(2\pi)\right].
$$

This uses the improper integration measure $dM$ in magnitudes. Its absolute evidence normalization is undefined. The determinant terms are constant for this fixed covariance, but are retained for consistent density accounting. No $H_0$ transformation Jacobian belongs in this magnitude integral; any choice of cosmological parameter prior remains separate.

The covariance has 121,659 nonzero entries among the 121,660 calibrator-to-noncalibrator pairs; their maximum absolute correlation is 0.29361. Splitting this into independent calibrator and cosmological factors would discard released information. The selected symmetric covariance has minimum eigenvalue $8.03857\times10^{-4}\,{\rm mag}^2$.

## Why replacing the Cepheid covariance is not certified

`STATONLY` also contains terms shared between different Cepheid hosts. Its off-diagonal blocks are constant across repeated observations of each host pair. Of 666 distinct host pairs, 630 are nonzero; the 36 zero pairs involve N1365. These terms differ from our reconstructed host covariance by RMS $9.86355\times10^{-5}\,{\rm mag}^2$, with maximum absolute difference $5.04160\times10^{-4}\,{\rm mag}^2$.

Four hosts contain genuinely distinct supernovae, exposing the following entries without confusing them with repeated observations of one explosion:

| Host | Physical supernovae | Released STATONLY shared entry (mag²) | Reconstructed Cepheid-only variance (mag²) | Ratio |
|---|---|---:|---:|---:|
| N1448 | 2001el, 2021pit | 0.00265025 | 0.00134361 | 1.97248 |
| N3147 | 1997bq, 2008fv, 2021hpr | 0.05463299 | 0.02734516 | 1.99790 |
| N5468 | 1999cp, 2002cr | 0.01095859 | 0.00549704 | 1.99354 |
| N5643 | 2013aa, 2017cbv | 0.00538959 | 0.00272228 | 1.97981 |

These are measured cross-release discrepancies, **not proof of a factor-of-two implementation error**. The exact embedding code, earlier host covariance and branch-specific assumptions have not been identified. The independent [primary-source review](calibration-source-review.md) checks the published covariance equations, pipeline configuration and release history without resolving the construction. For the other 33 hosts, the host variance is confounded with the common scatter of their one physical supernova. Off-diagonal patterns alone cannot recover those missing diagonal components. The matrix also contains statistical links to excluded `1992A`/`2007on` rows; none survive as calibrator-to-noncalibrator entries in the selected `STATONLY` block.

Subtracting our candidate host covariance from the released total leaves a positive-definite matrix, with minimum eigenvalue $7.29790\times10^{-4}\,{\rm mag}^2$. That only demonstrates mathematical feasibility; it does not identify the component that was actually added. Accordingly, the interface neither substitutes our host factor nor treats `STATONLY` as a verified supernova-only covariance. A certified replacement requires the original embedded host covariance, its row map and its construction or version history. The [release README](https://github.com/PantheonPlusSH0ES/DataRelease/blob/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/4_DISTANCES_AND_COVAR/README) explicitly includes Cepheid uncertainties; its short description of statistical components does not resolve the observed decomposition.

## Numerical qualification

The printed matrices differ from their transposes by at most $3\times10^{-8}\,{\rm mag}^2$, in 778 ordered entries. A Gaussian covariance must be symmetric. This interface explicitly uses $(C+C^T)/2$, retaining and hashing the original bytes.

Twelve fixed synthetic distance curves were evaluated against unmodified method bodies from the official CosmoSIS code and Cobaya 3.6.2. The mean predictions agree exactly. Fixed-$M$ quadratic values using the original raw matrix inverse differ by at most $1.37\times10^{-10}$; independently profiling that raw quadratic differs by at most $2.94\times10^{-9}$. Cobaya's flat-$M$ method supplied with symmetric covariance agrees within $3.27\times10^{-10}$.

The literal default Cobaya projection applied to the nonsymmetric printed matrix differs by $\Delta\chi^2=0.000300$–0.000679 on these curves: its rank-one formula assumes a symmetric inverse. This small discrepancy is preserved, not relabelled exact agreement. The initial strict symmetry check and the raw-projection check failed before the explicit symmetric convention was adopted; their logs remain in the ignored acquisition directory.

Sixteen independent synthetic Gaussian systems compare direct correlated-host-plus-$M$ integration with the summed-covariance expression, including determinant factors. Maximum discrepancy is $1.07\times10^{-14}$; adding a global 53-mag offset changes the result by at most $1.65\times10^{-12}$. These validate the algebra, not the empirical calibration assumptions.

## Reproduction and reuse

First reproduce the [distance-ladder acquisition and Cepheid-only factor](../code/distance_ladder/README.md). Then, from the repository root:

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/distance_ladder/calibration_interface.py --acquire
```

Public inputs are hash-pinned, and the script verifies them on reuse. Downloads and complete 77-row host, 37-host mean, and 666-pair covariance ledgers go to `.work/unified-cosmology/calibration-interface/`. The compact result binds their hashes, the source, official method bodies and input versions.

`ReleasedCalibration` exposes the ordered noncalibrator redshifts, selected data and covariance. Its `evaluate(D_A)` accepts 1,580 angular-diameter distances in Mpc in that order and returns the integrated log likelihood, projected quadratic, conditional $M$ moments and normalization terms. It performs no new light-curve or cosmological fit. Use it as a separate Pantheon+SH0ES target, retaining its full covariance and its released calibration assumptions. Do not multiply it by Dovekie, an overlapping Pantheon-only likelihood, our Cepheid factor, or a SH0ES-derived $H_0$ prior.
