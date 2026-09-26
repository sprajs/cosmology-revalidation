# Foreground dust structure and shared correction uncertainties

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../../README.md) for current execution status.

21 September 2026. New checks outside the separate SALT/SNM correction work. Original observations, releases and fitted results are preserved. Numerical outputs are under `runs/assumption_audit/`. This audit distinguishes an incompatible release product, a demonstrated sensitivity to covariance structure, and an unpropagated foreground-map alternative.

## Confirmed: Pantheon+ grouped systematics use a different statistical baseline

The supplied `sytematic_groupings/Pantheon+SH0ES_122221_MWEBV.cov` and `MWCOLORLAW.cov` are **not** the final `Pantheon+SH0ES_STATONLY.cov` plus the named dust systematic. This matters when interpreting a difference between the named group and final statistical covariance as the effect of dust alone.

For example, final STATONLY correlates observations 47 and 48 of SN 2007le by **0.05249753 mag²**. The entire MWEBV group matrix has only **0.00001110 mag²** in that entry. The difference cannot be an additional positive covariance contributed by a single dust-scale uncertainty. The release README explicitly says repeated light curves share intrinsic statistical scatter.

We recovered each one-mode dust covariance from an isolated high-redshift reference row and subtracted its outer product. The two independently recovered backgrounds agree to **1.76e−8 mag²**, at the files' rounding scale. Both backgrounds differ from final STATONLY in 17 diagonal and 376 off-diagonal entries above 1e−6 mag²; 209 selected cosmology rows are touched. There are also 12 off-diagonal entries between different literal CIDs, so the discrepancy must not be reduced to a simple duplicate-only repair without checking the released row-pair table. An outer product's sign is undetermined; this reconstruction identifies covariance, not the sign of a dust correction.

On the same 1,590 selected rows, simple flat-LCDM profile fits give:

| Input covariance | Ωm at profile optimum | Local curvature uncertainty | Minimum χ² |
|---|---:|---:|---:|
| Final STATONLY | 0.34784 | 0.01311 | 1484.11 |
| MWEBV group as supplied | 0.34104 | 0.01366 | 1512.52 |
| Recovered MWEBV systematic + final STATONLY | 0.34210 | 0.01393 | 1482.79 |
| MWCOLORLAW group as supplied | 0.34299 | 0.01331 | 1512.74 |
| Recovered MWCOLORLAW systematic + final STATONLY | 0.34332 | 0.01358 | 1482.70 |

These isolate a file-comparison trap. They are not new full-systematics cosmology results. **The project's main Pantheon+ likelihood uses the final total covariance and is unaffected.** This test does not establish which files were used in the authors' original fits, nor does it turn the covariance discrepancy into evidence for an erroneous published cosmological result. It establishes that these public group files must not be compared or subtracted as though they share the final statistical baseline.

An independent review also checked this without any reconstruction: the 47/48 principal submatrix of `MWEBV−STATONLY` has eigenvalue −0.0524751 mag², so it cannot be a positive-semidefinite additional systematic. Four different reference rows recover mutually consistent group backgrounds. See the [independent check](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/assumption_audit/bao/group-peer-review.json). The tabulated minima are quadratic discrepancies under different fixed covariances, not likelihood-ratio significance tests.

Evidence: [numerical record](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/assumption_audit/covariance-groups.json), [affected row pairs](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/assumption_audit/covariance-group-differences.csv), [implementation](../../../light_curve_fitting/code/assumption_audit/covariance_groups.py), and the [primary release](https://github.com/PantheonPlusSH0ES/DataRelease/tree/main/Pantheon%2B_Data/4_DISTANCES_AND_COVAR). Reconstruction is verified using both dust groups, positive-definite likelihood inputs and independently profiled magnitude offsets. The rank-one recovery threshold is recorded in code.

## Measured foreground alternative: spatial structure is not a global dust scale

Pantheon+ §3.2.3 uses SFD reddening scaled by 0.86 with F99 Rv=3.1, a global 5% reddening-scale uncertainty, and an alternative extinction law. These are explicit perturbations, not a general spatial error model. [Pantheon+ methods](https://arxiv.org/abs/2202.04077).

Chiang's CSFD map removes an estimated cosmic infrared background imprint from SFD; that imprint can correlate with extragalactic structure and redshift. It is a distinct mechanism from an overall Galactic reddening normalization. The public maps and reliability mask are available from [NASA LAMBDA](https://lambda.gsfc.nasa.gov/product/foreground/fg_csfd_reddening_map_info.html), with the reconstruction and limitations described in [Chiang 2023](https://arxiv.org/abs/2306.03926).

We downloaded the actual NSIDE=2048 RING maps, transformed the released SN coordinates to Galactic coordinates, and sampled both maps with the same bilinear interpolation. Differences are `0.86*(CSFD−SFD)` in **E(B−V)** units. Of 1,590 cosmology rows, 1,579 lie inside the correction footprint and 1,304 have their nearest pixel in the cosmology reliability mask. The stricter reliable subset also requires all interpolation pixels to be reliable and contains **1,297 rows**. No map value is silently treated as reliable outside this mask.

For all selected rows, the map difference has mean −0.000303 and standard deviation 0.001024 mag in E(B−V). Fitting a constant and one global MWEBV scaling leaves **99.1% of the centered variation unexplained**, or **98.36% within the reliable subset**. Thus the spatial difference cannot be represented by simply enlarging that one scale parameter. In the reliable subset, the row-weighted mean shift is about −0.000777 at 0.01≤z<0.05, −0.0000645 at 0.1≤z<0.3, and +0.0000425 at 0.3≤z<0.6. This is an observed map-pattern/survey-footprint difference, not a detection of SN luminosity evolution. Repeated SNe and nearby fields are correlated; no independent-row significance is attached. An [independent read-only review](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/assumption_audit/foreground-peer-review.json) reproduced the units, mask, interpolation and these summaries.

**These are not distance-modulus shifts.** Converting them into standardized SN distances requires observer-frame bandpasses, redshift, colour fitting, training, and selection response. Multiplying by a fixed B-band coefficient and applying another correction to released `m_b_corr` would skip those responses. Neither the present map calculation nor a fitted global scale establishes the missing covariance's size in cosmology. The reliable CSFD region also retains reconstruction uncertainty and some small-scale residual contamination. The concrete next experiment is a coherent foreground-map swap applied to the shared host and SN photometry, propagated through the separately maintained fitting pipeline, with map uncertainty and survey selection retained.

Evidence: [sampled rows](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/assumption_audit/foreground-map-samples.csv), [summary and masks](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/assumption_audit/foreground-map-summary.json), [input hashes and URLs](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/assumption_audit/foreground-map-manifest.json), [code](../../code/assumption_audit/map_dust.py). Nearest-pixel versus bilinear differences are recorded as an interpolation check; their RMS is 0.0000928 mag in E(B−V) for the full selected sample.

## Quantified: a shared age slope cannot be represented by independent row errors

For an age template `t` and one normal amplitude `a~N(1,s²)`, marginalizing that amplitude gives `C_new=C+s²ttᵀ`. Replacing the outer product with `diag(s²t²)` describes 1,590 different independent nuisance parameters. Those are different models even though each diagonal element is identical.

The current main analysis already uses a shared slope correctly. We integrated it analytically and compared with the original sampled posterior, then reweighted that marginal posterior to two counterfactuals. The original 0.030±0.004 mag/Gyr prior and frozen C14 median template remain fixed throughout:

| Covariance treatment | q0 mean ± SD | P(q0<0) |
|---|---:|---:|
| Fixed slope | +0.0625 ± 0.0742 | 0.198 |
| One shared uncertain slope | +0.0419 ± 0.0967 | 0.330 |
| Incorrect independent row uncertainties | +0.0625 ± 0.0743 | 0.198 |

The diagonal approximation loses almost all of this uncertainty's cosmological effect. Under the shared model the slope amplitude and q0 have correlation **0.636**. Weight ESS exceeds 65,900 out of 80,000 source draws in the counterfactuals; this assesses importance-weight concentration, not independent-chain convergence. The fixed-slope result also agrees with the separately run original fixed-slope chain.

The analytic projected covariance agrees with direct inversion to 1.12e−11. Explicit elimination of the Gaussian nuisance agrees with its marginalized quadratic form to 5.65e−9 in χ²; full versus interpolated likelihood differs by at most 7.20e−7 on tested draws. A determinant/normalization factor is constant here because the template is fixed. It generally would not remain constant if the template or its covariance changed with cosmology.

These constant normalization factors may be dropped for each normalized parameter posterior, but not for comparing absolute likelihoods or evidences between different covariance models. The output's `mode_chisq` fields are not evidence ratios. The original amplitude truncation is 22.5 prior standard deviations from the mean on both sides; the analytic untruncated normal has negligible missing tail mass here.

This validates the treatment of **one declared** uncertainty. It does not marginalize SFH/DTD choices, clock parameters, age–dust–metallicity posterior dependence, shared calibration, or selected-sample transport. Nor does a Gaussian prior derived from overlapping historical SNe automatically provide independent information when applied to Pantheon+. That dependence would require the joint data/model, not just a wider scalar prior.

Evidence: [numerical output](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/assumption_audit/shared-uncertainty.json), [manifest](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/assumption_audit/shared-uncertainty-manifest.json), [code](../../../light_curve_fitting/code/assumption_audit/shared_uncertainty.py).

## Reproduction

From the project root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/assumption_audit/covariance_groups.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/assumption_audit/shared_uncertainty.py
uv pip install --python .venv/bin/python --no-deps --target /tmp/supernova-assumption-deps astropy-healpix==1.1.2
PYTHONPATH=/tmp/supernova-assumption-deps OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/assumption_audit/map_dust.py
```

The last command requires `sfd_ebv.fits`, `csfd_ebv.fits`, `mask.fits` and `readme.txt` from `https://lambda.gsfc.nasa.gov/data/foregrounds/CSFD/` under `data/dust/csfd-v2/`. These total roughly 0.5 GB. The extra HEALPix package is installed separately; the project environment and lockfile are unchanged. Each manifest hashes inputs and outputs; the original fits and covariance files are never overwritten.
