# Absolute calibration from the released distance ladder

The released SH0ES linear system gives **H₀ = 73.0428 km s⁻¹ Mpc⁻¹**, with a propagated fit uncertainty of **1.0071 km s⁻¹ Mpc⁻¹**. This reproduces the published baseline **73.04 ± 1.01** fit. The published **±1.04** additionally includes the analysis-variant systematic allowance; we have not manufactured that allowance from the fixed matrix. This is an independent numerical solution of the released, selected and corrected observations, conditional on their physical assumptions—not a new reduction of HST images or an additional independent H₀ prior. The [primary paper](https://arxiv.org/abs/2112.04510v3) and [pinned official release](https://github.com/PantheonPlusSH0ES/DataRelease/tree/c447f0fea703fcd0fff57de5000947b5ca81286b/SH0ES_Data) define those inputs.

## What was recovered and checked

The public release supplies the complete data vector, design matrix and covariance: **3,492 rows, 47 full-rank parameters**. Its covariance download is a Git LFS object; the 133-byte pointer was resolved to the actual 48.8 MB FITS file and its declared SHA-256 was verified. Author reference code, photometry Table 2 and the paper are retained as pinned downloads under `.work`, not copied into the repository. Acquisition is recorded in [acquisition.json](../results/distance_ladder/acquisition.json).

The fitted relation is $y=L^Tq+e$, with the released full covariance for $e$. Cholesky whitening followed by SVD, QR and normal-equation solutions agree to $4.4\times10^{-11}$ in the fitted coordinates. A separate implementation using general LU solves, a scalar Schur complement and direct nuisance profiles agrees independently. The residual is **χ² = 3552.7593 for 3445 degrees of freedom**, or **1.03128 per degree of freedom**. Its conditional upper-tail probability is 0.098; this tests this selected Gaussian model, not every upstream calibration assumption. See [reconstruction](../results/distance_ladder/reconstruction.json) and [independent review](../results/distance_ladder/independent-review.json).

The last coordinate is $5\log_{10}H_0$. With flat measure in the native linear coordinates, its exact transformation is lognormal: median **73.04282**, mean **73.04976**, standard deviation **1.00729**, and 95% equal-tail interval **71.09528–75.04371**. These are not flat-H₀-prior results. The author initialization file is a broad-box reference rather than the final least-squares answer. Its zero-width auxiliary coordinate is algebraically isolated; removing that auxiliary leaves all physical fitted coordinates unchanged. The other author bounds are at least 7.95 marginal standard deviations away.

The native period-slope coordinate is an offset from −3.285; adding this back gives **−3.29852 ± 0.01492**. The metallicity slope is **−0.21658 ± 0.04521**. The paper prints a metallicity uncertainty of 0.046, and some other nuisance uncertainties differ slightly from the released matrix, so agreement of every printed uncertainty is not claimed. The printed paper's Equation 16 reverses the absolute-magnitude sign relative to Equation 3; the actual matrix uses the correct calibrator relation $m_B=\mu_{\rm host}+M_B$. The [visual equation check](../results/distance_ladder/equation-visual-review.json) preserves this as a printed algebra inconsistency, with no effect on the recovered H₀.

The baseline Hubble-flow ordinate incorporates the authors' low-redshift expansion with **q₀ = −0.55**. The resulting H₀ summary is therefore not independent of the adopted expansion history. Selected Cepheid photometry, Wesenheit/dust treatment, metallicities, crowding corrections, geometric anchors, compressed Milky Way constraints and SN standardization remain inherited assumptions. Existing Gaussian constraints are already matrix rows; none were added twice. Covariance was neither rescaled nor diagonalized, and no observations were clipped.

## A usable Cepheid-only calibration factor

Removing **all 354 SN rows** leaves **3,138 Cepheid/external-constraint rows and 45 parameters**. The two now-zero columns, SN absolute magnitude and $5\log_{10}H_0$, are removed. Marginalizing the remaining eight nuisance coordinates yields a **37-dimensional Gaussian host-distance likelihood with its full covariance**. It contains no supernova brightness, Hubble-flow redshift or standalone H₀ constraint. Its mean vector and covariance are in [cepheid-factor.json](../results/distance_ladder/cepheid-factor.json), with a generated NPZ interface under `.work`.

All 37 host labels were recovered uniquely from Table 2 by matching each column's ordered Cepheid periods and metallicities within their printed rounding intervals; fitted distances were not used to choose labels. Individual marginal uncertainties range from **0.03446 to 0.24984 mag**, with median **0.08891 mag**. The Schur-complement identity closes to $1.7\times10^{-15}$; twelve direct nuisance refits agree with the compressed quadratic to $1.3\times10^{-11}$ in Δχ². Relative to the full SN-informed solution, a host's mean can move by as much as **0.17795 mag**. Thus full-ladder fitted host distances cannot be reused unchanged as independent Cepheid measurements.

This factor is conditional on the original Cepheid/anchor model and flat nuisance measure. It supplies a concrete absolute-calibration component for a future consistently standardized SN model. It supplies neither a progenitor-age correction nor a license to reuse the full ladder's supernova information. The original release sets Cepheid-to-SN data covariance exactly to zero; that is a property of this supplied model, not a universal proof of independence of arbitrary reprocessed data.

## Overlap and the remaining integration requirement

A conservative exact-alias match, corroborated by separation ≤1 arcsec and |ΔzHEL|≤0.001, identifies **273 Dovekie events in 301 current Pantheon+ measurement rows**. Of these, **76 events in 95 rows** carry the released SH0ES Hubble-flow flag. No calibrator passes this match. Another 77 alias pairs fail the strict position/redshift check, and five sky/redshift-only pairs lack epoch corroboration; all remain explicit. Repeated light curves retain separate row identities. See [overlap.json](../results/distance_ladder/overlap.json).

These are current catalogue flags, not a certified one-to-one map of the earlier ladder matrix. The current Pantheon+ table contains 77 flagged calibrator measurements and 43 CID strings, while the paper counts 42 physical SNe; names are not automatically physical identities. Corrected brightness also differs between products: the first two SN 2011fe table values are 9.74571 and 9.80286, whereas the ladder ordinates are approximately 9.752 and 9.808. Matching row counts cannot establish calibration equivalence.

Within the released ladder, **all 21,329 calibrator–Hubble-flow covariance entries are nonzero**, with maximum absolute correlation **0.2783**. Those terms are available and retained here. The corresponding transformation to the Dovekie calibration and covariance is not supplied by these products. Actual fusion requires an identified calibrator/event mapping, consistent standardized SN observables and shared calibration-response or cross-covariance terms; no H₀ Gaussian has been multiplied into the current CMB+BAO+SN analysis. This is a precise missing interface, not evidence that unified absolute calibration is impossible.

## Reproduction

Run from the repository root in the existing scientific environment:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
.venv/bin/python studies/unified_cosmology/code/distance_ladder/acquire.py
.venv/bin/python studies/unified_cosmology/code/distance_ladder/solve.py
.venv/bin/python studies/unified_cosmology/code/distance_ladder/cepheid_factor.py
.venv/bin/python studies/unified_cosmology/code/distance_ladder/overlap.py
.venv/bin/python studies/unified_cosmology/code/distance_ladder/independent_review.py --output studies/unified_cosmology/results/distance_ladder/independent-review.json
.venv/bin/python studies/unified_cosmology/code/distance_ladder/independent_factor_review.py --output studies/unified_cosmology/results/distance_ladder/independent-factor-review.json
.venv/bin/python studies/unified_cosmology/code/distance_ladder/validate.py
```

The overlap step uses the pinned current [survey ledger](survey-selection.md) and [Pantheon interface](../results/inference/pantheon-interface.json); restore these with their own acquisition instructions first. Acquisition downloads missing source files and verifies existing ones; it fails on changed bytes instead of overwriting them. The solve and review commands regenerate their designated derived outputs. The visual-check record pins the two paper-page renders; on a fresh checkout, its recorded `pdftoppm` command restores those image fixtures before the final identity audit.

The [baseline design](../code/distance_ladder/design.json) preceded fitting; the [Cepheid compression addendum](../code/distance_ladder/cepheid-design.json) was declared after the baseline result and before fitting the SN-free factor. An [independent orthogonal-projection review](../results/distance_ladder/independent-factor-review.json) confirms both host-label recovery and the complete factor covariance. Compact authored results are retained; downloadable sources and regenerated arrays remain in `.work`.
