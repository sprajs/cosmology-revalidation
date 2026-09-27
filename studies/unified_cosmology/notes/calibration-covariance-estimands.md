# What the sibling covariance comparison tests

The near factor-of-two pattern is real, and comparing it with a covariance conditioned on fixed Cepheid parameters does not explain it. Our reconstructed 37-host factor already marginalizes all eight remaining Cepheid/anchor nuisance coordinates. Its four relevant errors also agree with the explicitly SN-free errors in Riess et al. Table 6(b). However, the released **SN residual covariance** and a **host-distance covariance** are different objects: identifying an erroneous duplicate requires establishing every term connecting them. The original covariance has not been changed.

The [independent calculation](../results/distance_ladder/covariance-estimands.json) reloads the original 3492-row, 47-parameter system, removes every row containing the SN absolute-magnitude coordinate, and verifies that all 354 SN rows are absent. The remaining 3138-row, 45-parameter system has zero covariance with the removed rows. A general-LU solve of its full covariance followed by Schur elimination of the eight nuisance coordinates reproduces the existing 37-host mean to $1.43\times10^{-13}$ mag and covariance to $1.22\times10^{-17}$ mag². This is separate from the producer's whitened least squares and the [earlier independent QR/SVD review](../results/distance_ladder/independent-factor-review.json).

Let $\widehat{\boldsymbol\mu}$ be the Cepheid-only host-distance estimate with covariance $H$, and let $J$ assign a light-curve row to its physical host. For independent SN and Cepheid measurement errors, the covariance of the calibrator residual vector $\boldsymbol r=\boldsymbol m-J\widehat{\boldsymbol\mu}-M\boldsymbol1$ is

$$
S=C_{\rm SN}+JHJ^T.
$$

For two distinct SNe in the same host, the Cepheid contribution is **one** $H_{hh}$. Repeating the host estimate for a sibling does not create a second independent host measurement. Shared anchor, Cepheid zero-point, period-slope and metallicity uncertainty is already contained in the marginalized $H$; it must not be added again under a second name. If SN and Cepheid errors share calibration, with $\Gamma=\operatorname{Cov}(\boldsymbol m,\widehat{\boldsymbol\mu})$, the more general expression is $S=C_{\rm SN}+JHJ^T-\Gamma J^T-J\Gamma^T$. Those cross terms require their own measured response/covariance.

The [published-table check](../results/distance_ladder/calibration-table6-review.json) and [source review](calibration-source-review.md) establish the following comparison:

| Host | SN-free host SD (mag) | Released sibling covariance / host variance | SD of an additional independent host-shared term that would explain the excess (mag) |
|---|---:|---:|---:|
| NGC 1448 | 0.036655 | 1.97248 | 0.036147 |
| NGC 3147 | 0.165364 | 1.99790 | 0.165190 |
| NGC 5468 | 0.074142 | 1.99354 | 0.073902 |
| NGC 5643 | 0.052175 | 1.97981 | 0.051646 |

The last column is a **conditional algebraic attribution**, not a measurement of new intrinsic SN scatter. An independent host-shared SN or measurement term with covariance $JDJ^T$, for diagonal positive $D$, would increase same-host entries and leave different-host entries unchanged. That is a valid counterexample to the claim that the pattern alone proves duplication. Its required size happens to track the Cepheid uncertainty closely, and the inspected paper prescription does not establish such an added sibling term. A different host-distance input with larger host-specific noise could also have similar off-host terms. Alternatively, adding an already included host-diagonal contribution a second time would make the same pattern. Without the actual construction inputs these possibilities remain distinguishable hypotheses, not identified causes.

This is not evidence for two copies of the **whole** matrix $H$: that would also double different-host entries. Their scale is only approximately that of our factor. Across the 630 pairs excluding NGC 1365, the median released-to-reconstructed ratio is 0.94583; a descriptive least-squares multiplier through zero is 0.93843. The RMS difference is $6.99\times10^{-5}$ mag² for one $H$, compared with $4.52\times10^{-4}$ mag² for two. All 36 NGC 1365 cross-host entries vanish in STATONLY while our factor has nonzero entries. These cross-release differences independently rule out assuming byte-identical covariance ancestry.

Fitting or marginalizing the common SN absolute magnitude also does not explain a positive host-specific increment in the raw matrix. With $a=\boldsymbol1^TS^{-1}\boldsymbol1$ and the GLS estimate $\widehat M$, the fitted-residual covariance is

$$
\operatorname{Cov}(\boldsymbol r-\boldsymbol1\widehat M)
=S-\frac{\boldsymbol1\boldsymbol1^T}{a},
$$

whereas integration over a flat $dM$ measure gives the singular projected precision

$$
P=S^{-1}-\frac{S^{-1}\boldsymbol1\boldsymbol1^TS^{-1}}{a}.
$$

These are not interchangeable with raw $S$ or with the host covariance $H$. For the actual selected 1657-row total covariance, the fitted-residual subtraction is the same $1.3910\times10^{-5}$ mag² for every entry. A common Cepheid zero point affects the calibrator subset, so it is not generally removed by fitting an offset common to all calibrator and Hubble-flow SNe. Using a full-ladder host posterior instead of the SN-free factor would be another change of estimand, and would reuse SN information.

Sibling magnitude differences cancel $JH J^T$ and any additional purely host-common $JDJ^T$. They therefore cannot determine the size or provenance of this common term by themselves. The executable audit verifies this cancellation, the positive-semidefinite extra-host counterexample, and the equivalence of joint host-plus-$M$ integration to the embedded-covariance calculation, with maximum numerical discrepancy $3.56\times10^{-15}$.

The bounded conclusion is a localized, independently verified covariance discrepancy relative to the reconstructed and printed SN-free host uncertainties. It is not yet a demonstrated double count, and supplies no justified covariance subtraction, revised $H_0$, or new cosmological result. The decisive evidence is the exact pre-embedding host covariance, SN-only covariance (including any sibling term), host/SN cross-calibration covariance, and the code assigning those components to light-curve rows. The [public-input search](calibration-construction-search.md) has not recovered this complete construction.

Reproduce the audit after restoring the pinned ladder and calibration-interface inputs:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python studies/unified_cosmology/code/distance_ladder/covariance_estimands.py
```

All producer inputs and prior review hashes are checked before use and again after calculation. No cosmological model, new data acquisition, or released-covariance mutation is involved.
