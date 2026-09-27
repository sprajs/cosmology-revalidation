# Independent late-time comparison of the supernova releases

Dovekie, Pantheon+ and the DES three-year supernova release give similar **best-fitting constant-w models** when each is combined separately with the same DESI DR2 BAO measurements: `w = −0.908`, `−0.913` and `−0.920`. Allowing time dependence adds only **1.38, 0.24 and 0.28** of improvement in minimum χ² beyond those constant-w fits. These are conditional maximum-likelihood results, not posterior credible intervals or a calibrated detection significance.

The comparison uses the unchanged [late-time geometry implementation](../code/inference/late_geometry.py): flat matter plus dark energy, no radiation over the fitted late-time interval, a free `H₀r_drag` scale, and one independently profiled supernova magnitude offset per alternative. It includes **no CMB, Cepheid absolute anchor, imposed brightness evolution, or early-time `w₀+wₐ` cut**. The same 13 DESI BAO measurements, spanning redshift 0.295–2.33, enter each separate fit. They are never duplicated in one likelihood. This is a sensitivity to supernova release choice within that model, not a replacement for the full CMB analysis.

## Best-fitting parameters

`H₀r_drag` is in km s⁻¹. ΛCDM fixes `w₀=−1,wₐ=0`; wCDM fixes `wₐ=0`; CPL uses `w(z)=w₀+wₐz/(1+z)`.

| Supernova release + DESI BAO | Model | Ωm | H₀r_drag | w₀ | wₐ | Minimum χ² |
|---|---|---:|---:|---:|---:|---:|
| Dovekie | ΛCDM | 0.30605 | 10087.21 | −1 | 0 | 1645.278 |
| Dovekie | wCDM | 0.29713 | 9970.61 | −0.90784 | 0 | 1639.497 |
| Dovekie | CPL | 0.31416 | 9943.39 | −0.84344 | −0.53571 | 1638.120 |
| Pantheon+ | ΛCDM | 0.30421 | 10101.30 | −1 | 0 | 1416.139 |
| Pantheon+ | wCDM | 0.29753 | 9978.72 | −0.91306 | 0 | 1411.525 |
| Pantheon+ | CPL | 0.30518 | 9971.01 | −0.89056 | −0.21477 | 1411.280 |
| DES3YR | ΛCDM | 0.29911 | 10140.98 | −1 | 0 | 27.348 |
| DES3YR | wCDM | 0.29764 | 9992.84 | −0.92024 | 0 | 25.544 |
| DES3YR | CPL | 0.30878 | 9943.90 | −0.86833 | −0.31745 | 25.265 |

Compare model improvements **within a row group**, not the absolute χ² between releases:

| Supernova release | ΛCDM minus wCDM χ² | ΛCDM minus CPL χ² | wCDM minus CPL χ² |
|---|---:|---:|---:|
| Dovekie | 5.781 | 7.158 | 1.377 |
| Pantheon+ | 4.614 | 4.859 | 0.245 |
| DES3YR | 1.804 | 2.083 | 0.279 |

The extra CPL freedom changes the best-fit parameters more than it improves the fit beyond constant w. The data releases, calibration assumptions and compression differ, and their supernova events overlap; all alternatives also share the identical BAO data. Their parameter differences therefore cannot be divided by independently combined uncertainties to produce a valid tension significance.

All nine best fits accelerate at the present epoch in their assumed background model (`q₀<0`); the CPL values are −0.368, −0.428 and −0.400, respectively. This statement about optimum locations is not a probability that acceleration is present and is not a test of unrestricted non-accelerating histories.

## Pantheon+ rows and covariance

The raw release has 1,701 light-curve measurement rows. The declared `zHD>0.01` cut retains **1,590 rows for 1,473 distinct released CID names**. There are 97 repeated objects: 78 with two measurements, 18 with three, and one with four. The resulting 117 extra measurement rows are retained with the **full released covariance**; they are not counted as independent new physical supernovae or averaged without covariance propagation.

The normalized adapter's CID ordering, survey identifiers, `zHD`, `zHEL`, `m_b_corr` values and selected covariance agree exactly with a reconstruction from the pinned raw release. Its tiny original covariance asymmetry, 3×10⁻⁸, is symmetrized before selection as documented by the adapter. The covariance–precision identity closes within 1.4×10⁻¹⁴. Same-CID covariance correlations range from −0.005 to 0.915, with median 0.354; repeated measurements must not be treated as independent scalar error bars.

Ten retained rows have the release's calibrator flag, but their Cepheid distances do not enter this fit. It uses corrected apparent magnitude `m_b_corr` with a free common offset, not `MU_SH0ES` or `CEPH_DIST`. The free BAO scale is `H₀r_drag`; neither `H₀` nor the sound horizon is separately measured in this late-time comparison. See the [pinned Pantheon+ adapter record](../results/inference/pantheon-interface.json).

The DES3YR alternative has **20 released Gaussian bin measurements representing 329 supernovae**, not a recovered 329×329 event covariance. Eighteen bins contain data; two retain the author's 999-mag empty-bin error sentinel and negligible weight. The model is evaluated at each released bin redshift. Its raw χ² is therefore especially unsuitable for comparison with an unbinned sample's χ². See the [DES3YR source and assembly record](../results/survey_selection/des3yr-interface.json).

## Numerical scope and reproduction

Each of the nine fits uses three independently seeded differential-evolution searches, each refined locally, plus eight unrelated local starts: **99 optimization runs in total**. The three global-start results agree within 3.2×10⁻¹² in χ² for every fit, and all reported optima are away from the declared parameter boundaries. Some unrelated local starts terminate at worse minima; their coordinates, objectives and termination status remain in the output. Successful termination alone is not taken as evidence of the best solution.

At every reported optimum, a separate calculation integrates the distance at every original supernova row using adaptive quadrature and explicitly profiles the scalar magnitude offset. BAO predictions are independently integrated as well. The largest discrepancy from the compressed implementation is 5.1×10⁻¹² in a component χ². This tests the numerical compression and normalization; it does not validate the release's physical bias corrections.

The reported objective is the sum of SN and BAO residual quadratics. Fixed Gaussian normalization, covariance determinants and offset-integration constants are omitted. Those terms cancel when comparing models on the same fixed data and covariance. They do **not** make raw χ² comparable as a relative likelihood between different releases, and no Bayesian evidence has been calculated. No posterior sampler was run and no uncertainty or significance is inferred merely from optimizer agreement.

After regenerating the three documented adapters and DESI BAO files, run from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/release_comparison.py
```

The [complete result](../results/inference/release-comparison.json) pins input and source hashes, all starts, parameter bounds, individual SN/BAO contributions and the independent direct checks. The frozen geometry source was checked for changes before and after execution and was not edited.
