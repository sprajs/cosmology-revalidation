# Background-model sensitivity of the calibrated supernova check

The tested neutrino and radiation conventions make very small changes to the calibrated supernova likelihood near its flat-ΛCDM solution. They do not explain the roughly 8% difference between our deterministic uncertainty and the paper's rounded uncertainty. This is a bounded background sensitivity check, not a replacement posterior or an exact historical reconstruction.

The [numerical record](../results/distance_ladder/calibration-background-review.json) compares nine fixed points, $\Omega_m=\{0.30,0.33245,0.36\}$ and $H_0=\{72,73.55,75\}$ km s⁻¹ Mpc⁻¹. Every case uses the same 1,657-row calibrated supernova likelihood, full covariance, Cepheid means and flat-$M$ integration. The baryon fraction is fixed to 0.048 across cases to isolate background effects; these points are not the complete modern CMB parameter model.

## The historical parameter names need interpretation

The released chain header declares `omnuh2=0.00083`, `massive_nu=3` and `massless_nu=0.046`. These labels do not by themselves establish the configuration actually executed.

Both the inspected [pre-chain CosmoSIS standard-library snapshot](https://github.com/joezuntz/cosmosis-standard-library/blob/5faec20ff633942aca2b21116d04aaf1f68b88e0/utility/consistency/consistency.py) and the [current snapshot](https://github.com/joezuntz/cosmosis-standard-library/blob/9e3dd611fc20bf39489d43bf4afeb503d26b4e79/utility/consistency/consistency.py) define $\Omega_m=\Omega_b+\Omega_c+\Omega_\nu$. The pre-chain conversion sets $\sum m_\nu=93.14\,\omega_\nu$ eV. However, its [CAMB interface](https://github.com/joezuntz/cosmosis-standard-library/blob/5faec20ff633942aca2b21116d04aaf1f68b88e0/boltzmann/camb/camb_interface.py) accepts `num_massive_neutrinos`, does not read the legacy `massive_nu` input, and explicitly ignores `massless_nu`. [CAMB 1.3.5](https://github.com/cmbant/CAMB/blob/1.3.5/camb/model.py) supplies one massive species by default, with $N_{\rm eff}=3.046$ in that version's constants. The chain does not identify its installed source revisions, so these snapshots establish a possible mapping rather than proving the historical runtime.

We therefore compare three explicitly distinct cases using the current pinned CAMB installation:

1. **Literal three-species interpretation:** $\omega_\nu=0.00083$, three degenerate massive species with total degeneracy 3, and 0.046 massless effective species. The detailed CAMB fields are set explicitly.
2. **Modern one-species assumptions:** $\sum m_\nu=0.06$ eV, one massive species, $N_{\rm eff}=3.044$. CDM is adjusted so the requested total $\Omega_m$ remains fixed.
3. **Contemporaneous-interface candidate:** the pre-chain conversion to $\sum m_\nu=0.0773062$ eV, one massive species and declared $N_{\rm eff}=3.046$. CDM retains the original 0.00083 subtraction. Current CAMB maps that mass to $\omega_\nu=0.0008312791$, shifting the actual density sum by $2.36\times10^{-6}$ in $\Omega_m$ at the central point. This mismatch is recorded, not silently corrected.

All cases set $T_{\rm CMB}=2.7255$ K, zero curvature and $w=-1$. The current solver's effective massive/massless splitting is recorded in every row; using historical parameter defaults does not reproduce historical solver internals.

## Quantitative comparison

Relative to the declared matter-plus-Λ calculation without radiation:

| Background case | Maximum fractional distance change | Maximum distance-modulus change | Range of Δlog likelihood over the nine points |
|---|---:|---:|---:|
| Literal three species | 0.00010342 | 0.00022458 mag | −0.00742 to +0.00935 |
| Modern one species | 0.00014720 | 0.00031967 mag | −0.01064 to +0.01340 |
| Interface candidate | 0.00014941 | 0.00032446 mag | −0.01085 to +0.01370 |

After projecting the shared magnitude mode out with the full covariance, the largest model displacement has Mahalanobis norm 0.00723. Thus the tested background differences are small compared with the released supernova noise in this region.

A quadratic fitted to the nine likelihood changes gives local relative changes in the $H_0$ standard deviation of $4.07\times10^{-6}$, $5.89\times10^{-6}$ and $6.08\times10^{-6}$, respectively—about **0.0004–0.0006%**. The local inferred $H_0$ shifts are at most 0.000106 km s⁻¹ Mpc⁻¹. The fit residual is at most $2.07\times10^{-5}$ in log likelihood. These are curvature diagnostics using the already computed posterior covariance, **not newly integrated posterior uncertainties or bounds on unsampled tails**. They supply no indication that the late-time background approximation causes a material uncertainty change near the solution. The separately inspected released-chain startup sensitivity is a much larger effect on the quoted sample standard deviation.

## Numerical checks and reproduction

The successful calculation used 33 background calls and no spectra or sampling. At the central point for each case, independent adaptive integration of CAMB's $H(z)$ recovers its angular-diameter distances to relative error at most $2.38\times10^{-11}$. Doubling `AccuracyBoost`, or computing the thermal background instead of the background-only mode, leaves the tested distances and likelihood values unchanged at stored precision. Spectral entry points are explicitly disabled.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/distance_ladder/calibration_background_review.py \
  --acquire --quadrature .work/unified-cosmology/calibration-lcdm-v1.json
```

The [design](../code/distance_ladder/calibration-background-design.json) fixes the grid and diagnostics; the [source manifest](../code/distance_ladder/calibration-background-sources.json) pins the inspected upstream code. The result also binds the actual CAMB Python files and native library. No existing target was modified. If a fully integrated background sensitivity becomes necessary, a small separately validated interpolation of the distance difference over $\Omega_m,H_0$ would support deterministic quadrature without a large CAMB calculation.
