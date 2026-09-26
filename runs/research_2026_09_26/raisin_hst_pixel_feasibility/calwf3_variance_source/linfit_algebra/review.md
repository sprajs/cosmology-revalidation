# CALWF3 3.7.3 IR-ramp slope-variance algebra

This is a synthetic source audit of the pinned `hstcal` commit `6a1147d7bdccb7e2a7b73af276f4597c6420fbb8`, whose version header declares CALWF3 3.7.3. Matching that source version to the public FLT header does not prove the exact archive executable or its runtime branches. No real image pixel or Stage B noise outcome entered this calculation. The seven-read, 50-second cadence and all six powers were frozen in `protocol.json` before the coefficients were computed.

For an unsplit accepted interval, [`fitsamps`](../source-3.7.3/cridcalc.c) calls the optimum-weight `linfit` at lines 916–942, combines intervals at 950–997, and returns `out_err = sigb` if there is only this one interval. `linfit` sets `rdns=21/gain` DN per correlated double sample at lines 1095–1097; uses the observed interval SNR to choose a power in {0, 0.4, 1, 3, 6, 10} at 1128–1159; forms dimensionless `u_i=|(i−3)/3|^p` and weight `W_i=u_i/rdns²` at 1162–1177; and sets `fit_uncert=√(S/denom)` at 1181–1203. At 1206–1218, the returned squared slope error adds read, dark-Poisson, and source-Poisson terms through `terma=(gain·fit_uncert·T)²`, `termb=gain·ddg`, `termc=gain·dy`, then divides the summed square root by `gain·T`. The source also has a denominator clamp and separate 0/1/2-point behavior; none is active in the declared seven-read synthetic cases.

Let `t=(0,50,…,300)` seconds, `T=300`, `Q=Σu_i(t_i−t̄_u)²`, and `h_i=u_i(t_i−t̄_u)/Q`. The fitted slope is `Σh_i y_i`, with `Σh_i=0` and `Σh_i t_i=1`. In electron units (`gain=1` as a unit convention), and for expected nonnegative source-plus-dark rate `r` electrons/s, the source-reported variance is

`Var_code = 21²/Q + r/T`.

For the declared comparison model, independent single reads have variance `σ_read²=21²/2` electrons², because the approximately 21-electron [WFC3 correlated-double-sampling read noise](https://hst-docs.stsci.edu/wfc3ihb/chapter-5-wfc3-detector-characteristics-and-performance/5-7-ir-detector-characteristics-and-performance) is measured from a **difference of two reads**. Subtracting a common zeroth read produces covariance `σ_read²(δ_ij+1−δ_i0−δ_j0)`. Since `Σh=0`, its common term cancels in the fitted slope, leaving `Var_read,actual=(21²/2)Σh_i²`. Independent cumulative Poisson increments give `Cov(P_i,P_j)=r·min(t_i,t_j)`, so `Var_P,actual=r·hᵀmin(t_i,t_j)h`. A separate six-increment formula, `r·Σ_{k=1}^6 50(Σ_{i≥k}h_i)²`, agrees numerically to `4.34×10⁻¹⁹` in its unit-rate coefficient.

| Fixed power | Reported read variance | Actual read variance | Reported/actual read | Reported Poisson coefficient | Actual Poisson coefficient | Reported/actual Poisson |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.006300 | 0.003150 | 2.000 | 0.003333 | 0.003571 | 0.9333 |
| 0.4 | 0.006761 | 0.003189 | 2.120 | 0.003333 | 0.003506 | 0.9508 |
| 1 | 0.007350 | 0.003335 | 2.204 | 0.003333 | 0.003441 | 0.9686 |
| 3 | 0.008628 | 0.003947 | 2.186 | 0.003333 | 0.003357 | 0.9929 |
| 6 | 0.009431 | 0.004553 | 2.071 | 0.003333 | 0.003336 | 0.9993 |
| 10 | 0.009725 | 0.004826 | 2.015 | 0.003333 | 0.003333 | 1.0000 |

Read variances are `(electrons/s)²` for `σ_CDS=21` electrons. Poisson coefficients multiply `r` in electrons/s to yield that same variance unit. The full comparison at any declared nonnegative `r` is the ratio `(reported_read + r·reported_P)/(actual_read + r·actual_P)`; no source rate was chosen from observed data. The code weights are in principle selected by **noisy** SNR, and the implementation may reject reads or split intervals, so the table is conditional on a fixed branch and one seven-point interval. It omits correlated or 1/f read noise, interpixel covariance, flat/dark-reference uncertainty, cosmic-ray selection, persistence, and detector nonlinearities. It neither estimates the actual public FLT ERR calibration nor explains the observed signed-repeat statistic. In particular, it does not justify editing or rescaling any ERR array.

The source files and Git blobs are pinned in `protocol.json` and `../source-3.7.3-manifest.json`; the executable freeze is `execution-freeze.json`. `compute.py` writes the six numeric coefficients to `result.json` and `coefficients.csv`, with a separate matrix-versus-increment covariance closure in every branch.
