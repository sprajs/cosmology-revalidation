# Fixed-cohort WFC3/IR dark RAW repeat variability

The root-released [frozen executor](dark_executor.py) completed all **seven primary disjoint pairs and the one prespecified NORMAL search sensitivity pair** in 4.28 seconds. The fixed union BPIXTAB mask and 247/256 spatial sites were unchanged. The maximum direct-versus-matrix contraction gap was `1.11×10⁻¹⁶`. The saved [Γ and mean arrays](read-matrix.npz), [all powers and both mask scopes](pair-quadrant-power.csv), [site contrasts](aperture-pair-power.csv), and [arithmetic gates](arithmetic-gates.json) preserve the outcome. This is total digitized RAW repeat variability in recorded nominal SPARS50 time coordinates; it is **not** isolated electronic read noise, a calibrated FLT error measurement, or an error rescaling.

The table reports the declared power-0 masked statistic `hᵀΓh = mean_pixels[(hᵀ(y_later−y_earlier))²]/2`, in DN² per nominal second². Its quadrant pixel counts are fixed across pairs: B 247,091; C 250,250; A 247,448; D 248,961. All rows and extra weighting powers are retained in the CSV.

| Pair (earlier→later) | B | C | A | D |
|---|---:|---:|---:|---:|
| Search 1 `idbx38naq→idbx39nsq` | 0.07158 | 0.08697 | 0.52084 | 0.05243 |
| Search 2 `idbx40o7q→idbx41onq` | 0.06721 | 0.12753 | 0.10610 | 0.03895 |
| Search 3 `idbx42opq→idbx43p7q` | 0.02426 | 0.21296 | 0.03974 | 0.03109 |
| Template 1 `idp240p7q→idp241psq` | 0.06669 | 0.17652 | 0.04726 | 0.10504 |
| Template 2 `idp242q0q→idp243slq` | 0.04187 | 0.03441 | 0.04500 | 0.06463 |
| Template 3 `idp244t1q→idp245t7q` | 0.02742 | 0.05535 | 0.02792 | 0.03551 |
| Template 4 `idp246tlq→idp247tnq` | 0.08175 | 0.10359 | 0.04166 | 0.03949 |
| **NORMAL search sensitivity** `idbx41onq→idbx43p7q` | **0.02319** | **0.05650** | **0.04595** | **0.02074** |

The NORMAL pair overlaps the primary search observations and has a different separation; it is not an independent replication or a controlled TDF comparison. The three NORMAL roots provide one search pair and no template pair. Eleven of fourteen fixed inputs carry `EXPFLAG=INDETERMINATE` and aggregate TDF-down telemetry warnings. The nominal seven read times are a coordinate convention, not proof of the physical read times in those records. Thus the variation among pairs or quadrants cannot be assigned to electronic noise, TDF condition, cosmic rays, changing dark current, or another single mechanism.

For every pair and quadrant, the [post-score decomposition](decomposition.csv) reprojects the saved seven-by-seven uncentered Γ and seven-element repeat-mean vector μ with the *unchanged* power-0 h:

`hᵀΓh = ½(hᵀμ)² + ½ mean_pixels[(hᵀδy − hᵀμ)²]`.

The recomputed totals agree with the score CSV within `1.11×10⁻¹⁶`. The squared coherent repeat-mean part is only **0.00030%–0.661%** of the total across all 32 pair/quadrant rows; the larger centered spatial part still contains any unflagged cosmic-ray steps, hot-pixel changes, shot noise, digitizer effects, and detector-state differences. It is not a pure white-read variance. [descriptive-summary.json](descriptive-summary.json) records all four NORMAL decomposition rows and full ranges.

The [original-flag/digital ledger](raw-dq-digital-ledger.csv) retains every read's RAW DQ and digital endpoint counts. Among the 993,750 fixed-mask good active pixels, **zero pixels in every pair** had an input endpoint value (0 or 65535) or nonzero RAW DQ in the seven nonzero reads; hence their separately saved [power-0 contribution fractions](raw-dq-digital-contributions.csv) are zero. At least one raw zero-DN pixel can appear outside the fixed mask, so this is a conditional statement about the scored region. RAW DQ is a constant zero HDU in these first-seven reads, and the absence of a flag does not exclude unflagged events or establish linear detector behavior.

The result is linked to source and archive processing only to this extent: CALWF3's `dqicorr.c`/`dodqi.c` coordinate rules fixed the external BPIXTAB mask; `noiscalc.c` and the CCDTAB establish the physical amplifier gain secondary. The separate [private native replay report](../../../../docs/research-2026-09-26/raisin-calwf3-native-replay.md) shows exact SCI/ERR/SAMP/TIME equality for two *current archive science FLTs*, with archive-only DQ bit 4096 differences. It does not turn these 2016–2017 dark RAW pair moments into the CALWF3 slope-error covariance or reproduce the historical RAISIN aperture catalog. No detector-noise correction or cosmological quantity is inferred.

The saved score release binds [execution-freeze.json](execution-freeze.json) SHA256 `1e0ed399d65a53dd5857a945693fecd20700cc86715b576a5b229d062e5bddd1`. The [post-score report protocol](report-protocol.json) specified the descriptive decomposition after primary scoring and before computing that derived table; it is not a new preregistered test. `report_results.py` reads only saved arrays and CSVs; it does not reopen RAW images. All intermediate and final primary outputs remain unchanged.
