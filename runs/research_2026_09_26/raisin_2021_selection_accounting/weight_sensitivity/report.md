# Weight-rule sensitivity on the fixed 2021 FITRES cohort

The [protocol](protocol.json) and [script](../../../../scripts/research_2026_09_26/raisin_2021_weight_sensitivity.py) were fixed before weighted outcomes (SHA-256 `5d08f3fe…` and `0600d000…`). Parent input hashes and the 25,606-object selected count were checked before calculation. The full 60-point curves, including actual count and Kish effective N for every point, are in [curves.csv](curves.csv); [summary.json](summary.json) and [manifest.json](manifest.json) give ranges and hashes. No new cuts or fits were made.

For each fixed region and each `s=0,.005,…,.295` mag, I weighted the archived NIR residual `DLMAG−SIM_DLMAG` with `(DLMAGERR²+s²)^(-p/2)`. `p=1` implements the Nov 2021 source's inverse inflated **error**; `p=2` implements the later Jan 2022 source's inverse inflated **variance**. The same `s` is used in both rules at each comparison point. The first bin is empty and retained as null.

| zHD region | N | p=1 mean range (mag) | p=2 mean range (mag) | matched-`s` p2−p1 range (mmag) |
|---|---:|---:|---:|---:|
| Global | 25,606 | .042637–.043630 | .041781–.043617 | −.857 to −.013 |
| [0,.2) | 0 | empty | empty | empty |
| [.2,.3) | 3,232 | .041860–.042825 | .040891–.042819 | −.970 to −.007 |
| [.3,.4) | 4,945 | .046713–.047837 | .046726–.048747 | +.013 to +.913 |
| [.4,.5) | 14,499 | .041874–.043269 | .040830–.043249 | −1.044 to −.020 |
| [.5,.6) | 2,381 | .037603–.040453 | .034785–.040424 | −2.818 to −.029 |
| [.6,1] | 549 | .043968–.045869 | .043979–.046955 | +.011 to +1.335 |

Globally at `s=0`, p=1 gives .042637 mag with effective N 24,005 and p=2 gives .041781 mag with effective N 17,467. At `s=.295`, both approach the unweighted mean: .043630 (effective N 25,606) and .043617 (effective N 25,606). The full-grid extrema are descriptive algorithmic ranges, not confidence intervals.

This is conditional on the already archived, already fitted, source-cut cohort. The historical `get_sigint` requires the unavailable imported `cosmo.mu`, and the actual sigma may differ between the 2021 and later executions because their underlying data and branches may differ. These curves therefore neither reconstruct either applied correction nor bound a physical bias. They show precisely how much the two documented weighting algorithms differ if forced to share a specified sigma on these same rows.
