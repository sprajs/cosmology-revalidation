# Numerical precision of weighted posterior quantiles

The [quantile helper](../code/inference/quantile_precision.py) measures how much the reported 2.5%, 16%, 50%, 84% and 97.5% marginal quantiles vary under a finite-chain block bootstrap. It is an **auxiliary Monte Carlo precision diagnostic**, not additional astrophysical uncertainty, a new acceptance gate, or evidence that an unseen mode is absent. No cosmological outcome was used to develop or validate it.

Each of four independent chains is divided into 10 or 20 contiguous blocks in its original chronological selection order. Five hundred fixed-seed replicates independently resample blocks within each chain. With the intended 500 selected rows per chain, the blocks contain 50 or 25 rows. The helper preserves every repeated row and its original untrimmed importance weight, then normalizes those weights again within each replicate. Explicit duplicate rows matter because combining their weights would change the existing weighted midpoint-CDF interpolation convention. Equal-sized complete blocks are required; the helper rejects an incompatible chain length rather than dropping rows.

For every marginal and partition, output includes the bootstrap standard deviation, the 2.5–97.5% range of bootstrap quantile estimates, and the mean bootstrap displacement from the original quantile. The four independent-chain quantile estimates are also reported separately, with their weight fractions and weight ESS, to expose tail disagreement. Four chains do not by themselves supply a calibrated combined precision bound. The bootstrap range describes numerical variation of an **estimated quantile**, not another posterior credible interval. Replicates with concentrated weights are retained. Their weight-ESS range is reported rather than clipped or used to tune a new gate.

## Known-target validation and limitation

The validation draws four chains of 500 points from N(0.4, 1.25²), with exact importance weights for a N(0, 1) target. It uses 32 independent cohorts for independent draws and 32 for stationary AR(1) draws with correlation 0.9. Each cohort receives 500 bootstrap replicates for both prescribed partitions. A single-row bootstrap is an additional validation comparator only. The independent-draw reference variance is calculated from the analytic self-normalized importance influence function by quadrature.

| Check | 10 blocks per chain | 20 blocks per chain |
|---|---:|---:|
| Independent draws: mean bootstrap SD / theoretical importance SD, range across five quantiles | 0.925–1.003 | 0.942–1.003 |
| Correlated draws: median-quantile bootstrap SD / observed SD across cohorts | 0.703 | 0.646 |

For correlated draws, the 10-block estimates are 2.31–2.75 times larger than single-row-bootstrap estimates. The block method detects dependence, but these partitions still underestimate the observed finite-chain variation in some quantiles. The 32-cohort comparison itself has simulation uncertainty. Thus the two partitions provide sensitivity evidence; they do **not** guarantee calibrated numerical confidence bounds when dependence is long relative to a block. Longer original chains, additional independent chains, or a separately validated blocking scheme would be needed to resolve such a case. The helper does not select a block size from these results.

Adding a constant of 137 to every logweight changes the bootstrap SD by at most 1.7 × 10⁻¹⁴. A constant marginal has exactly zero bootstrap variation. These are implementation checks, not proof of tail coverage. [Validation record](../results/inference/quantile-precision-validation.json).

## Reproduction and future qualified results

Run the synthetic validation from the repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/quantile_precision.py --validate \
  --output studies/unified_cosmology/results/inference/quantile-precision-validation.json
```

For an actual result, supply `--chain-folder`, `--correction-summary` and `--output` instead of `--validate`. The helper first invokes the qualified measurement consumer, which rechecks parent convergence, untrimmed-weight and chain-stability gates, source/data identities and every native correction record. It then verifies chronological ordering and exact agreement with the original marginal quantiles. The output retains the parent input hashes and the auxiliary helper's own source hash. These checks authorize a conditional numerical diagnostic, not a stronger scientific claim than its parent measurement.
