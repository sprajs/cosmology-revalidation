# DES revalidation: what is checked and what follows

The independent numerical checks pass for the frozen DES flux, selected-sample
prediction and finite calibration-mode calculations. The executable is
[`des_checks.py`](../../validation/des_checks.py); its inputs, output numbers and SHA-256
identities are in [`des.json`](../../validation/reports/des.json). Run it with:

```bash
.venv/bin/python validation/des_checks.py
```

This requires the three default `revalidation-20260926-v1` DES result directories,
the colour-only run `revalidation-20260926-des-colour`, and the independent
Tripp-seed run `revalidation-20260926-des-second-seed`. The colour run uses the
supplied configuration; the second Tripp run changes only the seed to 2026092602.
All original results and scientific workflow code remain unchanged.

## Calibrated flux reconstruction

All 24 full fits and 36 held-out fits are numerically valid. Independently
reconstructed native objectives agree within **4.98 × 10⁻¹⁴** in χ²; recomputed
fit objectives within **2.14 × 10⁻¹⁴**. Rebuilding the full conditional predictive
covariance and evaluating its density with SciPy agrees with the stored
held-out log densities within **1.43 × 10⁻¹⁴**.

The largest local covariance-scaled gradient norm is **1.01 × 10⁻⁵**. The two
optimizer starts differ by at most **7.35 × 10⁻¹¹** in χ²; halving Hessian steps
changes the matrix by at most **3.51 × 10⁻⁶** in relative norm. Refining wavelength
integration from 5 Å to 1 Å changes flux by at most **7.04 × 10⁻⁶ of that object's
peak model flux**; the separate public `sncosmo.bandflux` path differs by at most
**1.16 × 10⁻⁵** on that same normalization.

Compared with the frozen native parameter reference, the largest fitted
differences are **0.001111 mag** in amplitude, **0.02856** in width, **0.0005505**
in colour and **0.09517 day** in peak date. These agree with the research package's
stated tolerances. They are small implementation/fit differences, not a new
cosmological correction. The covariance, mask, calibrated observations and
model assets were supplied: neither detector reduction nor survey selection
has been independently reproduced in this calculation. Laplace predictive
uncertainty remains a local approximation.

## Does width add predictive information beyond colour?

The same 1,063 distinct high-probability DES objects are used: 850 for fitting
and 213 held out, spanning ten fields. All selected 3 × 3 measurement covariance
matrices are positive definite. A separate NumPy implementation recovers every
held-out density, posterior predictive mean and uncertainty from all 8,000
saved draws per run, with maximum discrepancies below **2 × 10⁻¹³**. The compiled
JAX likelihood agrees with NumPy for all nine supplied mean families and both
Gaussian/Student-t noise laws; coefficient gradients also pass finite-difference
checks. This directly exercises the covariance projection implicated in the
historical compiler defect.

| Predictor | Held-out log density sum, nats | RMSE, mag |
|---|---:|---:|
| Width + colour, default seed | 91.53737 | 0.158546 |
| Colour only | 51.43870 | 0.189919 |
| Width + colour, independent seed | 91.51785 | 0.158570 |

Width plus colour improves the summed score by **40.09866 nats** over colour
alone. The paired-object bootstrap gives a 95% interval of **0.11890 to 0.25872
nats per object** for the mean improvement. Resampling whole fields gives
**0.14355 to 0.22671**; every leave-one-field-out mean remains positive. These
intervals resample the fixed held-out predictions; they exclude calibration,
training-sample and model-selection uncertainty.

Both comparison runs and the additional seed pass the stored sampler gates.
Independent classical and rank-normalized/folded split R-hat checks also pass.
The second seed changes the total score by **−0.01952 nat**, consistent with the
approximate batch Monte Carlo uncertainty recorded in the JSON report. Effective
sample sizes larger than the 8,000 retained draws are not evidence of extra
draws: negative lag-one correlations are present and can produce antithetic
sampling gains.

This supports width as an additional predictor **within this frozen selected
sample and model**. The fitted width and colour slopes, −0.10836 ± 0.00662 and
2.4733 ± 0.0783, describe apparent-magnitude prediction. They are not an
independent dust-law measurement. Flexible distance offsets, noisy predictors,
fixed classifier membership and absent selection normalization prevent reading
them as population parameters or cosmology. The previously reported roughly
40.06-nat advantage is consistent with this independently sampled comparison.

## Calibration information and its limits

The 43 discovery and 1,020 validation objects are disjoint. Every supplied
object-level Gram matrix and residual projection is reconstructed from the
saved projected design rows, rather than merely rerunning the posterior formula.
The high-255 minus low-255 redshift contrast is recovered exactly using an
independent CID/redshift join. An augmented least-squares/QR calculation recovers
the posterior mean and covariance; a separate predictive Gaussian calculation
recovers the validation evidence increments.

| Mode/prior family | Combined response mean, mag | Conditional response SD, mag |
|---|---:|---:|
| Twelve released modes | 0.018543 | 0.009417 |
| Released modes + inherited observer prior | 0.050599 | 0.012067 |
| Released modes + isotropic observer prior | 0.049123 | 0.012364 |

For twelve modes, the combined conditional variance is **46.62%** of its prior
value. The corresponding validation log-evidence increments against zero
shared response are **50.0709**, **58.7522** and **58.6816**. These are conditional
finite-Gaussian-model comparisons; alternative observer terms describe residual
structure without identifying its physical source.

The mean response varies with the admitted family. Its uncertainty shrinks
only within the supplied design and prior. Unmodelled grey luminosity evolution
can remain invisible after the per-object amplitude direction is projected out.
The table therefore neither measures an actual calibration correction nor
provides an upper bound on all calibration/population bias. The original
[distance-response analysis](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/docs/research-2026-09-26/shared-distance-response.md)
and [selected-sample comparison](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/docs/research-2026-09-26/conditional-repair.md)
retain the larger identification context.

## Relation to published DES results

The frozen assets belong to the DES-SN5YR research chain associated with the
[original cosmology analysis](https://arxiv.org/abs/2401.02929). This revalidation
checks modern exported fixed objectives and the specified conditional models;
it does not reproduce that paper's complete historical execution or infer its
cosmological parameters. Differences between these selected-sample slopes and
published bias-corrected standardization parameters are not a matched test.

The later [DES-Dovekie reanalysis](https://arxiv.org/abs/2511.07517v3) changes
cross-calibration, retrains SALT3 and repairs a host colour-law approximation.
Its existence and scope were checked against the authors' abstract on
2026-09-26. Its assets have not been substituted here, and its numerical
cosmological results have not been revalidated by these checks. The present
audit finds no new numerical defect in the three DES workflows; the
remaining limitations are measurement provenance, selection, approximation
and physical identification.
