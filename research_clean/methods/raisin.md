# RAISIN distances, signed photometry and timing

## Paired released-distance accounting

`raisin` compares the same 79 physical SNe across optical, NIR and combined branches: 42 at `zHD<0.1` and 37 at `zHD>0.2`. Identity, ordering and redshifts must match before a contrast is calculated.

The primary descriptive contrast is the high-redshift mean of `(mu_optical-mu_NIR)` minus its low-redshift mean. Pairing cancels the common distance–redshift relation; the high-minus-low contrast cancels a branch-constant offset.

The variants use exact exported signs:

| Variant | Distance expression |
|---|---|
| Released | `DLMAG` |
| Without mass term | `DLMAG-MASS_CORR` |
| Without bias term | `DLMAG+DLMAG_biascor` |
| Without either term | `DLMAG+DLMAG_biascor-MASS_CORR` |

Expect approximately `+0.001236 mag` for released optical−NIR and `+0.07541 mag` after undoing both exported terms. This is accounting at fixed fitted photometry and selection. It is not a newly validated correction.

The exported covariance's systematic diagonal is restored from the difference between group-specific and statistical `lcparams` variance. Unknown cross-branch covariance is retained through sharp conditional uncertainty bounds, not guessed as zero. A separate calculation removes the common intercept and tries to reconstruct covariance shape using nonnegative outer products of released variant vectors. A residual indicates failure of that reconstruction; the code does not repair the published covariance or reinterpret optimized weights as physical priors.

The photometry inventory should contain 23,007 rows with no negative released fluxes after resolving the documented alternate file outside the distance sample. Peak-time/header identities are checked. This positive-only product is distinct from signed precursor observations.

## Signed pre-explosion baseline

`signed-baseline` reads ten frozen DES16 author precursor files directly. Rows are kept when finite, with positive quoted error, before peak-header date minus 180 observer days; a 365-day sensitivity is also supplied. It imposes no clipping, significance or PHOTFLAG cut.

Within each object/band it fits a weighted constant and calculates the residual squared norm relative to diagonal quoted errors. Expect 3,476 rows across 40 groups, 1,736 negative rows and aggregate `Q/dof=1.05553` for the 180-day cut. The 365-day cut gives 2,856 rows and `Q/dof≈1.02868`. Object-band and complete selected-row tables make the calculation inspectable.

These values describe a chosen signed baseline and a diagonal-error reference. They do not establish independent Gaussian noise, a valid empirical error renormalization, a constant baseline during the SN, or the historical conversion rule for every release object. No baseline is subtracted from the fitted SN data.

## Archived simulated timing

`timing` joins the coherent 2021 NIR and combined FITRES files by CID, verifies equal generating truth fields and retains five missing combined fits. The fixed population contains 30,000 NIR fits and 29,995 common fits. The explicit combined-fit selection is `AV<0.3*RV` and `0.75<STRETCH<1.185`, retaining 25,606 rows.

For the common cohort, expect population SD about `0.0101184 day` for NIR fitted peak minus truth and `0.548510 day` for the combined fit. These use `ddof=0`; a sample-SD report differs slightly. Redshift bins and cut failures are retained. Near-truth fixed NIR timing is a property of the archived simulation procedure, not a measurement of a real-data timing bias or its distance impact.

```bash
uv run --frozen python research.py run raisin --name raisin
uv run --frozen python research.py run signed-baseline --name signed-baseline
uv run --frozen python research.py run timing --name timing
```

Full historical/native refitting and profile-branch experiments remain outside this portable subset. Their prior numerical ambiguity prevents turning an omitted-epoch or timing sensitivity into a single certified correction. This edition retains the closed table-level and signed-photometry calculations without implying that those native inference questions are solved.
