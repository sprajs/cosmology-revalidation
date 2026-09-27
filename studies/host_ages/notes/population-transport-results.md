# Host age, progenitor delay and transport across redshift

**27 September 2026.** The available observations do not yet identify a residual age correction at high redshift. We can measure how the observed low-redshift host summaries change, test whether a single host-to-progenitor conversion is adequate in an author model, and quantify the precision required. Those are different conclusions from proving that the standard correction is sufficient or that an additional correction is necessary.

## An additional public author dataset

The availability search found [Murakami's public `age-of-titans` repository](https://github.com/SterlingYM/age-of-titans/tree/9b9ed9e5faf8b2f6267c0dce8cea97704d9b95ea), an older snapshot last pushed on 19 March 2026. It contains **8,610 transient-associated host summaries**, including stellar age, mass, dust attenuation, metallicity, specific star formation and inferred progenitor age. The earlier statement that no object-level TITAN age catalogue was available was too broad.

We exactly recover its 8,610-name selection from the original 9,361 host rows, its reliable-object list and `d_dlr <= 10`. All 57 original numerical host-property columns agree to `3.6e-15`, the CSV round-trip precision. The table is **not the final 6,983-host sample** in [Murakami et al.](https://arxiv.org/abs/2604.16597v1). The missing final mask cannot be reconstructed merely by matching the final count.

The notebook resolves an important definition: `prog_age_50` is the posterior median of the **mean delay calculated for each sampled SFH**, rather than the median of the progenitor-delay distribution. Its average across these hosts is 3.511 Gyr; their mean mass-weighted stellar-age summary is 6.514 Gyr. The numerical resemblance to the paper's approximate 3.5-Gyr result is not an exact final-sample reproduction.

The full SFH/age samples are absent from the public tree; the notebook reads an external `SN_age_samples.h5`. The original joint likelihood, priors, synthesis configuration, cosmic clock and final selection are therefore not independently reconstructed. The official [DR1 download page](https://titan-snia.github.io/dr1.html) still announces a future release; the [project's paper page](https://titan-snia.github.io/papers.html) describes a newer host catalogue but supplies no corresponding object-level download. This is a bounded search result, not proof that no additional public copy exists.

## A constant conversion does not describe these model-derived ages

A fixed split by transient-name hash assigns 6,898 entries to training and 1,712 to a held-out comparison. The target remains an author's model-derived posterior summary, with shared host-inference assumptions.

| Predictor of the progenitor mean-delay summary | Held-out RMSE |
|---|---:|
| Constant plus linear host age | 0.921 Gyr |
| Cubic host age | 0.689 Gyr |
| Cubic age plus mass, specific star formation, attenuation and metallicity | 0.455 Gyr |

Relative to the affine model, the cubic reduces RMSE by 0.232 Gyr, with a paired transient bootstrap interval of [0.200, 0.264]. Adding the other properties reduces it by 0.466 Gyr, interval [0.434, 0.497]. These intervals resample posterior centres; they omit joint age/property uncertainty and shared calibration. Distinct transient names cannot be guaranteed to represent distinct galaxies because host identifiers are missing.

Within this host-inference model, the mapping therefore contains useful nonlinear and population-dependent information. This does **not** establish independent progenitor ages or establish that the same mapping holds at high redshift. A universal affine conversion, needed for exact slope–contrast compensation, cannot simply be assumed from this table. Conversely, an imperfect affine prediction is not evidence for a new brightness effect.

The public notebook's UVJ quiescent and star-forming selections leave 1,461 entries outside both classes; we retain and report them rather than silently forcing a binary classification. Its high-redshift illustration enters the W26 age–redshift relation as fixed numbers, rather than measuring high-redshift host ages.

Literal name matching produces 16 Pantheon+ rows for **11 distinct SNe**, only six above `zHD=0.01`, all below `zHD=0.048`. Repeated distance rows do not add independent ages. There are no redshift or sky-coordinate columns in this host-summary table, so this linkage adds no observed high-redshift transport information. A separate environmental analysis can use the much larger ZTF overlap while preserving the older snapshot's limits.

## What the disputed 196-host sample actually supports

The original C25 summaries matched to Pantheon+ contain 100 unique SNe at `0.06<zHD<0.20` and 96 at `0.20<zHD<0.42`; the largest actual redshift is **0.37436**.

| Age choice for overlapping catalogues | Mean age, lower bin | Mean age, upper bin | Mean difference; summary bootstrap 95% interval |
|---|---:|---:|---:|
| G11 first | 5.094 Gyr | 4.443 Gyr | 0.651 [0.125, 1.214] Gyr |
| R19 first | 5.032 Gyr | 4.443 Gyr | 0.590 [0.029, 1.120] Gyr |

The resampling unit is the physical SN. These intervals describe variation of the published summaries, not uncertainty in true stellar ages or their shared SED priors. Multiplying the G11-first contrast by an imposed 0.030 mag/Gyr gives 0.01954 mag, **not a measured residual distance correction**.

With the declared 2-Gyr histogram bins, the two measured-age distributions overlap by 83.2% or 86.2%. Entropy-weighting the 100 lower-bin objects to match only the upper-bin mean retains an effective sample size of 92.0 or 93.2. Matching instead an **assumed** 3.1-Gyr mean reduces that to 49.8 or 51.2. Matching one mean does not match the full age/property distribution and does not supply observed `z~1` hosts.

Unknown selection matters even inside this observed support. If independent reweighting in each bin is allowed a maximum-to-minimum density ratio of two, the G11-first age contrast can range from −0.507 to +1.808 Gyr. At the imposed slope, that is −0.0152 to +0.0542 mag. This is a **conditional sensitivity range under an arbitrarily specified weight bound**, not a confidence interval or empirical bound on survey selection.

## Which compensation identities survive

For `tau=a+k*A` and a linear brightness relation, recalibrating `b_tau=b_host/k` preserves `b_tau*Delta tau=b_host*Delta A`. The experiment verifies this on the actual observed host distribution. A bounded, mean-preserving scatter around the affine map also preserves the linear mean identity. Scatter alone does not necessarily break compensation.

We also test constructed quadratic and redshift-dependent maps, always refitting the progenitor slope to recover the same imposed −0.030 mag/Gyr low-redshift host slope. Their between-bin brightness contrasts are 0.02333 and 0.02638 mag, compared with the affine 0.01954 mag. Bounded scatter combined with hypothetical delay-dependent selection gives 0.01041 mag. The discrepancies range over several millimagnitudes despite equal low-bin regression slopes. These scenarios demonstrate non-invariance; their coefficients are not inferred physical mappings.

The scatter construction respects delays between 0.04 Gyr and the age of the Universe in its declared clock. Selection reweights both the delay distribution within hosts and the host population, rather than changing only a median.

## Means, medians, clocks and population selection

A linear luminosity model must be averaged over the **mean** delay. Using the median is a different template. For the existing Son-like clock, B13 cosmic star-formation model and a fixed 0.030 mag/Gyr coefficient, the zero-to-one redshift contrasts are:

| DTD | Mean-delay template | Median-delay template |
|---|---:|---:|
| C14 smooth | 0.13745 mag | 0.16320 mag |
| W26 40-Myr cutoff | 0.10227 mag | 0.07428 mag |

Thus replacing means with medians can move the result in **opposite directions for different delay distributions**. A nonlinear luminosity law requires integrating the luminosity law itself, not evaluating it at either summary age.

Across the three declared clocks and B13/MD14 histories, unselected mean-delay templates span 0.08545–0.13745 mag; median templates span 0.05517–0.16320 mag. Each clock is recomputed inside its SFH convolution. These retain a fixed coefficient solely to compare imposed models; they do not jointly calibrate a physical progenitor slope.

A hypothetical selection factor `exp(±delay/4 Gyr)` changes the Son/B13 C14 mean template from 0.13745 to 0.05190 or 0.21239 mag. These deliberately broad stresses show why a cosmic parent population cannot substitute for a measured survey-selected population. They are not an allowed empirical cosmology-error interval.

## Precision and numerical checks

Before inspecting new outputs, the tolerance was set to a cosmological shift below 0.1 of the baseline flat-ΛCDM saved `sigma(q0)=0.0273248`, for the specific shape `B(z)=a*min(z,1)` with a free magnitude intercept. The full Pantheon+ covariance gives a linear response `Delta q0/a = −1.06747 per mag`, corresponding to **a≈2.560 mmag**. Independent nonlinear distance fits give `Delta q0=−0.002729` and `+0.002736` for the two signs of that amplitude. The small difference from the target 0.0027325 is the nonlinear response. This tolerance belongs to this shape and cosmological model; it is not a universal millimagnitude bound.

On the 196-object design, if the published host-age values were exact, the full-covariance slope uncertainty would be 0.00456 mag/Gyr with an intercept and 0.00562 after projecting out redshift, mass, colour and stretch. The latter retains 65.8% of the age information. Actual host-age uncertainty makes this an optimistic design calculation. The naive classical measurement-error correction gives negative latent-age variance, −0.349 Gyr², so a sample-size forecast based on that approximation is invalid. Better host likelihoods and independent age constraints are necessary, rather than assuming that more copies of the same noisy summaries resolve the problem.

All input bytes are checked against frozen hashes. Across 36 SFH/DTD/clock/selection scenarios, numerical integration is checked by doubled grids, independent logarithmic quadrature and Astropy clocks. Maximum checked differences are `7.63e-7` Gyr for grid medians, `5.14e-9` Gyr for independent mean-delay integration and `9.96e-8` Gyr for cosmic time. Separate fractional-optimization and GLS calculations verify the selection extrema and information formula. These checks validate arithmetic and identities, not physical identification.

## Group-based relative distances: executed observational inventory

[Blum et al. (2026)](https://arxiv.org/abs/2609.13525v1) examine 80 TITAN supernova cousins in 32 galaxy groups. Their group-averaged host-age split gives only a 0.71σ scatter difference; a spiral-fraction split gives 2.06σ. These are scatter comparisons, not a signed age–brightness slope. The complete paper/source archive supplies aggregate tables and figures but no per-SN group crosswalk or distance/age bundle. The [Tully group catalogue](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/AJ/149/171) does supply individual galaxy memberships, allowing a separate observed-cohort inventory.

We match the **host centroids** of the 401 selected ZTF–TITAN events to individual Tully galaxies, requiring one candidate within 3 arcsec and a heliocentric velocity difference below 300 km/s. Membership comes from that galaxy's published group assignment, never proximity to a group centre. Of 89 positional matches, three fail the velocity check and none are ambiguous within the radius. The 86 accepted matches include **56 events in 54 groups**. Only two groups contain repeated SNe; just **one contains two different host galaxies**. These repeated-group counts survive 1- and 5-arcsec matching. Unmatched hosts are not classified as field galaxies. The older TITAN table has no host coordinates, so its age-to-ZTF-host association remains conditional on the common SN identifier and the two DLR cuts.

The distinct-host pair is SN 2018hck and SN 2019ouf, in published group 200007, with host-age posterior centres 5.737 and 9.611 Gyr: a 3.874-Gyr contrast. Their actual host-centroid separations from the matched catalogue galaxies are 0.297 and 0.663 arcsec, and velocity differences are +26 and −160 km/s. The other pair, SN 2020sjo and SN 2020zhh, shares PGC 15090. Its differing fitted global-age summaries are **not two independent physical host ages**, so we exclude that apparent age leverage and pool its host age in the design check.

Even with exact ages and fixed standardization, the single distinct-host pair has an optimistic slope error of **0.0458 mag/Gyr**, under an imposed 0.12-mag intrinsic scatter and uniform group-depth model. That model places each distinct galaxy uniformly within ±1.5 times the catalogue turnaround radius, giving 0.0362 mag depth scatter per member and a 0.1254-mag endpoint pair span. A 0.030-mag/Gyr effect then has only **10.1% power** at a two-sided 5% threshold. These are conditional design calculations, not measured depth bounds or confidence intervals. The two sibling SNe share one host distance and depth offset.

A group intercept absorbs coherent group distance and velocity uncertainty. The catalogue velocity dispersion would correspond to 0.205 mag scatter if individual redshifts were used as distances; it must not be added as an independent error to the group-distance/depth calculation. More decisively, the four-event design with two group intercepts, width, colour and age has five columns but rank four; its rank is already four without age. The residual age information after those controls is numerically zero. Additional mass adjustment cannot restore rank. **No new age slope is fitted:** this cohort cannot identify one while refitting standardization, irrespective of the minimum-group rule. The [executed inventory](../results/population_transport/group-feasibility.json) and [availability audit](../results/population_transport/cousins-availability.json) preserve the counts, provenance and missing measurements.

## Reproduction and remaining measurements

From the repository root:

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/population_transport/run.py
.venv/bin/python studies/host_ages/code/population_transport/acquire_titan.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/population_transport/titan.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/population_transport/validate.py
```

After producing the 401-event join with `studies/host_ages/code/environment_validation/titan.py`, reproduce the group inventory with:

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/population_transport/groups.py
```

This checks pinned Tully source hashes, unique galaxy/group identifiers and an independent spherical-separation calculation (maximum disagreement `1.62e-10` arcsec). The declared matching and error scenarios are recorded before the inventory outcomes in `group-scenarios.json`.

The [compact records](../results/population_transport/) contain scenarios, input/code hashes, numerical summaries and validation. Downloaded author code/data and regenerable object tables remain under ignored `.work/population-transport/`; no author source is vendored into this repository.

The needed measurements are now more specific: the final TITAN mask and object redshifts; joint SFH–dust–metallicity likelihoods with original priors and calibration; independent local-age validation; selection/rejected-event information; and high-redshift hosts measured with compatible likelihoods. Those data would allow a jointly fitted luminosity slope and selected-population map. The existing observations and new public summary table do not yet authorize a physical age correction or a claim that every such correction is excluded.
