# Independent RAISIN, CSP and sign-selection revalidation

The six workflows reproduce their numerical summaries. Independent calculations agree with the published outputs under `results/baseline/`; no defect in these six default workflows was demonstrated. That result verifies the documented calculations, not every physical assumption or the full historical RAISIN cosmology analysis.

The executable audit is [raisin_checks.py](../../validation/raisin_checks.py), with all numerical results and 244 inspected-file hashes in [raisin.json](../../validation/reports/raisin.json). It reads the FITRES tables, original CSP archive, signed precursor rows and spectral tables directly. It does not import the workflow readers, numerical helpers or fitting routines. It records 46 explicit numerical comparison gates plus structural assertions. The methods below were chosen after the existing results were visible; these are validation and sensitivity checks, not a blinded discovery analysis.

```bash
.venv/bin/python validation/raisin_checks.py --results baseline
```

## 1. What changes when released distance corrections are removed?

The same 79 physical supernovae are present, in the same order with matching redshifts, in the infrared, optical and combined distance tables: 42 at low redshift and 37 at high redshift. For each object we first subtract its infrared distance from its optical or combined distance. We then subtract the low-redshift mean from the high-redshift mean. This pairing removes the common distance–redshift relation, while the second subtraction removes a constant difference between branches.

| Branch minus infrared | Released contrast (mag) | Remove both exported bias and mass terms (mag) |
|---|---:|---:|
| Optical | +0.001236 | +0.075412 |
| Optical + infrared | −0.027098 | +0.074821 |

The signs close algebraically: undoing the released bias term means **adding** `DLMAG_biascor`, while undoing the mass term means **subtracting** `MASS_CORR`. Calculations that reverse these conventions would create an artificial discrepancy. Removing these terms does not refit any light curve, change the selected objects, or establish that the terms should be absent.

A new descriptive sensitivity check resamples whole supernovae within each redshift group, preserving the paired branches. Across 10,000 fixed-seed resamples, the optical released contrast has a 95% percentile interval of approximately −0.0880 to +0.0847 mag; without both terms it is −0.0133 to +0.1583 mag. The corresponding combined-branch intervals are −0.0808 to +0.0263 mag and +0.0214 to +0.1298 mag. **These are conditional resampling intervals, not full uncertainty intervals for a physical correction:** common calibration/model systematics, uncertain cross-branch covariance and selection effects are excluded. In particular, the last interval must not be presented as evidence for a new cosmological signal.

All 23,007 inventoried released photometry rows have positive finite flux and positive quoted error. This inventory describes the release; it does not determine the original censoring/selection mechanism. Signed author precursor observations exist and are separately audited below.

## 2. Covariance reconstruction remains incomplete

The workflow restores each systematic diagonal from the difference between total and statistical `lcparams` error variances. The off-diagonal covariance export has zero diagonal. The independent audit checks the associated redshift ordering and recomputes the restoration.

The reconstruction asks a specific question: can centered outer products of the raw `DLMAG` shifts in FITOPT 1–28 represent the centered released systematic matrix? The original nonnegative fit leaves relative residuals of 13.4536% for infrared, 1.4377% for optical and 0.7353% for combined. A new unrestricted least-squares fit, allowing negative coefficients, leaves respectively 13.4535%, 1.4377% and 0.7353%. Thus positivity constraints are not the reason for the residual in this particular representation.

This does **not** prove that a published covariance is physically wrong. It diagnoses failure of the supplied vector representation, which is not a full execution of every historical source operation, including possible redshift-dependent transformations, centering conventions and variant definitions. Original covariance-generation provenance remains necessary.

Tiny negative eigenvalues of the reconstructed systematic-only matrices, about −2 to −3 × 10⁻⁷ mag², lie within conservative bounds propagated from their printed precision. The total statistical-plus-systematic matrices have positive minimum eigenvalues. Consequently, those tiny negative eigenvalues are not evidence of an invalid physical covariance, and no eigenvalue clipping or covariance repair was applied.

## 3. Archived simulation timing is an input convention, not measured precision

The coherent archived simulation tables contain 30,000 infrared fits and 29,995 combined fits. The five absent combined fits remain identified explicitly. Generating peak, distance, dust and stretch truth match exactly for the common objects. The documented combined-fit cuts retain 25,606 objects.

The infrared fitted peak equals `PKMJDINI` **exactly in all 29,995 common fits**. Its fitted-peak-minus-truth population standard deviation is 0.0101184 day; the combined fit gives 0.548510 day. The small infrared scatter therefore describes a near-truth fixed initializer, not timing precision recovered from noisy data. Changing from population to sample standard deviation is a separate, small reporting convention.

The new direct selection accounting also shows that the mean infrared distance residual changes from +0.03292 mag in all common fits to +0.04364 mag among selected fits; the rejected fits average −0.02963 mag. For combined fits the corresponding values are −0.00281, +0.00543 and −0.05091 mag. These are descriptive residuals in an archived simulation under its fitted model and cuts. They are not observed-data biases or a proposed correction. A causal timing-to-distance conclusion would require controlled simulations that keep the other generating/fitting assumptions fixed.

## 4. CSP differences: labels, source coverage and spectral response

The original archive and its SNooPy intermediates contain the same 5,491 near-infrared rows when the physical labels `Jdw → J` and `Hdw → H` are accounted for. Decimal-literal multiset comparisons preserve duplicates and match time, magnitude and uncertainty. No magnitude-based matching is used to assign instrument identity.

Of 76 LISTed release objects, 73 pass the complete converter identity. That count comprises **65 objects with near-infrared rows and eight empty identities**. The remaining three objects—2012fr, 2012ht and 2015F—supply 241 released rows but have no source-object rows in this frozen CSP-I archive. They also lie outside the selected 79-supernova cosmology sample. These are demonstrated source-coverage gaps, not demonstrated conversion failures. Of 520 raw `Jdw` rows, 303 carried into the compared release have an unambiguous physical WIRC-J label; comparing those two counts alone does not establish that the other rows were incorrectly converted.

A separate passband calculation integrates a chosen archived supernova spectrum through physical throughput curves, normalized to the BD17 reference spectrum. An independent analytic cubic integral on every union-knot interval agrees with all 140 grid results within **2.83 × 10⁻¹⁵ mag**. For the WIRC-minus-RC1 J comparison, the grid spans −0.06388 to +0.00462 mag and the largest absolute phase contrast is 0.06850 mag. These relatively large conditional spectral differences are numerically real for the supplied curves. They are not measured distance corrections: spectral realism, actual historical calibration, template training and a consistent light-curve/bias refit remain separate questions.

## 5. Signed precursor baselines and positive-only selection

For ten DES16 objects, the pre-peak-minus-180-day sample contains 3,476 signed rows in 40 object/band groups; 1,736 are negative. Independent QR projection of a weighted constant reproduces `Q/dof = 1.055531`. The stricter 365-day cut has 2,856 rows, 1,428 negative, and `Q/dof = 1.028682`.

Whole-object resampling preserves all within-object measurements and gives conditional 95% percentile intervals of approximately **0.930–1.192** and **0.935–1.113** for the two `Q/dof` values. With only ten objects these intervals are approximate and do not cover shared cross-object systematics. The 180-day aggregate falls to about 1.0023 if the high-scatter object DES16X3zd is omitted; the object is retained in the primary calculation. This sensitivity argues against treating the pooled excess as a universal error-rescaling factor. The cuts are nested, so they are not independent confirmations, and `Q/dof` close to one does not establish independent Gaussian noise.

Finally, the 128-replicate synthetic sign-selection experiment uses a true amplitude of 0.4. Mean fitted amplitudes are 0.40693 using all signed fluxes, 1.20234 when positive values alone are treated as uncensored Gaussian data, and 0.39595 using the full known eligible schedule and omitted signs. Independent closed-form dense-covariance GLS reproduces every all-signed and naive-positive fit within 4.42 × 10⁻⁸ amplitude units. Adaptive integration of the conditional latent factor, followed by refitting, agrees within 4.20 × 10⁻⁸ for 12 prespecified evenly spaced replicates.

Because each method uses identical random draws, paired differences give a sharper Monte Carlo comparison. The naive-positive minus all-signed mean difference is +0.79540 with Monte Carlo SE 0.02630; censored minus all-signed is −0.01098 with SE 0.00893. This validates a useful method demonstration: discarded signs and the eligible schedule can matter greatly. It does not establish that this Gaussian model or known schedule describes the real survey, and it supplies no observational cosmological correction.

## Sources and interpretation boundary

The frozen original products and retrieval identities are in [inputs.json](../../provenance/inputs.json) and [SOURCES.md](../sources.md). The direct primary inputs are the pinned RAISIN DataRelease commit `a383c4b`, author RAISIN_cosmo simulation commit `aaa709e`, signed precursor commit `b888214`, original CSP DR3 archive, and the archived KCOR SN/BD17 spectra and physical passbands. This audit does not reacquire or silently replace them.

The practical outcome is a reproducible chain of table accounting, lineage checks and numerical controls, with the unresolved steps left visible. A full new cosmology measurement still requires source-complete calibration and covariance provenance, realistic selection/noise validation, consistently refitted light curves and bias simulations, and a declared joint likelihood.
