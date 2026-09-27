# What survives a common-sample age comparison

**New observational calculations, 27 September 2026.** The additional host-age term is still unresolved. Its estimated size depends materially on the age-error and population model, while its incremental prediction benefit is small. These observations do not justify either declaring all age evolution removed or appending the entire historical correction to current distances.

The comparison uses 196 distinct SDSS supernovae at .06 < zHD < .42, matched to the same Pantheon+ release and full released STAT+SYS covariance. G11-first and R19-first age choices are both retained. The 33 objects shared by the original catalogues are not independent confirmations. Source rows and each signed correction are recorded in a regenerable object ledger.

## Current corrected distances

Holding the published age medians fixed, the marginal age slope is −.00495 ± .00456 mag/Gyr. Fitting redshift, host mass, SN colour and width alongside age gives **−.01101 ± .00562**. R19-first gives −.01086 ± .00559. These are conditional generalized least-squares errors; they omit age uncertainty and do not establish a physical progenitor-age slope.

With those same nuisance predictors, the slopes from corrected, bias-reversed and reconstructed Tripp residuals become −.01101, −.01098 and −.01097 mag/Gyr. Much of the difference between their marginal slopes is therefore associated with measured nuisance relationships. That is an accounting result on this selected sample, not proof that the released corrections remove the correct physical effect. Refitting width and colour while retaining the released covariance and bias model also does not reproduce a complete survey standardization refit.

A Gaussian errors-in-variables sensitivity model gives a substantially less restrictive result. It allows the latent mean host age to depend on redshift, mass, colour and width, integrates over individual Gaussian age uncertainties, and retains the full distance covariance:

| Residual definition | Age slope MLE, mag/Gyr | Nominal 95% profile range |
|---|---:|---:|
| Corrected distances | −.02440 | [−.056, +.008] |
| Exported bias reversed | −.02704 | [−.060, +.004] |

The ranges use Δχ² = 3.84 on a .004 mag/Gyr grid. All 62 profile optimizations converged; three starts agreed in negative log likelihood within .000003. The residual-scatter estimate approaches its lower boundary, and nominal interval coverage is not calibrated. These are neither Bayesian credible intervals nor measurements with the unavailable full host likelihoods. In particular, the input ages are posterior medians and interval summaries, not demonstrably classical unbiased measurements. The naive variance-minus-mean-error-variance calculation is negative in both catalogue choices. **The proposed −.03 scale is allowed in this model, but so is zero.**

## Held-out predictions and remaining age support

In 20 repetitions of five-fold prediction, whole physical supernovae are held out. The conditional Gaussian predictions include the full covariance between training and held-out objects and uncertainty in the fitted nuisance coefficients. This is validation within one observed survey sample, not independent-survey replication or retraining of the released corrections.

For corrected distances, adding age changes G11-first RMSE from **.154086 to .153416 mag**, a .000670 mag improvement. The paired mean-squared-error gain is .000206 mag²; a descriptive physical-SN bootstrap spans −.000542 to +.000995 mag². R19-first improves RMSE by .000499 mag.

A further, explicitly exploratory calibration simulates 5,000 correlated Gaussian error vectors with this covariance and the same repeated fits. It gives one-sided gain p-values .0482 and .0562 for G11-first and R19-first. This is borderline evidence under that conditional model and not stable to the age-catalogue choice. It does not include unknown age, population or correction-training uncertainty. The nominal high power for an injected −.03 relation assumes the *measured medians themselves* are the exact injection coordinate; it must not be transferred to an uncertain true-age slope.

After projecting out redshift, mass, colour and width, about 66% of the measured-age information remains. Redshift alone therefore does not erase most of this particular sample's age leverage. However, the .35–.42 bin contains only four objects, with ages 3.23–4.52 Gyr. This is poor support for the old-versus-young contrast at the upper end, and there are no matched high-redshift ages here to establish cosmological transport.

The previous 4.1% age-information loss from fitting a local Ωm derivative used tabulated plotting errors. On the identical G11-first design the loss is 4.062% with those errors, 3.921% with the released covariance diagonal, and **1.228% with the full covariance**. Thus the newer retained fraction .987719 and the older 4.1% figure are consistent calculations with different weighting. None bounds arbitrary population evolution.

## The young-host redshift adjustment

The original [Gupta et al. table](https://cdsarc.cds.unistra.fr/ftp/J/ApJ/740/92/) contains the light-curve parameters and errors needed for the quality cuts. Matching original IDs to the updated C25 ages reproduces **199 objects and 82 young hosts**. Applying the published Park cuts to the original columns reproduces **175 objects and 70 young hosts exactly**. The earlier statement that this original quality crosswalk was unavailable was too strong.

For C25's young subset, weighted least squares gives HR = .11745 − .46823z, close to the published .117 − .467z. Park's quality subset gives .14170 − .62738z, compared with its published .154 − .641z. That latter regression is still approximate: the precise author regression setup and full age posterior are not recovered by matching counts. Unweighted fits differ substantially, showing that the estimator and error choice matter.

The same data are used to estimate and apply this redshift adjustment. We propagated the resulting shared covariance as an explicit linear operator, and independently verified it with 20,000 Gaussian noise realizations. We also performed 2,000 paired whole-SN bootstrap refits, and a five-fold correction fit trained on other physical objects. The bootstrap results below condition on observed ages and the age<4 membership rule:

| Original 199-object G11 sample | Slope, mag/Gyr | Paired-bootstrap SD |
|---|---:|---:|
| Original residual, no preliminary trend subtraction | −.02666 | .00665 |
| Same-sample young-trend subtraction | −.02939 | .00736 |
| Cross-fitted young-trend subtraction | −.02844 | .00762 |
| Age and redshift fitted jointly | −.02796 | .00691 |

The old-data association survives cross-fitting, but the preliminary correction is not an independent discovery. Its paired slope increment is −.00276 ± .00205, with a 95% interval extending slightly above zero. Propagating the fitted trend increases the conditional age-slope standard error by about 4.8%, rather than an arbitrarily large penalty. The original diagonal errors do not explain all residual scatter (joint-fit χ² = 560.4 for 196 degrees of freedom); their much smaller formal fixed-age error cannot be treated as a reliable full uncertainty. Neither bootstrap nor linear propagation restores unreported cross-SN calibration covariance or age/dust covariance.

Changing the host-age estimator is also consequential. On the same 199 objects with the same original redshift, mass, colour and width controls, original G11 ages give an age coefficient −.00066, while C25-updated ages give −.02332 mag/Gyr. The corresponding quality-sample values are +.00065 and −.02457. These are exploratory fixed-age regressions, not a reason to choose one age estimate after seeing its slope. They identify the host-inference model as an essential part of the disputed measurement.

## What remains unavailable

The [C25 official supplement](https://academic.oup.com/mnras/article/538/4/3340/8098234#supplementary-data) contains two summary tables. Its downloaded arXiv source contains manuscript and figures, not a joint posterior bundle. The cited [MC-Age repository](https://github.com/benjaminrose/MC-Age) links [Zenodo 3875482](https://doi.org/10.5281/zenodo.3875482), explicitly the older Rose2019 chains. Those historical chains cannot silently replace C25's updated ages, SFH assumptions and priors. The [C26 arXiv source](https://arxiv.org/src/2605.21586v1) supplies manuscript and figures but no mock-generation code or original age-PDF arrays. Its published data-availability statement reports no new data. This bounded direct check does not prove that no unlinked author products exist; live journal requests returned 403 and that limitation is recorded.

The remaining exact-input needs are the updated joint host likelihood/posteriors and priors; author choices for overlap ages and regression details; original mock arrays and code; and the population/selection/training information needed to turn a local association into an additional distance correction. No new cosmological result follows from this branch alone.

## Reproduction and numerical checks

[Analysis code](../code/age_reconciliation/run.py), [frozen design and later amendment](../code/age_reconciliation/design.json), [compact result](../results/age_reconciliation/summary.json), [input/output hashes](../results/age_reconciliation/manifest.json), and [direct availability audit](../results/age_reconciliation/public-input-check.json) are retained. Downloaded inputs, object ledgers, profile arrays and bootstrap draws remain outside version control.

```bash
.venv/bin/python studies/host_ages/code/age_reconciliation/restore_inputs.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/age_reconciliation/run.py --archive .work/age-inputs
```

The restoration helper copies the two extracted C25 tables from the frozen main input bundle, copies or downloads the pinned Pantheon+ files, and downloads the original Gupta tables. It verifies all six hashes and refuses changed inputs. Restore the main age supplement first using the [data guide](../../../docs/methods/data.md). An existing read-only archive can instead be supplied with `--archive`. The separate `check_public_inputs.py --archive <restored-snapshot-directory>` availability audit also needs the archived C25 supplement ZIP and the two published paper text extractions named in its record; those are not inputs to the numerical comparison.

The numerical implementation independently checks full-covariance GLS, the integrated age likelihood against a 2N-dimensional block Gaussian, shared correction covariance, repeated-fold prediction operators, physical-SN uniqueness, and signed correction closure. The Gaussian block likelihood agrees within 3.6×10⁻¹⁵. Numerical agreement validates those calculations, not their physical assumptions.
