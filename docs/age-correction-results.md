# Does an additional age correction survive modern standardization?

**Research results, 27 September 2026.** The new calculations do not yet establish a universal additional correction, or establish that the existing treatment removes every cosmologically important age effect. They narrow the disagreement: some apparent differences arise from correction definitions and error models; fitting the usual nuisance parameters can absorb part of an age association; and uncertain host ages still limit the test of what remains.

This report follows the [experimental plan](age-correction-plan.md). The [evidence review](age-correction-evidence.md) separates the published arguments. The calculations below include new observational comparisons, controlled recovery tests and newly located public data. The full survey-wide physical inference remains incomplete for the specific reasons given below.

## What the observations say

### The original disputed sample

We recovered the original Gupta light-curve columns and reproduced the **199-object parent sample, 175-object quality sample and 70 young objects in that quality sample**. The earlier statement that these original quality variables were unavailable was too strong. The C25 young-host redshift regression is closely reproduced; the precise Park regression and updated full host-age posterior calculation remain approximate.

On the 196 matched Pantheon+ objects, jointly accounting for age, redshift, host mass, colour and width gives a corrected residual slope of **−0.01101 ± 0.00562 mag/Gyr if age is treated as exact**. A model that instead represents uncertain ages with a conditional Gaussian population gives **−0.0244 mag/Gyr**, with nominal 95% profile interval **[−0.056, +0.008]**. The latter accommodates both zero and the proposed −0.030 scale. These are different models of uncertain measurements, not independent confirmations. Intrinsic scatter approaches its nonnegative boundary, and the nominal profile interval has not been calibrated for coverage.

That Gaussian age approximation is not a replacement for the missing joint host likelihoods. In particular, the naive classical estimate of intrinsic age variance becomes negative when the quoted age-error variance is subtracted from the sample variance. Posterior age medians cannot safely be treated as ordinary unbiased noisy measurements without checking their priors and inference procedure.

Adding age improves held-out corrected-distance RMSE only from **0.15409 to 0.15342 mag**. A conditional correlated-Gaussian calibration gives a borderline tail fraction of **0.048** with G11-first ages and **0.056** with R19-first ages; the descriptive paired-SN bootstrap interval for the gain crosses zero. These measures answer different uncertainty questions. This is a small, age-choice-sensitive hint, not robust independent confirmation of an additional physical correction. [Observational methods and complete results](../studies/host_ages/notes/age-reconciliation-results.md).

### An independent local-age comparison

The DustPedia UV–IR host measurements allow a new comparison with Cepheid-calibrated supernova brightnesses, avoiding a redshift-derived distance for these very nearby objects. Nineteen light-curve rows correspond to **11 SNe in nine hosts**, with repeated events and Cepheid covariance retained as supplied by Pantheon+.

The age slope is **−0.0063 ± 0.0245 mag/Gyr**, or **+0.0034 ± 0.0252** after adding local attenuation and global host mass. These errors condition on measured age summaries. The known-age design has only **23% power** for a −0.030 slope, so its non-detection does not decide the dispute. Leaving out whole hosts gives no predictive improvement from age.

The host measurements themselves expose an assumption worth changing: global minus local 3-kpc stellar age has median **−0.312 Gyr**, with host-bootstrap interval [−0.480, −0.042] and robust scatter **1.087 Gyr**. A global host age is not an interchangeable measurement of the explosion environment. Overlapping photometry induces covariance that the marginal errors do not specify. [Independent-environment results](../studies/host_ages/notes/environment-validation-results.md).

### A newly accessible, larger comparison

A direct search of the TITAN author's repositories found an older public table of **8,610 host summaries**. This corrects our previous overly broad catalogue-availability statement. The final 6,983-host paper sample and full joint SFH posterior files remain unreproduced.

The table supplies a **401-object common sample with ZTF**, using exact unique transient aliases and the declared host-association cuts. It allows a new test of incremental age prediction after width, colour, redshift, host mass, local colour and other host summaries. The added global-age coefficient is **+0.0103 ± 0.0110 mag/Gyr**; its held-out log predictive score changes by **−0.23**, with a conditional host-bootstrap interval **[−2.10, +1.58]**. An inferred-delay predictor likewise gives no clear gain. These errors condition on the posterior age medians. Their median half-width of the 68% interval, **0.929 Gyr**, exceeds the **0.784-Gyr** age spread left after the controls. The brightness and host-model limitations require the qualification in the [full environmental analysis](../studies/host_ages/notes/environment-validation-results.md): these are global, model-inferred ages; their covariance with dust and metallicity is missing; and ZTF's released brightness calibration and blinding must be respected. This is not an independent measurement of individual progenitor ages or a cosmological distance release. A separate sensitivity analysis with a free brightness intercept and narrow redshift-bin intercepts protects the age comparison against common and bin-constant offsets; they do not verify unrestricted object-dependent blinding or calibration. In 300 paired null/nonzero working-model simulations, nominal interval coverage is 96%, null rejection 4%, and detection of a −0.030 coefficient 80%. That calibration does not include the missing host likelihood or survey selection.

### Additional distance and common-group checks

The newer Dovekie distance release supplies **12 usable TITAN matches**, all nearby Foundation SNe rather than the DES main survey. Using its full released covariance gives an age slope of **−0.0416 ± 0.0326 mag/Gyr**, changing to **−0.0091 ± 0.0384** with attenuation and host mass included. Even treating the age summaries as exact, detection power for the proposed slope is only 12–15%. No known TITAN or DustPedia age matches survive the name/sky/epoch searches for the **79-object RAISIN** optical/NIR sample. Those data therefore cannot yet supply the planned independent age–optical/NIR contrast. [Coverage and distance checks](../studies/host_ages/notes/environment-validation-results.md).

We also tested the newly identified galaxy-group route. Matching actual ZTF host positions and velocities to individual published Tully members finds **56 supernovae in 54 groups** within the 401-object common sample. Only one repeated group contains two distinct host galaxies; the other is a same-galaxy sibling pair, whose differing global-age summaries cannot count as independent physical global-age leverage. One cousin pair cannot separately fit age, colour and width, while uncertain group depth adds distance scatter. This is a concrete sample-size and identification limit, not a null age measurement. [Group membership and precision inventory](../studies/host_ages/notes/population-transport-results.md#group-based-relative-distances-a-further-observational-lead).

## Can the standard fit hide an age signal?

Yes, it can absorb part of a signal. But the size depends on the population and the fitted quantities; it is not established merely by centering residuals in redshift bins.

We injected age slopes of zero, −0.010 and −0.030 mag/Gyr into the actual 196-object redshift/predictor design and full covariance. Each amplitude receives the same **500 noise realizations**. Here the age coordinates are deliberately treated as known, so recovery can be checked against a specified truth.

| Analysis | Fraction of small injected slope left in residuals | Recovered residual slope for −0.030 injection |
|---|---:|---:|
| Fit intercept and actual nonlinear ΛCDM distance | 98.77% | −0.02969 mag/Gyr |
| Also refit width, colour and host-mass step | 64.03% | −0.01927 mag/Gyr |
| Jointly fit age with those quantities | The age coefficient is estimated directly | −0.03010 mag/Gyr |

The joint age intervals achieve approximately 94–96% coverage for nominal 95% intervals under this generator. This is a validated conditional calculation, not a validation of uncertain physical ages or actual survey selection.

Cosmology can still shift appreciably even when most of the age slope remains visible. Conversely, imposing the full −0.030 template on the simulated zero-age-effect population shifts mean fitted Ωₘ from about **0.30 to 0.40**. That is a concrete overcorrection counterexample. It shows why the full historical template cannot be assumed to equal the residual correction, without establishing which simulated truth describes the sky.

We also injected grey, chromatic and phase-dependent changes before refitting twelve DES light curves and recalculating a declared synthetic detection rule. Across **1,728 event–intervention–noise combinations**, a +0.090-mag grey shift is recovered as +0.09118 mag at native brightness. A complete F99 dust screen and a phase-dependent colour perturbation leave different standardized responses, and change which simulated events are selected. The calculation uses frozen covariance and a fixed trained SALT surface. It is not the missing complete survey retraining/BBC test. [Recovery results, scope and covariance records](../studies/host_ages/notes/age-recovery-results.md).

## Does changing the progenitor-age model leave the correction unchanged?

Only under restricted conditions. An exact affine host-to-progenitor mapping with a jointly rescaled linear slope preserves slope times mean-age contrast. Mean-preserving scatter can also preserve that identity. Nonlinear mappings, redshift dependence and selection need not.

The newly acquired TITAN summaries do not follow a sufficiently simple universal affine mapping: on 1,712 held-out transient identifiers, prediction of the inferred mean-delay summary improves from **0.921 Gyr RMSE** for an affine age mapping to **0.689 Gyr** for a cubic and **0.455 Gyr** when other host properties are included. These are predictions of another model's inferred quantities, not independent validation of true progenitor delays; shared hosts and joint posterior uncertainty remain unaccounted for.

Means and medians also matter. For the same existing Son-like clock and B13 star-formation model, a fixed 0.030 coefficient gives zero-to-one-redshift templates of **0.137 versus 0.163 mag** for mean versus median C14 delay, but **0.102 versus 0.074 mag** for W26 delay. A linear brightness law averages over the mean, not the median. Neither template is a measured correction to a survey-selected population.

The original measured-age sample stops at **z = 0.37436**. The newly linked low-redshift host tables do not supply the missing observed high-redshift age distribution. The current data therefore cannot certify transport of either a zero or nonzero local correction to z ≈ 1. [Population and transport results](../studies/host_ages/notes/population-transport-results.md).

## Which hypotheses survive?

| Hypothesis from the plan | Current assessment |
|---|---|
| Existing standardization is sufficient | Not established. Weak additional prediction is not a demonstrated millimagnitude bound on residual evolution. |
| An additional age correction is needed | Still possible. Its magnitude and selected-population redshift dependence are not identified by these tests. |
| Appending the full historical template overcounts | Demonstrated for declared simulated populations; a real-data correction of that size remains unvalidated. This is not proof that every age-aware model overcounts. |
| Available observations leave the correction unidentified | Supported as the present inference boundary: age likelihoods, physical mechanism, population transport and complete selection response are not jointly determined. |

The precision target was set before the recovery outcomes to limit a declared distance-bias template to **0.1 statistical standard deviations** in cosmology. For the declared distance-bias shapes, this requires roughly **2–3 millimagnitudes** over the specified redshift interval. This is a shape-specific design target, not a new error bar. The observed-age uncertainties and transport limitations have not met it.

The feasible branches have been executed and checked, but the full physical programme cannot truthfully be declared complete. The remaining requirements are concrete:

- Updated C25 joint host-age/dust/mass/metallicity/SFH likelihoods and original priors; the public historical R19 chains are not interchangeable with them.
- Exact final TITAN sample, joint SFH products and compatible observed high-redshift hosts, together with selection and host-association information.
- Complete, pinned survey training/classification/selection and simulation inputs. All **250** previously missing Dovekie simulation files still return HTTP 404 at the checked public endpoints; the released-equivalence classifier gate also remains unresolved.
- A full refit or joint selected-data likelihood that reruns every affected correction and passes coverage tests. Only then can a new age correction, or a defensible bound on it, be propagated into cosmology.

No new physical correction has been applied to the manuscript's cosmological results. The existing acceleration estimates therefore retain their stated assumptions; the uncertainty above is not converted into an invented deceleration measurement.
