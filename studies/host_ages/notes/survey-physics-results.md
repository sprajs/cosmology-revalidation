# Does a simulated age signal survive survey selection and correction?

An additional grey brightness term is **not automatically removed** in the physically supported mock population tested here. With an injected slope of −0.030 mag per Gyr, refitting the standardization and an independently trained correction leaves a differential age association of **−0.02120 ± 0.00224 mag/Gyr**, about 71% of the imposed signal. Holding the nominal correction fixed transmits almost the whole term, **−0.02988 ± 0.00131 mag/Gyr**. These are different estimands: the first includes refitting and partial absorption; the second measures transfer through a frozen nominal correction.

This result does not establish that real supernovae have that additional term. It demonstrates, within one explicit population and dust family, that “the usual correction must already remove any age signal” is not a mathematical necessity. Determining whether an extra term exists in observations still requires independently constrained ages, dust, measurement uncertainty, selection and a supported survey likelihood.

## Simulated experiment

We generated 60,000 attempts in nine native SNANA runs using the released DES five-year cadence, noise/flux-error model, detection trigger, host-redshift efficiency, Dovekie SALT3 surface and calibration. Detector selection and light-curve fitting therefore respond to the changed photons. The reconstructed SNNV19 classifier receives measured timing and redshift; its original-probability equivalence remains unverified. These are pure Type Ia mocks, so this exercises signal acceptance, not core-collapse contamination, classifier purity or BEAMS calibration.

The host population is the released Wiseman 2022 **mock** library. Its `SN_age` is an assumed progenitor delay in Gyr; `mean_age` is a different quantity in Myr. The width-age relation is the released Nicolas/W22 prescription. Neither quantity is a new age measurement. The supplied configuration assigns zero host-mass measurement error, so the correction sees noiseless mock mass. The brightness intervention is

$$
\Delta m=-0.030\,[t_{\rm mock}/{\rm Gyr}-3].
$$

It is applied before noise, detection and fitting. No age value is passed to the distance predictor. Independent correction-training and evaluation seeds share fixed survey cadence and a fixed host library. Common random numbers couple nominal and injected arms only after occurrence-level identity checks; a repeated simulation CID alone is not a unique occurrence, so matching uses CID plus LIBID.

The first 42,000 attempts use the released P23/W22 dust distribution and a separate released P23sys1 colour/dust alternative. That distribution generates negative F99 optical extinction for 25 of 386 primary-selected nominal objects. We therefore generated a separate **18,000-attempt population with fixed $R_V=3.1$** before selection, retaining the other declared distributions. It is a distinct physical-support hypothesis, not a post-selection deletion of unwanted draws. Every populated dust draw in its paired evaluation arms has nonnegative $A_V$ and positive extinction across 1000–30000 Å. Rejected attempts whose physics fields were not populated carry −9 sentinels and are not counted as negative dust.

| Physical-support arm | Attempts | Detection/host/minimum-cut accepted | Native fitted | Quality + reconstructed pIa>0.5 |
|---|---:|---:|---:|---:|
| Nominal correction training | 6,000 | 905 | 598 | 339 |
| Age-injected correction training | 6,000 | 891 | 590 | 335 |
| Nominal evaluation | 3,000 | 474 | 318 | 171 |
| Age-injected evaluation | 3,000 | 462 | 309 | 165 |

## Standardization, correction and remaining response

The independent predictor fits a weighted intercept, width, colour and mass step to simulated distance on the training pool, then learns its remaining residual with 40 nearest neighbours in measured redshift, width, colour and continuous host mass. It applies the fitted standardization and that residual correction once. It does not apply BBC on top. Native light-curve covariance supplies working weights with a declared 0.1-mag scatter floor; this is not the full published survey covariance.

Five prespecified redshift bins, a maximum neighbour radius, and minimum training/evaluation occupancy guard against unsupported extrapolation. The physically supported analysis retains 162 nominal evaluation objects, 158 injected objects under nominal training and 157 under age-injected training. Its z>0.9 bin fails support and is left unestimated. Reported errors come from 500 LIBID-block bootstrap realizations that independently resample training and evaluation pools while coupling paired arms. They quantify conditional Monte Carlo uncertainty, not uncertainty over real galaxy populations or alternative physical models.

| Injected age-response contrast, mag/Gyr | Released P23/W22 dust | Fixed $R_V=3.1$ population |
|---|---:|---:|
| Injected evaluation minus nominal evaluation, nominal correction held fixed | −0.02962 ± 0.00055 | −0.02988 ± 0.00131 |
| Age-refitted correction/evaluation minus nominal correction/evaluation | −0.02168 ± 0.00118 | −0.02120 ± 0.00224 |
| Change from refitting, evaluated on common injected-object support | +0.00794 ± 0.00108 | +0.00864 ± 0.00192 |

The two populations support the same limited conclusion: this correction family absorbs part, but not all, of the injected term. Approximate fractions should not be transported to observed cosmology. The remaining redshift-mean response after refitting is much less precisely determined. For fixed $R_V$, its four supported bins give −0.0091±0.0210, −0.0005±0.0132, +0.0032±0.0215 and +0.0229±0.0203 mag. There is no identified cosmological shift here.

Nominal correction is not proven unbiased. Its residual age slope is +0.00943±0.00352 mag/Gyr in the released-dust benchmark and +0.00561±0.00535 in the physical-support population. These closure diagnostics must be distinguished from the paired *change* caused by injection. The P23sys1 dust/colour alternative produces an age-slope change −0.00055±0.00388 mag/Gyr in the original benchmark, which neither identifies a true dust distribution nor excludes other dust-age degeneracies.

## Native BBC and recovered inputs

All 250 original missing Dovekie files were verified as valid LFS pointers against the actual current Git tree, including object hashes and byte counts. The real Git LFS batch endpoint returns HTTP 403, “Git LFS is disabled for this repository.” Published GitHub releases have no attached assets. The official [November 2025 SNANA data bundle](https://zenodo.org/records/17591282) restores the Dovekie model/calibration, cadence, host libraries, population PDFs, classifier weights and pipeline inputs. It does not restore the original mock realizations. The [July 2024 DES archive](https://zenodo.org/records/12720778) precedes Dovekie and is not substituted for it.

The released host-redshift map contains 12 probability nodes above one, up to 1.008. Modern SNANA correctly rejects that support. The primary experiment explicitly bounds grid nodes before interpolation; an independent valid alternative bounds the interpolated probability, reproducing the older uniform-threshold semantics. Across 6,000 matched attempts these alternatives retain exactly the same 1,054 objects, despite 123 runtime evaluations requiring legacy bounding. This empirical comparison is conditional on that random sample; it does not make the invalid source nodes physical measurements.

The production-dimensional BBC configuration is **not recovered at this Monte Carlo volume**. A native capacity guard first falsely rejects the final valid allocated high-S/N row; a preserved minimal bounds-check patch exposes the genuine next limitation: only 42 training events meet the released requirement for at least 50 with S/N>60. An explicitly separate S/N>50 diagnostic proceeds to interpolation and rejects all 386 evaluation objects for unsupported multidimensional cells. Those gates are retained.

A separate native **redshift-only** BBC calculation completes. In the positive-extinction population it retains 169 nominal and 162 injected evaluation objects, including a separate age-injected training correction. The released alpha/beta/mass-step fit runs and all three native 20×20 statistical distance-bin covariance matrices are positive definite. Their covariance is conditional on the training correction and is not the independently propagated bootstrap covariance above. This supported native route is a benchmark, not production BBC or SALT-surface retraining.

## Numerical checks and scope

Independent review verifies 6,000 occurrence pairs in the original benchmark and 3,000 in the physical-support campaign. For the latter, 29,576 noiseless paired epochs recover the imposed grey shift to 3.8×10⁻⁶ mag; generated distance is unchanged. Selection changes 474→462 objects and 382 epoch flags, showing that the intervention was propagated before selection. The true `SIM_mB` is unchanged because SNANA applies `gammaDM` separately; mistaking it for the final generated brightness would reverse the validation conclusion.

Two displaced optimizer starts for 24 quality-selected objects agree within 0.00027 mag in fitted brightness and 0.0037 in chi-square. Nonpositive fit covariances are excluded by the prespecified correlation-matrix check. Classifier single/batch predictions agree within 10⁻⁶. Supported response covariance matrices are positive definite, and first-half versus full-bootstrap slope standard errors agree within about 4.1%. Unsupported bootstrap bins are identified and covariance states how many jointly supported draws remain; repeated-survey frequentist coverage has not been established.

The nine scientific native jobs used 370.6 seconds for generation and 3107.5 seconds for fitting, cumulatively; classifier inference used 8.2 seconds on CPU. Additional sensitivity and engineering runs are separate. Native fitting, rather than classifier inference, dominates this scale. Larger independent simulations can run as parallel native shards with distinct seeds; exact production-dimensional coverage would require substantial additional training support, and original LFS mocks remain unavailable. No observations, shared SALT surfaces or classifiers were retrained, and no new cosmological correction is claimed.

Reproducible code: [survey_physics](../code/survey_physics/README.md). Compact numerical evidence: [response](../results/survey_physics/correction-response.json), [physical-support response](../results/survey_physics/positive-dust-response.json), [native BBC](../results/survey_physics/positive-dust-native-bbc.json), [validation](../results/survey_physics/validation.json), and [independent physics review](../results/galaxy_validation/survey-physics-review.json).
