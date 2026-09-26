# Calibration, extinction and population-bias investigation

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

There is no basis for certifying that all supernova biases are zero. The useful
scientific target is to identify which signed corrections are already applied,
which uncertain directions were propagated, which alternatives are actually
distinguishable, and what size of remaining bias the measurements can constrain.
This programme now makes that audit central while retaining the independent
expansion and classification work.

## Accounting before adding a correction

The existing [product ledger](../../light_curve_fitting/notes/correction-ledger.md) records verified release
algebra, signs, overlaps and covariance formats. The
[SALT/dust audit](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/docs/salt-dust-audit/HANDOFF.md) maps the underlying implementation
and all 47 released single-systematic products. The layers requiring separate
tests are:

| Layer | Already embedded | What remains unconfirmed |
|---|---|---|
| Detector/photometry | DES scene-model photometry, passbands and zero points | Raw-image closure; noise/correlations at source flux and host brightness; transfer from calibration stars to transient photometry |
| Foreground extinction | SFD-based MW extinction in observer-frame fitting; HEAD values already include .86 rescaling | Spatial map contamination, extinction-law shape and shared field errors; applying .86 again is wrong |
| Stellar cross-calibration and training | PS1/CALSPEC/stellar SED anchors, colour relations, passband shifts, empirical SALT surfaces | Exact executed Dovekie configuration, stellar extinction branch, stellar-to-SN spectral transport and joint training response |
| Light-curve fit | SALT amplitude/width/colour/time, model covariance, masks and quality cuts | Empirical colour is not physical E(B-V); covariance and clipping reconstruction; implementation equivalence |
| Population/selection | Simulated dust/intrinsic colour/width/host populations, detection and spectroscopic/classification selection | Parent populations and normalized selection likelihood; physical support of low-RV tails; exact classifier reproduction |
| Standardized distance | alpha*x1-beta*c, residual host step and subtraction of BBC bias | A small host step does not remove host-dependent BBC; changes need alpha/beta and bias surfaces refitted |
| Cosmology/covariance | Shared calibration/dust/systematic response modes and statistical uncertainty | Covariance is not a measured signed correction; alternatives outside its basis are not bounded by it |

Original DES supplies systematic covariance to which the final statistical
diagonal is added once. Pantheon+ supplies total covariance. Dovekie supplies
packed total **precision**; invert in full before marginal subsetting. Their
matrix operations are not interchangeable. Original DES and Dovekie are largely
the same observations, so their difference is not an independent experiment.

## Foreground-map alternative: measured response, weak identification

The externally specified primary alternative is .86*(CSFD-SFD), sampled at the
actual 64 pilot SN coordinates. Original-map resampling and direct replacement
are recorded separately. With fixed accepted epochs, trained SALT3 and
measurement-plus-model covariance, the pre-BBC standardized-distance shift has
RMS **.000263 mag**, range **-.000317 to +.001803 mag**. Forty-eight objects have
reliable interpolation footprints; their RMS is .000296 mag.

The fixed map alternative improves conditional residual log likelihood by only
.000901. Its projected residual signal-to-noise is .00987: this pilot has almost
no power to estimate its amplitude. The unconstrained amplitude 9.75±101.30
therefore is not a detected map error. Measurement-only weighting gives a very
different amplitude, reinforcing the identification warning. The map difference
is small in this pilot but cannot certify either map's truth, bound other dust
errors, or supply a survey-wide limit.

Six paired nonlinear refits (five fixed redshift ranks plus the largest map
change) verify the sign and local approximation. Maximum nonlinear-minus-linear
distance difference is 7.38e-8 mag for noiseless reference curves and 7.92e-6 mag
for observed flux. Both arms refit amplitude, width, colour and time with fixed
covariance. The original HEAD float is reread directly, because CSV float32
round-tripping failed an intentionally tight model identity check. No historical
result was overwritten. No selection, BBC or training change is included.

Reproduce with [foreground_response.py](../../dust/code/foreground_response.py)
and [foreground_refit_check.py](../../dust/code/foreground_refit_check.py).
Results: [map response](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/foreground_response/summary.json),
[paired refits](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/foreground_response/nonlinear_check/summary.json).

## Stellar extinction branch: a deliberate convention to test

The pinned Dovekie `get_all_obsdfs(redo=False)` code reads stored star magnitudes
then sets all `_AV` extinction columns to zero. The `redo=True` branch computes
and retains them. This is an explicit difference in the supplied source;
the actual published run's branch is not established by its presence.
The archival posterior/configuration audit now provides a more precise
[provenance report](calibration-provenance.md): the default normal-data MCMC
would use zero extinction **if run from this exact source**, but its current
13-survey configuration does not match either stored 12-survey posterior.
Neither posterior is linked to a specific released DES calibration execution.

A refreshed primary-source check adds relevant context: Grayling et al. explicitly
describe Dovekie's omission of stellar Milky Way corrections as intentional,
following a test that found them negligible. Thus the zero-AV branch is not, by
itself, evidence of an accidental omission. Our independent branch sensitivities
below quantify a conditional effect; exact stored-run provenance remains open.
See [BayeSN × Dovekie, section 2.2.2](https://arxiv.org/html/2606.19429v1#S2.SS2.SSS2).

The supplied DES/PS1 star file has 1,255 rows. The intersection of documented
cuts under the two conventions contains 1,178 stars; zero extinction changes
which stars pass. On the fixed cohort, independent synthetic-minus-observed
offset estimates change by the following millimagnitudes when restoring stored
extinction (linear synthetic colour relation, zero passband shift):

| Spectral library | g | r | i | z |
|---|---:|---:|---:|---:|
| NGSL, 211 standards | -.727 | -1.017 | -1.382 | -1.058 |
| CALSPEC, 23 standards | -.780 | -.567 | -1.508 | +.862 |

These are conditional branch sensitivities, not fitted published zero-point
errors. The z-band sign depends on the library. The full MCMC, white-dwarf
anchors, PS1 anchor offsets, passband shifts, spectral-library mixture and
training response are not reproduced by these descriptive fits. The supplied
AV file was already selected; the zero-extinction arm cannot recover stars
excluded before it was written. No per-star error covariance is supplied.

Linear/quadratic observed colour fits, spatial internal holdout residuals and
own-cut alternatives are saved in [calibration_stars](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/calibration_stars/summary.json),
with [executed code](../../dust/code/calibration_star_extinction.py)
and hashes. The holdout shares PS1/SED/dust assumptions and is not an external
validation of calibration already trained on those stars.

The supplied DES stars cover only RA 40.205–43.816 degrees and Dec -1.154–.999
degrees. The frozen sky split has 830 training stars and 348 test stars, with
only two occupied held-out cells. It is not a survey-wide field replication.
Restoring extinction changes held-out RMS by less than .0001 mag in each band;
this cannot select a dust convention from these residuals alone.

Applying each fixed stellar branch vector to the SN tangent responses, with no
amplitude retuning, gives high-minus-low-redshift quartile pre-BBC shifts of
.00034 mag (NGSL) or .00241 mag (CALSPEC) in the high-Ia subset under
measurement-plus-model covariance. Their fixed residual log-likelihood gains
are .198 and 1.055. This exercises the signed transport of a specific potential
calibration change; it neither establishes that the change is required nor
accounts for retraining and selection. The library-dependent response and
individual-object shifts are preserved in
[SN transfer](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/calibration_stars/sn_transfer/summary.json).

## Leading residual clue and decisive next control

An exploratory template analysis finds a filter-dependent pattern in 64 SN
light curves and in their 43 high-Ia-probability subset. Under native
measurement-plus-model covariance, leave-field-out gains have exact field-sign
probabilities .00781 and .01172 respectively. Removing the fitted template would
produce about .065 mag of low-to-high-redshift pre-BBC distance shape in the
high-Ia subset if alpha/beta were frozen. **This is not a measured correction.**
Measurement-only covariance changes the family preference and removes that
strength of evidence. Rest-frame and observer-frame templates are nearly
degenerate, and regularization matters.

The matched 43-object modern-SNANA control retains all 1,649 epochs. With exact
SNANA mean/covariance, maximum leave-field-out observer gain is 7.624 and the
nine-field sign-family fraction is .015625. Using the independent native mean
with that same exact covariance gives 7.635. Their projected mean mismatch has
chi-square only .002075; modern implementation mismatch is too small to explain
the observed pattern. The native shape/colour/time tangent is still used.

The frozen 43-object vector predicts residuals in all **1,020 other high-Ia DES
objects**, with 39,606 accepted epochs. Its fixed-vector log-score gain is
**27.233**, and all ten fields have positive matched products. Under the
conditional field-sign null the exact fraction is 1/1,024. Only seven fields
have positive gains at the *fixed amplitude*: the descriptive transfer
amplitude is .582±.055 times the pilot value, so the pilot overpredicts the
held-out amplitude. It was not retuned for the primary score. The frozen
rest/phase forecast gives -.692, rather than a comparable fixed-vector gain.
Integrating one shared discovery coefficient posterior gives joint gains
59.229 (observer) and 30.593 (rest); these are separate secondary predictions,
not sums of marginal object scores or a model-evidence claim.

A declared duplicate-epoch sensitivity finds 147 extra repeated epochs in 59
groups across 37 objects. Keeping one occurrence, taking marginal covariance
submatrices and repeating the nuisance projection gives gain 28.221 with the
same field-sign fraction. This neither removes data from the primary nor proves
every repeated record is erroneous. All 1,020 objects remain represented.
Results are preserved in
[frozen transfer](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/validation1020/analysis/result.json).
This is retrospective transfer on published selected objects, not a new survey.

The known classifier mismatch does not drive this particular gain. A declared
post-validation intersection with reconstructed official >.999 membership
retains 1,006 objects; its unchanged-vector gain rises to 30.285. None of the
43 discovery objects disagree. This is a membership sensitivity, not proof of
either classifier's correctness; [exact IDs and fields](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/classifier_residual_sensitivity/result.json)
are retained.

These residual tests deliberately project out each object's fitted amplitude,
width, colour and time. A gray DES zero-point change is absorbed into amplitude;
the tests therefore do not calibrate DES's absolute brightness relative to the
nearby surveys. Dust or SED changes absorbed into colour can also change
standardized distances without leaving a strong residual. Independent stellar,
spectral, population or geometric constraints remain necessary for those
directions.

Both the original and exact-covariance tests condition on fixed global SALT
surfaces and calibration. Common training/calibration uncertainty can violate
the field-sign symmetry assumption; per-object SNANA covariance does not
replace the released cross-SN systematic covariance. A surviving pattern is
not automatically outside the already modelled systematic budget.
See the evolving
[science review](measurement-identifiability.md) and its linked protocols/results.

### A broader spectral alternative changes the identification question

The old rest-frame comparator does not represent all plausible empirical
spectral variation. An outcome-blind calculation using actual passbands admits
a prespecified smooth wavelength/phase mean and smooth redshift drift. A
52-coefficient member with a declared .05-mag full-grid RMS budget, chosen
using discovery response geometry, reproduces **98.12%** of the frozen observer
vector's validation information after finite broadband integration. Its RMS on
actual wavelength/phase support is .04263 mag and maximum absolute change there
is .14762 mag. These budgets are sensitivity assumptions, not measured
population priors. No observed flux residual was fitted in this calculation.

This demonstrates that the observer pattern cannot be assigned uniquely to
photometric calibration by defeating the earlier small rest-frame model. The
constructed SED and observer changes also give substantially different
fixed-reference distance responses. The completed full nonlinear refit retains
all 1,020 objects and 39,606 epochs. On the unchanged high 255 minus low 255
redshift groups, the observed-flux responses are **+.060585 mag** for the
observer change and **-.186125 mag** for the SED change: a **.246710-mag**
separation. Their fitted flux-change vectors have cosine .99060 and squared
mismatch .01900 relative to the observer norm. Three initial numerical gradient
failures were retained and resolved with smaller difference steps under the
original gates; no object was removed. Covariance, spectral training, global
standardization and selection remain fixed. The
[SED report](sed-identification.md) preserves the whole-subspace test,
amplitude frontier and finite-response checks; the
[full nonlinear report](sed-nonlinear-validation.md) records numerical closure
and the remaining physical restrictions. Neither example is an inferred bias
or a corrected cosmology. Independent population or spectral measurements,
and their calibration dependencies, are now the decisive gate.

[Constructed nonlinear responses, not measured biases](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_ambiguity_figure/constructed-ambiguity.png)

Accounting for shared modes changes the comparison. In a separate
[fixed-offset predictive probe](sed-shared-probe.md), the same twelve-mode
discovery-conditioned null gives gains -1.3034 and -7.4799 nats for the full
observer and SED responses. Both full amplitudes worsen that null prediction;
the observer probe fares better by 6.18 nats. Shared calibration uncertainty
covers much of their common response but less of their difference. This is
partial discrimination within the specified conditional calculation, not a
physical-model posterior or a validation-tuned correction.

The calculation exposed ten early epochs with nonpositive independent SALT
means, across eight validation objects. Re-exporting native SNANA derivatives
on all eight changes the original gain from 27.233214 to **27.229805**;
halving the derivative step gives 27.230103. Removing exactly those ten rows
with marginal covariance and native reprojection gives 27.219293. Thus this
specific implementation edge does not explain the transfer, while native
derivative closure for the entire sample remains unclaimed. The exact native
means, masks and data reproduced in all targeted exports; see the
[targeted check](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/negmean/result.json).

The next executed control uses the actual released calibration shift lists and retrained
SALT surfaces. The [shared-uncertainty source audit](shared-calibration-uncertainty.md)
distinguishes what per-object covariance contains from the cross-SN distance
systematic budget. Upstream Fragilistic covariance is also present in the pinned
Pantheon+ calibration directory; its precise relation to the released variant
draws still needs explicit mapping. No new calibration term may be counted twice.

The native 43-object coupled-variation check has passed exact parameter,
epoch, observed-flux and quoted-error identity gates for all nine models.
Their span covers 95.43% of the frozen observer prediction's squared norm.
With the original .3 amplitude weights, a descriptive shared-covariance stress
retains 55.4% of its information. This is strong geometric overlap with existing
tested uncertainties, not an explanation of the amplitude or full-budget
posterior inference. The [result](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/shared43/result.json)
retains all modes, unchanged weights and the independent upstream covariance
comparison. The larger test conditions the same shared calibration variables
on discovery before predicting validation, as reported below.

That discovery-conditioned comparison is now complete for the nine coupled
modes. Calibration-only predicts the joint validation residuals with log-score
gain **49.803** over the fixed nominal model; calibration plus an extra
observer-band term gives **58.747**, a difference of **8.944**. A physically
isotropic observer prior with the same total squared prior magnitude gives
58.678, so this particular prior-orientation sensitivity has little effect.
Both models learn their shared variables from the 43 discovery objects only.
An independent Astra calculation from the saved projected arrays agrees within
7.7e-13; see [shared-mode review](shared-mode-review.md).

The originally frozen observer direction has conditional matched statistic
3.286 after centering on the calibration-only forecast and propagating its
shared uncertainty, compared with 10.586 under fixed nominal calibration.
Applying its *full original amplitude* now has gain **-.964** relative to that
calibration forecast. Thus the original amplitude is not justified as a new
correction. The extra empirical family still improves this limited Gaussian
predictor, but this is a retrospective nine-mode attribution test, not a
full-budget significance or physical calibration posterior.

The declared twelve-mode sensitivity now also includes CALSPEC, MW-map-scale
and MW-colour-law variations with their actual source weights. Calibration and
extinction alone give joint gain **50.071**; adding the empirical observer term
gives **58.752**, a difference of **8.681** (isotropic-prior difference **8.611**).
The conditional fixed-direction statistic is 3.234, while applying the original
full amplitude gives **-1.225**. These added released variations do not remove
the remaining predictive preference, but they also do not identify its cause.

CALSPEC changes the observed-flux convention. Its source-derived band scale
reproduces both exported fluxes and quoted errors exactly; model means are
mapped back to nominal flux units before differencing. The MW variant correctly
rescales the already processed E(B-V) by .95. All 43 discovery and 1,020
validation objects retain their coordinates and epochs. The
[expanded result](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/expanded12/result.json)
preserves the prior nine-mode comparison. The next control asks whether fitting,
selection and intrinsic scatter generate a similar residual in simulations;
these twelve modes are still not a complete systematic-error model.

The independent review also found an unresolved variance-allocation assumption:
nominal SNANA covariance already contains a within-object MW term based on
`MWEBV_ERR=.05*MWEBV`; the added released MW-scale mode is shared globally.
These are distinct covariance structures, but the code alone does not establish
that their physical sources are independent. The twelve-mode calculation keeps
the declared nominal covariance and released variation rather than silently
subtracting a component. It is not twelve independently validated error sources.

All 3,189 added object/variant pairs pass an independent source-mapping check.
Independent SVD score calculations agree within 9.3e-13. The expanded observer
model has 15 latent coordinates but only rank 14: the CALSPEC band response
overlaps exactly with the observer family after the gray tangent is removed.
The prior therefore determines one part of that attribution; data cannot
separately identify every named component. This does not invalidate the joint
predictive calculation, but it prevents interpreting its coefficients as
separately measured physical errors.

A source-weight discrepancy is now quantified as well. The methods paper
specifies colour-law covariance weight 1/3, whereas the saved configuration and
covariance code imply amplitude .3 and hence covariance weight .09. We preserve
the configuration result and repeat the calculation with amplitude sqrt(1/3).
The added observer-family gain becomes **8.301** (isotropic prior **8.226**),
and the conditional fixed-direction statistic becomes 3.173. Its full original
amplitude still worsens the forecast, with gain **-1.573**. This modestly weakens
the same finite-mode preference; it does not resolve which weight was actually
used in the published execution. The separate
[weight sensitivity](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/shared_mode_review/colorlaw-weight-sensitivity.json)
retains both conventions and independent numerical verification.

A complementary [Gaussian fitting control](../../light_curve_fitting/notes/gaussian-fitting-control.md) now
checks nonlinear nuisance fitting on six fixed cadences. All 768 fits and 48
second-start checks pass. Against paired linear projections, nonlinear fitting
changes the summed matched product by -.0235±.0045 Monte Carlo SE: a small
effect opposite the positive data pattern in this restricted generator. Fixed
covariance and masks exclude selection, clipping, physical scatter and shared
calibration, so this is not a survey-wide bound on fitting bias.

## Nonlinear distance response is not physical identification

The [shared-mode distance diagnostic](shared-distance-response.md) now measures
how much of the released-variation distance response escapes the residual test.
For a fixed high-minus-low-redshift quartile contrast in the 1,020 validation
objects, the twelve-mode prior SD is .01379 mag and the discovery-conditioned
SD is .01279 mag: **86.0% of the prior variance remains**. This is a finite-mode
contribution at fixed standardization, not total uncertainty or a measured bias.
Using all 1,063 objects' design information, without fitting residual means,
would reduce that same finite-mode SD to .00942 mag; 46.6% of prior variance
would remain. The 86% figure therefore describes the discovery prediction,
not an absolute information limit of the full data.
Adding the empirical observer family changes the conditional mean response to
.0583 mag with SD .0257, but does not make that response a physical correction.
Exact gray dimming remains invisible to the projected residual test.

Six fixed redshift ranks of the discovery sample were independently refit after
the hypothetical opposite signed filter correction. Observed flux and its full
covariance were transformed together: for diagonal scale S, `C_new=S C S` and
`L_new=S L`. All paired fits converged; the largest nonlinear-minus-local-linear
standardized-distance difference was .00176 mag. Alpha/beta and training remain
fixed and BBC/selection is not rerun.

For CID 1339232 at z=.969, that hypothetical correction shifts the standardized
distance by **+.1789 mag**, while chi-square changes from 23.4529 to 23.4449.
This illustrates weak identification through fitted colour: a nearly unchanged
light-curve fit can coexist with a substantial inferred-distance change. It
does not measure a .1789-mag bias in that object. The full signed response and
other five cases are in the [paired refits](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/calibration_pattern_refit/summary.json).
An independent second-start check agrees in distance shifts within 1.06e-7 mag
and verifies the covariance transformation; see
[numerical verification](calibration-verification.md).

[Hypothetical calibration response, not an estimated bias](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/calibration_response_figure_v2/calibration-response.png)

Any supported effect must then pass nonlinear refits, calibration/SALT training
response, regenerated detection/classification/BBC, shared uncertainty
propagation and injection/coverage tests before corrected cosmology is claimed.

## Testing the calibration and population assumptions with infrared data

A [paired RAISIN audit](../../infrared/notes/raisin-differential.md) now follows the corrections
through optical and infrared distances for the same 79 objects. Comparing each
object with itself cancels its common cosmological distance. The high-minus-low
redshift optical–NIR contrast is +.00124 mag after the released corrections,
with descriptive paired scatter SE .04497 mag. Removing the exported bias and
host-mass terms gives +.07541 mag. The differential contributions of those terms
are -.04455 and -.02962 mag. This is arithmetic on fixed released fits, not a
new physical correction or justification for removing existing corrections.

The comparison also shares optical timing, sample selection, calibration and
many physical objects with DES/Pantheon. Cross-branch fitted-distance covariance
is not supplied. Several systematic variants fail to reconstruct their named
covariance exports. A source-linked defect is now established: the mass-divide
variant fits its step at 10.44 dex but applies it using a hard-coded threshold
of 10. At fixed published step amplitudes, the corrected variant changes the
high-minus-low distance response by .003289/.006885/.009238 mag in NIR,
optical and combined branches. Its mass-group uncertainty gains a second
direction; optical and combined contrast SDs change .02800→.03168 and
.03534→.04344 mag. The nominal distances are unchanged. Independent raw-table
reconstruction verifies these arrays within 7.4e-15. Full covariance closure
still fails, preventing a supported total-likelihood update. Final optical–NIR
agreement alone therefore cannot confirm that
the pre-existing corrections are unbiased. Conversely, the contrast after
removing terms does not establish physical evolution. Forty independent raw-table
arithmetic checks reproduce the results within 1.64e-15.

The [intercept-independent covariance check](../../infrared/notes/raisin-covariance-quotient.md)
now rules out a harmless common-magnitude-offset explanation for the total
covariance mismatch. The mismatch remains in the distance-shape subspace;
neither published nor reconstructed covariance is silently substituted for
the other.

The [BayeSN audit](../../spectral_models/notes/bayesn-distance-identification.md) traces the actual 23
survey/filter conventions, verifies a calibrated forward operator on six
metadata-selected objects, and derives proper distance-amplitude integration.
It also shows why a fit with an external cosmological distance prior and fixed
trained intrinsic-colour distribution cannot independently establish that
distant-SN dust or luminosity evolution is absent. Local sensitivity calculations
are design diagnostics, not observed dust posteriors.

An upstream measurement-selection issue now takes priority. The
[RAISIN sign audit](../../infrared/notes/raisin-flux-sign-audit.md) finds that all 23,007 published
flux rows are positive. In the original DIFFIMG series of 17 matched DES
objects, the historical optical phase window contains 56 negative rows omitted
from the product and 811 positive rows all retained. Ten author DES16 precursors
corroborate the sign-dependent retention, but their positive-row error estimates
do not exactly match the release. The
[archived simulation cadence](../../infrared/notes/raisin-simulation-sign-accounting.md) repeats the
retained times. Drawing new noise on those times alone does not reproduce a
cut on each new realization's measured flux sign.

The [frozen paired-refit protocol](../../infrared/notes/raisin-signed-refit-design.md) therefore
separates release-to-author processing changes, positive-to-signed author
measurements, and a labelled release-plus-negative sensitivity. It preserves
the actual optical model/calibration and checks fixed versus floating peak
timing. This will quantify a conditional fitted-distance response; a corrected
cosmology additionally requires validated population/selection simulations,
the effect of optical timing on NIR distances, and covariance closure.
