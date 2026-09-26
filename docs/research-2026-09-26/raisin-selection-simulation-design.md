# RAISIN signed-flux selection: a paired generative simulation design

The archived positive-only cadence and a new positive-flux cut are different
statistical operations. Reusing the former with newly drawn signed noise does
not reproduce the latter. A fixed-truth Gaussian demonstration verifies that
distinction below. No observed light curve, dust posterior, native fit or
cosmology was fitted for this review. The proposed native experiment measures
a conditional estimator response before any claim about a published bias
correction.

The design and executable toy are in
[raisin_selection_simulation_design](../../runs/research_2026_09_26/raisin_selection_simulation_design/).
The [protocol](../../runs/research_2026_09_26/raisin_selection_simulation_design/protocol.json)
was frozen before synthetic outcomes, SHA256
`ae0fad5eea5876080d86d53222b67cb61f0e6c4a4ee4733de5a306601c48724c`.
It hashes the three input reviews, source configurations, cadence result,
native input manifest and source-row ledger. Observed flux values are never
used as the generating signal. Fixed error columns can be used only in the
explicit conditional engineering experiment described below.

## What the source evidence establishes

The [sign audit](raisin-flux-sign-audit.md) finds all 23,007 released rows
positive. For 17 DES aliases, the original DIFFIMG ledger has 56 negative
and 811 positive rows inside the header-defined optical phase window
[-7,+45] rest days: zero negatives and all 811 positives survive. The author
DES16 precursor permits a more closely matched ten-object comparison, but
its quoted errors differ slightly from the released product. The DES15
reduction bridge is materially worse. No exact script or execution record
creating the sign-selected optical product has been found. These facts
support a sign-omission counterfactual for the documented DES subset; they
do not establish an identical transformation for CSP or PS1.

The [simulation cadence audit](raisin-simulation-sign-accounting.md) finds
archived SIMLIB presence identical to released retention for every one of
9,843 original DES rows. `genSimlib.py` constructs a cadence from existing
object times/bands and imports sky/PSF/zero-point metadata from other survey
libraries; its private input is not execution-linked. The archived DES
generator sets `SMEARFLAG_FLUX: 1`, `SMEARFLAG_ZEROPT: 0`,
`GENMAG_SMEAR_MODELNAME: OIR.J19`, per-SIMLIB peak/redshift and a flux-error
model. Its search-efficiency and generation cuts are separate from an epoch
sign cut. `CUTWIN_SNRMAX: 3 griz 2 -20. 80.` conflicts with the adjacent
comment saying SNR>5; the executable numeric value, parsed by the pinned
binary, must govern a reproduction. Neither data nor simulation fit NML
specifies a positive epoch cut: `EPCUT_SNRMIN` is empty. The available
simulation optical window is [-15,+45], whereas the released data uses
[-7,+45]. These are available inputs, not proof of historical execution.

The [native refit design](raisin-signed-refit-design.md) and Astra's source
review add another pathway: iteration 1 can use a broader cadence before
iterations 2/3 enforce the phase window. Off-season rows can therefore change
initialization even if they never enter the final objective. The released
SNooPy.B18 grid is tabulated over phase [-20,+70] and stretch [.7,1.3]. Its
numerical extrapolation is not validated late-time or pre-explosion truth.
The initial simulation must keep its generating phases inside support and
separate this initialization pathway.

## Generative variables and the four minimal paired arms

Let X contain the full identified exposure schedule, band, zero point,
quoted-error convention and premeasurement quality metadata. Let theta
contain the true native SNooPy parameters and any declared population or
intrinsic-scatter latent variables. Generate a complete measurement vector
Y=f(theta,X)+L(theta,X)U, with U a standard-normal vector in the engineering
Gaussian case. More generally generate measurements, reported errors and
flags jointly from the documented extraction model. A known covariance C
has L L^T=C; this notation does not assert that the missing extraction
covariance is known.

Separate truly premeasurement flags from flags derived from the noisy
measurement or a detection. The latter must be regenerated with each draw
when reproducing the extraction pipeline. Holding an observed PHOTFLAG
fixed is only a conditional metadata experiment unless its noise
independence is established. Native cuts on fitted phase, SNR or residuals
are applied after branching, recorded per iteration and included in the
selection mechanism; they are not silently absorbed into Q.

Freeze a metadata quality mask Q and an archived mask M from the provenance
ledger. For realization r define T_r=1(Y_r>0). Using the **same full vector**
of model values, noise, reported errors and latent scatter in every arm:

| Arm | Retained exposure set before native fitting | Purpose |
|---|---|---|
| F | Q | Full signed eligible cadence |
| T | Q intersect T_r | Fresh sign omission for each new noise realization |
| M | Q intersect M | Archived observed-positive cadence, new signed noise |
| MT | Q intersect M intersect T_r | Interaction control: archive plus a fresh sign cut |

The primary paired effect T-F isolates fresh sign omission on a common
generating cadence. M-F measures fixed-cadence loss conditional on that
archive. MT-M measures truncation on the archived cadence. The interaction
(MT-M)-(T-F) need not vanish for a nonlinear native fit with heterogeneous
signals/errors. Equal row counts alone do not make arms equivalent.

Never regenerate until positive, replace negative flux by its absolute value,
replace it by zero, or resample noise/errors after branching. Common seed
labels are insufficient if changing cadence changes the generator's random
number consumption: generate once on the complete schedule, save truth and
noise by stable exposure identity, then subset. Paired objects also share
the same intrinsic-SED draw, calibration draw and timing initialization.
Rows repeated in source products are not automatically independent
exposures; preserve the multiplicity/identity ledger and stop an ambiguous
noise assignment rather than draw independent noise for copied data.

For the author-raw cadence M is the frozen author-positive membership; this
directly corresponds to the observed A/B source convention. Released R and
hybrid H have a different quoted-error bridge. Do not combine these error
conventions in the primary noise experiment or label it a release-exact
simulation. An additional release-versus-author error-convention sensitivity
requires a common physical exposure mapping and separately frozen errors.

An independent [error-bridge review](raisin-cadence-error-bridge-review.md),
added while this toy was running, identifies the exact numerical positive-row
mapping `sigma_product=round(sigma_author*1.086*0.4*ln(10),3)` on 3,133 unique
DES16 author/product matches. The executed optical conversion script remains
unidentified. This permits a precisely labelled error-convention sensitivity;
it does not identify the extraction covariance or establish how unexported
negative rows would have been processed. The primary author-error toy/native
proposal is unchanged. The final provenance check records this appended
source-review change; all frozen executable/configuration/data-ledger hashes
remain unchanged.

## Analytical checks and the appropriate likelihood

For one independent measurement Y~N(f,sigma^2), define s=f/sigma and
lambda(s)=phi(s)/Phi(s). Then

\[
P(Y>0)=\Phi(s),\qquad
E[Y\mid Y>0]=f+\sigma\lambda(s),
\]
\[
\operatorname{Var}(Y\mid Y>0)
=\sigma^2\{1-s\lambda(s)-\lambda(s)^2\}.
\]

For a fixed mask generated from an **independent** old noise realization,
E[Y_new|M_old,theta,X]=f. Loss of rows changes information and can change
nonlinear estimator bias; this equality does not establish unbiased distance
moduli. Nor does it remove the correlation between an archived mask and the
population's true theta. Drawing new population parameters independently of
a cadence selected using the old object's flux is not generally a sample
from the survey's joint population/cadence distribution.

Independence of old and new noise is a declared toy assumption. If standardized
old and new errors have correlation rho, a one-epoch calculation instead gives
E[Y_new|Y_old>0]=f_new+rho sigma_new lambda(s_old). Persistent template or
calibration errors can therefore make the frozen-mask problem different
again. No value of rho is inferred from the release.

There are three distinct observation likelihoods:

1. **Signed measurements available:** use their joint signed measurement
   density. Restoring measured negative values is not censoring.
2. **Positive values with membership conditioned upon:** for independent,
   known-sigma epochs the retained-value density is
   `N(y;f,sigma^2)/Phi(f/sigma)`, y>0. The normalization depends on f. Its
   score is `(y-f)/sigma^2 - lambda(s)/sigma`; its expectation is zero.
   An ordinary Gaussian score on positives has expectation lambda(s)/sigma.
3. **Full eligible schedule and each missing sign known:** if the only
   omission rule is y<=0, use the joint censored likelihood
   \[
   L=\prod_{j:R_j=1} N(y_j;f_j,\sigma_j^2)
     \prod_{j:R_j=0}\Phi(-f_j/\sigma_j).
   \]
   This retains information in the missing counts/identities. It is not
   multiplied by a second conditional-positive normalization. Conditioning
   on the whole retention pattern instead discards that pattern's
   information and yields the appropriate truncated likelihood.

The third formula is unjustified for epochs absent for unknown weather,
quality or processing reasons. Missing positives outside the fitted phase
window are an example requiring a separate mask. Correlated errors require
multivariate Gaussian orthant/conditional probabilities: the missing-vector
CDF is conditional on retained observed values, and a product of univariate
Phi factors is generally wrong. If quoted sigma is itself a noisy,
flux-dependent estimate, the joint density of flux, error and retention is
needed; treating selected estimated errors as fixed truth can bias weighting.
Native model/K-correction/MW error terms likewise must not be mistaken for
independent measurement draws or duplicated in the generating scatter.

This is a selection-likelihood diagnostic, not a request to replace SNANA's
objective silently. The first native experiment deliberately keeps the same
native objective in all arms and measures its behavior under declared
generators. A separately implemented truncation-aware estimator would need
its own derivative/normalization and recovery gates.

For a linear amplitude fit at fixed design and covariance,
`a_hat=(h^T C^-1 y)/(h^T C^-1 h)` is unbiased for an unconstrained amplitude.
Positivity constraints and the distance transformation
`mu=constant-2.5 log10(a)` introduce different boundary/Jensen effects even
in the full signed arm. A native SNooPy parameter named DLMAG must be read
and exported in its own units; the SALT x0 conversion is not applicable.

## Fixed-signal demonstration, executed

The [executable](../../runs/research_2026_09_26/raisin_selection_simulation_design/gaussian_toy.py)
uses Python 3.12.14, NumPy 2.2.6 and SciPy 1.15.3 in `phase2/env-official`.
It draws 300,000 measurements at each fixed SNR 0,.5,1,2,5, then 50,000
paired 64-epoch realizations at f=.5, sigma=1. The archived mask is drawn
once from a separate fixed seed; it retains 39 of 64 epochs. No measured
SN flux sets these values. The estimated quantity is a retained arithmetic
flux mean, without constraints or logarithmic conversion.

| Arm | Mean estimate | Monte Carlo SE of mean | SD across realizations | Analytic mean |
|---|---:|---:|---:|---:|
| F, full signed | .499481 | .000559 | .124949 | .500000 |
| T, fresh positives | 1.008743 | .000471 | .105384 | 1.009160 |
| M, archived mask, signed | .499402 | .000716 | .160129 | .500000 |
| MT, archived mask, fresh positives | 1.009177 | .000607 | .135687 | 1.009160 |

Paired T-F is +.509263 ± .000480 Monte Carlo SE; M-F is -.000079 ±
.000447. The smaller scatter of the truncated mean is not accuracy:
its expectation is displaced by .509160 flux units. Every empty arm would
be reported; none occurs here. The script also computes the exact finite-N
variance using the binomial retained-count distribution conditional on at
least one retained row. The four arm means and all scalar mean/retention
checks pass the frozen six-Monte-Carlo-SE numerical gate (largest scalar
deviation 1.64 SE). This is a simulation validation tolerance, not a
scientific significance threshold.

Numerical integrals verify retained-density normalization and zero expected
truncated/censored scores to 8.9e-16 maximum absolute discrepancy. The naive
positive Gaussian score has exactly the predicted positive expectation.
Execution took .43 seconds. The [result](../../runs/research_2026_09_26/raisin_selection_simulation_design/result.json)
and compressed paired estimates/counts preserve seeds, software versions,
all outcomes and Monte Carlo uncertainties. This constant-signal toy does
not estimate the sign or magnitude of a native distance shift.

## Native pilot to freeze after source/forward gates

No native simulation is run by this review. The proposed pilot selects the
minimum and maximum zHEL of the already frozen ten-DES16 provenance cohort,
with CID tie-breaking: DES16E2clk (.367) and DES16X3cry (.612). Use the
released native SNooPy.B18/KCOR and the native refit settings. Freeze three
truths `(stretch,AV)=(1,0),(.85,.3),(1.15,.3)`, host RV=1.518, shared
header peak, and a reference distance generated with H0=70, OmegaM=.3,
flat LCDM. This distance is an explicit signal-generation convention, not
external validation of that cosmology or an inferred distance. All truth
parameters are independent of the paired observed-fit outcomes.

The first generator uses **author-raw quoted errors fixed as a diagonal
engineering covariance**, with no intrinsic scatter added. This isolates
selection conditional on the available error estimates; it is not the
source SNANA photostatistical/error-fudge model and is not claimed to
reproduce extraction noise. Run 64 full-vector Gaussian draws per truth
and cadence, using all four arms. That makes 1,536 primary fits. Retain
deterministic second starts for all primary fits in the inference-quality
pilot, giving at most 3,072 fits plus noiseless closure cases. At the
collaborator's current native timing of roughly .46 seconds/object after
initialization, or about .96 seconds/object including small-batch overhead,
this is approximately 24–49 CPU-minutes, 12–25 wall minutes at two workers,
before additional I/O. Time the first 32 fits and revise this forecast;
do not launch an open-ended expansion.

At handoff the numerical prerequisite is **not yet closed**. The native
refit team's [240-fit stability check](../../runs/research_2026_09_26/astra_design/raisin_signed_refit/stability-result.json)
finds 20 of 80 cases move DLMAG by more than .01 mag across starts; the
maximum spread is .224064 mag with an unchanged accepted mask, despite all
return flags being zero and all reported covariances positive. Four cases
also change masks. The team is checking grid profiles and row ordering. The proposed
Monte Carlo must wait for that audit: an optimizer's return flag alone
cannot certify a stable distance estimator. The root's independent
off-season signed-noise audit is a useful additional error-scale check;
it does not identify the complete near-peak covariance or source error
generation process. Neither finding changes the frozen toy truths or
membership rule.

The order of gates is:

1. **Freeze inputs and supported exposure sets.** Hash the native executable,
   sources, loaded libraries, model/KCOR, every generated input, seed and
   exposure mapping. Use the author-raw cadence Q within the fixed metadata
   phase interval [-7,+45] for this first pilot. Keep source flags in the
   ledger; apply no new flag rule without a documented native mapping.
   Export reader acceptance and iteration-specific phase/filter masks.
   Do not generate unsupported off-season SN flux or call this restricted
   experiment the complete raw-cadence restoration.
2. **Forward and noiseless closure.** Source-gate the generic native
   `PHASE2_FLUX`/`USRFUN` export for SNooPy coordinates. Verify that the
   forward generator and fitted model agree in physical flux, DLMAG units,
   zero points, MW correction and model errors. A zero-noise nominal fit
   must recover generating DLMAG/AV/stretch within 1e-4 and fixed peak
   exactly; every input row must have an explained acceptance outcome.
   If the data do not identify these parameters, preserve the degeneracy
   instead of declaring optimizer success sufficient. A copied-input
   rerun must meet the tighter existing native zero-change gates.
3. **Paired measurement closure.** Before fitting, require bit-identical
   surviving signal/noise/error values across arms, exact threshold logic,
   independent archived/new noise, and successful negative-flux ingestion
   in F/M. Record complete truth, draws, masks, orders and cut reasons.
4. **Optimization and native-objective closure.** Save all Minuit flags,
   native errors/model terms, minima, boundary/support warnings and accepted
   rows at every iteration. Use the full native FCN including priors;
   exported FITCHI2 excludes priors and cannot alone rank minima. Two starts
   must agree within 1e-4 in DLMAG/AV/
   stretch and 1e-6 relative objective, or the realization is unresolved.
   Never choose a minimum by its shift direction. Flag model-grid boundaries
   and failure rates by arm; no dropping a failed arm from the headline
   counts. Use matched input order; any ordering sensitivity is a separate
   deterministic control.
5. **Timing/initialization sensitivity.** Only after the fixed-peak pilot,
   repeat the same saved draws with floating peak and the same shared
   header initialization and native iteration rule. The subsequent
   [source review](raisin-profile-solver-review.md) confirms the native
   peak prior recenters on the preceding fitted value each iteration;
   do not call it an independent fixed-header Gaussian prior. Export each
   prior/covariance state and moving phase masks. Separately test the
   historical [-15,+45] simulation window using valid native support.
   Restoring broad ITER1 initialization requires a defensible full-season
   signal model, an explicit no-SN/off-season treatment, or a separately
   declared preprocessing restriction. It cannot be achieved by extrapolating
   the SNooPy grid without a gate.

The generating diagonal errors, frozen phase mask and fixed header peak are
conditional restrictions. The real header peak comes from optical/NIR fits;
it is not independent timing information. A survey-level simulation must
reproduce that timing pipeline and pass each arm's optical timing into its
NIR fit. Source photostatistics, intrinsic OIR.J19 scatter and actual
selection can then be added as separate, paired stages after their assets
and implementation are pinned. Passing the conditional pilot does not
authorize interpreting it as those later stages.

## Estimands, selection normalization and the correction boundary

For each fixed truth/cadence report arm-specific success/acceptance
probabilities and the bias and scatter of fitted DLMAG, AV, stretch and peak.
Report paired T-F, M-F, MT-M and the interaction with Monte Carlo SE based
on the paired draws. Native marginal fit errors added in quadrature are not
an uncertainty on these differences. With 64 draws the precision is
`SD(paired difference)/8`; an inconclusive pilot leads to a separately
frozen size calculation, not outcome-based truth or cut tuning. Also report
truth recovery before any distance-bias correction and failure counts in
each arm.

If fit/detection failures occur, two different estimands must remain
explicit: each arm's selected-population mean, and the paired mean conditional
on both arms passing. Their difference includes selection changes, and the
intersection estimate does not recover the unconditional population bias.
Record all simulated attempts, including objects with no usable epochs;
never give failures a made-up distance. Bootstrap/resample at the shared
object/draw level when population randomness is added, preserving within-SN
noise and calibration correlations.

An event-selection correction requires the probability of **the complete
pipeline** passing, including detection, spectroscopic/HST follow-up,
preprocessing, sign cuts, fitted quality/phase cuts and timing feedback. For
population parameters eta and selection S, the selected-data density is
proportional to
\[
 \frac{\int p(d,S=1\mid\theta,X)\,p(\theta\mid\eta,X)\,d\theta}
 {\int P(S=1\mid\theta,X)\,p(\theta\mid\eta,X)\,d\theta}.
\]
If cadence/host variables are themselves selected, their population
distribution belongs in these integrals too. Event and epoch selection
probabilities cannot generally be multiplied as independent factors.
An end-to-end bias simulation can integrate the selection empirically, but
must apply the same complete selection in data and simulations and retain
rejected attempts. Do not additionally divide an already selected empirical
distribution by the same selection probability.

The source DES generator references `DES3YR_SIM_ERRORFUDGES.DAT`,
`SEARCHEFF_SPEC_DES_Moller_G10_v7.DAT`, `SEARCHEFF_PIPELINE_RAISIN.DAT`,
the external pipeline-logic file and an original survey SIMLIB. They are
not present among the pinned files inspected in `raisin_sign_source` or
the current `sources` tree; resolve exact linked versions before source
simulation reproduction. OIR.J19 implementation/version, downstream fit
cuts and historical command overrides also require execution linkage.
Do not replace them with a generic efficiency curve or infer them by matching
the preferred distance answer. Missing original exposure metadata/errors
and unknown template covariance limit a full-cadence physical simulation.

**Acquisition update:** the [subsequent source-asset review](raisin-simulation-assets.md)
has acquired and Git-blob-verified the three named DES error/efficiency files
and author OIR model. The OIR file exactly matches the bundled 2024 file;
available DES trigger logic is also recorded. Their historical execution
linkage remains unclosed. This updates the local acquisition gap above,
without altering the frozen simulation design or substituting a new model.

A cosmology correction therefore needs more than an observed A-to-B/H-to-R
refit or a successful toy. It requires a justified signed measurement layer,
source-matched generating errors/scatter and population, full event/epoch
selection and optical-to-NIR timing propagation, then a newly trained bias
correction tested on independent simulated draws and prespecified plausible
population/noise variants. Its calibration, host-step and covariance
propagation must be recomputed consistently. The ten DES16 objects neither
identify low-redshift CSP/PS1 sign behavior nor a cross-survey redshift
correction. Until those gates close, the defensible output is a conditional
fitted-distance sensitivity and a demonstrated simulation-selection mismatch
in the available configurations, not a new dust or cosmology correction.
