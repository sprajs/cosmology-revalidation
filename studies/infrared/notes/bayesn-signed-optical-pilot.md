# Signed optical BayeSN fits with held-out NIR: a feasible conditional pilot

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

Proceed with a new, explicitly conditional ten-object experiment. The earlier
BayeSN stop due to unknown optical signs no longer applies to these ten DES16
objects. Their complete signed author DIFFIMG ancestors are available, the
positive-row lineage is established, and all ten belong to the released
79-object nominal cosmology tables. This permits a measurement-level question
without pretending that the selected population or historical analysis is fully
reproduced: **does the external distance prior change optical dust/SED inference
in ways that predict the measured NIR differently?**

This note freezes a design, metadata masks and algebra checks. It does **not**
report a new measured-flux likelihood, posterior, NIR prediction score or
cosmology fit. The [design protocol](../specifications/experiments/bayesn_signed_pilot/pilot-design-protocol.json)
has SHA256 `544131bf73dcd7de5069a1e25ca10f9cc3fa00e36ebabc0644a610429395d433`.
Existing [BayeSN forward work](../../spectral_models/notes/bayesn-distance-identification.md),
[signed-cohort results](raisin-signed-cohort.md),
[baseline checks](raisin-signed-baseline.md) and failed attempts remain unchanged.
The old statements that sign provenance prevents *every* observed BayeSN fit are
superseded for this source-defined subset. Exact reproduction of the 2024 paper
is still a different task.

## A timing construction that does not reuse NIR information

The released `PEAKMJD` cannot be an independent timing prior: historical timing
can incorporate NIR, and the available native optical refits inherit that
initializer and a moving peak penalty. Instead, every signed optical ancestor
contains `PRIVATE(DES_mjd_trigger)`. Use that recorded optical trigger, not the
released peak date, to define

`tau = (t_peak - t_trigger)/(1+zHEL) ~ Uniform(-10,20)` rest days.

Infer tau jointly with the optical flux likelihood. Retain a fixed row only if
`10 < (MJD-t_trigger)/(1+zHEL) < 30`. Every retained measurement then remains
strictly inside M20's phase interval `(-10,40)` for **every** allowed peak. There
is no parameter-dependent epoch removal, extrapolation, boundary clamping, or
NIR-derived initialization. Four numerical starts at tau = -5, 0, 10 and 15
are fixed independently of the released header. The trigger is a detection
statistic: this is an explicit timing assumption conditional on the selected
sample, not an uninformative physical distribution of discovery-to-peak times.

The [metadata-only inventory](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/bayesn_signed_pilot/frozen-cohort.csv)
and [row masks](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/bayesn_signed_pilot/frozen-mask-metadata.csv)
retain 175 signed optical rows and 36 held-out NIR rows:

| Object | zHEL | Optical r/i/z | Held-out NIR | NIR bands |
|---|---:|---:|---:|---|
| DES16C1cim | .5310 | 24 | 2 | H |
| DES16C3cmy | .5564 | 15 | 3 | H |
| DES16E1dcx | .4530 | 30 | 6 | J,H |
| DES16E2clk | .3670 | 15 | 6 | J,H |
| DES16E2cqq | .4260 | 18 | 2 | H |
| DES16S1agd | .5040 | 15 | 3 | H |
| DES16S1bno | .4700 | 18 | 3 | H |
| DES16S2afz | .4830 | 14 | 3 | H |
| DES16X3cry | .6120 | 14 | 5 | J,H |
| DES16X3zd | .4950 | 12 | 3 | H |

All have three optical bands and at least five distinct optical nights. The
engineering pair is **DES16E2clk and DES16X3cry**, selected by the minimum and
maximum redshift, not dust, SNR or residuals. The ten-object cohort is immutable;
an unsuccessful object stays in the denominator and failure ledger.

The price of this clean timing construction is loss of optical rise data and
some late NIR epochs. Timing and shape may consequently remain weak. A declared
wider-window sensitivity uses tau Uniform(-10,30) and the smaller common mask
`20 < trigger-rest-time < 30`; refit *both* time-prior choices on that identical
smaller mask under each distance arm. Comparing different data masks alone
would not isolate the time-prior effect. If more than 1% posterior probability
lies within one rest day of a timing boundary, report the inference as
timing-limited; do not silently expand the prior or remove the SN. Released
header dates were used only for the preliminary feasibility table, never for
these frozen masks or fitting priors.

## Inputs and physical operator

Use the actual raw signed optical values and quoted errors, including duplicate
occurrences, with no sign, SNR or residual cut. Use the released NIR values and
errors only after optical posterior outputs and diagnostics are frozen. The raw
optical headers contain obsolete host photometric redshifts and `MWEBV=0`;
replace those *metadata* with the spectroscopic zHEL, corrected zHD and MW
reddening already established in the released cohort. Keep the raw measurement
values. The known positive-error roundtrip factor 1.0002429643966138 and decimal
rounding provide a tiny deterministic coordinate sensitivity, not a fitted
noise correction. Do not mix the different DES SMP reduction into this input.

All ten r/i/z and J/H filters have their entire positive tabulated transmission
inside M20's rest 3000–18500 Å range. Exclude g prospectively: its out-of-domain
positive photon weight ranges from 3.93% to 44.23% across these redshifts.
No blue extrapolation or guessed missing UV spectrum is supplied. The exact
released DES KCOR curves, primary spectra and reference magnitudes define the
flux operator. The existing six-SN bridge has **no DES J adapter**, since its two
DES pilot objects had H only. The first metadata attempt stopped on that missing
file; the preserved failure did not evaluate flux scores. Add J directly from
`kcor_DES_NIR.fits`, column `WFC3_IR_F125W-J`, and repeat the reference-integral
and forward-resolution gates. Do not silently borrow PS1 J. The legacy column
name `WFC3_IR_F125W-H` must continue to use its actual F160W transmission.

M20's trained W0, W1, intrinsic residual covariance, M0=-19.5 and gray scatter
sigma0=.088 mag stay fixed. None of these ten names appears in the existing
M20 training-name overlap ledger, but this does not make their calibration or
population independent of M20. The earlier Fisher calculation showed that free
intrinsic SED coordinates can absorb high-redshift dust directions; the trained
prior is a substantive source of identification.

For the small pilot use independent per-object priors theta~N(0,1), the full 42
unit-normal whitened residual-SED coordinates, AV~Exponential(scale .329 mag)
and RV~Uniform(1.2,6). The AV scale is an inherited training-population assumption,
not learned from ten selected SNe. A scale-1-mag sensitivity is declared. Do not
start with a weakly identified selected-sample population hierarchy and call
its hyperparameters a survey dust measurement. This is deliberately **not** the
2024 paper's RV>=.5 population prior or an exact paper reproduction.

The RV range is fixed before posterior outcomes and avoids the previously
identified nonpassive F99 tail. A [pure source-law check](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/bayesn_signed_pilot/design-algebra-result.json)
tests every cubic wavelength stationary point at 4,097 RV values over this
range; minimum A_lambda/AV is .133948. This checks the chosen mathematical
operator, not empirical correctness of F99 or a universal physical RV limit.
Any negative extinction found in implementation is a hard support failure;
never clamp attenuation. At fixed AV=0 the RV derivative must vanish. That does
not require a marginalized finite-noise RV posterior to equal its prior:
integrating uncertain AV can itself induce prior-volume effects.

## The two distance-prior arms and the prediction

Everything except the distance prior is shared. The first arm uses the paper's
flat LCDM relation H0=73.24, Omega_m=.28 at zHD and a declared redshift/150 km/s
uncertainty formula. Exact vectors are saved in `frozen-cohort.csv`; these are
paper-like inputs, not a recovered paper execution. The second arm assigns
mu~Uniform(20,50) mag independently of redshift, with the already declared
Uniform(15,55) and Uniform(25,50) sensitivities. Both retain
D=mu+deltaM, deltaM~N(0,.088²).

For fixed other latents and fixed measurement C, integrate D with its **proper
normalized prior**, including the convolution of uniform mu and Gaussian gray
scatter. Do not substitute a profile amplitude, flat-amplitude prior, improper
flat distance prior or posterior mode. The saved synthetic probability check
agrees between an independent dense D grid and adaptive quadrature to
2.37e-12 in predictive log density, including a negative optical flux value.
This is algebra, not real posterior recovery.

The primary per-SN quantity is

`log p(y_NIR | y_opt, LCDM prior) - log p(y_NIR | y_opt, broad prior)`.

Compute the density of the **joint NIR vector**, retaining the same posterior
D, time, AV, RV, shape and intrinsic SED draw across all its rows. Draw D from its
conditional optical posterior or use the equivalent joint/optical marginal
integral. NIR fluxes cannot enter preprocessing, chains, starts, convergence
choices, masks or hyperparameter choices. A synthetic payload-swap test must
leave the optical inputs and inference exactly unchanged.

Report all ten contributions, posterior predictive distributions and MC errors,
not only a sum or a favourable band. Seven objects have H alone; their tests
mostly concern amplitude and phase evolution. Only the three J/H objects
provide a NIR colour contrast. Do not fit a corrective NIR offset to improve the
same held-out score and then call that prediction validated.

A predictive preference distinguishes observable consequences of the two
**model-plus-prior** hypotheses. It cannot uniquely attribute a discrepancy to
dust rather than M20 intrinsic colour, calibration, timing, noise or selection.
If both arms predict nearly the same NIR despite different AV/RV, that is an
informative demonstration of nonidentification. Neither outcome identifies
arbitrary gray luminosity evolution separately from distance. These objects
were selected using existing observations, so no survey-population significance
or calibrated population Bayes factor follows from the conditional score.

## Noise and calibration controls without manufacturing a correction

Start with fixed diagonal quoted measurement errors and explicit latent M20
SED variation; add no SALT model covariance or released distance covariance.
The signed pre-explosion check supports an approximately correct *aggregate*
noise scale, but its heterogeneity and same-night structure prevent a claim of
independent Gaussian near-peak noise. Quoted errors may depend on measured
flux. A normalized Gaussian is a declared working likelihood here.

Two fixed sensitivities are implementable. First, derive a per-object/band
constant-offset distribution from the disjoint optical baseline
`MJD<trigger-180`, then integrate its uncertainty in the SN fit. This tests a
stationary host/template-offset hypothesis; it is not permission to subtract
a known bias from the data. Second, use rho=.25 within same-object, same-band,
same-integer-MJD blocks, preserving the quoted marginal variances. That is a
labelled covariance stress, not a measured universal rho. No error rescaling to
force chi-square near one is allowed. The main likelihood retains all duplicate
occurrences; an identity-defined duplicate sensitivity must precede scoring if
needed.

Nominal calibration is fixed in the first prediction. An independently declared
sensitivity can integrate the released HST_CAL wavelength mode, amplitude
0.00714×lambda_um mag and weight1, with **one common latent amplitude for the
whole cohort**. Its nominal-unit transformation must be proved; errors/covariance
must transform consistently if data coordinates are changed. Do not independently
redraw the shared offset per SN or add its distance covariance again. This one
mode is not the entire calibration budget. In particular, historical HST CRNL
inclusion remains [unresolved](raisin-crnl-inclusion.md); neither a generic .04-mag
correction nor current 2026 pixel-noise measurements may be imposed on these
historical NIR fluxes.

## Minimal next implementation and gates

1. Build a new standalone signed-pilot adapter in its own directory, reusing the
   isolated BayeSN environment and pinned M20 assets. Separate the optical input
   payload from a sealed NIR-flux payload. Implement the missing J filter and
   dynamic phase, spline and Hsiao interpolation. Preserve all old files.
2. Gate the two metadata endpoint cadences first: eager/JIT predictions, phase
   endpoints and knot crossings, units, gray scaling, RV null at AV=0,
   quadrature, derivatives and likelihood normalization. Reuse the existing
   numerical tolerances in the protocol, including <.1 quoted error for the
   last wavelength-resolution refinement and <1e-6 for log-D integration.
3. Generate the predeclared two-cadence synthetic set: zero dust, moderate dust
   with a residual-SED draw, and a dust/gray-offset case deliberately inconsistent
   with the LCDM prior; two signed-noise seeds per case. Use the exact fixed masks
   and joint optical/NIR latent identity. Noiseless/operator recovery comes first.
   These twelve datasets are an implementation and sensitivity screen, not a
   frequentist coverage certificate. Do not demand that a false narrow distance
   prior recover a deliberately shifted truth.
4. Benchmark 100 gradients and short warmup under two CPUs before launching full
   chains. The forward/algebra cap is 15 minutes, benchmark cap five minutes;
   the initial synthetic and observed stages each have a proposed 60-minute
   two-CPU cap, subject to the measured forecast **before** outcomes. Four chains,
   Rhat<1.01, bulk/tail ESS>=400/200, no divergences, E-BFMI>.3, and predictive
   log-score MCSE<.05 remain required. Do not call capped or unconverged chains
   successful. If runtime is larger, report the benchmark and required budget;
   do not weaken inference gates to reach real data.
5. After review/release and synthetic checks, fit the two engineering objects'
   optical data under both distance priors. Freeze chains and diagnostics before
   opening their NIR payloads. Then extend unchanged to all ten, retaining every
   failure and prior dependence. No population or cosmology fitting is part of
   this experiment.

The current concrete blockers are finite implementation tasks—new J calibration
closure, dynamic timing, proper D integration and recovery—not missing all-survey
selection or the lost exact 2024 executable. Those larger gaps limit interpretation
rather than prohibiting this clearly labelled selected-sample experiment.
