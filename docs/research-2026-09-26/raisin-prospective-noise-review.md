# Native noise and support review for the prospective timing experiment

This is a source and design review of exact SNANA `v11_04d`, commit
`10ec91297e4482d593cb5d3d055b10d4aa915071`. It executes no native generator or
fitter and reads no new timing-response outcomes. It preserves the failed
archived precision gate. The proposed eight-attempt experiment can test a
conditional timing mechanism; it cannot by itself estimate the published
survey bias or validate the physical population/noise distribution.

## Error-map applicability is narrower than the map's declared defaults

The nominal author input names `DES3YR_SIM_ERRORFUDGES.DAT`. Its eight maps
cover griz in shallow `E1+E2+S1+S2+C1+C2+X1+X2` and deep `C3+X3` groups.
However, the archived `DES_RAISIN.simlib` supplies 23 complete object names
as `FIELD`, for example `DES16E1dcx`. All 120 present band/field combinations,
covering 6,877 rows, fail the native map test. There are no JH maps either.

This follows the actual substring direction in
`sntools_fluxErrModels.c:853–902`: the full observed field must occur within
the map's expanded field list. `snlc_sim.c` copies the SIMLIB field without
extracting `E1`. `get_FLUXERRMODEL` first returns both errors equal to the
input, then returns immediately for a missing map. The result is a no-op,
not a fatal error. No REDCOV is declared. A configured modern error-map file
also suppresses legacy SIMLIB error correction; this particular SIMLIB has
no legacy FLUXERR keys, so that fact adds no demonstrated lost correction.

The [metadata/source audit](../../runs/research_2026_09_26/astra_design/raisin_timing_assets/prospective_noise_review/map-result.json)
and complete [band/field ledger](../../runs/research_2026_09_26/astra_design/raisin_timing_assets/prospective_noise_review/map-matching.csv)
record the source, configuration and hashes. The earlier statement that map
defaults affect both true and reported errors remains correct **conditional
on a match**. It does not describe these fields. This is not evidence that
observed errors need inflation, nor proof that every historical run executed
this exact configuration. The engineering test should retain this no-op,
not silently rename fields to force a correction.

## What the SIMLIB noise inputs actually represent

The exact pinned author `raisin_cosmo/genSimlib.py` matches Git blob
`3099613d50dae997be0d22380fa5547a6690f647`, SHA-256
`264c9aa6109726b8200a1f2ce7fe58ba4638b6a488280659a57c8153173560b3`.
This refines, rather than replaces, the existing
[sign/cadence source review](raisin-simulation-sign-accounting.md) and
[cadence/error bridge](../../runs/research_2026_09_26/raisin_cadence_error_bridge_review/independent-cadence-result.json).

For optical epochs, `mkdessimlib` uses the existing photometry's MJD and
filter, finds the **first time-only** match within 0.01 day in a private
`DES_DIFFIMG.SIMLIB`, and copies sky sigma, PSF, zero point and zero-point
error. It writes gain=1 and read noise=0. The visible matching function does
not additionally require source band or field. It does not infer sky sigma
from the photometry's quoted flux error. This describes available code,
not a proven executed match failure.

For JH it reads the first matching private `*.e01/*_sub_masked.fits`, takes
sigma-clipped image scatter, multiplies it by `sqrt(0.13/0.25)*EXPTIME`,
then adds 0.01 in quadrature. It uses filter-specific fixed PSF and zero
point constants, adding `2.5 log10(EXPTIME)` to the zero point, again with
gain=1/read noise=0. The script therefore reuses one selected image's
background estimate for an object's epochs. The actual private image,
upstream DES simulation library, glob expansion and executed command are
not linked by the inspected public products; those exact private filenames
were not found in the local `sources/` and current-run asset inventory.

Consequently, the available source does not establish whether upstream
noise inputs already contained an empirical correction, or whether its
background/extraction approximation reproduces actual photometric
covariance. Applying the eight maps after changing fields would be a new
model intervention, potentially double counting some noise, not an
established repair. Minimum closure inputs are the executed constructor,
its original optical SIMLIB and selected NIR image/EXPTIME records, plus a
row-level mapping and component variance ledger. The existing cadence
review verifies membership; it does not certify that noise operator.

## Generated scatter and quoted errors are coupled but distinct

`snlc_sim.c:22741–22895` constructs variances in photoelectron units. Base
search variance contains source photon signal, optional host photon signal,
sky variance integrated over the noise-equivalent area, and read variance.
Optional zero-point and reference-image terms have separate controls.
Search and extra-fudge Gaussian draws are per epoch; reference-image noise
can share a Gaussian across the same field and band. These are Gaussian
photostatistical draws, not discrete Poisson draws.

Let `p=GAIN*10^[0.4(ZP−27.5)]` be photoelectrons per FLUXCAL, `mu` the true
FLUXCAL signal and `y` its noisy realization. After drawing noise using the
true variance, `gen_fluxNoise_apply` changes the reported variance by

`Delta C_reported,ii = [max(y_i,0) − mu_i] / p_i`.

This is the exact conversion of the source update at lines 23411–23426,
for the ordinary `SMEARFLAG_FLUX=1` branch before other special overrides.
It replaces the true source-shot term by a nonnegative measured-source
term. Therefore quoted errors are random and correlated with their fluxes,
even when the error-map branch is a no-op. All downstream peak and NIR fits
must share the **same generated flux/error rows**, not merely the same RNG
seed or nominal error curve. Different generator options can consume
random draws differently. Saturation, nonfinite values and any special
noise options must be recorded and gated rather than filtered silently.

`SMEARFLAG_FLUX=0` sets noise shifts to zero but retains calculated errors.
This is useful for mean/recovery closure, not for claiming a zero-variance
likelihood. `GENMODEL_ERRSCALE=0` disables a separate draw from template
errors; removing OIR spectral scatter likewise does not automatically
remove the SNooPy model uncertainty used by the fitter. A reduced fit
chi-square below one can therefore arise in a deliberately measurement-
noise-only generator fitted with additional model uncertainty. Neither
forcing chi-square to one nor adding OIR again is justified by that alone.

The engineering ledger should preserve true R8 means, each noise shift,
pre/post reported variance, true component variances, unit factors, the
actual FITS R4 data/error values and generated/kept attempt identities.
The proposed zero MW scatter controls must include both absolute and
fractional sigmas: generated extinction uses `MWEBV_SMEAR`, while the
reported header uses the map value. Matching MW law 94, RV and calibration
files is necessary but not sufficient for mean closure because generator
and fitter have distinct rest-template/K-correction call paths.

## Why the fixed-weight timing expansion needs a qualification

The [previous derivation](raisin-timing-quantization-design.md) remains a
valid conditional fixed-weight calculation. It is not the full expectation
of this native estimator. If `y=A h+n`, `W0=C0^-1` and
`delta W=−W0 delta C W0`, the usual timing slope
`a1=(h1^T W0 h)/(h^T W0 h)` itself changes at first order by

`delta a1 = [h1^T delta W h − a1 h^T delta W h] / (h^T W0 h)`.

For small timing error `eta`, positive-flux/interior solutions and a
covariance perturbation linear in noise, the paired distance expansion
acquires the additional second-order term `K eta delta a1`, where
`K=2.5/ln(10)`. Its expectation contains another covariance between timing
error and the realized weight perturbation. An independent Gaussian timing
injection matched only to the observed timing standard deviation misses
this term as well as the previous mean, curvature and same-photon noise
correlations. A deterministic synthetic check gives cubic truncation
remainders with halving ratios 7.86, 7.94 and 8.06; no SN outcomes enter it.

This extension still freezes the other covariance rules. Native template
variance and per-iteration covariance feedback need the actual state
ledger; they are not proved by this Taylor calculation. For each realized
fixed covariance, exact amplitude derivative formulae can be used as a
local diagnostic, while the paired native test measures the declared
algorithm's total response.

## Phase support must cover trial evaluations

The available fit template uses a broad epoch window to keep the same rows
when peak changes. That does not constrain the model evaluation domain.
`snlc_fit.car:8060–8071` initially gives peak parameter bounds of ±100 days;
the ±4-day grid only initializes the fit. The 5-day peak prior is soft and
recentered between iterations. Setting its width to zero divides by zero
in `FITINI_PKMJDPRIOR`, and is not a free-peak switch.

`INDEX_GRIDGEN` clamps a cell index, but SNooPy then extrapolates the end
segment with an unrestricted phase ratio. Its table here spans −20 to
+70 days; the experiment's true-epoch range −15 to +45 is a separate
choice. Final `|eta|<=4` or final callback states do not show that the
optimizer, initial grid, or covariance construction avoided extrapolation.

A minimal output-only guard is a count/min/max phase record immediately
after `Trest=Tobs/(1+ZSN)` in `USRFUN`, before template/K-correction calls,
with CID, iteration, observer band and call context. It should distinguish
fit-relevant calls from later plotting, retain any out-of-domain calls,
and stop scientific scoring if the declared support gate fails. The C
`genmag_snoopy` entry can independently cover actual rest-template calls,
but by itself loses the observer/CID context and does not guard K-correction
support. An unchanged ordinary/instrumented replay gate is still required.
No new constrained optimizer is required just to diagnose the first eight.

## Minimum interpretation after engineering gates

If exact noiseless means, reader units/R4 fields, attempts, masks, support,
multiple starts and per-iteration state gates pass, the proposed fixed
population/cadence experiment can isolate how a joint estimated peak
changes the NIR estimator using the same photons. Joint timing includes
NIR noise and is not an independent optical clock. A fixed-eight attempted
cohort must retain failed/retried draws; success-conditioned replacement
would alter the estimand. No selection cut, population distribution,
intrinsic-scatter assignment or calibration correction is established by
this conditional engineering pass.

All quoted source excerpts, construction blob checks and synthetic algebra
are retained in [source evidence](../../runs/research_2026_09_26/astra_design/raisin_timing_assets/prospective_noise_review/source-evidence.json),
[executable review](../../runs/research_2026_09_26/astra_design/raisin_timing_assets/prospective_noise_review/source_audit.py)
and [algebra output](../../runs/research_2026_09_26/astra_design/raisin_timing_assets/prospective_noise_review/weight-timing-algebra.json).
