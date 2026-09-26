# Bounded P21/G10 epoch residual mechanism diagnostic

This follow-up is selected after the real-data residual outcome. It asks whether
existing zero-injected-coherent-calibration simulations can generate a comparable
observer-shaped residual through scatter, noise, fitting or measured selection.
It is not an exact real high-Ia null, a population-model ranking, or a cosmology
test. Use P21 and G10 only; no new population tuning or nine-arm scan.

## Fixed metadata cohort and inputs

Use raw HEAD/PHOT/DUMP under phase2/literature/simulations/outputs/
PH2_pilot02_{P21,G10}; saved measured fits/flags under phase2/hierarchy/forward/
{arm}-fitted.csv.gz and {arm}-attempts.csv.gz; original fit NML under
phase2/checkpoint/official-inputs/snana_forward_{p21,g10}.nml. Exact source,
binary, model/KCOR/noise/population/input hashes must be recorded. Generated
attempt index owns the26,518-row denominator; written CID is unique only after
its verified written-stage join. Do not join failed generated attempts by CID.
Counts currently recorded are written2678/3085, fitted1825/2164 and basic
quality1314/1631 for P21/G10. Preserve and independently verify these counts.

Freeze up to256 archived fitted objects per arm before inspecting their epoch
residual scores. Require finite measured x0>0,x1,c,t0 and zHEL in [.05,1.2].
Use first linked PHOT-row FIELD for both simulated coverage and the real target;
do not silently use compound FITRES field strings. Six z bins have boundaries
[.05,.2,.35,.5,.65,.8,1.2]. In each nonempty field/z cell take the first four
objects in SHA256('20260926-sim-residual-v1|arm|CID') order, then fill to256
by global hash order among remaining eligible rows. No quality, classifier,
truth-dust or residual cut beyond this declared selection. Retain archived
basic-quality flags for stage comparisons. Selection code, IDs, metadata,
expected-stage counts and hashes are written under astra_design/simulation_design.

The only real-derived flux direction is the unchanged exact43 coefficient mean
in validation1020/frozen-discovery-coefficients.npz. Its zero-sum griz convention
and sign are preserved; no template amplitude, sign, band shape or prior is
learned from simulation residuals. The real comparison uses the already saved
validation1020 per-object a,I,G statistics, on common transport cells only.

## Native generation/fitting implementation control

Verify the actual generator binary from its manifest. Historical option99 was
the approximate F99 law, whereas current option99 is exact splineF99. Run the
bounded selected cohorts through the output-only native fitter with both:
(1) generation-consistent current option-99 at RV3.1, once source mapping is
verified; (2) current exact99/RV3.1 as the real-analysis implementation control.
Keep full SALT fit parameters free; initialize from archived measured fits,
not truth. Keep original fitted epoch/quality/clipping options, add output and
CID initialization only, and save accepted epoch membership. Do not use
LFIXPAR_ALL for the primary fitter/selection diagnostic. Preserve prior/peak
initialization semantics and compare new fitted parameters, PKMJDINI and
membership with archived fits; no exact original mask claim is possible because
old simulated fits did not save epoch masks. Do not tune cuts to reproduce IDs.
If a native fit fails, retain it in the accounting; do not silently replace its
CID. Require at least95% selected-object export success for an interpretable
pilot comparison; otherwise diagnose the implementation gate first. All score
comparisons use explicitly reported common successful IDs when pairing laws.

Use the existing parameterized exact objective exporter, requiring objective
closure and finite positive covariance. At each exported fit and its actual
accepted epochs compute nominal independent local shape/colour/time derivatives,
exact native amplitude derivative and native-mean observer contrasts. Use the
exported per-object C, rank4 projection and the same known data/model flux-unit
checks as the real analysis. Report per-object epoch counts, clipping/coverage,
chi-square/dimension, a=r dot(Tc), I=||Tc||², G=a-I/2. No refit ofc.
If the generation-consistent independent derivative is not available directly,
use the exact native mean and the established native tangent approximation,
quantify it on a few finite-difference native checks, and label the boundary.
Never interpret a law/implementation mismatch as a population or selection effect.

## Stages, transport and uncertainty

Report full generated→written→archived fit→archived basic-quality counts before
any pilot weighting. Flux residuals do not exist for unwritten attempts; do not
invent them. Within the frozen pilot report all common refit successes and the
frozen archived-basic-quality subset separately; report changes under newly
computed basic quality only as a distinct fit-closure sensitivity. These stage
contrasts do not isolate selection causally, since their populations differ.

Primary transported descriptive amplitude is sum(w*a)/sum(w*I), with fixed
field/z-cell weights proportional to the real validation target cell counts
per pilot object in that cell. Use only cells represented by at least two
successful simulated objects in BOTH P21/G10 for the same stage/law; report
retained real/sim fractions and do not extrapolate unsupported cells. Normalize
weights to sum1 for amplitude and descriptive mean gain. The raw unweighted
score and field contributions remain available. Compute the real comparator
with the same supported cells and cell weighting; do not contrast unmatched
cohorts. Add one declared transport sensitivity splitting SNRMAX1 at15 using
measured archived FITRES values; require the same at-least-two support rule and
report increased support loss. No bandwidth or bins are chosen after scores.

Use1000 LIBID-block bootstrap replicates, resampling complete simulation cadence
blocks, to describe Monte Carlo sampling variability; preserve arm pairing by
shared generated design where possible and state where writes/fits differ.
Freeze transport weights/support in the bootstrap and retain empty-information
replicates as failures, not silently dropped favourable draws. These intervals
are not calibrated physical-null probabilities or independent-seed coverage.
No independent-field sign significance is asserted for shared systematic modes.

## Classifier and physical-support boundaries

The existing simulations are truth-Ia and currently lack reconstructed SNN
probabilities. This primary diagnostic ends at basic measured quality; do not
call it a matched PROB_SNNV19>.999 null. If classifier processing is independently
made executable, freeze a separate amendment before inspecting classifier-
selected residuals. Use the same real measured preprocessing, including
PKMJDINI/data-derived epoch alignment; never substitute SIM_PEAKMJD or other
truth coordinates. Retain the existing real classifier disagreement and missing
membership chain as limitations.

Audit P21's actual dust/scatter-generation path, latent AV/RV support and noise
recipe before interpreting its result physically. Tag known passive-screen
support failures if that law is actually invoked; do not clamp or delete them
to make the generator look physical. A stated support-valid subset may be shown
only as an additional diagnostic with all counts retained. G10 is a contrasting
empirical scatter mechanism, not known physical truth. Existing parameters were
trained on overlapping DES data, and one generated realization is not repeated
coverage. A simulated match proves that the declared mechanism can mimic a
calibration-shaped residual; it does not establish that mechanism or a correction
for the real survey. A simulated mismatch leaves model/selection/noise failures
and source differences open.

## Concrete allocation

Sol executes under runs/research_2026_09_26/simulation_residual_control/ with
scripts/research_2026_09_26/simulation_residual_control.py, preserving all
existing simulations/fits.
Start with8 frozen pilot IDs per arm to validate exporter, source law and runtime,
then execute the fixed256 cohorts and both laws if gates pass. This8-object
engineering gate cannot change cohort/scoring choices. Astra reviews masks,
unit/parameter/noise assumptions, score algebra and conditional interpretation;
parent integrates results and decides any larger classifier/selection rerun.

## Declared P21 support stratification and handoff refinements

Before simulated residual scoring, tabulate SIM_AV/SIM_RV for every frozenP21
pilot object and evaluate the actual historical extinction routine on the fixed
2800–8000 Angstrom grid in10-Angstrom steps, including A8000. Use the already
built, source-verified SnanaExtinction library. Flag negative attenuation below
-1e-12, invalid AV/RV, and values outside the wrapper's validated input range
separately. Save a contribution ledger by those truth-support labels; do not
clamp or mutate truth. Report the nonnegative-on-declared-grid subset as a
secondary sensitivity, with transport support recalculated transparently and
the full cohort result retained. This grid is a necessary optical support check,
not proof of validity over every generated/extrapolated wavelength. Confirm
whether P21 actually invokes the corresponding passive-screen generator path
before attributing a mismatch to that support failure. All frozen G10 AV/RV entries are -9 sentinels. Historical SALT source
(genmag_SALT2.c, host attenuation branch) sets HOSTXT_FRAC=1 unless both
RV and AV are positive. Tag these as no_explicit_host_screen_sentinel, not
negative physical extinction. Do not evaluate an undefined sentinel RV through
the screen formula. The initial generic AV<0 tag ledger is retained with an
explicit initial-sentinel-tag filename; it is superseded before scoring.

The engineering8 IDs are fixed evenly spaced ranks across each sorted256 cohort,
not chosen after scores. Sol's actual owned paths are
scripts/research_2026_09_26/simulation_residual_control.py,
runs/research_2026_09_26/simulation_residual_control/ and
docs/research-2026-09-26/simulation-residual-control.md. These are the exclusive execution output locations. Original generated artifacts and
previous fits remain read-only. A fixed-coordinate exporter may be used after
the new free-parameter fit to validate its accepted-mask/parameter closure;
it must not replace the primary free-fit mechanism control.
