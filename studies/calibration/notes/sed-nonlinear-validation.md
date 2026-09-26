# Full validation nonlinear constructed-ambiguity gate

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

**The full nonlinear gate confirms the constructed ambiguity across all
1,020 validation objects.** On the exact existing high 255/low 255 membership,
the noiseless observer perturbation changes the fixed-reference standardized
contrast by **+.060383 mag**, while the fixed discovery-designed SED
alternative changes it by **-.186184 mag**. Their fitted flux changes remain
close: relative squared mismatch 1.893% and cosine .99064. The observed-flux
diagnostic gives nearly the same contrast separation, **-.246710 mag**.
These are responses to specified hypothetical mean changes with frozen C,
not measured corrections, physical population evidence or cosmology results.

## Frozen protocol before full-cohort outcomes

Expand the completed six-object gate to all 1,020 previously frozen validation
objects. Keep the original exact accepted rows, full frozen SNANA covariance,
four nuisance coordinates, SALT3 DES 5YR model and passbands, observer vector,
52-column discovery-designed SED coefficients with .05-mag full-grid RMS,
two starts, bounds and convergence gates. Do not tune coefficients or remove
objects after inspection. Preserve the six-object artifacts.

Run native-noiseless and actual-observed targets separately. For noiseless
nominal closure, verify the known exact zero solution and rank-four tangent
instead of spending optimization runs recovering it. Fit the two perturbed
noiseless models and all three observed-target models from both original
starts. The actual SED multiplier stays inside the photon integral, and its
phase follows fitted t 0. Signed spectral flux is retained. Official/native
row scaling is not interpreted as a spectrum, and no mean or spectrum is
clamped.

Run the official-mean target only for the six previous metadata-selected
objects plus all eight unique objects containing the ten previously flagged
nonpositive native mean epochs. This sensitivity retains all their rows.
The exact eight-object edge list and the prior six are frozen before execution.
This is a declared implementation sensitivity, not selection for a physical fit.

Freeze the exact high 255/low 255 membership from the existing twelve-mode
distance analysis. Aggregate nonlinear and both tangent responses on this
unchanged membership only after all relevant fit gates pass. Report failures
without selecting them away. The two-start agreement, gradient, half-step
Jacobian, rank and boundary thresholds are identical to the six-object gate.
Also require finite model predictions and complete nonzero passband support
inside the native wavelength and phase domains at each fitted solution.
Record maximum |Delta_m|, fitted phase range, negative spectral support and
nonpositive integrated predictions as flags; negative empirical spectra are
not silently repaired or used as an automatic deletion rule.

Cap execution at two worker processes, each with BLAS/OpenMP threads=1. Write
atomic per-object results, predictions, source/input hashes and timing records.
Support resume only with unchanged frozen source/protocol/inputs. Initial
resource estimate is approximately 3–6 minutes from the six-object timing,
with a longer allowance for difficult objects; update this from completed
objects without changing scientific rules.

Primary products are the nonlinear-versus-tangent fixed-reference responses
and covariance-weighted similarity of the fitted observer/SED mean changes.
After gates pass, the observed-target objective differences can be reported
as **descriptive local fitted-flux scores** only. They are not integrated
model evidence, a measured physical correction or an inferred SED population.
All covariance matrices remain fixed at their exported references. Native
parameter-dependent covariance, retraining, classifier/selection and BBC are
separate gates; no cosmological fit is performed here.

Results follow after execution.

## Full-cohort result and unchanged comparison

All 1,020 objects and 39,606 accepted epochs are retained. The original 255
lowest and 255 highest zHEL membership is copied and independently verified.
Every perturbed noiseless fit and every observed-target fit uses the same
original two starts. The fixed SED coefficients were selected from discovery
response geometry before this run; they were not learned from these observed
fluxes. Its spectral multiplier remains inside the actual passband integral
and follows fitted phase. The observer multiplier uses the same frozen griz
vector with its finite exponential flux response.

| Target | Observer nonlinear contrast, mag | SED nonlinear contrast, mag | SED minus observer, mag |
|---|---:|---:|---:|
| Native noiseless | +.060382741 | -.186183947 | -.246566688 |
| Observed-flux diagnostic | +.060585163 | -.186125094 | -.246710257 |

The target is the same fixed-alpha/beta pre-BBC coordinate as before:
`mu_ref=mB+.16087*x 1-3.1178*c`. Changes compensate a hypothetical change to
the fitted mean while holding the target flux and covariance fixed. These
are not changes to a cosmological posterior or independently inferred SN
luminosities. Achromatic luminosity evolution remains outside the projected
test, even though the implemented spectral basis excludes its explicit
phase-independent gray coordinate.

After nonlinear nuisance fitting, the noiseless observer and SED fitted mean
changes have whitened norms 18.2453 and 17.7627. Their difference has norm 2.51032,
cosine .990637 and squared mismatch relative to the observer response .0189302.
The observed-target equivalents are 18.2593,17.7776,2.51673,.990601 and.0189978.
Thus the close projected resemblance survives actual nonlinear nuisance
refitting, together with the opposite-sign standardized responses. This is
an equivalence example within the specified finite family, not proof of
actual SED evolution or an independently supported prior over that family.

## Size of nonlinear effects

The noiseless observer high-minus-low finite-mean tangent prediction was
+.060888 mag, versus the actual nonlinear +.060383 mag. For SED it was
-.179104 mag, versus actual -.186184 mag. Its infinitesimal tangent contrast,
-.186308 mag, happens to be closer because finite SED and nuisance nonlinear
effects partially cancel; that cancellation is not a general linearity result.

| Target and perturbation | Median absolute nonlinear-minus-finite-tangent, mag | 95th percentile | Maximum |
|---|---:|---:|---:|
| Noiseless observer | .000436 | .000661 | .001162 |
| Noiseless SED | .000565 | .007665 | .028101 |
| Observed observer | .000780 | .003326 | .012104 |
| Observed SED | .001194 | .009437 | .028181 |

The tail is why the earlier six-object gate could not establish the full
contrast. The full result now uses actual nonlinear fits, with the same
fixed covariance and empirical spectral model. Refitting alpha/beta, BBC,
selection or spectral training was not part of this gate.

## Numerical closure and preserved initial failures

The original execution performed 10,284 optimizer runs, plus exact zero-shift
closure/rank checks for all 1,020 nominal noiseless models. It completed in
232.84 seconds using two single-thread CPU workers. Median per-object worker
time was .448 seconds and maximum .849 seconds. Each object has an atomic
checkpoint, full optimizer records, predictions, input hashes and timing.

Three objects initially missed the unchanged 1e-4 Fisher-gradient criterion:
CIDs 1250357,1620051 and 1304183. Their parameters, objectives and failed gates
remain in the primary results; the primary aggregate correctly withheld
headline contrasts. A separately documented numerical refinement reran these
three objects from the original starts with tenfold smaller finite-difference
steps, identical mean models, targets, bounds and covariance. It performed
30 additional optimizer runs, with no relaxed threshold or scientific change.

All three passed the original gates and additional refined/half/quarter-step
gradient checks. Independently differencing the actual objective along all
four Fisher-normalized parameter directions also passed the original 1e-4
threshold at both final step scales. Their maximum objective-direction
gradient norms are 1.05e-5,3.96e-6 and 1.73e-6. The largest change to any
standardized response under refinement is 1.40e-5 mag. Initial failures:3;
resolved:3; remaining:0. No object was excluded.

The resolved full set has maximum Fisher gradient 8.63e-5, Jacobian half-step
relative discrepancy 6.51e-6, two-start standardized-response disagreement
8.15e-6 mag and objective disagreement 1.18e-8. Source phase domains and
nonzero passband wavelength support pass for every saved solution. Independent
array arithmetic reconstructs objectives within 1.78e-14, all standardized
shifts and contrasts exactly, and all prior-six responses exactly. These
are numerical closure checks, not proof of global minima or interval coverage.

The resolved view points to the unchanged primary checkpoints for 1,017
objects and separately retained refined checkpoints for 3 objects. It does
not overwrite the primary failure record. Use the resolved result for the
numerical findings above.

## Model support and bounded official-mean sensitivity

Accepted/fitted rest phases remain inside native SALT's[-20,50]-day domain;
the extreme fitted range is approximately [-17.8874,45.1783]. All nonzero
filter response lies inside the native [2000,11000]-Angstrom rest-wavelength
domain. The constructed SED maximum absolute change over evaluated support
remains approximately .14762 mag.

Native signed spectral flux is preserved. The reference noiseless nominal
means include the ten previously identified nonpositive integrated epochs.
After fitting, observed-target nominal, observer and SED means have 9,8 and 8
nonpositive epochs respectively. Negative spectral cells occur more widely:
782 reference epochs have negative contributions exceeding 1e-6 of absolute
photon weight. No negative values were clamped or rows removed. Complete
support means interpolation stays inside the model domain; it does not make
negative empirical spectral cells physical.

The explicit mean integration matches the previous independent native
archive to relative 9.21e-14 across the full cohort. Native-minus-official
mean discrepancy has median exact-C norm .00517 and maximum .63164; its
largest quoted-error discrepancy is .63163 at an already flagged early epoch.
The official mean is therefore not silently treated as identical to the
independent mean on all objects.

The separate official-mean target was run on 14 declared objects: the prior
six plus all eight objects containing the ten early edge epochs, retaining
every epoch. Relative to each corresponding noiseless-native response, the
largest change is .000791 mag for observer and .000204 mag for SED, both
at CID 1314317. This is a bounded implementation sensitivity, not a third
official-mean fit on all 1,020 objects or certification of an exact SNANA
spectral implementation.

## Descriptive observed-flux scores and remaining boundary

After all numerical gates pass, the summed change in maximized fixed-C
Gaussian log score is+27.17931 for the observer perturbation and+26.03757
for the constructed SED alternative, relative to the nominal model. The
difference is 1.14174. Baseline chi-square is 32,153.22975; observer and SED
values are 32,098.87113 and 32,101.15461. The comparison includes nonlinear
per-object nuisance optimization, but no shared coefficient integration or
calibrated prior probability for either physical interpretation.

These are descriptive local fitted-flux scores. They are not model evidence,
a replacement for the earlier shared twelve-mode predictive calculation,
or evidence that the flexible SED family describes the true population.
The construction intentionally imitates a frozen empirical observer vector;
its modest score difference demonstrates the limitations of assigning a
physical label using this flux pattern alone.

This gate is complete. Parameter-dependent **native covariance feedback**
remains separate: C here is the exact exported reference covariance held
fixed throughout. A new covariance treatment needs an explicit prescription
for the hypothetical SED and its model scatter, without silently replacing
SNANA covariance by a different implementation. Independent spectral/population
constraints, joint training/calibration uncertainty, actual selection and BBC
remain necessary before interpreting a luminosity or cosmology consequence.

Artifacts: [main script](../code/sed_nonlinear_validation.py),
[frozen protocol/configuration](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_nonlinear_validation/freeze.json),
[preserved primary result](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_nonlinear_validation/result.json),
[resolved full result](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_nonlinear_validation/resolved/result.json),
[resolved paired responses](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_nonlinear_validation/resolved/paired-responses.csv),
[support ledger](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_nonlinear_validation/resolved/support-ledger.csv),
[numerical refinement checks](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_nonlinear_validation/refined_steps/refinement-verification.json),
and [independent verification](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_nonlinear_validation/resolved/independent-verification.json).
