# Full validation nonlinear constructed-ambiguity gate

## Frozen protocol before full-cohort outcomes

Expand the completed six-object gate to all 1,020 previously frozen validation
objects. Keep the original exact accepted rows, full frozen SNANA covariance,
four nuisance coordinates, SALT3 DES5YR model and passbands, observer vector,
52-column discovery-designed SED coefficients with .05-mag full-grid RMS,
two starts, bounds and convergence gates. Do not tune coefficients or remove
objects after inspection. Preserve the six-object artifacts.

Run native-noiseless and actual-observed targets separately. For noiseless
nominal closure, verify the known exact zero solution and rank-four tangent
instead of spending optimization runs recovering it. Fit the two perturbed
noiseless models and all three observed-target models from both original
starts. The actual SED multiplier stays inside the photon integral, and its
phase follows fitted t0. Signed spectral flux is retained. Official/native
row scaling is not interpreted as a spectrum, and no mean or spectrum is
clamped.

Run the official-mean target only for the six previous metadata-selected
objects plus all eight unique objects containing the ten previously flagged
nonpositive native mean epochs. This sensitivity retains all their rows.
The exact eight-object edge list and the prior six are frozen before execution.
This is a declared implementation sensitivity, not selection for a physical fit.

Freeze the exact high255/low255 membership from the existing twelve-mode
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
