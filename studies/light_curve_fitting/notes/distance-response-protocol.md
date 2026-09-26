# Local distance response of shared calibration modes

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

This retrospective diagnostic is specified after the nine-mode residual scores
and before computing the quantities below. It is not a significance test or a
new cosmology fit. Preserve the earlier residual scores and their limitations.

Use the exact nominal epoch covariance and saved nuisance Jacobian for all 43
discovery and 1,020 validation objects. The nuisance coordinates are
`p=(mB,x1,c,t0)`; replace the amplitude derivative by its exact nominal value
`dF/dmB=-0.4 ln(10) F`. For each already released finite mode, retain the
original amplitude weight and nominal flux convention. With `Jw=L^-1 J`, the
first-order compensation at fixed data is

`delta p = -pinv(Jw) L^-1 delta F`.

The standardized response is `[1,.16087,-3.1178,0] delta p`, in magnitudes.
These fixed alpha/beta are the existing reference convention. No host step,
selection, classifier, BBC, covariance response or training likelihood is
refitted. Coupled calibration/SALT finite mean variants already contain their
released retrained surfaces; do not add that mean response twice.

Run the nine-mode calculation first. Run the expanded twelve-mode sensitivity
only after its source mapping and native output gates pass. For CALSPEC divide
its model flux by the independently derived band scale before differencing.
Keep both outputs. The observer extension uses the existing three contrast
columns and unchanged .02 coefficient scale.

The primary descriptive target is the unweighted mean standardized response
of the highest 255 validation redshifts minus the lowest 255, with deterministic
CID tie breaking. Do not optimize bins, weights or target after seeing answers.
Save all per-object, per-mode responses as well. For each calibration-only and
calibration-plus-observer approximation, propagate the specified independent
standard-normal mode prior and the discovery-only shared posterior to this
contrast. Report mean, standard deviation and remaining variance fraction.
This finite-mode conditional uncertainty is not an empirical bound on all
calibration errors or luminosity evolution. Do not use validation residuals to
refit shared coefficients.

Check the arithmetic independently through weighted normal equations and
verify reconstruction into tangent plus orthogonal residual parts. A gray
one-magnitude model-dimming tangent must project to zero residual and produce
`delta mB=-1`, with other coordinate responses zero. This checks sign and
demonstrates the lost information; its large unit amplitude is a derivative,
not a proposed finite one-magnitude change.

Hash all consumed numerical inputs, script, protocol and outputs. Retain
failed checks rather than silently changing tolerances or dropping objects.
