# 64-object retrospective residual transfer extension

Frozen after the 12-object null/power result, before inspecting the 64-object
residual scores. This is exploratory and conditional on a previously chosen
redshift-ranked sample, published accepted epochs and reference SALT fits.

Use all 64 objects and 2581 epochs already saved in the flux-response NPZ.
Project whitened (observed minus reference model) flux and each candidate
response off the complete local amplitude, stretch, colour and peak-time SALT
Jacobian using SVD. No magnitude distance or cosmology model is fit. This is
a linearized refit; do not call it a nonlinear flux-model comparison.

Observer templates: three signed griz contrasts g-r, i-r, z-r. Rest/phase
templates: archived host-dust R_V=3.1 and R_V=2.0 responses plus the colour
Jacobian multiplied by tanh(rest phase/20). These finite empirical directions
are not complete intrinsic/dust physical families; the two dust columns can
be nearly collinear. Report ranks and response information explicitly.

Compute Gaussian shared-template posterior prediction with identical zero
mean ridge scales .01, .02, .05 (in each template's declared units). The .02
scale is primary; none is an empirical prior. Leave one whole object out when
learning amplitudes, and also leave one whole field out in a separate
predictive variant. This addition was requested before any 64-object score
was calculated. Use both measurement-plus-model C (primary) and
measurement-only C (noise sensitivity), so there are 24 candidate scores in
the declared maximum test. Retain every object's score and no outcome-driven
object exclusions. Field identity comes from the frozen hierarchy row table.

Calibration: 10000 Gaussian samples under the fixed projected design; 10000
whole-object random sign flips of the observed residual vectors; and complete
field sign enumeration if at most ten fields. The same random Gaussian
epoch noise is transformed into each covariance arm's projected coordinates
to preserve paired dependence. If implementation instead cannot pair noise
arms coherently, calibrate each twelve-candidate covariance arm separately and use a
Bonferroni factor two, rather than an invented cross-arm correlation.

Sign flips assume independent centrally symmetric object (or field) residual
blocks; selection/clipping, shared training/calibration, nonlinear fitted
coordinates and correlated field errors can violate that assumption. They
are robustness diagnostics, not assumption-free significance. All null
simulations retain the nuisance projection and template fitting rule.

Report baseline residual adequacy, rest/observer feature overlap and injection
power. Do not convert a preferred residual basis into a physical correction,
luminosity-evolution estimate or new cosmology. A stronger follow-up requires
nonlinear refits, independent covariance validation, matched classifier and
selection, and transfer to genuinely separate wavelength/field information.
