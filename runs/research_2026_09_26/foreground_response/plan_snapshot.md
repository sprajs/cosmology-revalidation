# Assumption-first continuation, 26 September 2026

This continues the existing programme with the user's authorization for sustained
local research. Earlier inspected data and outcomes remain known; this is a
registered reanalysis, not pristine blinding. No external communication is part
of this work. Historical outputs and dirty working-tree files are preserved.

## Order of questions

1. **Identification:** which observable combinations can distinguish distance,
   intrinsic luminosity evolution, dust and selection at all? An arbitrary
   achromatic luminosity function remains an exact SN-distance degeneracy.
2. **Independent observational leverage:** test the expansion shape using the
   compressed BAO observations without a supernova luminosity law, CPL family,
   H0 calibration or present-q extrapolation. Separately identify a phase/band
   discriminator in flux data that can reject a population explanation.
3. **Selection and measurement:** reproduce classification and finish the
   generated-to-selected chain; diagnose fitter/noise and covariance failures.
4. **Numerical and empirical calibration:** require chain movement, convergence,
   synthetic coverage, calibrated null tests and power before interpretation.
5. **Correction attribution and transport:** estimate the required correction
   versus independently constrained responses, test age/dust alternatives on
   common observations, and test held-out wavelength/survey transfer.
6. **Cosmological interpretation:** only then compare adequate correction models,
   with geometry, ruler, propagation and remaining evolution assumptions explicit.

Completed work may support standard cosmology within the tested alternatives,
reveal a reproducible discrepancy, or demonstrate a remaining identification
limit. None proves that every possible approach is issue-free. A numerical or
selection failure is not new physics.

## First work allocation

- Principal investigator: a BAO-only shape-constrained acceleration test and
  conditional reconciliation budget; scientific integration and verification.
- Astra scientific investigator: independent derivation/design challenge and a
  feasible additional observation-level discriminator.
- Sol numerical researcher: diagnose conditional Tripp failure, add a gate before
  scoring, and verify a repaired run without changing scientific priors.
- Sol selection researcher: preregister and reproduce SNNV19 inference, then apply
  matched classification to the nine forward variants only if validation passes.

All new runs capture executed source and input hashes. New diagnostics do not
overwrite historical experiments or retroactively mark E00–E17 complete.

## User's additional calibration requirement

The user explicitly asks to investigate extinction, miscalibration and selection
that can alter the distant-SN luminosity distribution, including all corrections
already embedded in training, fitting, bias simulations and published products.
The correction ledger must distinguish applied signed changes, marginalized
uncertainty and unvalidated population assumptions. A supported new correction
requires matched flux refits and regenerated selection/BBC before production
cosmology. A stress response is not an estimated bias; a null residual test only
constrains directions in which it has demonstrated power.

An immediate additional, externally specified alternative is the CSFD-v2 minus
SFD foreground-map pattern at actual SN coordinates. Propagate 0.86 times the
map difference through the existing 64-object observer-frame dust derivative,
not a constant R_B, and compare the projected residual likelihood at fixed map
amplitudes zero and one. Map reliability masks, field dependence, baseline map
conventions and the sign of refitting response must be explicit. This is a
conditional alternative-map pilot, not a declaration that CSFD is ground truth
or that its errors are independent of the SN survey.

## BAO shape experiment: protocol frozen before new fit outcomes

Use the six anisotropic DESI DR2 pairs at z=0.51, 0.706, 0.934, 1.321,
1.484, 2.33 from the existing pinned 13-row release. Preserve the within-pair
covariance and Ly-alpha ordering by reading quantity labels. Remove the BGS DV
row using a marginal covariance submatrix, not a submatrix of precision. Hold
BGS out of this first linear experiment because DV is nonlinear in DM and DH.
The prior audit has already established the release's block-diagonal assumption.

Assume flat homogeneous/isotropic metric geometry, positive expansion and one
constant comoving BAO ruler. Do not assume a matter density, gravity field
equation, dark-energy equation of state, calibrated ruler, or stellar luminosity.
The released Gaussian compressed likelihood and its upstream extraction remain
inherited assumptions, not raw-galaxy remeasurement.

Let t=ln(1+z), d=DM/rd, and g=(1+z)DH/rd. Then d'=g and
g'=-q*g (derivatives in t). A history with q>=0 has positive nonincreasing g.
With observed endpoints (t_j,d_j,g_j), its closure obeys:

    g_j >= 0,
    g_{j-1} >= g_j,
    (t_j-t_{j-1}) g_j <= d_j-d_{j-1}
        <= (t_j-t_{j-1}) g_{j-1}.

For the first endpoint d_0=0 at t=0 and unmeasured g(0) is unrestricted above
g_1, so only d_1>=t_1*g_1 is imposed. There is no low-z H extrapolation.
Monotone step histories attain interval means within these brackets; continuous
histories can approach their boundary. Thus this is a conservative closure of
the no-acceleration family, allowing arbitrarily rapid transitions. It cannot
identify q(0), a recent narrow-bin q, or a unique expansion history.

Primary statistic: minimum full-covariance quadratic distance from the 12-vector
to that convex cone. Also report radial-only isotonic distance and the six signed
AP contrasts d_j-t_j*g_j with their covariance. The latter are diagnostic, not six
independent uncorrected detections. Report physical fitted endpoints and residuals.

Validation: verify label/order mapping; analytic constant-q limits; feasible
constructed monotone histories; whitened versus original objective; alternative
QP implementation/KKT checks; invariance to common positive rescaling of data
and covariance. Independent Astra review must challenge sufficiency and inference.

Null calibration: for a closed convex cone C, mu in C implies
dist(mu+epsilon,C)<=dist(epsilon,C), pointwise, because C+mu is a subset of C.
Accordingly Gaussian draws at the zero-mean cone vertex give a conservative
composite-null upper-tail calibration (the zero-distance vertex itself is a
mathematical bound, not a physical universe). Also simulate an optimally scaled
coasting boundary and fitted null, clearly separating their pointwise calibration
from the globally conservative bound. Begin with 1,000 draws; increase to 10,000
if a scientific threshold could depend on Monte Carlo resolution. Report binomial
intervals, never zero tail probability or an automatic Wilks sigma conversion.

Test power using independently generated flat LCDM histories at Omega_m=0.2,0.3,
0.4 and a common ruler scale set without selecting a preferred sign. Report the
whole grid. Failure to reject the cone is not evidence for nonacceleration.
Curvature/ruler evolution would change this test's physical interpretation; they
require separate measured-observable equations and cannot be declared excluded.

The E05 SN correction calculation, if pursued, will use this test only as a
conditional expansion constraint. A covariance penalty is not a measured signed
physical correction and cannot bound arbitrary gray luminosity evolution.
