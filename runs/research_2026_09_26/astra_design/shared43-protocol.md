# Coupled released calibration/SALT variations: bounded discovery43 audit

This is a retrospective attribution/budget diagnostic after the observer residual
and its validation outcome, not a new discovery test. Use the unchanged43
high-Ia discovery objects, nominal exact SNANA coordinates, accepted epochs,
nominal frozen covariance, and nominal independent shape/colour/time tangent
with exact amplitude derivative. Evaluate native MODEL000 and MODEL001–009
with LFIXPAR_ALL enabled and exact43 parameter seeds. Verify all four exported
parameters and every epoch are unchanged, and nominal MODEL000 predictions
close to the original nominal release. No residual exclusions.

FITMODEL_NAME selects the full surface plus the native SALT3.INFO calibration
shifts: MAGSHIFT changes primary magnitudes, WAVESHIFT changes filter grids.
Keep automatic shifts enabled, model-specific MAG_OFFSET and wavelength
support, the original KCOR and modern exact-F99 configuration. The signed lists
match all nine alternative INFO files per independent Sol source audit.
Native output must show DES shift application. Do not separately apply a
second copy of those shifts. Nominal per-object C is held fixed for this
attribution experiment, rather than mixing variant noise changes with means.

Project each finite alternative-minus-nominal flux column through the nominal
nuisance projector. These are realized coupled shifts, not derivatives against
independently measured latent amplitudes. Quantify rank, singular values,
projection overlap with the frozen observer direction, and signed correlation
with observed residual. Retain all9 alternatives. No choice based on match.

As an explicitly labelled second-moment stress, form W_j=.3*projected_delta_j
and C_residual=I+W W^T, applying this shared covariance once across all43
objects. Source covariance code scales differences before their outer product;
thus .3 is an amplitude scale and each squared contribution is .09. Do not
renormalize nine weights to one, centre the realized draws, add the released
CALIBplusSALT3 distance covariance, or treat columns as9 independent observed
calibration measurements. Compute the fixed observer matched filter and gain
under this covariance, its information reduction, and conditional nuisance
posterior scales only descriptively. The Gaussian continuous-mode interpretation
is a new approximation, not the original release's full latent prior or executed
epoch likelihood. Report separately the unweighted span geometry.

This does not include CALSPEC, all other systematics, joint upstream covariance
conditioning, calibration/SALT likelihood regeneration, BBC or selection. The
local Fragilistic102x102 zero-point covariance exists; its seed/draw mapping and
wavelength-shift covariance remain additional provenance questions. A failure
of these nine finite directions to cover the residual cannot establish a new
calibration bias or complete systematic-budget failure.
