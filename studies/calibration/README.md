# Calibration and spectral degeneracy

The chosen DES distance contrast has conditional standard deviation 0.00942 mag with systematics-only modes, 0.01207 mag with inherited observer modes, and 0.01236 mag with an isotropic observer prior. Its mean changes too. A precise contrast therefore does not uniquely identify instrumental calibration. Shared spectral, colour and luminosity directions limit transfer from stars or training samples into cosmology.

## Code and scope

Native derivatives, shared-mode design matrices, observed and simulated predictive scores, calibration-star comparisons, spectral response and independent matrix checks.

The [source directory](code/) retains the original calculations. See the [execution guide](../README.md#execution) before running them: historical imports and relative paths are assembled into an isolated workspace, and some studies require additional data or native software. Source preservation does not certify that an incomplete experiment succeeded.

## Supporting findings

- [Calibration, extinction and population-bias investigation](notes/calibration-investigation.md)
- [Dovekie calibration artifact provenance and release linkage](notes/calibration-provenance.md)
- [Read-only verification of calibration and foreground pilots](notes/calibration-verification.md)
- [A low-redshift training prior is not exactly independent of expansion shape](notes/lowz-training-prior.md)
- [Scientific design and independent challenge, 26 September 2026](notes/measurement-identifiability.md)
- [Smooth rest-frame SED identification geometry](notes/sed-identification.md)
- [Six-object nonlinear nuisance-refit gate](notes/sed-nonlinear-refit.md)
- [Full validation nonlinear constructed-ambiguity gate](notes/sed-nonlinear-validation.md)
- [DES colour/shape cut sensitivity in the frozen 1,020-object cohort](notes/sed-quality-cut-probe.md)
- [Frozen DES colour/shape cut probe](notes/sed-quality-cut-protocol.md)
- [Constructed alternatives under the existing shared calibration model](notes/sed-shared-probe.md)
- [Shared-systematic discrimination of constructed mean changes](notes/sed-shared-protocol.md)
- [Shared calibration and SALT3 uncertainty in the DES release](notes/shared-calibration-uncertainty.md)
- [What the residual test leaves uncertain about distance](notes/shared-distance-response.md)
- [Independent review: conditioning shared calibration modes before validation](notes/shared-mode-review.md)
- [Public training-input overlap diagnostic](notes/training-overlap-protocol.md)
- [Public SALT training-input overlap](notes/training-overlap.md)
