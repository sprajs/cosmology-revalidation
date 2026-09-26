# Calibrated light curves and native fits

The twelve-supernova DES comparison closely recovers mean fluxes and fitted brightness. For CID1896213, the native local peak-time Hessian gives 0.123924 day while the smooth Hessian gives 0.564431 day; native MINOS gives 0.564311 day and profiling supports the larger uncertainty. This is a local Hessian failure, not a failure of the published scalar timing error. Alternative F99-shaped and phase-dependent colour laws give no decisive predictive advantage in the pilot.

## Code and scope

Independent photon integration, accepted-epoch extraction, native SNANA instrumentation patches, objective and covariance exports, profiles, synthetic recovery and held-out predictions.

The [source directory](code/) retains the original calculations. See the [execution guide](../README.md#execution) before running them: historical imports and relative paths are assembled into an isolated workspace, and some studies require additional data or native software. Source preservation does not certify that an incomplete experiment succeeded.

## Supporting findings

- [DES/Dovekie uncertainty and correlation audit](notes/assumption-audit/des.md)
- [Pantheon and host-age assumption audit](notes/assumption-audit/pantheon.md)
- [Generative model and identifiability of supernova population evolution](notes/causal-model.md)
- [Claim-and-evidence ledger](notes/claim-evidence-ledger.md)
- [Distance-product correction and dependence ledger](notes/correction-ledger.md)
- [Using known omitted signs without discarding correlated errors](notes/correlated-sign-likelihood.md)
- [Observation-level DES data audit](notes/data-audit.md)
- [Local distance response of shared calibration modes](notes/distance-response-protocol.md)
- [Independently specified spectral alternatives: acquisition audit](notes/external-spectral-models.md)
- [Conditional Gaussian check of nonlinear nuisance fitting](notes/gaussian-fitting-control.md)
- [Next experiment: independent chromatic prediction, with the gray mode explicit](notes/identification-next-experiment.md)
- [Independent falsification audit](notes/independent-audit.md)
- [Independent engine: validation and deliberate failures](notes/independent-engine-validation.md)
- [Independent calibrated-flux likelihood: 12-SN adversarial pilot](notes/independent-flux.md)
- [Where calibration and bias corrections enter the inference](notes/measurement-to-cosmology-ledger.md)
- [Original DES-SN5YR: executable reconstruction from calibrated fluxes](notes/official-reproduction.md)
- [Pre-explosion SMP noise control: support audit and frozen design](notes/offseason-noise-protocol.md)
- [Access gaps and requests that would make replication exact](notes/open-gaps.md)
- [Experimental preparation, without running the experiment](notes/replication-boundary.md)
- [Observation and measurement assumptions beneath the DES SNANA fit](notes/snana-assumptions.md)
- [Independent arithmetic check of the frozen 1,020-object validation](notes/validation1020-verification.md)
- [ZTF BayeSN environmental analysis: source update (20 September 2026)](notes/ztf-source-update.md)
