# Host ages and brightness residuals

Changing the residual definition changes the fitted age slope. Corrected Pantheon+ residuals, bias-reversed residuals and the original authors’ residuals are different outcomes. In the original R19 age code/archive, an additive versus fractional star-formation parameter convention differs; this is not automatically a defect in the later Chung analysis. Missing original age posteriors and estimator choices still limit exact replication.

The [27 September results](../../docs/age-correction-results.md) add common-sample reconciliation, paired recovery tests, independent environmental comparisons and population-transport checks. They locate new public host data while retaining the unresolved distinction between an age association and a validated additional distance correction.

## Code and scope

Supplement extraction, object matching, covariance regression, latent-age likelihoods, age-model semantics and mock recovery.

The [source directory](code/) retains the original calculations. See the [execution guide](../README.md#execution) before running them: historical imports and relative paths are assembled into an isolated workspace, and some studies require additional data or native software. Source preservation does not certify that an incomplete experiment succeeded.

## Supporting findings

- [Independent age–Hubble-residual investigation](notes/age-signal.md)

## New executed studies

- [Same-sample age and correction reconciliation](notes/age-reconciliation-results.md)
- [Cosmology and calibrated-flux injection recovery](notes/age-recovery-results.md)
- [Independent local environments and TITAN–ZTF comparison](notes/environment-validation-results.md)
- [Host-to-progenitor mapping and transport](notes/population-transport-results.md)

These new entrypoints run directly from the repository and document their input restoration commands; the historical source below `code/age_signal/` retains its original layout.
