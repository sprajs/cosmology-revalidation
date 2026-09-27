# Qualification of the four cosmology results

The final audit distinguishes a completed calculation from a supported scientific result. It reads the four explicitly pinned continuation plans for ΛCDM with fixed SN brightness, CPL with fixed brightness, CPL with a sampled brightness drift, and CPL with a smooth 0.1-mag luminosity prior. It neither evaluates a cosmology nor launches, resumes or stops any process. The [current audit](../results/final-audit.json) is a timestamped snapshot; pending results remain pending.

Before exposing a parent parameter table, the audit calls the existing full `measurement_summary.summarize_run` qualifier. That rechecks current scientific sources, actual likelihood assets and versions, all four chains, convergence diagnostics, the selected native records, target identity, raw importance weights and every declared weight/stability gate. The saved measurement must exactly equal the newly qualified row. A result file without a matching pipeline receipt remains pending, even if the file already says “passed.” A completed pipeline or zero command exit code confers no qualification.

Every child result is checked against its own evidence graph: producer source, output receipt, parent inputs, target/settings, lineage, cache manifests and sealed records. A sibling's valid lineage cannot supply a missing parent binding. Failed gates and numerical exceptions remain separate, named entries. Unknown status strings fail closed. The expansion and luminosity histories become available to the renderer only after their respective qualifications and hash checks pass.

Several successful calculations remain diagnostics:

- Quantile precision describes Monte Carlo sensitivity; it cannot exclude an unvisited mode.
- The SN quadratic check uses **N−1** degrees of freedom after the free magnitude offset. Its tail averages are conditional, using the same observed data; they are not calibrated frequentist probabilities or Gaussian significances. Smooth-luminosity replication redraws coefficients from their stated prior rather than predicting another realization with the same fitted luminosity curve.
- The native precision screen uses **32 predetermined chronological points** and compares two accuracy settings. The final audit independently recomputes its centered total/component variation flags. Passing means that the screen found no large variation at those points, not that the full posterior has been requalified at higher accuracy. A flagged screen remains visible and the conditional lensing continuation is recorded as not run.

Each probe omission and each joint-lensing variant must satisfy its own unchanged overlap gates. The joint-lensing container's `separate_lensing_sensitivities` status is not a qualification. One failed variant cannot disappear behind a passing sibling. These comparisons retain the original untrimmed weights and inherited physical domain. They do not estimate missing cross-probe covariance, provide model evidence, or establish independent-fit tension significance.

Interpretation remains model dependent. **B(z)=0** is fixed in the baseline, so its zero-width curve is imposed rather than measured. The sampled and smooth luminosity alternatives are sensitivity models, not observationally identified progenitor-age corrections. Positive B denotes dimmer standardized SNe at fixed distance. H and H₀ are in km/s/Mpc, the sound horizon in Mpc, luminosity drift in magnitudes, and cosmic age in Gyr. q and j are dimensionless: q<0 means scale-factor acceleration, while j describes the third time derivative of the scale factor, not dq/dt. Pointwise intervals are not simultaneous bands. Flat GR, fixed neutrino physics, the declared analytic CAMB CPL domain, released SN corrections and probe factorization remain explicit assumptions.

Run the read-only audit in the same pinned modern environment as the cosmological targets:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/final_audit.py --output studies/unified_cosmology/results/final-audit.json
.venv/bin/python studies/unified_cosmology/code/final_audit_validate.py
```

The [validation](../results/final-audit-validation.json) contains 24 temporary synthetic pipeline cases, four independent NumPy precision-classification controls and four manually specified report-schema cases. Synthetic parent qualification is explicitly mocked in those tests only; the command-line audit has no bypass. Tests include failed and pending outcomes, wrong targets, missing parent bindings, altered source/cache bytes, incorrect SN degrees of freedom, false precision passes and mixed lensing-variant outcomes. They validate reporting and identity handling, not the still-pending observational posterior.

For downstream use, `final_audit.audit_all()` returns one entry per cohort. A `verified_children` entry for `expansion_history` or `luminosity_history` contains the qualified output path, SHA-256, target identity and parent-correction binding. Consumers must rehash the file before reading it, inspect the cohort's integrity errors and preserve its precision and scientific qualifications. The [design](../code/final-audit-design.json) and [schema](../code/final-audit-schema.json) specify this contract.
