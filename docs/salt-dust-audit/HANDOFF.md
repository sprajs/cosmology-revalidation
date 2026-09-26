# SALT / SNANA dust audit: resume handoff

Saved on 2026-09-21 at the user's request to finalize this audit and preserve progress for future agents. The independent audit deliverable is complete: evidence, uncertainty channels, correction/response matrices and collinearities are recorded. **Bias-free analysis is not established. A finite global physical bias bound remains unidentified.** Continuing the scientific work below is separate from accepting the completed audit record.

## Start here

Read [the report](README.md), [the uncertainty register](uncertainty-register.json), then [physical alternatives and follow-up design](physical-dust-alternatives.md). The register deliberately has `bias_free_established=false` and `global_bias_bound=null`. Do not replace unknown uncertainty with zero or combine alternative models in quadrature.

Workspace: `/home/szymon/Documents/ChatGPT/supernova`. This audit owns only:

- `docs/salt-dust-audit/`: interpretation, reviews, source registries and this handoff.
- `scripts/salt_dust_audit/`: reproducible calculations and checks.
- `runs/salt_dust_audit/`: numerical outputs, row keys, contracts and hash manifests.

All work is saved on disk. It has not been committed, pushed or published. The three owned directories were untracked at handoff; other tracked modifications and untracked directories belong to concurrent work. Do not stage, discard or overwrite them indiscriminately. Existing source data/build assets are ignored, and this audit additionally ignores extracted `.so` files and cloned stress FITS directories. A Git-only copy is not a reproducible input bundle; preserve the local `sources/` and `phase2/official/` assets or reacquire their manifest-pinned inputs.

The other active thread, **Read Codex goal objective**, ID `01a0c0d4-847d-7c51-878b-41b9bc8dec7d`, has been working on the independent calibrated-flux baseline, flags/repeated observations, hierarchy and selection likelihood. Read its latest state before duplicating work; this snapshot is not a live status claim. This audit thread is `01a0c124-5c7e-7bd0-81ed-17afa791de7f`. No message or scientific result was sent externally.

## Stable results worth retaining

1. **Occupied negative-extinction support.** All 25 original public Ia mock HEAD files contain 71,946 written simulated objects. Their stored RV/AV yield historical A(8000 A)<0 for 4,601 objects; substituting exact F99 yields 4,791. Source tracing confirms the multiplier enters the simulated SED without clipping. These are written mock counts, not real-object or final BBC-selected frequencies. Negative attenuation is incompatible with the instantaneous passive-screen interpretation; it does not establish a cosmological bias.
2. **Actual-SNANA bounded response.** A selected 12 low-RV plus 12 matched-control fixture floors negative host extinction as an explicitly artificial stress. Four low-RV fits change; eight low-RV and all controls do not. Maximum fixed-alpha/beta pre-BBC delta-mu is +0.0006548857674 mag. All three fitter runs exit successfully. Baseline/no-op/prior-reference agree exactly in 104 numeric columns. Independent global epoch assignment confirms 929 accepted epochs, 51 affected accepted epochs and zero acceptance-mask changes. Small response in this fixture is not a survey/BBC bound.
3. **Dust can be absorbed by SALT fits.** The 64-object, 2,581-epoch local flux experiment at host RV=3.1 gives median standardized response +0.01117 mag for added E=.01, with median absorbed squared weighted flux norm 99.9898%. This is a response under declared perturbations and fixed observing designs, not an inferred real dust bias. Nonlinear and independent QR/photon-integration checks corroborate it.
4. **Versioned law differences are real.** Historical SNANA option 99 equals current option -99 exactly on the tested grid. Current 99 agrees with an independent exact F99 implementation to 5.81e-11 mag per E. Host extinction also uses the globally named MW law setting; changing it can inadvertently change both host and foreground.
5. **All 47 single-systematic products audited.** Original DES supplies systematic covariance; Dovekie supplies packed total precision. Both total covariances are positive definite. Statistical whitening is necessary because rejected rows have enormous errors. CALSPEC is already contained in the calibration/SALT group: removing that duplicate closes Dovekie's sum to 6.46e-7 relative whitened error and leaves nine surface modes. Original DES retains an unassigned positive rank-one remainder already included in its total. SIGINT_MODEL is a provenance lead, not an identification.
6. **Physical decomposition is not identified by fit quality.** SALT c is empirical, beta is not automatically RV+1, and intrinsic color/dust/environment transformations can leave summary observables unchanged. Unrestricted luminosity drift can mimic distance-shape changes. Recorded covariance modes and geometric overlaps cannot uniquely recover signed physical perturbations or a missing population prior.

## Artifacts and exact numerical checks

| Work | Main outputs / documentation |
|---|---|
| Source and dust-law numerics | `runs/salt_dust_audit/snana_extinction/`; [implementation report](snana-implementation.md) |
| Population and simulated support | `runs/salt_dust_audit/snana_population/`, `snana_mock_support/`, `snana_config/` |
| Full released matrices and cosmological operators | `runs/salt_dust_audit/matrix_audit/`, `matrix_grouping/`; [matrix report](correction-matrices.md) |
| Per-object signed flux/parameter response | `runs/salt_dust_audit/flux_response/response.csv`, `matrices.npz`, `summary.json`; [experiment](flux-response.md), [independent review](flux-response-review.md) |
| Actual SNANA paired refits | `runs/salt_dust_audit/snana_stress/results.json`, `fit_response.csv`, `epoch_acceptance.csv`, `contract.json`, FITRES/LCPLOT outputs; [independent review](snana-stress-review.md) |
| SALT training/model assumptions | [model report](salt-model-assumptions.md), `theory-model-grid-audit.json`, primary-source registries |
| Cross-reviews | [synthesis](synthesis-review.md), [mock/matrix review](matrix-cross-review.md), `runs/salt_dust_audit/matrix_mock_crossreview/results.json` |
| Whole-record integrity | [verification JSON](../../runs/salt_dust_audit/record-verification.json) |

Use the existing environments rather than installing a new one. `phase2/env-official/bin/python` has the pinned SNANA-analysis dependencies including sncosmo 2.12.1 and extinction 0.4.9. Matrix calculations use `.venv/bin/python`. All three subagents completed their reviews; no simulation or numerical process from this audit is left running.

First run the inexpensive integrity check:

```sh
.venv/bin/python scripts/salt_dust_audit/verify_record.py
```

It checks recorded SHA256 values and local documentation links, not physical validity. The final observed result is stored in the linked verification JSON: nine manifests, 277 file hashes and 60 local links, with zero errors. Do not repeat expensive runs just to reproduce an already verified result. [The report's reproduction section](README.md#reproduction) lists every main computation when fresh outputs are needed.

The actual fitter executable is `phase2/official/build/SNANA-current/bin/snlc_fit.exe`. The pinned modern source commit is `886408a4e171896db5eaa97e735a655f50cec2db`; historical source is `2fe0f564361a873860661ff61080b4db9c607edf`. Original DES release tag 1.3 is at `e3493cb3b9505fc1f3f392d887364b50dae21439`; inspected Dovekie release is `c9a4fcafc4cbd19bd750dee47fc76194a45c181f`. Refresh external provenance only when the next question depends on newer material.

## Next scientific work, in order

1. **Coordinate with the existing independent-flux/hierarchy work.** Establish its current baseline and selection normalization before integrating new dust models. Avoid starting a competing likelihood with inconsistent cuts, row identities or accepted epochs.
2. **Run matched, physically admissible host-law alternatives.** Begin with positive power-law screens over a declared wavelength interval, matching integrated B and V attenuation rather than equating nominal RV labels across families. Use the released SED and passbands, shared random draws, explicit host-only changes and a no-op control. Treat extreme slopes as stress cases. Arbitrary RV>=2 or clipped extinction is not established truth; low RV around 1.4 has empirical support.
3. **Refit the population and test held-out optical/NIR information.** Vary intrinsic color, dust, host dependence and evolution jointly with the same selection normalization. Public CSP/RAISIN flux inputs are concrete next resources documented in the alternatives note, but were not fetched/refitted here. The preserved BayeSN archive uses F99 and prior cutoffs; it cannot independently rule out a tail that its prior forbids. Use an NIR-trained model for NIR tests.
4. **Regenerate selection and BBC before claiming cosmological shifts.** Track generated, detected, fitted and selected truth; propagate host-redshift/classification selection, finite-simulation noise and support. Compare matched generators and deliberately mismatched generators. The fixed-flux intervention does not redraw Poisson errors or detection and cannot substitute for this stage.
5. **Close remaining covariance provenance.** Obtain high-precision original components and signed systematic realizations if available. Trace the unassigned original DES rank-one remainder against the exact executed configuration/output inventory. Do not add it again to the total, identify it from an eigenvector alone, or mistake component serialization indefiniteness for an invalid total likelihood.

For the signed-matrix continuation, retain per-variation MU, MUREF, CID+IDSURVEY selections, alpha/beta/gamma, COVOPT group membership, systematic weights, executed dust/calibration/BBC configuration and seeds, or equivalent matched regenerated outputs. These are the missing inputs needed to connect eigenmodes to physical parameters. The matrix reviewer confirmed at handoff that current code hashes match the successful audit/verification/grouping records; no matrix artifact requires regeneration merely to resume.

For a justified joint nuisance covariance S, use the stored signed response J in `C_mu = J S J^T`, retaining cross-object/shared modes. An unknown physical S must not be invented. Released systematics already include some effects; prevent double counting. Report any eventual bound conditional on the explicitly admitted model family and independent constraints.

## Known implementation traps already resolved

- Original HEAD MWEBV values are already scaled by .86; applying it again was a pilot error corrected before authoritative outputs.
- SALT3.DES5YR calibration requires MAG_OFFSET=.27 plus the KCOR primary-magnitude convention in the independent flux calculation.
- Join release rows on CID plus IDSURVEY. To subset a released precision matrix, invert the full matrix, subset covariance, then reinvert.
- Preserve signed spectral modes for materially indefinite serialized components; do not silently project them to positive semidefinite matrices.
- Fitter timestamps are rounded. Repeated-time/filter epochs require measured flux and one-to-one matching; nearest-time alone is ambiguous.
- Stress FITS rewriting changes string padding. Logical metadata are unchanged, baseline/no-op bytes agree, and noise residuals are preserved only up to the documented float32 effects (below 0.000118 sigma).
- Some grid SED nodes are negative, but the bounded integrated-flux/model-variance checks did not fail. Do not turn defensive-code branches into an observed failure claim.
- Twenty public Dovekie PIPPIN configs are byte-identical to original files. This limits execution provenance; it does not prove that final released distances used stale settings.

The appropriate continuation is model discrimination and propagation, not an automatic dust correction. Keep the uncertainty register honest as additional evidence arrives.
