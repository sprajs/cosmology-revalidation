# SALT / SNANA dust audit

**The analysis cannot currently be certified as bias-free.** This audit records verified numerical effects, the available uncertainty and response matrices, and what remains unidentified. It does not claim that a particular alternative dust correction is true or that the cosmological acceleration result has been overturned. The observation-level reanalysis in `phase2/` remains a separate active investigation.

A concrete new finding is that the released simulation population occupies an extinction-law extrapolation where the nominal host-dust screen can **amplify red light**. This is a physical-support problem in the model component, beyond the already documented historical F99 approximation change. A first actual-SNANA paired stress test gives a small direct effect on its selected fixture; the final selection/BBC/cosmology impact remains unquantified.

## Findings supported by executed checks

| Finding | Evidence and numerical scale | Interpretation |
|---|---|---|
| Some simulated dust lies in a negative-extinction regime | All 25 original public Ia mock HEAD files contain 71,946 written simulated objects. For 4,601 (6.395%), the historical law applied to their stored RV/AV gives A(8000 Å)<0. Substituting the exact F99 formula at the same truth gives 4,791 (6.659%). | Counts describe stored simulations after generation/write selection, not real dust or final BBC-selected objects. The source passes the resulting multiplier greater than one into the SED. A full distance bias is not yet established. |
| Broadband amplification is possible within relevant optical coverage | Independent integration of released SALT3.DES5YR and DES-i at z=.1, peak phase, x1=c=0, E=.1, RV=.4 gives A_i=−.021541 mag. | A constructed, reproducible consequence of the extrapolated law. It is not an observed excess in a real SN and does not assign that dust to the sample. |
| A first paired actual-SNANA stress test has a small direct effect | Flooring negative host extinction in a fixed selected fixture changes accepted epochs for 4 of 12 low-RV mocks; maximum fixed-alpha/beta pre-BBC shift is +.000655 mag. The other 8 and all 12 matched controls are unchanged. Baseline/no-op/reference agree exactly; no accepted-epoch mask changes occur. | This ad hoc floor is a stress intervention, not a proposed true law. Small direct changes in this 24-object fixture do not measure population-wide or regenerated-BBC effects. |
| Historical and exact F99 are numerically different | Extracted historical source and modern legacy option agree exactly. Modern exact source and independent `extinction` implementation agree within 5.8e−11 mag per E(B−V). At E=.1, the maximum optical difference is .00409 mag for RV=3.1, .01155 for RV=2, and .01592 for RV=1.5. | Establishes implementation response, not a universal extra distance correction. Dovekie's overall revision also changes calibration, training and selection. |
| Dust perturbations are nearly collinear with fitted SALT parameters | On 64 selected DES observing designs with measurement-plus-model weighting, additional E=.01 at RV=3.1 gives a median local pre-BBC standardized shift +.01117 mag. Median absorbed squared weighted flux norm is 99.9898%. | Good flux fits alone cannot identify the physical dust correction. Nonlinear injections and independent QR/integration checks confirm the local calculation. Population information can constrain more, conditional on its assumptions. |
| Released correction uncertainty overlaps cosmological and standardization directions | After statistical whitening and removing a common intercept, leading MW-law uncertainty-mode overlap with dmu/dOmega_m is .586 in original DES and .515 in Dovekie. Colour-coefficient and bias-scale columns have correlations −.664 and −.688. | These are local geometric overlaps, not posterior parameter correlations or proof that dust caused a measured cosmological shift. |
| Covariance format, grouping and precision matter | Original DES releases systematic-only covariance; Dovekie releases packed total precision. Both verified total covariances are positive definite. All 47 single-systematic products (23 original, 24 Dovekie) were audited. | CALSPEC is also contained in the calibration/SALT group. Removing that overlap closes Dovekie's component sum to 6.46e−7 relative whitened error. Original DES retains an unassigned positive remainder; coarse single-component serialization also matters. |
| SALT residual maps pass a bounded numerical check | All 63,971 tabulated M0/M1 covariance matrices are positive definite. Independent review finds no negative pre-floor model variance or nonpositive integrated flux at the 2,581 response-test epochs. | The presence of defensive variance/flux code is not evidence of an observed failure. This does not establish completeness of the two-surface model or physical validity of dust priors. |
| Several ambiguities are structural | Intrinsic colour, reddening and environmental luminosity admit exact reparameterizations in effective summary models. Arbitrary luminosity drift can offset arbitrary distance-shape changes. | A unique physical correction or finite global bias bound requires extra information or disclosed restrictions. A restrictive prior is not a new observation. |

## Where the evidence and matrices live

- [Resume handoff for future agents](HANDOFF.md): completed work, exact artifacts, verification commands, ownership boundaries and the next scientific tests.
- [SALT model assumptions and primary literature](salt-model-assumptions.md): training, colour-law conventions, covariance, spectral dimensionality, intrinsic/dust/age degeneracy, and discriminating tests.
- [SNANA implementation and population audit](snana-implementation.md): historical/current source traces, MW versus host extinction, PDF interpolation and actual simulated support.
- [Correction matrices and collinearities](correction-matrices.md): full row identities, covariance semantics, signed spectral decompositions, shared modes and local cosmological response operators.
- [Broadband response experiment](flux-response.md) and [independent review](flux-response-review.md): per-object parameter response, covariance and residual-information matrices with exact signs, units and limitations.
- [Physically admissible dust alternatives](physical-dust-alternatives.md): observational support for low RV, limits of arbitrary cutoffs, public optical/NIR discriminators and a selection-aware comparison design.
- [Machine-readable uncertainty register](uncertainty-register.json): twelve assumption channels, their evidence, what is quantified and what remains missing. `global_bias_bound=null` is intentional; unknown uncertainty is not zero uncertainty.
- [Actual simulated support ledger](../../runs/salt_dust_audit/snana_mock_support/results.json) and [broadband support example](../../runs/salt_dust_audit/f99_support/summary.json).
- [Actual-SNANA paired stress results](../../runs/salt_dust_audit/snana_stress/results.json), [per-object shifts](../../runs/salt_dust_audit/snana_stress/fit_response.csv), and [intervention contract](../../runs/salt_dust_audit/snana_stress/contract.json).
- [Independent actual-SNANA stress review](snana-stress-review.md): raw FITRES/FITS verification, independent epoch matching, exact controls and serialization limits.
- [Study plan and implementation amendments](plan.md): ownership, design choices and corrected diagnostic mistakes.
- [Independent synthesis review](synthesis-review.md) and [mock/source cross-review](matrix-cross-review.md): corroboration and limits of the main findings.

There are three different mathematical objects here. A signed response J maps an explicit small perturbation into fitted quantities. A covariance C describes assumed uncertainty and generally does not determine the signs or physical causes of its underlying perturbations. A BBC bias table is a simulator-conditioned mean correction, not a universal inverse matrix separating dust from intrinsic colour. They must not be substituted for each other.

Given a justified joint nuisance covariance S, local propagation is `C_mu=J S Jᵀ`. Shared modes create covariance between different SNe. The released systematic modes already include some such effects; adding them again as independent errors double-counts uncertainty. Neither this propagation nor an adequate chi-square proves that the mean correction is unbiased under omitted populations.

## Work still needed before a physical bias bound

1. Extend the bounded paired-refit test of low-RV support through population fitting and survey selection. Compare explicitly defined admissible dust laws/populations, re-estimating their nuisance distributions instead of presenting an arbitrary RV floor as truth. Preserve shared random draws and quantify simulation noise; regenerate BBC for a cosmological claim.
2. Complete the active independent calibrated-flux baseline, including clipping, flags, repeated rows and covariance issues. The local response calculation fixes the published accepted epochs and cannot settle those problems.
3. Compare intrinsic-colour/dust/host/evolution alternatives with the same selection normalization and held-out observations. Require flux/time/band residual, colour-tail, host, field-depth and redshift-transfer checks, with mis-specified-generator coverage tests.
4. Resolve provenance of signed systematic realizations and the original release's unassigned covariance remainder. The apparent nine-versus-ten calibration-mode discrepancy is explained numerically by nine surface modes plus the grouped CALSPEC contribution. Matrix eigenmodes suffice for uncertainty geometry but not for recovering missing physical perturbation directions.
5. Report inference conditional on tested population/evolution restrictions. Additional independent dust/spectral/NIR/environment information is needed to break the remaining physical degeneracies; flexible distance curves alone do not do so.

These are required scientific follow-ups, not conclusions from a software test suite. The new audit has made the alternative uncertainty record requested by the user concrete, but **end-to-end bias closure remains open**.

## Reproduction

Run from the repository root. The existing official environment is pinned in `phase2/official/inputs/requirements.lock.txt`; source data and toolchain assets are local and excluded from Git by the existing project policy. A Git-only checkout does not contain the input bundle.

```sh
phase2/env-official/bin/python scripts/salt_dust_audit/snana_extinction.py
phase2/env-official/bin/python scripts/salt_dust_audit/snana_population.py
phase2/env-official/bin/python scripts/salt_dust_audit/snana_mock_support.py
phase2/env-official/bin/python scripts/salt_dust_audit/snana_config_audit.py
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/salt_dust_audit/matrix_audit.py
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/salt_dust_audit/matrix_grouping.py
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/salt_dust_audit/matrix_verify.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/salt_dust_audit/flux_response.py
phase2/env-official/bin/python scripts/salt_dust_audit/f99_support.py
OPENBLAS_NUM_THREADS=1 phase2/env-official/bin/python scripts/salt_dust_audit/snana_stress.py prepare
OPENBLAS_NUM_THREADS=1 phase2/env-official/bin/python scripts/salt_dust_audit/snana_stress.py run
phase2/env-official/bin/python scripts/salt_dust_audit/theory_model_grid_audit.py
OPENBLAS_NUM_THREADS=1 phase2/env-official/bin/python scripts/salt_dust_audit/theory_flux_response_review.py
.venv/bin/python scripts/salt_dust_audit/verify_record.py
```

Manifests record input/code/output hashes and matrix row keys. Independent checks cover numerical reconstruction and the stated model calculations; their passing status does not validate every astrophysical assumption. All work remains local, with no input release changed or externally published result.
