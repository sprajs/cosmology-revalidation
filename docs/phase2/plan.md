# Phase 2: adversarial comparison against public DES observations

Opened 20 September 2026 before new population/cosmology fits. User explicitly requests a deeper independent reanalysis, not another comparison of corrections applied to released distances. Phase 1 at Git `2f19e34` remains immutable evidence and a baseline, not proof that standard corrections are best.

## Scope and order

1. Read primary DES five-year cosmology/methods/photometry/classification/simulation/calibration papers and relevant Pan-STARRS, Pantheon+, population/dust/age critiques, DESI and CMB analyses. Produce a versioned machine-readable map of assumptions, data levels, overlap, critiques and available products.
2. Audit every available DES HEAD/PHOT row, identifiers, pointers, flags, errors, covariances, host information, selection/classification, simulations and calibration assets. Preserve originals and all new provenance. Calibrated scene-model fluxes are lower-level observations, not unprocessed detector pixels.
3. Independently reproduce official DES distances/cosmology and a calibrated-flux SALT3 baseline, attempting the public executable pipeline. Distinguish exact pipeline, equivalent equations and approximate independent implementation. Record numerical discrepancies before fitting alternatives to actual data.
4. Specify generative models; implement an independent hierarchical likelihood, selected-sample normalization and synthetic recovery tests. Only open real-data comparisons after the baseline gate is assessed.
5. Compare explicit correction families using held-out predictive evidence and diagnostics. Infer a flexible expansion history with transparent resolution/curvature/population assumptions. Do not choose a conclusion or call the best-scoring empirical model a uniquely established physical mechanism.
6. Independent Astra adversarial audit attempts to falsify code, selection normalization, evidential ranking and acceleration interpretation. Publish no external material.

Three Astra Ultra child investigators run concurrently alongside the principal investigator (four active slots); later waves handle independent hypotheses and audit. The first wave owns exhaustive data quality, official reproduction, and literature/selection theory respectively. No global agent configuration changes.

## Data layers and safeguards

The primary cosmology product is original DES tag 1.3; Dovekie is a separate calibration/model branch. Existing distance-level agreement is necessary but not sufficient. Independent flux fits must retain the actual bandpasses, zero-point convention, Milky Way law, redshift, phase/wavelength cuts, covariance convention and trained SED model. Summaries of independently fitted fluxes may serve as a computational likelihood approximation, but must be validated against flux-level evaluations; a fit to published SALT summaries alone will be labelled as such.

The original metadata amplitude cross-covariances are printed to few decimals despite small x0, so full-precision covariance cannot be assumed. Reconstruct from independently fitted light curves or a higher-precision public release. Preserve correlations and detect nonpositive/ill-conditioned matrices; no silent diagonal replacement.

For simulations distinguish generated, triggered, measured, successfully fitted and final-selected objects. A selected mock catalogue alone is not a detection-efficiency denominator. It may support *relative* selection normalization if the exact generated latent density and adequate support are known: Z(theta)/Z(reference) = E_selected,reference[p(theta)/p(reference)]. Validate this identity with simulations having a known denominator before use, condition on the same redshift/field/host quantities, inspect effective weights, and do not call it an absolute efficiency reconstruction.

Photometric Ia classification is inferred from the same fluxes. Avoid treating its probability as independent data; use either a declared conditional mixture approximation or a jointly modelled Ia/non-Ia generator. Report contamination sensitivity, high-purity subsets and information lost by those cuts.

## Registered model families

The population engine will model the pre-BBC observables (and validate against calibrated flux), not subtract an extra historical age coefficient from corrected MU. The minimum comparison set is:

- G: latent stretch/colour and absolute-magnitude population with a Tripp relation, intrinsic scatter and host term; explicit measurement covariance.
- G-drift: redshift-dependent stretch/colour distributions and, separately, stretch/colour luminosity coefficients. The two effects must not be conflated.
- D: intrinsic colour plus nonnegative dust reddening, separate intrinsic-colour luminosity coefficient and dust coefficient, host-dependent dust scale; competing host-dependent dust coefficient versus intrinsic host luminosity term.
- D-drift: selected combinations of dust/population drift and intrinsic luminosity drift, not every possible high-dimensional model chosen after outcomes.
- Flexible empirical colour/host model: a mixture or nonlinear colour relation that can challenge the dust decomposition without asserting a unique dust mechanism.
- Age/evolution: where true age data are absent, a stated host proxy or a specified redshift evolution template with free amplitude. Such a fit is not a measured progenitor-age relation. Any externally fixed Son-type coefficient is a separate restricted hypothesis.

Core latent example: x~N(mu_x(z,H),sigma_x), c_intr~N(mu_c(z,H),sigma_c), E>=0 with exponential or stated flexible law; m=mu_distance(z)+M-alpha(z)x+beta_intr c_intr+R_B E+gamma H+delta_M(z)+epsilon; observed colour=c_intr+E. R_B is an effective band coefficient, not automatically a measured monochromatic R_V. Integrate latent Gaussian variables analytically where exact, validate reddening integration against independent quadrature, and retain selected-sample normalization. Host masses/colours carry uncertainty; inferred age, metallicity and host attenuation are not interchangeable.

No model may use the released BBC correction as both a training target and independent evidence of its validity. Official simulation assumptions supply hypotheses and selection machinery, not observations proving those assumptions.

## Registered tests and interpretation criteria

P2-SYN: Recover known latent parameters, population drift and weak q histories from synthetic measured samples, with and without selection and contamination. Check calibrated uncertainty/coverage and test wrong-model injections. Validate exact likelihood integrals and derivative/distance limits independently. An analytically nonidentified transformation must remain nonidentified in the numerical model.

P2-PRED: Use fixed deterministic object-level train/test splits, stratified by survey depth/redshift/host where possible; group duplicate observations. Primary comparison is paired held-out log predictive density on identical observations with uncertainty across objects/field groups. Added parameters must earn out-of-sample predictive improvement. Report magnitude, colour and stretch components and calibration, so a model cannot win only by fitting a nuisance marginal while failing luminosity prediction. Record fold IDs before outcome inspection. Check field/depth-held-out and redshift-held-out transfer separately from interpolation. Do not merge correlated folds as independent experiments.

P2-PPC: Inspect predicted versus observed flux/time/band residuals, colour/stretch distributions, magnitude-colour-host relations, scatter/tails, and trends with redshift, SNR, host surface brightness, phase and detector/field where available. Report failures for every model, including the preferred one. A model's dust/age label does not make it physically true.

P2-SEL: Verify candidate selection likelihoods against generated-plus-selected mocks when accessible. Examine relative-normalization importance-weight ESS and support. Compare complete/high-efficiency regimes, shallow/deep fields and plausible selection perturbations. If candidate support fails, do not report decisive Bayes factors or extrapolated weights.

P2-Q: Compare smooth finite-resolution q(z) or positive-H reconstructions that do not hard-code LCDM. Flat geometry is the primary disclosed assumption; curvature sensitivity where identifiable. Include free intercepts and calibration/population nuisance. Separate q(0), finite low-z interval and existence of an accelerating epoch. BAO/CMB are later labelled additions, not priors used to determine the SN correction's validity. If luminosity drift and distance remain degenerate, report their joint allowed region and which data actually break it.

P2-ROB: Alternate reasonable population priors, dust distributions, calibration modes, covariance treatment, host-error model and contamination prescriptions. Change one source of information at a time on an identical sample. Prior changes are registered before opening their results. Marginal evidence is optional only if prior measure and evidence estimation are credible; held-out predictive evidence is primary.

The baseline gate will state which official quantities reproduce, tolerances and unresolved differences. A failed flux baseline prevents claiming exact raw-flux model discrimination; it does not justify quietly retreating to corrected distances. Continue resolving public-input/toolchain issues and independently test lower-level likelihoods that are feasible.

## Completion standard

Deliver code, pinned environments, source/derived manifests, complete quality audit, official reproduction record, validated independent population/cosmology pipeline, correction-by-correction evidence table and machine-readable literature map. Another researcher must be able to run the documented commands against the recorded public inputs. A ranking must say which observable supports it, how selection/calibration affect it and what remains confounded. The project is not complete after merely showing that different corrections move q; it must attempt empirical model discrimination on the public measurements and establish any remaining identification limit through actual likelihood/simulation tests.
