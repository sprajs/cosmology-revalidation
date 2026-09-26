# Independent host-to-progenitor mapping investigation

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

## Assessment

The correction is not determined by an observed host-age slope alone. It requires a transported conditional luminosity model, the population distribution to which that model applies, compatible weighting and age definitions, and the corrections already present in the distances. Both the assertion that a shorter progenitor-age range makes a correction negligible and the assertion that a steeper converted slope necessarily cancels that reduction are incomplete without these conditions.

There is an exact invariance result for a common affine change of age variable, and a broader exact invariance for complete invertible reparameterization of the response and distribution. There is no general invariance of a *single fitted slope times a mean/median evolution* under nonlinear, noisy, population-dependent or selected mappings. Broad delay distributions alone do not destroy invariance when the conditional mean map remains affine and stable. This distinction avoids overclaiming against either W26 or C26.

The independent cosmic-SFH calculation reproduces the scale of Son's correction under a disclosed convention, but not its exact authors' configuration. W26's quoted cosmic-SFH-only 1.5 Gyr median evolution is not reproduced by the explicit B13 empirical formula used here. TITAN substantially improves the wavelength information needed for host SFH inference, but its headline cosmological-bias number is a model-dependent extrapolation using W26's existing slope, not a new TITAN Hubble-residual fit. Its source paper explicitly acknowledges this. Nothing in this branch estimates a posterior probability for acceleration.

## Source baseline and what was read

Full relevant methods and appendices read: published Son Sections 2–3 and correction/cosmology conventions in Section 4; published Wiseman Sections 2–3 and Appendix A; published Chung reply Sections 2–4 including the added cumulative-cut mock test; Murakami v1 Sections 2–6 and Appendices A–C; Park v1 Sections 2–4 and Appendix A; C14 SFH/DTD sections and Appendix A; W22 galaxy/SN/selection methods and bias-correction appendix. References were used for source tracing, not treated as independently reproduced evidence.

Live primary metadata snapshots are preserved with URLs/times/hashes in `runs/mapping/sources-2026-09-20/acquisition.jsonl`. As checked on 20 September: W26 arXiv remains v2 (8 May); C26, Murakami, Park and Son remain arXiv v1. The published C26 article is therefore the newer substantive reply, despite the unchanged preprint. W22 repository remote HEAD equals acquired commit `94665399e8cfc9a7dfc0f92533c80ca52dc30406` (1 October 2025), and the GitHub releases API returned an empty list. The TITAN DR1 page still has no data download; its papers page describes forthcoming host/overview products without supplying the required SFH posterior catalogue. This is a checked access limitation, not proof that no other private or unindexed product exists.

Park's DOI was verified through primary Crossref metadata and yields a publisher version-of-record PDF URL. Retrieval returned HTTP 403, and the INSPIRE arXiv query supplied no record. This branch initially evaluated the full arXiv v1. The age investigator subsequently recovered published HTML and confirmed the relevant methods; its evidence and exact Park reproduction barrier are recorded in `docs/investigations/age-signal.md`. No author was contacted.

## What each mapping actually assumes

| Source | Inference and assumptions | Boundary |
|---|---|---|
| Son 2025 | Combines slopes measured with host population ages; assumes they describe progenitor evolution. Uses cosmic-SFH×DTD median delay evolution and multiplies by 0.030 mag/Gyr, subtracting the resulting positive high-z correction from distances. Figure 2 distinguishes the cosmic stellar mean from progenitor median. | A cosmic-volume mixture is not automatically the selected SN-host population. All galaxies ever hosting SNe does not establish equal contemporary SN rate or survey inclusion. The mean correction of a linear luminosity law uses a mean delay, not generally a median. |
| Wiseman 2026 | Builds W22 galaxy SFHs, host abundance weighting, rate weighting and selection; adopts hard-cut power-law DTD with minimum 40 Myr and slope -1.13. Explicitly distinguishes galaxy, SN-host and SN-progenitor ages. | The SFH/DTD/selection joint model is predictive, not measured progenitor ages. The local slope must also be transformed consistently. W26's mass-step prediction assigns luminosity using *host mass-weighted age* (Section 3.3), which does not test every latent-progenitor model. |
| Chung published 2026 | Uses C14 mass–SFH models to map host ages 1.5–8 Gyr to inferred progenitor ages; finds about 25% steeper slope for C14 and about 4 times steeper for W26 DTD. Reports about 85% of Son's correction after selected W26 evolution is inserted. | A local cross-host compression factor need not equal the redshift-population compression factor. Full mappings, weighting, errors, configurations and covariance are needed to verify the 85% number; none located. Selection is not a common change of units. |
| Park v1 | Uses C14 SFHs and SN-rate×double-Schechter host weighting; draws 100,000 mocks; adopts C14 smooth DTD with alpha=20, tp=.3 Gyr, s=-1; imposes a linear progenitor luminosity law -0.030 mag/Gyr and .05-mag Gaussian noise. | Step recovery is a consistency demonstration conditional on that causal law, not independent discovery of it. Footnote 3 states low-mass SFHs underrepresent old populations and may overestimate the age contrast. Does not include explicit binary population synthesis. |
| Murakami TITAN v1 | BAGPIPES fits UV–MIR host photometry, convolves posterior SFHs with DTD, and averages per-observed-event SPADs equally. Distinguishes per-event delay mixture from distribution of per-host expected delays. Baseline slope -1.07, cutoff40 Myr. | 6,983 selected host histories, no finalized TITAN HR relation. Selection forward model deferred. High-z evolution is extrapolated from low-z hosts or an external CSFH. Section 6 estimates are explicitly scale estimates, not direct replacements for cosmological inference. |

S25/C26 median and W26/TITAN mean curves are not interchangeable labels for one age. Nor are formed-mass and surviving-mass host ages identical. W26's statement that progenitor age is always younger than host age must be read as a population tendency under its model: an individual event can come from the old tail of a younger-mean host; a strict universal inequality also depends on the cutoff and stellar-mass weighting.

## M1: controlled transport calculations

Command: `uv run --frozen python scripts/mapping/controlled.py`. Outputs and manifest: `runs/mapping/controlled-mapping.csv`, `controlled-mapping-manifest.json`, `mean-median-counterexample.json`, `causal-asymmetry-counterexample.json`. These are seeded synthetic populations, 200,000 per redshift, with no observational likelihood.

| Controlled case | Prediction minus true selected mean shift (mag) | Interpretation |
|---|---:|---|
| Deterministic common affine mapping, host-slope product | <2e-17 | Exact invariance verified numerically. |
| Stable affine conditional mean, broad latent scatter | -0.0000038 ± 0.0000231 Monte Carlo SE | Latent scatter alone need not spoil transport. |
| Classical host age noise, naive host-slope product | -0.008013 ± 0.0000322 Monte Carlo SE | An affine relabelling preserves this error; it does not repair measurement attenuation. |
| Independent delay imputation, fit on imputed delays | -0.002722 | A random age draw is not an observed progenitor. |
| Nonlinear true delay response, host-slope product | +0.015104 ± 0.0000796 Monte Carlo SE | Fitting the known actual delay recovers the shift, but a linear host proxy does not transport. |
| Host-linear response, nonlinear axis, fitted single delay slope | -0.014808 | The exact nonlinear pushforward remains invariant; replacing it by a single slope does not. |
| Mapping intercept changes with population | +0.024000 | Local slope conversion misses the changing conditional intercept. |
| High-z selection on brightness | +0.084381 ± 0.000426 Monte Carlo SE | Selection changes the residual mean as well as the age mixture; both naive products fail. |

The chosen scales illustrate possible failures; they estimate neither the prevalence nor size of those effects in real surveys. A separate finite-mixture calculation yields unchanged median delay but +0.072 mag mean-luminosity evolution. Park-style correction asymmetry also arises under a common latent driver with no age->luminosity arrow: age correction reduces a -0.07413-mag mass step to -0.00115 mag, while mass-step correction leaves the age slope at -0.01481 mag/Gyr. This is a constructive limit on causal interpretation of the asymmetry.

## M2: independent CSFH×DTD reconstruction

Command: `uv run --frozen python scripts/mapping/csfh_dtd.py`. The code implements B13 Eq. F1/Table 6 from the newly preserved primary paper, and MD14 Eq. 15 as a separate sensitivity. It computes the radiation-free flat-LCDM/CPL cosmic clock, then integrates `CSFH[t(z)-T]*DTD(T)`. Smooth C14 and hard-cut DTDs remain separate. Cosmic-volume means, medians, formed-mass stellar means, surviving-mass stellar means, and their corrections are output at z=0..2.5 in steps of .01. There is no survey selection model.

At z=1, with B13's empirical fit:

| Cosmology | C14 mean evolution (Gyr) | C14 median evolution (Gyr) | C14 median correction (mag) | W26 DTD median evolution (Gyr) |
|---|---:|---:|---:|---:|
| LCDM, Om=.3, H0=70 | 4.47354 | 5.29037 | 0.158711 | 2.42080 |
| Son CPL, Om=.353,w0=-.42,wa=-1.75,H0=70 | 4.18231 | 4.96563 | 0.148969 | 2.36213 |
| Same CPL, H0=63.6 from DESI Table 5 | 4.58154 | 5.44003 | 0.163201 | 2.47584 |

Thus Son's stated 5.3-Gyr/0.16-mag scale is approximately reproduced, especially under the conventional LCDM70 clock; this numerical agreement does not establish which clock/configuration generated Son Figure 2. The cited CPL63.6 reconstruction produces 5.44 Gyr. The mean-based correction for that case is 0.137446 mag instead of 0.163201 mag. A median may be scientifically intended, but a mean-luminosity correction cannot use it without justification.

W26 Appendix A quotes approximately 1.5 Gyr from its DTD×B13 cosmic-SFH curve. The identified B13 empirical fit instead gives 2.36–2.48 Gyr here. MD14 sensitivity gives 1.84–1.92 Gyr in the two CPL clocks. The reason has not been resolved: B13 tabulated/modelled history versus its empirical fitting formula, normalization/clock conventions, or other unpublished configurations are possible. It is not a numerical integration error at the tested resolution. This is an **unresolved approximate-reproduction discrepancy**, not evidence of fraud or a demonstrated W26 arithmetic error.

The smooth versus .3-Gyr hard-cut DTD distinction changes this particular mean/median evolution negligibly (about 0.0002 Gyr for means and <0.0001 Gyr for medians). It does not explain the W26 discrepancy. B13's Table 7 gives substantial CSFH systematic scatter, so the tiny numerical integration errors are not astrophysical uncertainties.

Dimensional check: changing CPL H0 from 70 to 63.6 changes the z=1 median correction by 0.014232 mag even though the SN magnitude intercept can absorb the ordinary distance zero point. Central finite differences give `d ln(Delta median delay)/d ln H0 = -0.95397`, and -0.95167 for the mean; fixed DTD time scales prevent exact inverse scaling. A self-consistent likelihood must either recompute the clock/correction or explicitly hold it as an external template with propagated assumptions. This does not imply Son necessarily omitted the dependence; its exact run bundle is unavailable.

Validation: Astropy independent age evaluations agree within 1e-7 Gyr; the flat-LCDM analytic age agrees within the enforced 2e-6-Gyr bound. Doubling delay and cosmic-time grids changes mean ages by <3e-8 Gyr and median ages by <1e-6 Gyr. The generated figure was visually inspected. Manifest: `runs/mapping/csfh-dtd-manifest.json`; tabular curves: `runs/mapping/csfh-dtd-curves.csv`. The empirical CSFH-versus-redshift fitting function is held fixed across these clock choices; the source galaxy luminosity/volume measurements have not been re-inferred in each alternative cosmology, so this is not a complete self-consistent remeasurement of the CSFH.

[Equation-level age curves and fixed-slope corrections](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/mapping/csfh-dtd.png)

The figure uses CPL63.6/B13 and fixed 0.030 mag/Gyr. Smooth C14 and hard-cut300 curves almost overlap. It is a conditional calculation, not a fitted cosmology or empirical confidence band.

## M3: genuine author-product audit

Command: `uv run --frozen python scripts/mapping/hostlib_audit.py`. The supplied compressed W22 library contains 1,510,000 simulated rows with distinct GALIDs, z=0.0005–1.24, host mean ages in Myr and sampled SN delays in Gyr. The inspected generator samples galaxy mass with rate×mass-function probabilities, then draws a delay from its SPAD. Equal row weights are used here; adding `pred_rate_total` as a second weight would be wrong for this product. The exact figure-generation lineage remains unverified.

At supplied z=.0005, the mean host age is 4.38697 Gyr and mean/median delay 3.20091/.86425 Gyr; at z=1 these are 1.36734 and .79617/.31250 Gyr. These are real summaries of the author's simulated product, not a synthetic replacement or a TITAN observation. They illustrate that mean and median evolution can differ greatly. They do **not** reproduce W26's selected curves, since no exact configuration/mask ties this library to those curves.

Data diagnostics found 173 exact floating-point redshift groups (166 at 12-decimal precision) and 39 groups with counts other than 10,000, including fractional low-z values. These remain visible in the output; no silent redshift regrouping or reweighting was used. The full output is `runs/mapping/author-hostlib-summary.csv`; audit and hashes are adjacent JSON files.

Exact execution barriers are concrete. Current public `aura.py` requires a separate HDF host grid and per-host `SN_ages/*` distribution files; configurations point to machine-specific `/media/data3/...` paths. OzDES efficiency files are separately read from those paths and were not located in the supplied tree. The supplied final CSV cannot reconstruct the full per-host posterior distribution or alternate DTDs. Worse, current `aura.py` catches missing delay-file reads and can fall back to a uniform age grid if no distribution remains; running it unmodified with missing inputs would silently change the model. It was therefore inspected, not executed as a supposed reproduction. The current repository precedes W26 and contains no identified W26-specific figure pipeline.

W22 Section 3.4.1 explicitly approximates host redshift success and matches an ad-hoc redshift distribution while omitting additional detection efficiencies; it does not simulate full light curves. W26 may use extensions, but exact extensions/configuration are not supplied. W26's wording about full selection treatment cannot be established solely by pointing to W22. Requesting those products would resolve this limitation; none was requested externally.

## Counterevidence and what would change the assessment

Against a large universal uncorrected age effect: W26's same-object correction comparison and high-mass morphology contrast are empirically relevant; a genuine local progenitor-age contrast together with a small luminosity contrast could strongly constrain a universal linear age law. TITAN's posterior-SFH estimates and alternative SFH tests support shorter delays within its selected low-z population. These are substantive counterevidence, though exact pipeline and sample-transfer validation remain needed.

Against an assertion that standard corrections settle every age issue: observed environmental correlations survive some correction choices; physical intrinsic-population effects and dust effects can coexist; W22 itself explores intrinsic steps plus dust. C26's slope-conversion objection is mathematically valid when its transport conditions hold. Null residual correlation with one noisy proxy or a stable mass contrast does not constrain every common-mode redshift evolution. Disagreement in the explicit B13/W26 reconstruction prevents claiming exact simulation replication.

Strongest missing discriminating products: (1) W26 selected host/age joint distributions by redshift, exact DTD/SFH/mass weights, efficiency grids, masks, seeds and figure scripts; (2) C26 local host-to-progenitor mapping tables/posteriors and the distribution/weighting used for the 85% correction; (3) the 6,983 TITAN host catalogue and SFH posterior draws with metallicity/dust priors, host association probabilities, selected-sample masks, and finalized calibrated HRs; (4) independent UV–MIR/spectroscopic ages and matched dust constraints over the high-z sample; (5) Son's exact correction array, cosmological clock, CSFH implementation, uncertainty covariance and runtime configuration. Release of these could promote approximate/blocked tests to empirical reproductions.

A measured transported conditional relation that predicts held-out residuals across surveys, after matched correction and selection treatment, would support an additional age correction. Conversely, sufficiently precise matched host/progenitor contrasts inconsistent with that relation would exclude it. Neither agent agreement nor success of an injection that assumes the desired age law supplies this evidence.

## Reproduction status

| Target | Status |
|---|---|
| Affine transport identity / controlled violations | Exact analytic identities and numerical synthetic checks completed. |
| Son CSFH×DTD correction scale | Approximate equation-level reconstruction; exact author clock/configuration not recovered. |
| W26 cosmic-SFH-only quoted 1.5-Gyr evolution | Approximate attempt completed; explicit discrepancy retained. |
| W26 selected curves and C26 85% mapped correction | Blocked exact reproduction: required run-specific distributions/configuration missing. |
| W22 supplied host-library summaries | Exact extraction/derived summaries of acquired author product; not W26 figure reproduction. |
| TITAN 6,983-host SFH/age distributions and residual constraint | Blocked independent empirical reproduction: catalogue/posteriors and final HR release absent at checked endpoints. |
| Park causal claim | Source audit and explicit noncausal counterexample completed; exact empirical regression reproduction assigned elsewhere. |

Primary links: [Son](https://doi.org/10.1093/mnras/staf1685), [Wiseman](https://doi.org/10.1093/mnras/stag797), [published Chung reply](https://doi.org/10.1093/mnras/stag1513), [Park preprint](https://arxiv.org/abs/2605.12596), [Murakami](https://arxiv.org/abs/2604.16597), [W22 repository](https://github.com/wisemanp/des_sn_hosts), [TITAN DR1](https://titan-snia.github.io/dr1.html), [B13](https://arxiv.org/abs/1207.6105), [MD14](https://arxiv.org/abs/1403.0007). Original acquisitions were preserved; new sources and outputs are separate. No external communication or publication occurred.
