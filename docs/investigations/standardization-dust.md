# Independent standardization and dust investigation

Initial independent report recorded 2026-09-20, before reading other investigator conclusions. Experiments STD-01/02/03 were registered in `docs/experiments/standardization-plan.md` before inspecting their numerical outcomes. This is an analysis of released corrected distances plus explicit analytic dust counterexamples, not a raw-photometry reanalysis.

## Findings

The claim that existing corrections are perfectly adequate is too strong: released analyses contain a real dust-law implementation correction, a calibration-systematic weighting correction, and a newly identified low-redshift host-mass revision. Dust-only models have documented failures to reproduce observed residual scatter. Conversely, the evidence examined here does not justify applying the full historical -0.030 mag/Gyr age slope again to already corrected distances. Host effects are present in simulation bias corrections even when the separately fitted residual mass step is nearly zero. A disagreement between galaxy attenuation RV and SN sightline extinction RV is not, by itself, a contradiction.

The most consequential new source acquisition is **Hoyt et al. (2026) Union3.1 Appendix F and Roy Choudhury 2607.24443v2, dated 15 September 2026**. They were absent from the initial dossier. An independent cubic reconstruction matches the public updated Pantheon+ table to 7.1e-15 mag. This update changes 114 rows (102 unique SNe), including 18 SDSS rows; it changes no calibrator rows. The changed rows become fainter in standardized apparent magnitude by a mean +0.07220 mag, or +0.01154 averaged over all 713 low-z rows. This is a standard host/bias correction revision, not Son's age law. Its cosmological consequences require a separate named sensitivity branch with uncertainty caveats.

## Controlled original-DES versus Dovekie calculation

The native uncalibrated flat-LCDM profile likelihood reproduces the published SN-only matter-density estimates:

| Vector and total covariance | N | Best Omega_m | Delta-chi-square=1 interval |
|---|---:|---:|---:|
| Original DES tag 1.3 |1829|0.351954|[0.335305,0.369003]|
| DES-Dovekie |1820|0.330317|[0.315229,0.345757]|

Dovekie Table 8 quotes 0.330±0.015. Our profiled chi-square is 1631.4205; including its released magnitude-marginalization term `log(1^T P1/2pi)` gives 1640.2732, matching published 1640.3. This explains the apparent discrepancy rather than treating unlike chi-square definitions as a failed reproduction. The likelihood uses an analytic magnitude intercept and 64-point Gauss-Legendre distance integration, checked against independent adaptive quadrature to numerical precision.

Covariance inspection confirms the dangerous format difference: original files supply systematic covariance only; Dovekie files supply packed total inverse covariance. Both reconstructed total matrices pass Cholesky positivity; Dovekie inverse reconstruction residual max 7.2e-15. Its stat-only diagonal reproduces released MUERR to 9.9e-6 mag, consistent with rounded table output. Extremely large maximum variances come from BEAMS downweighted objects and are not automatically corrupt entries.

There are 1718 exact CID/survey matches; their zHD and zHEL agree exactly. The matched mean Dovekie-minus-original distance change is +0.007253 mag for DES, -0.001381 for Foundation and -0.007855 for other low-z. These are descriptive means; the intercept is arbitrary and the inter-release covariance is unavailable. No sigma significance is assigned to the differences.

The matched-sample crossed experiment isolates consequences of changing distance vectors and weighting, conditional on each released covariance:

| Matched distances | Matched covariance | Omega_m |
|---|---|---:|
|original|original|0.352783|
|Dovekie|original|0.338187|
|original|Dovekie|0.367022|
|Dovekie|Dovekie|0.335186|

Thus the change is not a simple calibration offset or a covariance-only rescaling. On the same observations, changing distances under the old covariance shifts Omega_m by -0.01460; changing weights alone on old distances shifts it by +0.01424; under the new covariance the distance shift is -0.03184. The operations interact because correlated uncertainty changes which residual modes matter. Crossed covariances are controlled diagnostics, not four independently calibrated physical releases.

The documentation-based decomposition for matched DES SNe has mean terms -0.01545 from amplitude, -0.01285 from alpha*x1, +0.04174 from -beta*c, -0.000745 from the host step, -0.009395 from subtracting biasCor, +0.003890 from the absolute zero point, and +0.000061 closure remainder. Recalibration and retraining cause correlated changes among these fitted terms. A colour-term shift is not proof of an astrophysical colour change; a bias-term shift is not uniquely dust. The Dovekie paper's own sequential dust-law-only intermediate result goes in the **opposite** direction to total recalibration: F99 correction before Dovekie yields Omega_m=0.369±0.017 (Appendix B), whereas final Dovekie gives 0.330. We have not independently regenerated that intermediate pipeline product.

Rounded README coefficients do not close individual rows exactly. An algebraic fit to Dovekie metadata recovers alpha 0.168643908, beta 3.144043763, gamma 0.033151520 with RMS 0.0000292 mag. Original metadata has five host masses rounded to exactly 10.00, losing location within its narrow smooth step and leaving closure errors as large as0.0171 mag. Excluding the preregistered diagnostic neighborhood |logM-10|<=.02 gives original closure RMS 0.000101 mag. Released MU, not reconstructed MU, is used for all fits.

Figure and machine-readable tables: `runs/standardization/des_comparison/matched_distance_changes.png`, `flat_lcdm_fits.csv`, `results.json`; per-object derived table `data/derived/standardization/des_matched_comparison.csv`.

## Physical distinctions and identifiability

For a stellar point source behind a dust screen, extinction is `A_lambda=1.086 tau_lambda`, and `RV=A_V/(A_B-A_V)` characterizes that sightline. For an unresolved extended galaxy, attenuation includes the relative geometry of stars and dust, scattering into/out of the aperture, and luminosity-weighted populations. Host stellar age is an inferred property of the integrated stellar population; SN progenitor delay is the elapsed time from formation to explosion. SALT colour c mixes intrinsic spectrum and dust reddening; beta is not automatically RV or RV+1. Metallicity changes stellar populations and potentially SN physics; it is neither colour nor age.

STD-02 uses an explicit absorption-only uniform mixed slab. Its integrated transmission is `T(tau)=(1-exp(-tau))/tau`; point-source extinction still has the microscopic RV. Numerical integration over source depth agrees with the transmission formula within 1e-12, and the thin limit recovers microscopic RV. With RV=3.1,tauV=.1, galaxy effective RV=3.1348; with microscopic RV=2,tauV=3, galaxy effective RV=3.1449. The **integrated attenuation** trend can therefore increase even as microscopic extinction RV decreases, without selecting unusual dust grains. These are constructed populations, not measured galaxy fits. Real scattering and population geometry add degrees of freedom; they do not justify assuming equality of the two quantities.

An illustrative A_B<.5-mag detection requirement admits only 10.23% of uniformly located SNe in the latter slab, reducing mean selected E(B-V) from0.8143 to0.08333 mag, while preserving sightline RV=2. This exposes selection of transparent paths as a physically possible mechanism. It does not establish the actual DES selection function or its relation to host morphology.

For centred latent quantities, `c=c_int+E` and `m=beta_int*c_int+R_B*E+b_A*A+noise`. If intrinsic colour and reddening are independent, the regression colour coefficient is

`beta_eff=[beta_int Var(c_int)+R_B Var(E)]/[Var(c_int)+Var(E)]`.

Two distinct decompositions in STD-02 give exactly the same c/m covariance and beta_eff=3. Population A has beta_int=2,R_B=4 and equal colour variances .005; population B has beta_int=2.5,R_B=4.5 and variances .0075,.0025, with adjusted grey variance. Equal first two moments are not equal full distributions: skewness, IR information and justified latent distributions can distinguish them.

There is also an exact invariance if intrinsic-colour and dust means are allowed to depend on age: set `E'=E+kA`, `c_int'=c_int-kA`, `b_A'=b_A-(R_B-beta_int)k`. Then observed c and m are unchanged. The explicit bounded example has nonnegative E and shifts b_A from -.03 to -.04 mag/Gyr while leaving c/m unchanged at 1e-16. Identifiability therefore requires constraints beyond a marginal age–residual plot: multi-epoch optical+NIR SN spectra/photometry with an intrinsic-SED model, independently measured local dust/geometry, and selection-calibrated within-mass/within-redshift age variation. Strong parametric priors can identify model parameters conditionally; they are not new observations.

## Strongest evidence against both simple accounts

Against complete adequacy of existing standard corrections:

- Dovekie §5 documents the F99 polynomial approximation error and calibration weights summing to 0.81 rather than 1. These were actual implementation errors, not hypothetical astrophysical alternatives.
- Hoyt's homogeneous host analysis identifies a substantial low-z Pantheon+ mass offset; the public correction changes 114 rows. Consequently an older corrected-distance null test should not automatically be presented as the final word.
- Popovic 2024 §6/Table 4 finds all tested dust-only models underestimate residual scatter, particularly in high-mass hosts. An intrinsic age step improves fit, while still failing to describe everything. It discusses limited independent evidence for the required mass-dependent RV, and acknowledges neither its dust variants nor the direct Salim mapping fully explains the observed mass step. Model flexibility and posterior predictive shortcomings matter.
- Dovekie §7.2 says provenance of modified BS21 systematic dust parameters was lost. Its discussion tests the exact-F99 change against old scatter parameters but does not establish that all latent dust/age physics is identified.

Against automatically applying the full extra age correction:

- The near-zero explicit Pantheon+ mass-step coefficient coexists with a substantial host-mass-dependent bias correction. Son's test deleting DES gamma does not remove that second channel. The new low-z correction itself empirically illustrates its scale: low/high bias branches differ by 0.072 mag on affected objects.
- W26's quoted low-z morphology contrast 0.037±0.020 mag, nearly unchanged without mass/bias corrections, is an empirical challenge to the unqualified 0.09–0.18-mag prediction obtained by multiplying global-host age differences by .030. It is not yet independently reproduced here because exact morphology membership/labels are missing.
- Salim 2018 explicitly distinguishes attenuation from extinction; the slab calculation demonstrates why C26's direct RV inconsistency is not logically decisive. Failure of equality cannot establish which SN dust model is correct.
- Corrected residual trends can be small even when physical luminosity depends on age; cosmological validity instead requires that standardization transports correctly across the selected redshift-dependent population. Conversely, no evolving mass step cannot exclude an age evolution component common to both mass groups.

## Reproduction classification and missing inputs

**Exact to release precision:** DES/Dovekie row/covariance parsing and matched comparison; Dovekie published flat-LCDM central value and normalized chi-square; public Pantheon+ low-z cubic branch-transfer table; analytic dust identities.

**Approximate:** README-coefficient individual-distance reconstruction (rounded parameters/masses; measured closure documented). Profile intervals are frequentist delta-chi-square intervals, compared with published posterior summaries; they need not exactly coincide. Dust slab and latent examples are constructive tests only, not reproduction of galaxy observations.

**Blocked as an exact end-to-end reproduction:** F99-only intermediate cosmology without the precise intermediate distance/covariance/configuration; W26 morphology contrast without its exact crossmatch table and labels; full dust/selection refit without paper-specific simulation config/host likelihoods (250 Dovekie simulation lightcurve objects remain unrecovered); private Hoyt corrected vector comparison; a physical age/dust separation without additional observables or defensible restrictive latent assumptions. Required author requests would name these exact products and cuts; no authors were contacted.

New sources preserved with hashes under `sources/updates/2026-09-20-standardization/`. W26 arXiv remains v2, 8 May 2026; Chung reply remains arXiv v1 but published August version has extra tests; Dovekie remains v3, 27 March 2026. Newly acquired BayeSN×Dovekie 2606.19429v1 reports 12% scatter reduction on a likely-Ia, z<.7 sample **without bias corrections** and explicitly calls itself a step towards an end-to-end analysis, not a completed independent cosmological result. The source map must distinguish this promising discriminator from an acceleration measurement.

More precisely, [BayeSN×Dovekie](https://arxiv.org/abs/2606.19429) trains G26 on 1024 SNe, using Dovekie calibration priors and a hybrid distance prior: informative at z<0.08, flat above, to avoid conditioning high-redshift training on cosmology. Section 5 calls DES-SN5YR a validation set independent of its training data. That is **training/test independence as claimed by the authors**, not a new observing sample independent of the DES cosmology discussed here; we have not independently audited all training IDs. The comparison imposes P_Ia>=0.5 and z<0.7, and removes bias corrections from both SALT3 and BayeSN. NMAD scatter is 0.164 versus 0.185 mag. Its inferred selected-training-population mean RV=1.84±0.20 is further evidence that low SN-inferred RV cannot simply be declared impossible from galaxy attenuation, while the paper explicitly warns its selected population is not representative. Model code is public at [bayesn/bayesn](https://github.com/bayesn/bayesn), and photometry is available in the already acquired Pantheon+/DES repositories; a complete BayeSN selection-corrected cosmological likelihood is not released by this paper. The cited ZTF environmental analysis 2605.06799 is a separate lead for a more independent survey test, not reproduced in this bounded branch. Subsequent source/method/revision review is recorded in `docs/investigations/ztf-source-update.md` and the final report; no local ZTF BayeSN fit is claimed.

The primary physical reading supporting the distinctions and adversarial assessment is [Salim et al. 2018](https://arxiv.org/abs/1804.05850), especially its introduction and geometry discussion; [Popovic et al. 2024](https://academic.oup.com/mnras/article/534/3/2263/7778269), especially Sections 3–6 and Tables 4–6; [Pantheon+ methods](https://arxiv.org/abs/2202.04077); [DES-Dovekie methods and appendices](https://arxiv.org/abs/2511.07517); [Hoyt Appendix F](https://arxiv.org/abs/2601.19424); and the [public September reconstruction](https://github.com/shouvikrc/pantheonplus-hostmass-correction-reconstruction). Calculated quantities above are our outputs; author-reported statistics remain labelled.

Run all calculations with the pinned `.venv`:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/standardization/compare_des.py
.venv/bin/python scripts/standardization/dust_identifiability.py
.venv/bin/python scripts/standardization/pantheon_mass_update.py
```

Each run writes a manifest with code and source hashes. Plan amendments retain failed closure diagnostics; source acquisitions remain untouched. An independent audit should particularly challenge matrix ordering/subsetting, correction signs, boundary masses, the newly found revision's uncertainty, and interpretation of covariance-swapping results.
