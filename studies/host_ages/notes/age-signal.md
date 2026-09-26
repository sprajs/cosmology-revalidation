# Independent age–Hubble-residual investigation

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

Investigator age_signal; started 2026-09-20. This initial record was written before reading other investigators' conclusions. The registered protocol is `docs/experiments/age_signal-plan.md`. This branch analyzes public host-age summaries and distance/light-curve-parameter tables; it does not reconstruct photometry or host SED fitting, or independently fit cosmology.

## Initial source and input findings

Read the published W26 paper (doi:10.1093/mnras/stag797), published C26 reply (doi:10.1093/mnras/stag1513), Chung I including its correction, and relevant Park III methods from the available arXiv text. W26's stated Pantheon release is commit 7fc6805; the preserved snapshot is used. The C26 published reply adds the Figure 2 cumulative-cut and mock test missing from arXiv v1. Live checks on 20 September confirm W26 arXiv v2 remains the latest arXiv version and C26 remains v1; journal text is consequently the relevant C26 source. No later substantive response was identified in the targeted primary-source search; absence from the search is not proof of absence.

Public age-table rows: G11 199, R19 102. **33 SN IDs occur in both age tables**, often with different age estimates from their different photometry. The samples cannot be counted as statistically independent measurements of different explosions. Among Pantheon-matched rows, 30 overlaps produce 226 age rows but only 196 unique SNe. Exact numeric ID joining restricted to SDSS (`IDSURVEY=1`) matches uniquely in the distance table; this avoids collisions with other survey identifiers.

For unique SN IDs and **Hubble-diagram redshift zHD**, .06<z<zcut gives **100,144,176,192,196** for zcut=.20,.25,.30,.35,.42, precisely the published C26 counts. Heliocentric cuts instead give 99,143,174,192,196. Thus the count/cut reconstruction is exact. The choice of G11 or R19 age for overlapping IDs is not recovered by the counts: both are retained as explicit sensitivity variants. The full crosswalk includes nonmatches, original distance-row indices and age provenance.

The release's `biasCor_m_b` is a combined simulation correction including selection and astrophysical/host effects. Reversing it is **not a clean removal of only a scalar mass step**. The main corrected residual is `MU_SH0ES - mu_flatLCDM(zHD,zHEL; Omega_m=.3)` with a fitted intercept. The corresponding reversal adds `biasCor_m_b`; a separate reconstruction from mB+.148*x1-3.112*c is labelled an approximate Tripp variant because author coefficients/settings are not fully provided. The almost-constant difference between those variants is absorbed by the intercept, but small colour/stretch differences remain. Source reconstruction gives coefficients about .147537 and 3.092219 for the released corrected magnitude after reversing bias, which explains that difference.

The published C26 mock description explicitly forces mean HR to agree across redshift bins after assigning luminosity from age. This operation is not implied by the mathematical definition HR=mu_obs−mu_LCDM. It corresponds to removing an imposed redshift-dependent component, which may approximate cosmological fitting under further assumptions but must be justified separately. A controlled paired mock and analytic calculation are registered to test this issue.

## Preliminary calculations, before final numerical verification

Seeded pilot LINMIX Gaussian-age regressions recover Chung I's original-summary results: R19 roughly −.061±.0145 mag/Gyr and G11 (with its young-host redshift correction) roughly −.046±.015. G11 without that redshift correction but with the same revised ages is roughly −.036±.015. These are approximate reproductions because full author configurations and posterior ages remain unavailable. The same modern Pantheon-matched sample gives a much weaker corrected age association; numerical results are being rerun with tighter retained-chain diagnostics. These facts support a real dependence of the reported slope on distance/correction definition rather than a dispute over arithmetic alone.

The 5000-iteration pilot disclosed a reproducibility issue in public LINMIX: its serial constructor ignores the passed seed, and silent extra iteration mode does not increment its max-iteration counter. Source remains unmodified. Subsequent runs use seeded parallel chains and unsilenced output redirected to logs. Unseeded and early pilot products are preserved under `runs/age_signal/unsuccessful/`; they are not the final quantitative record. Independent retained split-chain Rhat is calculated because the upstream diagnostic uses the first halves, while its reported posterior uses the second halves.

## Interpretation limits and strongest counterevidence

Against the working hypothesis that an extra age term of −.03 mag/Gyr necessarily survives modern corrections: the same measured host ages, on modern corrected distances, show a weak residual association. Reintroducing the entire old slope as an independent correction risks counting host/dust/selection effects twice. C26's mock forces out the between-redshift mean, so reproducing its dilution does not by itself establish this happens in the real cosmology pipeline.

Against the hypothesis that modern corrections are completely sufficient: corrected non-detection on 196 objects is not proof of zero age evolution; its slope uncertainty is substantial, ages are noisy and model-dependent, host-age and progenitor delay are different quantities, and empirical standardization could remove a local association while leaving a redshift-dependent component. A full selection-aware joint age–colour–mass–redshift model, with host-age posteriors and independent validation, is not supplied by either these simple regressions or our approximate reconstruction.

What would materially change the conclusion: an author crosswalk and per-SN error ledger fixing the remaining exact-reproduction ambiguity; host-age posterior samples and covariance with mass/dust; a released forward pipeline showing that injected evolution of a known amplitude is hidden by fitted cosmology and survives held-out tests; or an independent age-stratified sample whose within-redshift, correction-controlled residuals predict out-of-sample distance shifts. This branch alone cannot infer present acceleration or overturn a combined cosmological likelihood.

## Verified final numerical results

The final seeded LINMIX fits use four chains, at least 20,000 iterations per chain, retain their second halves (40,000 draws total), and have independently calculated retained split-Rhat <1.05 for alpha, slope and scatter variance. They use Gaussian age-summary errors and the default three-component latent-age mixture. These diagnostics establish adequate agreement of those chains, not correctness of the age/selection model. No posterior event with zero sampled occurrences is assigned literally zero probability.

| Reproduction target | Published slope (mag/Gyr) | Public-table reconstruction | Status |
|---|---:|---:|---|
| C25 R19 Gaussian-summary relation | −.062 ± .014 | −.06175 ± .01466 | Approximate numerical reproduction; author configs unavailable |
| C25 G11, after its young-host redshift correction | −.047 ± .0155 | −.04577 ± .01454 | Approximate numerical reproduction |
| G11 same revised ages, original HR | Not a headline target | −.03635 ± .01505 | Controlled correction comparison |
| W26/C26 full corrected Pantheon+ | −.007, uncertainty about .012 | −.00550 ± .01203 (G11-first); −.00649 ± .01330 (R19-first) | Approximate numerical reproduction |
| Full sample, approximate Tripp residual with tabulated errors | C26 −.022 ± .012 | −.01954 ± .01374; −.01987 ± .01300 | Close approximate reproduction, not identical coefficients/errors |
| Same, COVADD removed from tabulated variance | C26 −.024 ± .012 | −.02199 ± .01267; −.02157 ± .01277 | Approximate; free residual scatter remains in likelihood |
| Same, .06<zHD<.20 | C26 −.034 ± .012 | −.03189 ± .01183; −.03096 ± .01238 | Approximate numerical reproduction |

Errors here are posterior standard deviations, not Bayes factors or probabilities of cosmic acceleration. The two entries in a cell use G11-first and R19-first, respectively. Neither age choice was selected to maximize agreement.

Published cumulative-cut slopes −.034,−.027,−.026,−.025,−.024 are reproduced approximately by −.03189,−.02438,−.02372,−.02361,−.02199 (G11-first); each uncertainty is approximately .012. All sample counts match exactly. See `runs/age_signal/published-cut-comparison.csv` and `matched-age-summary.png`. Exact point-for-point reproduction remains blocked by the original age crosswalk, standardization coefficients, error reconstruction and LINMIX configuration.

A paired 2,000-resample unique-SN WLS bootstrap gives narrow-minus-full Tripp slopes of −.00516±.00357 (G11-first), with a 95% percentile interval [−.01258,+.00126]. This is a **diagnostic that ignores age errors**, not a definitive latent-age test or an author-LINMIX difference posterior. It demonstrates why overlapping cuts cannot be compared as independent detections. The observed cut progression alone does not identify its cause: sampling, changes of age leverage, corrections, baseline and population mix remain possible.

## Covariance and uncertainty audit

An important release-level discrepancy emerged: on these 196 rows, `MU_SH0ES_ERR_DIAG` does not equal sqrt(diag(`Pantheon+SH0ES_STAT+SYS.cov`)), contrary to the column-description implication. Their RMS difference is .06037 mag; the STATONLY comparison is .06903 mag. The original row-index subset was used for both matrices; both are symmetric and positive definite on this sample. For example SDSS 1371 has tabulated error .127058 mag, whereas sqrt(STATONLY diagonal) is .0855 mag. Full corrected fixed-age WLS gives chi²≈94.1 for 194 degrees of freedom with table errors; full-covariance GLS gives ≈197.6. This identifies a substantive input distinction; it does not establish which exact errors W26/C26 used or why the release differs. The release itself warns against fitting cosmology with the tabulated plotting errors.

With the **actual STAT+SYS diagonal**, Gaussian-summary LINMIX gives corrected slopes −.00938±.01132 (G11-first) and −.00723±.01088 (R19-first). Reversing bias gives −.02702±.01172 and −.02652±.01224; restricting these reversed-bias samples to z<.20 gives −.03547±.01147 and −.03610±.01178. Corrected narrow samples are −.01888±.01176 and −.01852±.01276. Thus uncertainty choices materially affect the uncorrected significance, but the corrected residual association remains weak. These LINMIX calculations use a covariance **diagonal**; they are not full-covariance fits.

A distinct errors-in-variables calculation includes the **full 196×196 STAT+SYS covariance**, analytically integrates a one-normal latent host-age population, and leaves residual intrinsic scatter free. For G11-first:

| Gaussian population model | Corrected slope MLE | Approximate profile interval Δχ²≤3.84 | Reversed-bias slope MLE | Profile interval |
|---|---:|---:|---:|---:|
| Age alone | −.01199 | [−.032,+.010] | −.02933 | [−.050,−.010] |
| Joint age, redshift, host mass and SN colour | −.01498 | [−.038,+.008] | −.02097 | [−.044,+.002] |

These are likelihood intervals on a .002 mag/Gyr grid, **not Bayesian credible intervals**; regular likelihood-ratio coverage can be imperfect near the nonnegative intrinsic-scatter boundary. The joint model treats redshift/mass/colour as fixed measured predictors; unknown age–mass–dust error covariance is not modelled. It is an exploratory sensitivity model, not a causal decomposition. Diagonal versus full covariance in the same age-only Gaussian model shifts the corrected MLE by about .002 mag/Gyr. An independent direct 2N block-Gaussian calculation agrees with the integrated conditional likelihood to 1.8×10⁻¹⁵; the zero-age-error limit agrees with GLS parameters to 4.2×10⁻¹⁷ (`algebra-validation.json`).

The summary-age uncertainty interpretation is itself problematic: the naive classical-error moment estimate Var(age)−mean(sigma_age²) is negative (−.35 or −.50 Gyr²). This is a failed approximation, not evidence of negative physical variance. Reported ages are posterior medians with broad, sometimes multimodal intervals, not necessarily unbiased classical noisy ages. Full posterior samples, stellar-population priors and correlated host-property posteriors are needed for a stronger latent inference.

Changing the reference flat-LCDM Omega_m from .2 to .4 changes the full-covariance fixed-age corrected slope only from −.00453 to −.00533 mag/Gyr, much smaller than its uncertainty. On observed age summaries, the age–z correlation is only −.20 (G11-first) or −.17 (R19-first). A local projection onto one fitted Omega_m derivative removes 4.1% or 2.6% of weighted age variance. These are conditional measured-age calculations; they do not limit arbitrary luminosity evolution or resolve the uncertain true-age distribution.

## Published mock: analytic and controlled numerical audit

Let A be host age, Z a redshift bin, and Y=bA+epsilon with zero-mean independent noise. A fixed correct cosmological baseline gives HR=bA+epsilon plus an intercept. It does **not** force E[HR|Z] to vanish. If instead Y is explicitly bin-centered while A is not, the population OLS slope is

`b_centered = b * E[Var(A|Z)] / Var(A)`.

This follows from the law of total covariance. For finite samples, replacing each conditional mean by the sample bin mean gives the same identity in the noiseless case using within-bin and total sums of squares. Wider redshift spans with more between-bin age variation therefore flatten the slope by construction. Centering both age and residual recovers the within-bin slope under the stated noise assumptions.

The controlled surrogate uses 200 synthetic galaxies in each of eight .05-wide redshift bins, 1000 realizations, seed 20260920, age mean 7−10z Gyr, age scatter1.5 Gyr, injected b=−.034 mag/Gyr, .12 mag Gaussian noise, and an explicitly correct luminosity-distance baseline. At zcut=.20/.42:

- Fixed true cosmological baseline: recovered −.03410/−.03406, with realization SD .00265/.00162.
- Residual-only bin centering: recovered −.02978/−.02141, with SD .00254/.00140.
- Center both variables: recovered −.03407/−.03404.
- Null injection: two-sided 5% rejection rates .051/.049 (Monte Carlo uncertainty about .007).

The noiseless attenuation identity is verified across all four age-evolution gradients (0,5,10,15 Gyr per redshift) and every cut to floating-point precision. See `mock-results.csv`, `mock-identities.json` and `mock-audit.png`. The mock is **not an exact reproduction of C26's numerical Figure2**: their precise age PDFs, noise, seed and conversion code are not public here. It reproduces and tests the stated operation. Its agreement with dilution cannot independently establish that a fitted ΛCDM baseline produced that operation in actual data. Conversely, cosmological fitting can project out a component of real age evolution; an appropriately specified joint inference/forward injection is needed to quantify that, rather than declaring it impossible.

## Park and C25 dependency verification

The [published Park journal HTML](https://academic.oup.com/mnras/article/549/2/stag935/8679169) was recovered separately and saved with a hash. It confirms the relevant arXiv methods: updated C25 ages, original G11/R19 HR inputs, |x1|≤3, sigma_x1≤1, |c|≤.3, sigma_c≤.1, leaving175 of199 G11 objects; its70 young hosts (age≤4 Gyr) supply HR=(-.641±.284)z+(.154±.063), a different correction from C25's199-object analysis. This exact Park175-object result is **blocked** without the exact original light-curve-parameter/quality-cut crosswalk and posterior bundle. Applying these cuts to newer Pantheon parameters would silently substitute another analysis.

C25's tabulated corrected-minus-original G11 HR is based on a redshift regression of82 age<4 Gyr hosts, HR=(-.467±.256)z+(.117±.052). A correction fitted from the same sample has uncertain coefficients shared by all corrected objects; exact propagation requires induced covariance and the covariance with the age fit. The released table gives original-style per-object HR errors, and does not release that full propagated correction covariance. Our numerical C25 reproduction follows those published columns; it does not make their omitted shared uncertainty disappear.

Park simulations adopt an empirical age slope, a model for host SFHs/DTD, and a host-mass distribution before projecting magnitude steps. They show consistency of those chosen inputs with steps; they do not uniquely select age as the cause over every correlated dust/metallicity alternative. Park explicitly acknowledges that its low-mass SFHs underrepresent old populations and can overestimate low/high-mass age contrast, and that explicit binary population synthesis is absent. The claimed independence of G11/R19 in a combined-significance argument is also contradicted by the33 shared SDSS IDs found above.

## Contemporary host-mass revision sensitivity

After this independent initial report was recorded, the standardization investigator supplied a verified reconstruction based on [Hoyt et al. Union3.1 Appendix F](https://arxiv.org/abs/2601.19424) and [Roy Choudhury's September repository](https://github.com/shouvikrc/pantheonplus-hostmass-correction-reconstruction), commit c5583379eb06f9a4b045353bb39a0991f8c15e4c (arXiv2607.24443v2,15 September2026). Source verification is detailed in `docs/investigations/standardization-dust.md`. This evidence was absent from the acquisition dossier and from my initial narrow exchange search.

It changes **12 of the196 age-matched SNe**, with mean delta .004893 mag across the full matched sample. Add the exact delta `new.m_b_corr-old.m_b_corr` to HR; `MU_SH0ES` in this reconstructed file is unchanged. Using actual covariance-diagonal LINMIX, the revised corrected slopes are −.00888±.01100 (G11-first) and −.00795±.01165 (R19-first), compared with original −.00938±.01132 and −.00723±.01088. Full-covariance fixed-age GLS changes to about−.0056 mag/Gyr. This bounded sensitivity leaves the weak corrected age association intact. The covariance and broader correction pipeline were held fixed, so it is not an official complete release reanalysis or proof against a cosmological effect of that revision.

## What this branch establishes

The numerical claims on both sides are substantially reproducible from public summary data once corrections and sample definitions are distinguished. A host-age association in older or bias-reversed residuals is supported. A comparably strong **additional** linear association in fully corrected matched-sample residuals is not established, but slopes around−.03 remain near or within some95% uncertainty bounds. These196 objects cannot certify that all cosmologically relevant evolution is removed. The C26 mock does not by itself supply the missing identification because its crucial bin-centering operation imposes the dilution under test.

Exact barriers are now narrow and actionable: original matched age choices, exact error variants and standardization coefficients; full host-age/property posteriors; the G11 young-host correction covariance; C26 original mock age-PDF arrays/seeds/conversion code; and a selection-aware injection/recovery or independent held-out sample. No author contact, publication or raw-source modification was performed.

## Reproduction commands and evidence

Use the root pinned environment (`uv.lock`, Python3.12.14). LINMIX source is the immutable acquired archive at `sources/external/age_signal/linmix-933dbb1359dcb5404cf881bfdd5cf433b0195152.tar.gz`; its URL/hash and Park HTML are in `runs/age_signal/acquisition.json`. It is imported locally rather than installed globally.

```bash
.venv/bin/python scripts/age_signal/analyse.py prepare
.venv/bin/python scripts/age_signal/analyse.py diagnostics
.venv/bin/python scripts/age_signal/analyse.py linmix
.venv/bin/python scripts/age_signal/mock_audit.py
.venv/bin/python scripts/age_signal/sensitivity.py
.venv/bin/python scripts/age_signal/covariance_linmix.py
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/age_signal/fullcov_eiv.py
.venv/bin/python scripts/age_signal/validate_algebra.py
.venv/bin/python scripts/age_signal/mass_revision.py
.venv/bin/python scripts/age_signal/summarize.py
.venv/bin/python scripts/age_signal/freeze.py
```

LINMIX commands reuse already saved individual-fit files; move the relevant output subdirectory aside for a fresh rerun. Sources, original row indices and masks are preserved in derived CSVs. Individual experiment manifests name input/output hashes, source revision, seeds and settings. `runs/age_signal/provenance/` contains code and dated plan snapshots whose bytes match every recorded digest, including historical setup variants. `final-inventory.json` inventories lightweight evidence; numerical chains and early unsuccessful/pilot results remain local. The first independent findings and each later amendment remain above for the decision trail.
