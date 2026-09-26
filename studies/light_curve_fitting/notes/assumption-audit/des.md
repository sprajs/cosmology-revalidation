# DES/Dovekie uncertainty and correlation audit

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../../README.md) for current execution status.

Audit date: 2026-09-21. Scope: released distances/covariances, calibration priors, host-population and age-model coverage. No SALT/SNM correction work or `phase2` products were modified. This extends, rather than repeats, the existing `docs/investigations/standardization-dust.md` investigation. Source release: `des-science/DES-SN5YR` commit `c9a4fcafc4cbd19bd750dee47fc76194a45c181f`; input hashes are in the two manifests under `runs/assumption_audit/des/`.

## Main findings

| Finding | Classification | Consequence |
|---|---|---|
| `CAL_SALT3` already includes `CALSPEC` | Confirmed release bookkeeping fact | Summing all individual systematic matrices double counts CALSPEC. The existing comparison uses the published total and avoids this trap. |
| The released VPEC increment is not positive semidefinite | Confirmed algebraic inconsistency with an independent additive covariance interpretation; production provenance unresolved | It cannot be treated as the covariance of an independent random nuisance. The total remains positive definite. An illustrative repair has a small SN-only effect. |
| PS1MD and Foundation wavelength uncertainties are drawn independently for identical named filter transmissions | Confirmed independence assumption in generator and released inputs; physical adequacy unresolved | Shared instrument/filter uncertainty is not represented as a shared wavelength draw. Independent survey effects may exist, but do not justify assuming every component is independent. |
| Nine calibration modes and three dust-posterior directions are treated as fixed uncertainty products | Confirmed finite representation; magnitude of missing convergence uncertainty unresolved | Normalizing weights fixes their mean normalization, not convergence of a covariance estimated from a small random sample. |
| “Model SN age” is one alternative population-model direction | Confirmed coverage limitation | It does not establish marginalization over general luminosity-age-redshift relations. |
| Shipped likelihood defaults and table reader do not read its own current data release | Confirmed implementation/reproduction issue | Use the independently verified local reader or repair the upstream adapter; this is not evidence that the published fit used that broken adapter. |

## 1. Covariance nesting, finite precision, and an invalid additive interpretation

Let `S` denote the inverse of `STATONLY.npz`, `T` the inverse of `STAT+SYS.npz`, and `K_i = inverse(single_i.npz) - S`. All inversions use Cholesky solves, and all products use the HD row ordering, with metadata aligned by `(CID, IDSURVEY)`. The statistical matrix is diagonal. Its standard deviations agree with table `MUERR` at release rounding precision.

The 24 individual increments obey

`T - S = sum(K_i) - K_CALSPEC`

to **6.46e-7 relative Frobenius error after statistical whitening**. The unwhitened maximum remainder, 0.0170 mag², is dominated by very large BEAMS-downweighted variances; it is inappropriate to apply an absolute tolerance to all entries without considering their scale. The packed precision arrays are float32. Typical spurious whitened eigenvalues after inverse/subtraction are approximately ±1e-7.

`K_CAL_SALT3` has ten robust positive eigenmodes, but `K_CAL_SALT3 - K_CALSPEC` has **nine**. This nesting is also consistent with the legacy `+cal` substring selector at `7_PIPPIN_FILES/D5yr_analysis.yml:360–361` and the case-insensitive substring logic at `sources/repos/RickKessler__SNANA/util/create_covariance.py:1473–1479`. The standalone filenames are therefore not a disjoint error budget. Do not add the standalone CALSPEC matrix to the calibration aggregate, and do not interpret ten aggregate modes as ten independent calibration draws. The prose/table discrepancy about nine versus ten surfaces remains a documentation/provenance issue: ten surface directories exist, whereas the released pure calibration covariance has nine significant modes.

The **VPEC increment** is different from the other named uncertainty components. Its smallest eigenvalue is **−0.0007541 mag²**, or **−0.04411 after statistical whitening**, with 47 negative whitened eigenvalues below −2e-6. These are much larger than float32 artifacts. One principal submatrix, for `(2008ar,65)` and `(1343871,10)`, is approximately

```
[[ 0.00223, -0.00021],
 [-0.00021,  0.00000]] mag².
```

Its minimum eigenvalue using the unrounded reconstructed values is −1.9602e-5 mag². Nonzero covariance with zero variance violates the Cauchy–Schwarz condition for an independent random error. Every entry of the reconstructed VPEC increment lies within 3.8e-9 mag² of a five-decimal grid. This strongly localizes a precision/provenance concern to the supplied VPEC contribution, but does not identify its original generation step.

The total remains positive definite: its minimum statistically whitened eigenvalue is 0.9620. The matrices alone do not prove whether VPEC was intended as an independent covariance or as a modification of an existing statistical model; the latter need not be positive semidefinite separately. The paper describes two velocity alternatives, and the release does not provide a generation ledger resolving the discrepancy. Request that ledger and the full-precision velocity contribution before a physical correction.

For scale only, replacing the negative eigenvalues of the VPEC increment by zero changes the SN-only flat-LCDM best fit from **0.3303168 to 0.3304927** and the local Fisher error from **0.0151996 to 0.0152072**. This arbitrary PSD completion is not an astrophysical velocity model, but it shows that this particular diagnostic does not explain the much larger original-to-Dovekie shift.

## 2. A concrete shared-passband independence assumption

The [Dovekie calibration paper](https://arxiv.org/abs/2506.05471), Figure 18, uses PS1SN filter characterization for Foundation as well. The corresponding release templates use the same PS1 `FILTPATH` and the same `g_filt_revised.txt`, `r_filt_tonry.txt`, `i_filt_tonry.txt`, and `z_filt_tonry.txt` files:

- `sources/repos/bap37__Dovekie/templates/new_kcor_templates/PS1SN.input:6–13`
- `sources/repos/bap37__Dovekie/templates/new_kcor_templates/Foundation.input:6–13`

Nevertheless, `sources/repos/bap37__Dovekie/surfaces-dovekie.py:213–237` samples a new Gaussian wavelength shift for each survey/filter. The zero-point covariance is used separately at lines 172–175; it does not couple these wavelength draws. This is present in the actual released prior realizations, not just an unused code path: `SALT3.DOV0UNC/SALT3.INFO:65` gives PS1MD g = **−8.347 Å**, whereas line 105 gives Foundation g = **+1.449 Å**. All ten released realizations differ in every common band. The RMS PS1MD–Foundation differences are **6.42, 16.42, 17.63, 17.34 Å** in griz.

These facts establish the independence assumption. They do not establish how much of an effective passband error ought to be common: atmospheric conditions, focal-plane sampling and reductions can differ. A physically motivated covariance should separate common hardware/filter characterization from survey-specific effective transmission. The useful next test is a rerun with that shared component and the same survey-specific remainder, followed through calibration, distance biases and covariance. Merely changing final error bars cannot reproduce it. No such end-to-end correction is claimed here.

## 3. Finite uncertainty directions and conditional population assumptions

After accounting for CALSPEC nesting, calibration contributes nine modes. `P24a`, `P24b`, `P24c`, `W22`, `MWEBV`, `COLORLAW`, `MASSLOC`, `ALPHAEVOL`, `BETAEVOL` and `GAMMAEVOL` are each rank one. The three dust posterior variations thus provide three possible directions; the W22 age/population variant provides one. Low rank is normal for a coherent nuisance, but limits what “including a systematic” establishes.

For independent zero-mean Gaussian draws and linear distance response, an estimate `q^T C_hat q = mean[(q^T delta_mu)^2]` from N draws has fractional standard deviation `sqrt(2/N)`: **47.1% for nine** and **81.6% for three**, in that systematic variance component. These are analytic convergence diagnostics under stated assumptions, not measured uncertainty in the published cosmological error bar. They neither prove undercoverage nor justify inflating the total by those percentages. Repeated seeds, more draws or derivative propagation are needed to demonstrate covariance convergence. The released matrices have no information about omitted draws.

The [DES-Dovekie methods](https://arxiv.org/abs/2511.07517), §6.2 and §11.5, retain the older Dust2Dust population after agreement checks, while reporting slightly worse metrics after recalibration. This is conditional model reuse, rather than a joint recalibration–dust-population posterior. Varying redshift dependence of alpha, beta and gamma separately also does not span all joint dependencies among dust, age, metallicity and selection. The existing report already explains why a colour change cannot be assigned uniquely to dust.

The paper's §3.2 footnote explicitly distinguishes alpha/beta from marginalized Bayesian nuisance parameters. However, BBC uncertainty renormalization and systematic variants already propagate some effects of fitting them. This audit therefore **does not claim their uncertainty is wholly missing**, and does not add an independent alpha/beta covariance on top of the release. A joint likelihood comparison would be needed to establish an omission.

## 4. What the age and host-mass tests actually cover

The [original DES methods](https://arxiv.org/abs/2401.02945), §6.2.3, describe the adopted W22 alternative as galaxy star-formation histories plus delay-time distributions generating progenitor ages and hence stretch–host correlations. They report too few red hosts and insufficiently steep correlations in that alternative. Local evidence is `runs/assumption_audit/des/des-v2.txt:1194–1226`. A single W22 covariance direction cannot establish marginalization over the intrinsic luminosity-age steps explored elsewhere in [Wiseman et al. 2022](https://arxiv.org/abs/2207.05583), or arbitrary smooth age evolution. “Model SN age” must not be read as a complete age uncertainty model. Conversely, this does not justify reapplying a full empirical age slope to corrected distances.

The explicit mass step is effectively hard: the shipped BBC input has location 10 and width 0.001 dex at `7_PIPPIN_FILES/base_files/bbc/BBC_des5yr.input:60–61`, and the prior independent release-closure test recovers that narrow form. The released mass errors are much larger. Under a **diagnostic Gaussian approximation** to those errors, 55/1820 objects have 0.1 < P(logM > 10) < 0.9: 33 DES and 22 low-z. Averaging only the explicit 0.033-mag step over that approximation instead of using a hard assignment produces a mean absolute per-object difference of 0.000295 mag; it is 0.000213 for DES and 0.000971 for low-z.

These small averages do not demonstrate a major new bias. Nor do they test mass-dependent bias-correction branch transfer, correlated mass-model errors or age/dust degeneracy in host SED fitting. The metadata errors are not full posterior samples, and moving a global threshold by 0.3 dex is not equivalent to marginalizing each object's correlated host posterior. A full check must propagate the same host latent variables through both the explicit step and the simulation bias corrections.

## 5. Release adapter failures and verification

`4_DISTANCES_COVMAT/DES-Dovekie-SN_Likelihood.py:14–15` defaults to two absent old filenames. Even supplying current paths leaves line 47 using `Table.read(..., format='ascii.csv')`; the actual HD is a commented whitespace SNANA table. The installed Astropy reader produces one comment-header column, and line 52's `data['zHD']` fails. The local `compare_des.py` reader uses whitespace plus comment handling and successfully reproduced the fit, so earlier local cosmology results are unaffected. This is an adapter/reproduction defect, not proof of a published inference error.

Another documentation inconsistency is `4_DISTANCES_COVMAT/README.md:26`, which lists alpha uncertainty 0.0003, whereas Dovekie Tables 3 and 8 give 0.003. Do not use the README number as a statistical prior without resolving the typo.

Reproduction commands:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/assumption_audit/des_covariance.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/assumption_audit/des_release_checks.py
.venv/bin/python scripts/assumption_audit/des_validate.py
```

Outputs: `covariance_results.json`, `single_systematics.csv`, `release_checks.json`, `calibration_prior_draws.csv`, `host_mass_boundary_diagnostic.csv`, plus source manifests under `runs/assumption_audit/des/`. Validation checks the matrix-nesting identity at the precision available, robust ranks, explicit PSD counterexample, total positivity, reproduced Omega_m, reader failure, published shared-filter draws and source hashes. It does not substitute those narrow checks for an end-to-end astrophysical refit.
