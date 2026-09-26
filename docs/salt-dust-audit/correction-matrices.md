# Released dust covariance, correction responses, and collinearities

This audit establishes numerical properties of the released original DES-SN5YR and DES-Dovekie distance products. **It does not establish that the light-curve, dust, selection, or cosmology analysis is bias-free.** It supplies reproducible uncertainty modes and response operators, identifies directions that are difficult to distinguish, and records information that the public covariance products cannot recover.

The code is [matrix_audit.py](../../scripts/salt_dust_audit/matrix_audit.py), with a separate [verification program](../../scripts/salt_dust_audit/matrix_verify.py). Results, complete input hashes, and output hashes are under [runs/salt_dust_audit/matrix_audit](../../runs/salt_dust_audit/matrix_audit). The audit reads the pinned source repositories and never edits phase2 inputs, photometry, fits, or active hierarchy outputs.

## Release and matrix semantics

| Product | Exact local release | Rows | Stored quantity | Correct operation |
|---|---|---:|---|---|
| Original DES-SN5YR | tag 1.3, `e3493cb3b9505fc1f3f392d887364b50dae21439` | 1,829 | `STATONLY.txt.gz` is zero; `STAT+SYS.txt.gz` and individual files are systematic covariance, in mag² | Add `diag(MUERR_FINAL²)` once |
| DES-Dovekie | `c9a4fcafc4cbd19bd750dee47fc76194a45c181f` | 1,820 | Packed upper triangle of **inverse** statistical-plus-selected-systematic covariance, float32, in mag⁻² | Unpack symmetrically and invert; do not add another statistical diagonal |

These conventions follow the [original release README](https://github.com/des-science/DES-SN5YR/blob/e3493cb3b9505fc1f3f392d887364b50dae21439/4_DISTANCES_COVMAT/README.md), [Dovekie README](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/4_DISTANCES_COVMAT/README.md), and [Dovekie loader](https://github.com/des-science/DES-SN5YR/blob/c9a4fcafc4cbd19bd750dee47fc76194a45c181f/4_DISTANCES_COVMAT/DES-Dovekie-SN_Likelihood.py#L67-L88). A Dovekie single contribution is recovered as

\[
K_s=P_{\mathrm{stat}+s}^{-1}-P_{\mathrm{stat}}^{-1}.
\]

Subtracting precision matrices does **not** yield a systematic covariance. Neither this difference of covariances nor a square root of its diagonal is a signed shift of distances.

Every exported row contains zero-based row index, `CID`, and `IDSURVEY`. Matrix order is the Hubble-diagram order; metadata is joined by the pair `(CID, IDSURVEY)`. Dovekie metadata order differs from matrix order. Both published total matrices pass positive-definiteness checks. Their smallest eigenvalues are 0.002230 and 0.003679 mag², respectively. The statistical diagonal agrees with published errors to rounding precision. The existing [release comparison](../../scripts/standardization/compare_des.py) uses the same correct format and subsetting semantics; a principal block of a precision matrix would condition on omitted SNe rather than marginalize over them.

The [Dovekie paper, Eq. 9 and §7](https://arxiv.org/abs/2511.07517v3) defines systematic contributions using weighted outer products of distance differences. The inspected [SNANA covariance implementation](https://github.com/RickKessler/SNANA/blob/886408a4e171896db5eaa97e735a655f50cec2db/util/create_covariance.py#L1200-L1220) constructs a scaled difference of `MU−MUREF`, then its outer product; the reference-distance subtraction matters for redshift-changing systematics. It adds the statistical diagonal before inversion. The implementation is evidence for these operations, not proof that this exact SNANA revision produced every released file.

## What can and cannot be recovered

For a rank-one component, `K_s = b_s b_sᵀ` determines `b_s` only up to a global sign. For several components, rotations of the factors also leave their total covariance unchanged. A covariance cannot reveal which realization corresponds to increasing physical `R_V`, which exact parameter displacement generated it, or the probability that the alternative dust model is true.

The outputs therefore store **spectral loadings with arbitrary orientation**, never a supposedly measured dust correction. The largest absolute loading in each column is set positive for reproducibility. The complete `eigenvalue_sign` column is essential: reconstruction is `B diag(sign) Bᵀ`. A negative signed mode is a matrix diagnostic, not a real-valued Gaussian uncertainty draw. For components whose retained signs are all positive, `B Bᵀ` is a controlled approximation to the released component, but choosing independent unit Gaussian amplitudes remains an additional interpretation of the released covariance model. The factors are not independent empirical measurements of dust physics.

All dust/MW components listed below have one resolved positive mode after statistical whitening. Their labels are preserved verbatim. The original README identifies `P21SYS1/2/3` as dust-population alternatives, `P21_HOSTCOLOR` as an alternative division by host colour, and `W22_AGE` as the Wiseman model. The current README retains the old names while the files use `P24a/b/c`, `BS21`, `INTRSC_COLOR`, and `W22`. The paper discusses the corresponding classes, but the full current file-to-configuration map and signed alternative distance tables are needed before assigning numerical derivatives with respect to physical dust parameters.

## Numerical precision and nonadditivity

Statistical whitening uses `S = diag(sigma_stat)` and `K_white = S⁻¹ K S⁻¹`. Modes are retained if `|lambda| > max(10⁻¹², 10⁻⁵ lambda_max)`. This is a stated compression threshold, not a significance threshold. Full raw and whitened eigenspectra are also retained.

This scaling is consequential. A few nearly rejected BEAMS events have effective statistical errors as large as 450–468 mag; these numbers implement near-zero cosmological weight and are not ordinary light-curve measurement errors. The statistical covariance condition numbers are approximately 9.3×10⁷ and 1.3×10⁸. In Dovekie, converting packed float32 precision to covariance and then subtracting can generate conspicuous raw negative eigenvalues at extremely low statistical weight. For MWEBV, the negative part is only 2.35×10⁻⁶ of the whitened Frobenius norm. The independently reconstructed rank-one MW components agree with the released matrices to 3.42×10⁻⁶ and 7.13×10⁻⁶ in that metric. Raw unweighted ranks would misleadingly promote numerical artifacts into extra physical modes.

VPEC is a separate issue: it retains negative modes after whitening, with negative-spectrum Frobenius fractions of 0.01715 and 0.01769. Original VPEC lies on a coarse 10⁻⁵ mag² lattice. A direct witness is `(CID 1343871, CID 2007nq)`: `Cii=0`, `Cij=0.00021`, `Cjj=0.00515` mag², violating the covariance Cauchy–Schwarz bound. This is compatible with serialization damage and does not establish an astrophysical error. It prevents using this particular single-component file as an exact standalone Gaussian prior without resolving its provenance/precision. The supplied total covariance remains positive definite. Negative modes are preserved, not silently repaired or clipped out of the source matrix.

The naive sum of every released individual file differs from the released total systematic covariance because the file groups are not all disjoint:

| Release | Individual files | Relative Frobenius mismatch, raw | Relative Frobenius mismatch, statistically whitened |
|---|---:|---:|---:|
| Original | 23 | 0.003103 | 0.003783 |
| Dovekie | 24 | 0.005736 | 0.003527 |

**The calibration overlap is resolved.** Subtracting the standalone CALSPEC matrix from `CALIBplusSALT3` or `CAL_SALT3` leaves nine positive whitened modes and no resolved negative ones. The saved PIPPIN grouping rule `[CAL_SALT2] [+cal,=DEFAULT]` includes CALSPEC: SNANA's `apply_filter` uppercases the labels and performs substring matching for `+`. Thus nine SALT surface directions plus CALSPEC explain the ten modes in the combined calibration group. The Dovekie paper's §7.1/ Table 6 nine-versus-ten wording is not evidence of an additional SALT realization or an incorrect current weight.

After removing the separately listed, already grouped CALSPEC contribution, Dovekie's component sum matches the released total systematic covariance to **6.4624×10⁻⁷ relative whitened Frobenius norm**, consistent with finite precision. It is not an unexplained systematic-budget gap. Original DES retains one resolved positive, unassigned remainder: whitened eigenvalue 0.52949909, trace 0.01306775 mag², and relative norm 0.0034653. The original README lists `SIGINT_MODEL`, but no individual file by that name is present. A missing scatter-model component is therefore a provenance lead, not an established identity; the remainder is explicitly exported as `unassigned_grouping_remainder_modes.csv`. It is already included in the supplied total covariance and must not be added again.

The exact group calculation is independently reproducible with [matrix_grouping.py](../../scripts/salt_dust_audit/matrix_grouping.py), with [results](../../runs/salt_dust_audit/matrix_grouping/results.json); the main summary includes the corrected grouping calculation too. The total release is used for the total-weight response operator throughout. None of these grouping diagnostics changes distances, changes the supplied total covariance, or assigns a physical dust bias.

## Correction Jacobian and local cosmology response

The released distance convention is

\[
\mu=-2.5\log_{10}x_0+\alpha x_1-\beta c+\gamma h-b_{\rm BBC}-M_0.
\]

At fixed fitted photometry and fixed BBC correction, `dmu/dalpha=x1`, `dmu/dbeta=−c`, `dmu/dgamma=h`, and `dmu/dbias_scale=−biasCor_mu`, where the last column varies the multiplicative scale of the *already published* bias term. `h=logistic((logmass−10)/0.001)−1/2` is evaluated from rounded public masses and should not be regarded as exact near the mass boundary. These are algebraic response columns, not an end-to-end rerun: changing dust, beta, or calibration normally changes fitted light curves, simulated populations, bias cells, selection, scatter, and fitted nuisance parameters together.

The row table also records the actual published `−beta*c`, `alpha*x1`, and `−biasCor_mu` terms, using release-specific documented coefficients. Coefficient values and masses are rounded, so these diagnostic terms do not supersede the released `MU` vector. Per-event light-curve derivatives are `dmu/dx0=−2.5/(ln(10)*x0)`, `dmu/dx1=alpha`, `dmu/dc=−beta`, and `dmu/db_BBC=−1`; the full fit covariance and the dependence of BBC on these same quantities are required for their uncertainty propagation.

For a specified local model derivative matrix `H` and fixed released precision `P`, the exported operator is

\[
A=(H^T P H)^{-1}H^T P,\qquad
\delta\theta=A\,\delta\mu,\qquad
C_{\theta,s}=A K_s A^T.
\]

**The sign is positive:** `delta_mu` means a positive addition to the observed released distance; `delta_theta` is the corresponding change in the fitted model. A physical correction intended to remove a known positive observational bias would enter as its negative. Operator units are fitted-parameter units per magnitude, including a magnitude intercept. No physical bias has been inferred here.

Operators are supplied for flat ΛCDM (`intercept, Omega_m`), flat wCDM (also `w0`), and flat w0waCDM (also `wa`), linearized at `Omega_m=0.3,w0=−1,wa=0`, using the release's `zHD` and observer factor `1+zHEL`. Each model has a free intercept and no CMB/BAO/external prior. Both statistical and total weighting are provided. The four-parameter SN-only response is particularly ill-conditioned and must not be substituted for the paper's combined-probe constraints. Responses are local, fixed-covariance quantities; nonlinear refitting, a changed sample, and changed covariance are outside this calculation.

For illustration, under **total-weight flat ΛCDM** a `+0.01 mag*z` data perturbation maps to `delta_Omega_m=−0.00723` original and `−0.00728` Dovekie. A `+0.01 mag` shift of every DES event relative to low-z maps to `−0.00937` and `−0.00804`. These are controlled sensitivity directions, not measured dust errors or uncertainty bounds. The nuisance response changes with metric and release; simply multiplying a dust-law change by a universal beta correction is not justified.

The following selected component summaries use the **released component scale**, not a newly inferred physical one-sigma dust prior. The last column is `sqrt((A K_s Aᵀ)_Omega_m,Omega_m)` with total weighting; it is a conditional propagation diagnostic, not the increase in fitted error from an add-one-in analysis and not the bias of the alternative fit.

| Release | Component | Median diagonal amplitude (mmag) | Local propagated Omega_m scale |
|---|---|---:|---:|
| Original | BS20 | 6.950 | 0.001025 |
| Original | P21SYS1 / 2 / 3 | 5.330 / 5.520 / 9.260 | 0.000683 / 0.001959 / 0.000151 |
| Original | P21_HOSTCOLOR | 20.410 | 0.000089 |
| Original | W22_AGE | 13.240 | 0.000585 |
| Original | MWEBV / MWCOLORLAW | 0.870 / 0.363 | 0.000861 / 0.001985 |
| Dovekie | BS21 | 6.570 | 0.000662 |
| Dovekie | P24a / b / c | 4.810 / 4.914 / 9.125 | 0.000570 / 0.001993 / 0.000551 |
| Dovekie | INTRSC_COLOR | 20.265 | 0.000031 |
| Dovekie | W22 | 6.855 | 0.000194 |
| Dovekie | MWEBV / COLORLAW | 1.610 / 0.460 | 0.000623 / 0.001710 |

A large diagonal amplitude can project weakly onto one cosmological parameter after the other released uncertainties are included. That fact does not demonstrate that the corresponding physical model is correct or harmless under another likelihood.

## Collinearity results

All cosines below use the statistical inverse covariance and remove the arbitrary common magnitude intercept. Thus they compare coherent distance patterns at their actual statistical weight. They are geometric degeneracies, not posterior probabilities or evidence that the latent physical parameters are correlated by exactly these amounts.

| Compared directions | Original | Dovekie | Meaning |
|---|---:|---:|---|
| MW-law mode and `dmu/dOmega_m`, absolute cosine | 0.586 | 0.515 | A portion of the released foreground-law response resembles cosmological redshift curvature |
| P21SYS2 / P24b and `dmu/dOmega_m`, absolute cosine | 0.417 | 0.402 | These dust-population alternatives have non-negligible cosmological alignment |
| `−c` and `−biasCor_mu`, signed cosine | −0.664 | −0.688 | Colour standardization and simulated bias correction have strongly overlapping patterns |
| `dmu/dOmega_m` and linear `z`, signed cosine | −0.9963 | −0.9962 | Smooth redshift evolution can closely mimic the matter-density response over this sample |
| `dmu/dw0` and `dmu/dwa`, signed cosine | 0.9526 | 0.9528 | SN-only dark-energy responses are difficult to separate |

The normalized Gram matrix of all ten nonconstant cosmology/nuisance directions has condition number 1.13×10⁶ original and 1.05×10⁶ Dovekie. This extreme value is driven in part by deliberately included near-redundant smooth cosmological/redshift columns. Restricting to the seven nuisance/host/redshift/survey columns gives approximately 12.93 and 13.20; it would be misleading to describe every dust parameter as individually unidentifiable based on the ten-column condition number. Matrices and conditional variance-inflation factors are saved so this distinction can be inspected.

For comparisons between whole covariance components, `component_covariance_overlap.csv` records the Frobenius cosine of the statistically whitened, intercept-projected matrices. For rank-one positive components it equals the square of the absolute mode cosine. These overlaps do not authorize treating physical model alternatives as statistically independent, assigning signs to their shifts, or inferring an empirical dust prior.

## Artifact schema and reproduction

For each release, `rows_and_correction_jacobian.csv` fixes the ordering and correction-column signs; `cosmology_response_operator.csv` maps a row-ordered vector in magnitudes to parameter shifts. `propagated_parameter_covariances.json` contains full small matrices, including cross-parameter covariances. `systematic_diagonal_variances.csv` records the diagonal of every recovered component. `systematic_eigenspectra.csv.gz` retains raw and whitened eigenvalues, while `systematic_spectral_modes.csv.gz` stores compressed, signed factors with explicit row IDs. Original DES also has the unassigned remainder factor described above. The mode/design cosines, component overlaps, stat/total design Gram matrices, numerical checks and source hashes are separate human-readable files. Compression affects the factor representation only; covariance propagation uses the recovered full input matrices.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/salt_dust_audit/matrix_audit.py
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/salt_dust_audit/matrix_grouping.py
OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/salt_dust_audit/matrix_verify.py
```

[Verification](../../runs/salt_dust_audit/matrix_audit/verification.json) checks every manifest hash and row order, compares the separate grouping calculation, reconstructs MW components from exported factors, preserves negative VPEC modes, compares cosmological distances and derivatives to independent Astropy implementations, and tests the operator's sign and scale by actual nonlinear flat-ΛCDM refits under a small colour-response perturbation. Astropy distance differences are below 10⁻¹⁴ mag; the operator/refit beta-response differences are below 4×10⁻⁹ in Omega_m per unit beta. These checks verify the numerical calculation, not the completeness of the underlying astrophysical uncertainty budget.

## Inputs still required for a stronger bias claim

The present files cannot supply signed systematic shifts, complete dust-parameter response derivatives, or a joint latent dust/colour/host prior. Those require the exact alternative `MU`, `MUREF`, fitted nuisance parameters, per-variation row selections, simulation seeds/configurations and systematic weights. Exact file-to-configuration mappings and more precise covariance products are needed to identify the original unassigned remainder and resolve the VPEC serialization diagnostic; Dovekie's group sum is resolved above. Cross-release covariance is needed to give an uncertainty to a Dovekie-minus-original distance difference. The published matrices do not bound the effect of a dust law or population evolution outside their tested alternatives.

A demonstration of negligible cosmological bias additionally requires end-to-end closure and coverage tests across physically defensible dust, intrinsic-colour, host, calibration and selection alternatives, with refitting and retraining where the change requires them. These response matrices help choose and diagnose those tests; they cannot replace them. No event distances or published uncertainties have been automatically corrected by this audit.
