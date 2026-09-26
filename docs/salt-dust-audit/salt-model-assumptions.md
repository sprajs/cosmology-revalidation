# SALT2/SALT3 assumptions, dust non-identifiability, and falsifiable bias tests

Audit date: 2026-09-21. This document audits the model layer, independently of the ongoing phase-2 analysis. It does not alter that analysis or turn a proposed systematic variation into a measured correction. Repository paths below are relative to the repository root. Primary sources and source-file hashes are recorded in `theory-source-registry.json`; the executed model-grid diagnostic is `theory-model-grid-audit.json`.

**The present evidence does not establish a bias-free supernova analysis.** SALT is an empirical spectral model, not a unique separation of dust, intrinsic color, progenitor evolution, and cosmological dimming. This is a structural limitation, not evidence that an arbitrarily large alternative correction is real. The useful target is a stated tolerance on distance/cosmology bias, tested against explicit alternative populations and survey selection, with correlated uncertainties and unresolved directions retained.

## 1. Keep five different objects separate

1. **SALT spectral/color model:** maps fitted amplitude, shape, and color to observer-frame flux. It absorbs observed spectral diversity into a small number of coordinates.
2. **Milky Way extinction:** observer-frame foreground attenuation, imposed using a map, extinction law, and normalization. It also enters training.
3. **Host extinction in a simulator or physical model:** rest-frame attenuation, potentially with a distribution of reddening and extinction-law parameters. Standard SALT light-curve fitting does not separately measure these latent quantities.
4. **Luminosity standardization:** a relation between fitted SALT coordinates and standardized magnitude, commonly involving alpha, beta, and a host term.
5. **Selection/bias correction:** a conditional expectation learned from a simulated population and an observing/selection pipeline. It is not an extinction law or universal linear correction matrix.

SNANA is software that can implement these layers; running it correctly does not prove that their physical assumptions hold. Conversely, discovering a poor astrophysical assumption does not establish a software bug.

## 2. The actual local SALT model

The run log `phase2/official/diagnostics/snana-des-hessian0.log` identifies the model used by the refits as `sources/repos/des-science__DES-SN5YR@1.3/2_LCFIT_MODEL/SALT3.DES5YR`, not simply a generic model called SALT3. The archived SNDATA distribution contains a parallel copy; the executable log is the authority for the actual path. Its source revision is the original DES release tag 1.3, commit `e3493cb3b9505fc1f3f392d887364b50dae21439`.

The model has the usual two-surface form

\[
F(p,\lambda)=x_0\,[M_0(p,\lambda)+x_1M_1(p,\lambda)]
                   \exp[c\,CL(\lambda)].
\]

Here p is rest-frame phase, lambda is rest-frame wavelength, and `t0` is also fitted. This is the defining empirical form introduced in [Guy et al. (2007)](https://arxiv.org/abs/astro-ph/0701828). The color term is independent of phase; shape-dependent color evolution resides in the spectral component. The fitted c describes color relative to a trained fiducial spectrum. Neither c=0 nor c<0 establishes absence of dust.

In the pinned SNANA source, `genmag_SALT2.c::SALT2colorlaw1` writes the identical multiplier as `exp[(ln(10)/2.5)*c*C(lambda)]`. Thus its C has the opposite sign to conventions that write `10**(-0.4*c*CL)`. The actual code's reference wavelengths are 4302.57 and 5428.55 Angstrom, with C=0 and C=1 respectively. Comparing raw color-law arrays from different packages without converting conventions can reverse a conclusion. A positive c leaves the reference blue flux fixed and increases reference red flux at fixed x0; it does not act as a pure physical screen at fixed luminosity. The fitted amplitude absorbs the missing grey/blue dimming.

The actual `SALT3.INFO` uses five fitted polynomial coefficients over 2800–8000 Angstrom and tangent-line continuation of C outside that interval. Its native surface grid covers −20 to +50 days and 2000–11000 Angstrom. The SNANA interpolation table is finer than the original data grid; finer sampling is not additional training evidence. The refit namelist restricts phases to −15 to +45 days and **mean filter rest wavelength** to 3500–8000 Angstrom. Filter wings can extend outside this mean-wavelength interval. Neither the file extent nor the broader generic SALT3 advertised coverage certifies accuracy there.

The log reports two loaded surfaces. Although `salt3_template_host.dat.gz` is present in the directory, `NSURFACE_SALT2` loads numbered templates, and this run loads only templates 0 and 1. File presence is not evidence that host-dependent spectra were fitted. The `SIGMA_INT=0.106` entry is marked for simulation; it is not an all-inclusive physical error budget. The INFO also specifies a 0.005-mag error floor, a color-dispersion cap of 1, spline mean-surface interpolation, and linear error-map interpolation.

## 3. What training identifies, and what it fixes by convention

SALT training can leave each SN amplitude free; it does **not** need a cosmological distance–redshift relation to determine relative spectral/time behavior. Normalization and coordinate conventions remove otherwise redundant descriptions. SALT3 uses training-sample centering/scaling and a shape/color decorrelation convention; this does not prove physical independence of dust and progenitor properties. Spectral continuum recalibration, regularization of sparsely sampled regions, and fitted residual covariance are additional modeling choices. These features, and the distinction between population residual variation and finite-training uncertainty, are documented in [Kenworthy et al. (2021), sections 2.1–2.5](https://arxiv.org/abs/2104.07795).

For example, rescaling `M0,M1 → a*M0,a*M1` and `x0 → x0/a` leaves flux exactly unchanged. It shifts `mB=-2.5 log10(x0)+constant` by `2.5 log10(a)`, which the fitted magnitude intercept absorbs. Similarly, a color-origin shift can be traded against the mean spectral surface. A coordinate choice cannot measure absolute intrinsic color. A shape/color rotation can redistribute apparent correlation while preserving essentially the same predicted light curves.

This has a direct practical consequence: before comparing beta, population mean color, or host steps between retrained models, align their coordinate conventions and refit the luminosity relation. A reduced host step is not automatically removal of a physical host effect. [Taylor et al. (2024)](https://arxiv.org/html/2401.07304v1) explicitly found that a mass-step reduction in their split-model analysis arose from a changed fiducial normalization; their low-mass training also displayed UV ringing. Neither result says that every host-specific model fails.

Training and evaluating on overlapping SNe introduces shared information. A held-out object split must keep duplicate observations of a SN together; a calibration perturbation must still be shared across every object that uses that calibration, including training. A random held-out split drawn from the same selected population tests interpolation within that population. Survey, redshift, host, color-tail, and wavelength transfer require their own tests.

The training stage's free amplitudes also do not make the **whole analysis** cosmology independent. Population fitting from Hubble residuals, simulator luminosities/rates, redshift-bin assumptions, bias correction, and the final standardization/cosmology likelihood may condition on a fiducial distance relation. The dependence must be traced at those stages instead of inaccurately attributing it to an obligatory cosmology inside SALT spectral training.

## 4. Dust physics does not collapse to c or beta

For a foreground screen, write

\[
F_{\rm extincted}(p,\lambda)=F_{\rm intrinsic}(p,\lambda)
     10^{-0.4 E k(\lambda;R_V)},
\qquad E=E(B-V),\quad R_V=A_V/E.
\]

The idealized identity `R_B=R_V+1` follows only when `E=A_B-A_V` uses the same B,V definitions. A coefficient multiplying an effective SALT coordinate is not automatically this physical R_B. A broad-band attenuation is a ratio of wavelength integrals of `F*T*10**(-0.4*A)` and `F*T`, including the same photon weighting in both integrals; it depends on phase, SED, filter, and reddening. Consequently a phase-independent monochromatic dust screen can produce phase-dependent broad-band color changes. Observer-frame MW attenuation must be evaluated at observer wavelength; host attenuation at rest wavelength. Their operations cannot be replaced by an indiscriminate single color subtraction. The monochromatic-law versus photometric-system distinction is explicit in [Fitzpatrick (1999)](https://arxiv.org/abs/astro-ph/9809387).

Galaxy integrated attenuation additionally depends on source/dust geometry and scattering, whereas a SN probes a particular sightline. Equality between an integrated-host attenuation slope and SN sightline extinction is an extra physical hypothesis. A discrepancy between them cannot by itself falsify either estimator.

The often-used approximation `c=c_int+E` describes an effective latent population model. Its validity in fitted SALT coordinates must be measured by applying the physical extinction law to a spectral time series, simulating photometry, then refitting with the actual cadence/bands/cuts. Directly adding E to the simulation's intrinsic SALT c assumes the answer. The existing `docs/phase2/generative-model.md` correctly labels its decomposition as effective and preserves this limitation.

[Brout & Scolnic (2021)](https://arxiv.org/abs/2004.10206) motivate an intrinsic-color plus dust population with variable extinction slopes and host dependence. Their result is evidence that such a model can explain particular summary trends; it does not establish unique physical identification of its latent decomposition, the universality of its population priors, or the absence of additional luminosity evolution.

### Exact and approximate degeneracies

The following derivations are audit calculations, not claims of observed population values. Remove shape and write a deliberately simplified relation

\[
c=u+E,\qquad m=\mu+M+b u+R E+h a+\epsilon.
\]

Here u is intrinsic effective color, a an environmental variable, and R an effective dust luminosity coefficient. If u and E are independent, R is constant, and measurement error/selection are absent, the fitted linear slope is

\[
\beta_{\rm eff}={b\operatorname{Var}(u)+R\operatorname{Var}(E)
                    \over\operatorname{Var}(u)+\operatorname{Var}(E)}.
\]

With `s=Cov(u,E)`, the numerator gains `(b+R)s` and the denominator gains `2s`. Selection can induce this covariance even if the parent variables are independent. With random R the general numerator is `b*Cov(u,u+E)+Cov(RE,u+E)`; substituting the mean R without testing its correlations is unjustified. Measurement errors correlated between m and c add another term if they are ignored.

The population mean residual after subtracting a fixed beta color relation is

\[
\overline r(z,H)=(b-\beta)\overline u(z,H)
                   +\mathbb E[(R-\beta)E\mid z,H,S=1]
                   +h\overline a(z,H)+\cdots .
\]

Even a constant physical law can therefore yield evolving mean residuals as the **selected** population changes. Conversely, a nonzero intrinsic/dust effect need not bias cosmology if the correct selected conditional mean is modeled. Small global mean residual, scatter, or host step is insufficient to establish this transport across redshift.

There is an exact reparameterization for constant b,R:

\[
E'=E+k a,\quad u'=u-k a,\quad
h'=h-(R-b)k.
\]

It preserves both observed c and m. Nonnegative dust restricts allowed k and support; independent/exponential dust priors may break this symmetry *conditionally on those priors*. They do not create a new physical observation. Host ages, masses, and colors can be mutually correlated and measured with correlated errors; a marginal environmental slope does not identify h after this ambiguity.

A separate exact ambiguity is `m_std=mu(z)+M+DeltaM(z)`: replacing `mu→mu+g(z)` and `DeltaM→DeltaM-g(z)` preserves all SN magnitudes for any admitted g. Arbitrary luminosity evolution and a distance curve cannot both be measured from these magnitudes alone. This does not show that a physically plausible dust/progenitor population realizes an arbitrary g. The constant part is already degenerate with M/H0; independent observables and explicit restricted evolution models are necessary for the remaining shape.

## 5. Uncertainties and correction matrices that are meaningful

For a coherent fitted covariance `V_i` of `(mB,x1,c)`, a fixed-coefficient Tripp variance is

\[
\operatorname{Var}(\mu_i\mid\alpha,\beta)=
\begin{bmatrix}1&\alpha&-\beta\end{bmatrix}
V_i
\begin{bmatrix}1\\\alpha\\-\beta\end{bmatrix}.
\]

Off-diagonal terms are essential. A Hessian covariance is conditional on the fitted model/cuts and a local Gaussian approximation. It is not a posterior for unknown dust physics or a remedy for multimodality, epoch rejection, or a wrong mean model. The repository's output-only Hessian export avoids combining MINOS scalar errors with Hessian cross-covariances; that is an implementation improvement, not physical validation.

For shared physical/calibration parameters eta, define a response matrix `J_ik = d(mu_i)/d(eta_k)` by **rerunning the affected analysis stages**. If their joint prior/posterior covariance is Sigma_eta, then to first order

\[
C_{\rm shared}=J\Sigma_\eta J^T.
\]

A universal uncertain slope on f(z) yields the rank-one matrix `sigma_s^2 f f^T`. Treating it as diagonal noise wrongly makes it average down across SNe. If calibration, training, dust hyperparameters, and bias corrections are correlated, use one joint Sigma or joint realizations; adding independently derived outer products can double count or omit cross terms. Numerical finite differences are responses, not calibrated uncertainty amplitudes: a ±0.1 R_V perturbation is not a one-sigma mode without evidence for that scale.

To diagnose collinearity, let W be the inverse data covariance and N the nuisance design matrix (at least an intercept). Residualize with

\[
P=W-WN(N^TWN)^{-1}N^TW.
\]

For two response columns j,k, use the nuisance-projected overlap

\[
\rho_{jk}={j^TPk\over\sqrt{(j^TPj)(k^TPk)}}.
\]

A zero denominator means an unobservable mode after nuisance projection, not zero correlation. Whiten and scale declared physical columns before reporting singular values/condition numbers; otherwise units can manufacture apparent ill conditioning. With cosmology derivative matrix D, the local systematic displacement is `(D^T P D)^(-1) D^T P J delta_eta`, when the inverse exists. Rank deficiency must be retained or profiled; a pseudoinverse does not constitute physical identification. These expressions assume differentiability, fixed covariance, fixed membership, and local linearity. Selection boundaries, rejected epochs, or nonlinear dust corrections require paired simulations/refits, with membership changes recorded separately.

A BBC correction is instead `b(q)=E[mu_fit-mu_true | q,S=1, simulation model]`, for its particular conditioning coordinates q and weights. Different dust populations can share similar fitted color histograms and give different conditional mean biases. A universal matrix that deconvolves unknown dust from `(mB,x1,c)` does not exist without further assumptions.

### Which covariance is actually in the local fitter?

`genmag_SALT2.c::gencovar_SALT2` adds a color-dispersion covariance across epochs **in the same passband**, evaluated at the filter's mean rest wavelength; different-passband entries of this term are zero. Its diagonal also uses local M0/M1 variance maps and an error floor. Those maps combine model-error descriptions, but this per-SN matrix cannot encode correlations between different SNe caused by their common trained surfaces. A full training/calibration uncertainty treatment needs shared variations or the relevant joint parameter covariance. It is also not the arbitrary temporal/wavelength covariance of physical dust-law diversity.

SALT2 and SALT3 error-map units differ. The source has separate relative-error and absolute-surface-error paths; swapping map filenames or interpreting both as relative variances is invalid. The presence of `salt3_lc_model_variance_*` files does not prove their use: the actual log reads `salt3_lc_variance_0`, `_1`, and `_covariance_01`, plus the color-dispersion file.

## 6. Executed grid diagnostic and its limits

Run:

```bash
phase2/env-official/bin/python scripts/salt_dust_audit/theory_model_grid_audit.py
```

The script hashes its inputs and output-producing code; validates coordinate alignment and finiteness; computes eigenvalues of every local 2x2 surface covariance; minimizes `v00+2*x1*v01+x1^2*v11` analytically over declared intervals; and checks mean-surface signs. It reads 1,828 successful recovered refits, whose envelope is `x1=[-2.979277,4.016594]`, `c=[-0.272958,0.456116]`. These are all successful refits, not necessarily final BBC-quality objects.

| Original source grid | Nodes | Negative covariance eigenvalue | Minimum eigenvalue | Negative variance anywhere over recovered x1 envelope |
|---|---:|---:|---:|---:|
| Whole tabulated domain | 63,971 | 0 | 2.41168e-11 | 0 |
| −15:45 days, 3500:8000 Angstrom | 27,511 | 0 | 9.10696e-8 | 0 |

The code contains a branch replacing a negative SALT3 quadratic variance with its absolute value, and similarly treats a negative flux normalization. This audit finds **no negative covariance or quadratic variance on the original grids**, so the branch's existence is not evidence it corrupted the released fits. Convex linear interpolation of positive-semidefinite matrices also preserves positive semidefiniteness, subject to the actual interpolation staying within valid cells.

The mean linear surface does have support limitations: within the fit rectangle, 703 original nodes permit negative `M0+x1*M1` somewhere in `x1=[-3,3]`; 689 do so within the observed refit envelope. Positive color multipliers cannot change that sign. This is a mathematical property of the linear model, not a count of affected observations. The script does not reproduce SNANA's spline resampling, actual observation phases, or passband integrals. A monochromatic negative value can integrate to a positive band flux, and a problematic phase can be unobserved. An instrumented or independently reproduced passband/epoch check is required before assigning an actual fitted-distance effect. Covariance positivity establishes neither calibration nor correctness of the mean SED.

## 7. Published validation: substantial evidence, bounded claims

Primary-source results must not be converted into a universal tolerance for this dataset:

- [Dai et al. (2023; revised 2024), appendix B](https://arxiv.org/html/2212.06879v2) report a baseline simulation bias in w of `−0.024±0.006`, changing to a recovered `w=−1.005±0.009` when the low-z training population was changed to match their training compilation. Their study also examines training spectra, calibration, and model variations. These numbers are conditional simulation results, not a correction to DES. The baseline finding directly motivates training-population sensitivity tests rather than assuming matched software removes all bias.
- [Taylor et al. (2023)](https://arxiv.org/abs/2301.10644) compare SALT2 and SALT3 trained on the same 1,083 SNe and report a negligible framework-choice systematic in their SN+CMB comparison, `Delta w=0.001±0.005`. This is strong evidence for that controlled comparison, not for the shared single-color, population, calibration, or selection assumptions.
- [Kenworthy et al. (2025), SALT3+](https://arxiv.org/abs/2502.09713) find additional phase-dependent color variation that ordinary SALT3 partly absorbs into c. They report an omitted-component Hubble-residual trend of `0.039±0.005` mag, while finding **no evidence of a bias in current cosmological measurements**. Both findings matter: incomplete spectral dimensionality is observed, but an environmental/redshift distribution is required to convert its trend into a distance-shape bias. This audit has not refitted DES with SALT3+.
- [DES-Dovekie](https://arxiv.org/abs/2511.07517), already preserved as `papers/text/2511.07517v3.txt`, documents a host-law numerical approximation correction and calibration-systematic weights summing to 0.81. Its appendix's F99-only cosmological shift is not the same as the combined recalibration/retraining shift. Thus neither an obsolete implementation nor the net release difference should be labeled a direct dust measurement. The older/newer law distinction must be kept explicit in simulations.

## 8. Discriminating tests and completion boundary

These are concrete test requirements, not tests claimed to have passed here.

| Bias channel | Distinguishing test | Record required | Limitation after passing |
|---|---|---|---|
| Physical dust projected onto SALT c | Redden independent intrinsic SEDs with physical laws and variable R_V; simulate actual photometry and refit | Responses of mB,x1,c and fit residuals vs E,R_V,phase,bands,z; covariance and accepted observations | Valid only for the tested SED/dust family |
| Intrinsic color evolution absorbed by dust | Paired dust-only, intrinsic-only, and joint perturbations with matched observed-color distributions | Weak response directions, parameter recovery, predictive residuals on held-out optical/NIR data | Restrictive priors can remain decisive |
| Missing spectral dimension | Inject empirical spectral/phase variation, or refit an independently trained richer model | Conditional distance changes and flux prediction by phase, wavelength, environment, redshift | Model agreement can share calibration/population assumptions |
| MW/host/calibration confusion | Joint observer-frame MW and rest-frame host perturbations plus zero-point/filter changes | Correlated response columns and nuisance-projected singular spectrum | Same sightline/systematics may remain weakly identified |
| Training demographics and regularization | Retrain under declared calibration, color-support, host, and regularization alternatives | Trained assets, training IDs, out-of-sample predictions, matching bias simulations | Finite alternatives do not span arbitrary evolution |
| Dust-dependent selection | Generate undetected as well as detected objects; apply detection, redshift, quality and classification selection | Recovery/coverage vs z,host,color; full normalization and simulation truth | Injecting and fitting the same family tests implementation only |
| Mask/cut sensitivity | Freeze nominal masks, rerun controlled alternatives and the resulting bias calibration | Fixed-membership and membership-changing results separately | A hand-picked cut can remove rather than solve a bias |
| Weak Gaussian measurement approximation | Compare Hessian likelihood with a grid/profile/posterior for tails or rejected-epoch cases | Asymmetry/modes and induced hierarchical shifts | Does not validate the parent population |

At the likelihood level, inspect residual **means and covariance**, not just chi-square per degree of freedom. Inflate-noise operations can improve chi-square while leaving a coherent redshift drift unchanged. Coverage must be evaluated for intentionally misspecified alternatives as well as the nominal generator. Independent NIR/spectral/environmental data can add leverage, but each requires its own calibration and selection treatment; it is not automatically bias-free.

The report should therefore retain three separate statuses: (i) reproduced implementation under pinned assumptions; (ii) quantified robustness within a declared alternative family and tolerance; (iii) unresolved physical/model freedom. The current source-grid check supports a narrow numerical statement in (i). The exact identifiability limits and unexecuted discriminating tests keep a universal bias-free claim unproven. No extra deterministic dust or age correction is justified by this audit alone.
