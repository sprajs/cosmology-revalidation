# BayeSN RAISIN: distance-prior identification audit

Status: design frozen and six-object calibration/forward/Fisher gate completed on 2026-09-26; no new dust posterior or cosmology fit. Observed-data likelihood work is blocked by the separately discovered sign-dependent photometry omission. All new artifacts are under `runs/research_2026_09_26/bayesn_distance_identification/`.

The proposed paired experiment can measure **how dust inference changes when external distance information is removed, conditional on the same trained M20 SED and intrinsic-colour priors**. It cannot identify grey luminosity evolution separately from distance, and it cannot make the trained model cosmology independent. The public assets are sufficient to attempt a controlled forward calculation; an exact reproduction of the 2024 hierarchical analysis is presently blocked by missing execution-linked population code, preprocessing, calibration and distance-input vectors.

## Source and execution mapping

[Thorp et al. 2024, §§2–3](https://arxiv.org/html/2402.18624v2) analyse 42 CSP and 37 RAISIN SNe with fixed M20 intrinsic hyperparameters. Their flux likelihood uses external Gaussian distances from flat ΛCDM, H0=73.24 km/s/Mpc and ΩM=.28, with redshift and 150 km/s peculiar-velocity uncertainty. Individual dust follows an exponential AV distribution and, for the population model, a normal RV distribution normalized above .5. Hyperpriors are τA~HalfCauchy(1 mag), μR~Uniform(1,5), σR~HalfNormal(2). A separate common-RV model uses Uniform(1,6). Stan supplies posterior sampling. Selection correction is deferred; these are selected-sample inferences. The paper links the light-curve release and a newer NumPyro implementation, not its exact Stan population execution.

| Asset | Pinned evidence and purpose |
|---|---|
| 2024 paper | `../raisin_differential/2402.18624v2.pdf`; SHA256 `c850cdc495d85c2db177475447ef46add8b34c2bd259aacb588e9868b7bd5db4` |
| Author NumPyro source | [bayesn/bayesn@08eef9e54188f601506ef9f7f47fc65f9f52a946](https://github.com/bayesn/bayesn/tree/08eef9e54188f601506ef9f7f47fc65f9f52a946); 177 extracted regular files individually Git-blob verified in `official-code-acquisition.json` |
| M20 model | [bayesn-model-files@1d77e0be8203f09fff0a6b725abee15011a4d7b9/BAYESN.M20](https://github.com/bayesn/bayesn-model-files/tree/1d77e0be8203f09fff0a6b725abee15011a4d7b9/BAYESN.M20); ten files Git-blob verified. Its YAML is byte-identical to the NumPyro M20 YAML, SHA256 `bcfd03a3716babbfb55b400f8c880223d802d73d0f49d5edd3ebab45b18c33d6` |
| Legacy Stan source | [bayesn-public@4805269494dce00a78456b205d05e4413e9b35db](https://github.com/bayesn/bayesn-public/tree/4805269494dce00a78456b205d05e4413e9b35db); source tree and eleven selected files verified. Two single-SN photometric-distance Stan models are present; neither is the 2024 dust-population model |
| Photometry/calibration | Pinned RAISIN_DataRelease@a383c4b; prior 490-file acquisition manifest and `raisin_differential/frozen-membership.csv`. Same physical 79-object membership; no derived FITRES distances used as photometry |
| M20 training paper | [Mandel et al., arXiv:2008.07538](https://arxiv.org/abs/2008.07538); complete PDF SHA256 `a6a64f84964512bf12d97463c3f405eade5041339c553609b1746ecd8b83620d` |

The [M20 assets](https://github.com/bayesn/bayesn/blob/08eef9e54188f601506ef9f7f47fc65f9f52a946/bayesn/model_files/M20_model/BAYESN.YAML) have nine wavelength knots (3000,4300,5400,6200,7700,10400,12400,16500,18500 Å), six phase knots (−10,0,10,20,30,40 rest-frame days), 42 correlated residual-SED coordinates, M0=−19.5 mag and σ0=.088 mag. Their trained RV=2.886 and τA=.329 mag must **not** replace the dust hyperparameters inferred in the proposed hierarchy. The numerical 3000 Å endpoint is not a validation of u/U: the model INFO explicitly excludes those bands. T21, W22, G25 and G26 are different trained models and cannot be substituted.

### Public wrapper traps

The following are code facts at the pinned NumPyro commit, not inferred paper settings:

* `SEDmodel.__init__` defaults to **T21**, four host devices, and RV truncation 1.2. Explicit M20 and resource settings are necessary.
* `dust_model`, lines1267–1325, uses μR~Uniform(1.2,6), RV≥1.2 and **inferred σ0**. Its dust likelihood is flux Gaussian; it uses precomputed phases and does not sample peak time. It therefore needs a small, explicit statistical adapter to implement the paper equations.
* `fit_model_uniformRV`, lines829–892, samples RV independently inside the SN plate, fixes τA, and uses a five-magnitude distance prior. It is not the paper's shared-RV model.
* `train_model_*` also change the SED and use magnitude-space likelihoods. They are unsuitable for this distance-only comparison.
* The legacy wrapper defaults to S/N>3 and distance σ=100 mag; the modern flux arrays retain negative flux. Neither wrapper establishes the actual 2024 masks. Fixed versus inferred peak time also needs an explicit execution mapping.

See [modern likelihood source](https://github.com/bayesn/bayesn/blob/08eef9e54188f601506ef9f7f47fc65f9f52a946/bayesn/bayesn_model.py#L1267) and [legacy Stan source](https://github.com/bayesn/bayesn-public/blob/4805269494dce00a78456b205d05e4413e9b35db/BayeSNmodel/stan_files/photometric_distance_model.stan#L64).

## What remains inherited from training

M20 training itself used external distances: mostly a low-redshift fiducial Hubble relation, plus eight redshift-independent distances. Fixing the trained mean, shape component and residual covariance retains that training information in **both** proposed arms. The original low redshifts reduce sensitivity to expansion-history shape; this is not the same as demonstrating zero dependence. Intrinsic colour versus dust separation still uses the trained distribution. [Mandel et al., §§2.6–2.8](https://arxiv.org/pdf/2008.07538).

An outcome-free extraction of unique official SN designations from that paper's Table1 finds **29 of the 42 RAISIN CSP objects**, and none of its PS1/DES objects. These are two different 79-object samples. The name join is in `M20-training-name-overlap.csv`; the table lacks sky positions and an execution-linked training file roster was not obtained. This is a paper-input overlap proxy, not verified acceptance by an exact training run. Both pilot CSP objects are such matches. Their pilot use is an internal numerical check, not out-of-training prediction. This overlap is distinct from the K21/SALT overlap documented in the earlier RAISIN audit.

## Grey amplitude can be handled exactly in one dimension

This section is an algebraic derivation from the pinned flux operator, not a new posterior. In `get_flux_batch`, lines656–719, the integrated spectrum is multiplied by 10^(−.4[M0+D]), where D=μ+δM. For fixed θ, ε, AV, RV, peak time, redshift, MW attenuation and bandpasses, write

\[
y\sim N(a f,C),\qquad a=\exp[-K(D-D_*)],\quad K=.4\ln10,
\]

with fixed reference D*=35 mag. Here y and f are FLUXCAL vectors, f is the actual broadband model at D*, C is the measurement covariance in squared flux units, and a is dimensionless. Intrinsic ε remains an explicit 42-dimensional latent variable; **do not add SALT model covariance or count ε twice**. The source likelihood uses quoted independent flux errors, so its baseline C is diagonal. Any correlated measurement extension would be a separately declared model.

For fixed positive-definite C define q=fᵀC⁻¹f, b=fᵀC⁻¹y and c=yᵀC⁻¹y. Then

\[
\log L(a)=\mathrm{const}-\tfrac12(c-2ba+qa^2),\qquad
\hat a=\max(0,b/q),\qquad
\chi^2_{\rm prof}=c-\max(b,0)^2/q.
\]

At b≤0, a=0 is a boundary supremum for a strictly positive amplitude. This profile is a useful likelihood-identification diagnostic. It is not marginal model evidence and does not include the dust/SED integration volume. For a proper flat-a interval [a0,a1], the integral is analytic using normal CDF differences. That **is not a flat-distance prior**: uniform μ has measure da/(Ka).

An unbounded flat μ prior is improper in this flux model: as D→∞, a→0 while L→exp(−c/2)>0. A finite normalization is required even when the tail is numerically very small.

For the original Gaussian distance and independent δM~N(0,σ0²), convolution gives exactly

\[
D\sim N(\mu_{\rm ext},\sigma_{\rm ext}^2+\sigma_0^2).
\]

For the prespecified broad arm μ~Uniform(L,U), unchanged δM implies

\[
p(D)=\frac{\Phi((D-L)/\sigma_0)-\Phi((D-U)/\sigma_0)}{U-L}.
\]

Thus the remaining amplitude integral is only

\[
\mathcal L(\eta)=\int L(\exp[-K(D-D_*)]f(\eta))p(D)\,dD.
\]

The Gaussian-D and convolved-uniform-D integrals are not generally analytic under a Gaussian **flux** likelihood, but require no further SED integrations: after f, q,b,c are computed, each quadrature point is scalar arithmetic. A stable log-CDF difference, likelihood-centered subdivisions plus full prior-tail coverage, doubled quadrature order and independent adaptive reference are required. Replacing this by Gaussian magnitudes or flux ratios would alter the likelihood. If measurement C depends on amplitude, q,b,c are no longer sufficient; retain C(D), its determinant and the full one-dimensional likelihood instead.

`amplitude-algebra-check.json` records artificial-vector checks, including a negative datum: analytic flat-a integration agrees with numerical quadrature to1.53e−16; convolved-uniform-D integration agrees with explicit μ,δM integration to3.70e−14; its prior density integrates to1+2.22e−16. These validate the algebra only, not BayeSN integration, inference or SN recovery.

### Colour information and its limitation

Whiten f and its derivatives with C: u=C^−1/2 f and G=C^−1/2 ∂f/∂η. Profiling amplitude leaves local information

\[
I_{\eta\eta}^{\rm colour}=a^2G^T\left(I-\frac{uu^T}{u^Tu}\right)G.
\]

A purely grey perturbation vanishes under this projector. Wavelength-dependent dust need not vanish, so a broad-distance arm can constrain dust through relative optical/NIR fluxes. However, ∂f/∂RV is proportional to AV near AV=0; RV information then tends to zero. Dust derivatives can also overlap the intrinsic ε and θ directions. Keeping their trained priors is consequential, not an independent colour measurement. No posterior-variance ratio will be described as a unique percentage of colour information. Report D as the light-curve amplitude; the separate μ and δM decomposition remains prior-dependent.

## Frozen paired experiment and measurable gates

`protocol.json` SHA256 **`12a49f9c4da298f2f6fc87c3049dd20bd24712735fecb7a4764217f1f52f0e42`** predates new BayeSN likelihood evaluations. It fixes:

* Original Gaussian distance arm A, versus independent **μ~Uniform(20,50) mag** for every SN in arm B. Width/centre sensitivities are Uniform(15,55) and Uniform(25,50). These are proper stress priors, not measured distance distributions. Keep the same .088-mag grey prior and exact convolution.
* The same physical objects, photometric rows/errors, calibration, MW attenuation, intrinsic model/priors, dust hierarchy, phase convention and selection conditioning. Redshift remains in wavelength and time transformations. Only its amplitude-distance information changes.
* Full-sample target summaries: posterior ΔμR=μR,high−μR,low and ΔτA=τA,high−τA,low (mag), using the fixed42 CSP/37 RAISIN split. Report both arms' means, SDs,68/95% intervals and paired changes with Monte Carlo errors. These are selected-sample dust-population summaries, not cosmic evolution estimates. Common-RV inference is a separate prespecified sensitivity.

The metadata-only six-object numerical pilot sorts within survey by (zHEL,CID), choosing zero-based floor((N−1)/3) and floor(2(N−1)/3):

| Survey | Object | zHEL | One-based rank |
|---|---|---:|---:|
| CSP |2005ki|.0192|14/42|
| CSP |2009ad|.0284|28/42|
| PS1MD |PScD500301|.325|7/19|
| PS1MD |PScJ550202|.422|13/19|
| DES |DES15E2mhy|.4391|6/18|
| DES |DES16S1agd|.504|12/18|

This pilot is too small for a scientific population contrast; it is for calibration, likelihood, recovery and resource gates. No replacement object may be selected for an attractive dust outcome or easier convergence.

The μR estimand follows the paper's parameterization: it is the location of the underlying normal before truncation, not exactly the mean of the truncated distribution. For lower bound r0=.5, that latter mean is μR+σR φ(α)/(1−Φ(α)), α=(r0−μR)/σR. A future analysis should display this implied mean as a labelled secondary quantity while retaining the frozen ΔμR target.

**Calibration gate.** The active RAISIN KCOR input gives CSP B a BD17 reference magnitude9.882; current BayeSN `B_CSP` uses9.896 with the same named standard and filter basename. This fourteen-mmag difference already rules out simply feeding unchanged RAISIN FLUXCAL into modern defaults. The full six-object bridge must trace every bandpass, standard, zero point, epoch/peak and flux error. For identical passbands, a documented coordinate change m_new=m_old+Δ produces F_new=10^(−.4Δ)F_old and C_new=S C_old Sᵀ, accompanied by the same model-reference transformation. Different transmission shapes require SED transport, not a scalar zero-point fix. No empirical per-row rescaling to fitted data is allowed. The original paper's band/epoch bridge is not yet located, so any reconstructed bridge must be labelled conditional.

**Forward and likelihood gates.** Separate zHEL from the corrected external-distance redshift. Export exact μext/σext rather than trusting wrapper defaults. Freeze support and phase masks before likelihoods; preserve negative flux and signed spectral/model values. Require reference-spectrum calibration closure<1mmag and legacy/modern/refined broadband differences<0.1 quoted flux error per included row. No silent S/N cuts, amplitude-dependent error floor, or tmax default changes. Require marginal log-likelihood agreement<1e−6 and scaled gradient agreement<1e−5 between independent integration methods, including b≤0 and remote-boundary cases.

**Recovery and convergence gates.** Use the exact M20 forward operator with the frozen observing metadata. Prespecified injection coordinates are AV=0,.3,.8 mag; RV=2,3.1; θ=0; zero ε initially, then draws from trained ε. Grey stresses0 and+.15 mag test whether arm A reallocates brightness into dust; they do not posit measured grey evolution. Check zero-noise closure, amplitude-profile invariance under grey rescaling and RV non-identification at AV=0. Ten fixed-seed six-SN synthetic datasets screen gross inference failure; they cannot certify coverage. If sampling is warranted, use four chains, initially1000 warmup+1000 draws each, Rhat<1.01, bulk/tail ESS≥400/200, zero divergences, E-BFMI>.3, maximum-depth hits<1%, mean MCSE<.05 posterior SD, and independent-seed agreement within combined MCSE. Report distance-boundary probability and all prior-width sensitivity. Failed numerical gates cannot be repaired by tuning dust priors toward desirable outcomes.

The paper-original RV lower support and the physical requirement that a passive screen have Aλ≥0 throughout the actual integration domain are separate gates. Record any violation; never clamp attenuation or silently change the original prior. An independently checked law-support addendum will accompany the numerical pilot.

## Runtime and resource decision

The live inventory found Python3.12.14, NumPy2.2.6/SciPy1.15.3 in both principal environments. `phase2/hierarchy/.venv` has JAX/JAXlib0.6.2 and NumPyro0.19.0, but lacks sncosmo, extinction, ArviZ and ruamel.yaml. `phase2/env-official` has sncosmo2.12.1/extinction0.4.9, but lacks JAX/NumPyro and other imports. No installed BayeSN or CmdStan was found in the three checked environments. No package or posterior fitting was performed for this design freeze.

The next authorized stage will build an isolated, pinned local environment if required, without changing existing environments. Cap CPU concurrency at two. First budget15 minutes for the bridge/preflight and10 minutes for compilation plus100 gradient evaluations. Only after measured timings and numerical gates, budget at most60 wall minutes for the paired six-object sampling pilot, preserving checkpoints if this cap is exceeded. Extrapolate actual gradient cost, tree depth and ESS/hour before any79-object run. A 300-wavelength default grid and42 residual coordinates per object are known from source; no speed claim substitutes for that benchmark.

The decisive near-term product is an executable, calibrated M20 flux likelihood with an exact proper-amplitude treatment. If the calibration/phase execution link remains unavailable, a declared controlled reconstruction can test distance-prior sensitivity, but cannot be advertised as reproduction of the published dust posterior. Even successful paired recovery leaves trained intrinsic-colour assumptions, sample selection and arbitrary grey luminosity evolution unresolved.

## Executed six-object extension

The design freeze above is retained. Parent-authorized follow-up built an isolated `.venv` containing38 packages from `requirements.lock`, verified by `uv pip sync --require-hashes`; `runtime-isolated.json` records exact versions and lock hash. Existing scientific environments and all pinned source files remain unchanged. The M20 model initialized on CPU in14.7 seconds. No GPU runtime is required or claimed. The additional forward protocol is `forward-protocol.json`; `protocol-addendum.json` clarifies grey invariance, physical-support accounting and the absence of observed fitting.

### All pilot filters traced

`calibration_bridge.py` links all23 survey/band combinations actually present in the six files to KCOR `FilterTrans`, `PrimarySED`, reference magnitude and optional SN-photometry offset. All such optional SN-photometry offsets are zero. `release-filter-config.yaml` and `release-filters/` export these exact released tables into the public model interface, retaining signed transmission tails. Their primary spectra are flux density **per Å**; the obsolete FITS primary-header wording suggesting “per10Å” is contradicted by the array values and the source writer's explicit per-Å assignment (`SNANA/src/kcor.c`, `wr_fits_PRIMARY`). The metadata ledger retains duplicate consistent PEAKMJD headers and the reported flux errors.

These custom filters define a **conditional observation operator in the released coordinate system**, not a recalibration of measured photometry. The SED remains conditioned on M20's original training calibration. No derived residual was used to choose a zero point or filter.

The first `calibration-bridge.csv` follows the public implementation's Simpson quadrature on its native filter knots and compares shapes at released10Å nodes. The refined `calibration-bridge-exact.csv` is the appropriate table for physical reference/shape comparison: it integrates piecewise-linear spectra and transmissions on their **merged knots**, splitting absolute-difference roots for the photon-weight metric. This exposes native subgrid structure hidden by the coarse first comparison. Representative results are:

| Convention or shape check | Result |
|---|---:|
| CSP B, native-minus-release reference magnitude, common precise integral | +14.11436 mmag |
| CSP H, J, Y, j, y native-minus-release reference integral shifts | −23.514,−19.550,−15.331,−19.069,−15.018 mmag |
| PS1 g,r,i,z reference integral shifts | +26.634,+39.534,+28.235,+15.935 mmag |
| CSP B/g/i/n/r and PS1 optical normalized photon-weight half-L1 difference | approximately1e−8 |
| CSP H/Y/j/y half-L1 difference on merged knots | .00316/.00126/.00643/.00118 |
| DES g/r/i/z half-L1 difference | .000988/.003607/.003388/.008393 |
| HST F125W/F160W half-L1 difference | .000236/.000281 |

A reference-integral shift in a **different physical passband** is not a complete photometric transformation; its SN-dependent spectral term remains. Nor are any of these numbers measured photometry biases or evidence that the 2024 authors used mismatched defaults. CSP B illustrates an additional purely numerical issue: the first public-grid Simpson comparison gives+12.2645mmag, whereas the common precise integral gives+14.1144mmag. The independently calculated value agrees with the latter. The released custom-filter zero-normalization Simpson error is at most**.752625mmag among supported bands**, below the frozen1mmag gate. The excluded u band differs by4.805mmag; it is not used to claim closure.

### Forward closure, a failed adapter, and its correction

`forward_geometry.py` uses384 retained epochs, with fixed published PEAKMJD and no S/N criterion:126/116 CSP,41/21 PS1,36/44 DES rows. It excludes u and phases outside the open(−10,40)-day interval, and requires the1%-peak transmission support within3000–18500Å in the rest frame. `supported-band-wings.csv` additionally checks the full tabulated throughput: no retained band's absolute photon weight lies outside this numerical range; at most.000400 lies below3500Å. This is support accounting, not validation of a broader intrinsic population.

The **first adapter failed**: closing the fixed batched numerical arrays inside a JIT-compiled function produced grossly incorrect predictions, although its grey and derivative checks passed. Those failed outputs, code and protocol are preserved under `failed-closed-constant-jit/`; the independent reviewer also saved a snapshot. A one-object source simulator and eager batched calculation agreed. Passing **all numerical inputs explicitly** to a JIT of the source `get_flux_batch` restored agreement. The phase interpolation, spectra, zero ε and band weights were independently compared. `layout-debug.json` localizes the discrepancy, and `adapter-resolution.json` records the replacement. The underlying closure/compiler cause is not fully diagnosed; this is not a claim about the 2024 Stan execution or about physical SED behaviour. Fisher results from the failed adapter are withdrawn.

The corrected adapter adds an eager-versus-compiled check at every resolution. It passes without changing the model, data, priors or thresholds:

| Gate | Observed value |
|---|---:|
| Eager versus full-argument compiled flux, maximum quoted-error units |7.35e−13|
|300→600 wavelength points, maximum change/error|.12272|
|600→1200 points|.033681|
|1200→2400 points|.014039|
| Full versus half finite-difference step, relative whitened-J norm difference|1.54e−9|
| RV derivative at AV=0|exactly zero|
| D→D+.15 mag, relative error against10^(−.4×.15) grey scaling|3.24e−15|
| Noiseless positive-amplitude recovery|floating-point closure|
| Four resolutions plus derivatives, two-CPU affinity|55.7 seconds|

The300-point default fails the first .1-error resolution criterion, so the information calculation uses2400 points. Numerical resolution refinement changes only the in-memory integration-grid setting; the pinned model source is untouched. The grey test rescales a **generated signal**, or equivalently model reference with inverse amplitude. It does not claim that multiplying observed data alone while fixing C leaves a noisy likelihood unchanged.

### Local dust information after amplitude and intrinsic variation

The following are **Fisher design scales, not posterior constraints**, evaluated at prescribed AV=.3mag, RV=3.1, θ=0, ε=0 and fiducial D(zHD). Only the quoted errors, observing metadata and generated flux enter. Neither measured flux nor fitted dust enters the objective. The external-amplitude precision is1/(σext²+.088²); the alternative sets that precision to zero. Both columns retain unit-normal θ and42 whitened ε priors, fixed peak time and the model's zero residual-SED endpoint constraints. The absence of an external distance prior therefore does not mean absence of luminosity/colour assumptions.

| Object | Local σ(AV), no distance / external (mag) | Local σ(RV), no distance / external |
|---|---:|---:|
|2005ki|.05463/.05120|.9955/.9590|
|2009ad|.05859/.05317|1.205/1.135|
|PScD500301|.25665/.12067|3.446/2.331|
|PScJ550202|.22484/.13896|5.437/4.536|
|DES15E2mhy|.21917/.11615|3.088/2.176|
|DES16S1agd|.29855/.15020|4.096/2.893|

With **unconstrained θ and ε**, machine-epsilon rank thresholds give apparent high-z dust-information eigenvalues around10^−31–10^−27. Finite-difference noise can inflate nuisance rank, so `free-intrinsic-rank-sensitivity.json` also reports relative SVD thresholds1e−12,1e−10,1e−8 and1e−6: the largest residual high-z dust eigenvalue is at most**3.35e−7** across this grid (at1e−8, at most1.26e−8). The defensible conclusion is a strong local near-degeneracy under unrestricted intrinsic variation, not an exact physical identity. These objects have only21–44 retained measurements against44 free nuisance coordinates, so this is a flexible, nearly or fully saturated stress. For the two CSP objects the eigenvalues are approximately(5.80e−5,.00912) and(1.47e−5,.00404), stable except for a small change at the coarsest first-object rank threshold. Ranks, singular values and projected derivative norms are saved in `fisher-rank-details.json`. This stress can include physically implausible intrinsic changes; it establishes a dependence on the trained restrictions, not a competing realistic population model. No ratio in this table is a fraction of colour information. Nonlinear recovery with latent ε and population posteriors remains unexecuted.

The independent [next-experiment review](identification-next-experiment.md) checks the profile/marginalization distinction, the calibration sign, Fisher algebra and training/selection limitations. Its `bayesn-dust-support.json` also verifies the source F99 spline against an independent natural cubic interpolation: at RV=.5, minimum Aλ/AV=−.299924 near7214Å, with negative attenuation from approximately6180–9367Å. For this implemented law on3000–18500Å, nonnegative attenuation requires RV≳.683448. This is a necessary support bound for this law, **not a universal physical RV prior**. No posterior mass is inferred and no curve is clamped.

### Observed likelihood is presently blocked

The separate signed-measurement audit finds all23,007 rows across117 RAISIN photometry files positive. Matching17 DES aliases to released DIFFIMG identifies negative original epochs omitted from RAISIN, including56 within the historical optical fit phase range while811 corresponding positive epochs are retained. In the(−10,40)-day phase window the corresponding counts are47/47 negatives absent and836/836 positives retained, before this pilot's filter-support restrictions. See `../raisin_flux_sign_audit/retention-result.json` and `lineage-handoff.json`. The precursors and retained data also differ in some quoted errors; adding missing rows is not automatically an exact restoration of the executed measurement version. The current forward/Fisher gate therefore uses an **available-row design conditional on quoted errors**, not a validated full measurement schedule. An ordinary untruncated Gaussian likelihood on these released positive-only observations is not justified by the successful model-only checks. The upstream sign-selection mechanism, flags, complete epoch list and within-model-phase effects must be resolved before any observed dust posterior, noise recovery or held-out-NIR likelihood. This is a measurement-validity boundary, not a request for permission.

A next bounded step can verify the corrected explicit-argument operator on synthetic parameter recovery and proper-D quadrature, using complete signed observing metadata when available. Optical→held-out-NIR discrimination should then account for any NIR-dependent sample/epoch selection and treat the29 M20-overlap CSP objects as internal checks. The56-second forward benchmark supports further numerical checks within the planned budget, but provides no HMC ESS/hour forecast. **No population sampling is recommended or started while the measurement gate is unresolved.**
