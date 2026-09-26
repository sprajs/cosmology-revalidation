# From emitted radiation to calibrated observations and inferred distance

This chapter derives the observation equations used to connect the [luminosity](luminosity.md), [extinction](extinction.md), and [cosmology](cosmology-populations.md) chapters. Equation labels O01–O23 are stable references for subsequent assumption tests. The derivations concern physical observables first; the final sections identify the empirical and statistical approximations in the current implementation.

## 1. Quantities, units and the observation equation

Let \(L_{\lambda_e}(t_e)\) be emitted spectral luminosity, in erg s\(^{-1}\) Å\(^{-1}\), and \(f_{\lambda_o}(t_o)\) observed spectral energy flux, in erg s\(^{-1}\) cm\(^{-2}\) Å\(^{-1}\). Luminosity is integrated over outward directions; for an anisotropic explosion viewed in direction \(\boldsymbol n\), use the isotropic-equivalent \(4\pi\,dL/d\Omega\) in that direction. A single viewing direction does not determine the angle-integrated luminosity of an asymmetric explosion.

The observer/source mappings in an unperturbed expanding spacetime are

\[
\lambda_o=(1+z)\lambda_e,\quad \nu_e=(1+z)\nu_o,\quad
 dt_o=(1+z)dt_e,\quad E_{\gamma,o}=E_{\gamma,e}/(1+z). \tag{O01}
\]

An emitted energy packet \(L_{\lambda_e}dt_e d\lambda_e\) loses one factor \(1+z\) in photon energy and is spread over one factor in arrival time and one in wavelength interval. The sphere's present area is \(4\pi D_M^2\). Using \(D_L=(1+z)D_M\) therefore gives

\[
f_{\lambda_o}(t_o)=\frac{L_{\lambda_e}(t_e)}{4\pi D_L^2(1+z)},\qquad
f_{\nu_o}(t_o)=\frac{(1+z)L_{\nu_e}(t_e)}{4\pi D_L^2}. \tag{O02}
\]

The opposite spectral prefactors are required by \(f_\nu |d\nu|=f_\lambda d\lambda\), not competing conventions. Integrating O02 recovers \(F_{\rm bol}=L_{\rm bol}/(4\pi D_L^2)\). A further \((1+z)^{-2}\) applied to O02 would double count energy loss and time dilation already contained in luminosity distance. [Hogg's distance definitions](https://arxiv.org/abs/astro-ph/9905116) provide the reference convention; these factors follow directly from the packet calculation above.

For a foreground-screen model, and a single unresolved achromatic lens magnification \(\mathcal A\), the useful full equation is

\[
\boxed{f_{\lambda_o}(t_o)=
\frac{\mathcal A\,L^{\rm int}_{\lambda_o/(1+z)}((t_o-t_{\rm exp,o})/(1+z))}
 {4\pi D_L^2(1+z)}
\exp[-\tau_{\rm host}(\lambda_o/(1+z))-\tau_{\rm MW}(\lambda_o)-\tau_{\rm other}(\lambda_o,z)] .} \tag{O03}
\]

This is the physical forward relation that a phenomenological light-curve model replaces in part. It assumes an unresolved source, geometric optics, negligible wavelength-dependent lensing and delay structure, and no appreciable scattered light returned to the aperture. General scattering requires the transport equation in the extinction chapter. Multiple unresolved lensed images require a sum with individual \(\mathcal A_j\) and time delays, not one factor. Intrinsic angle dependence and microlensing of a wavelength-dependent source can invalidate the achromatic approximation.

For intergalactic absorbers, each slice sees \(\lambda_o/(1+z')\), giving, for proper density and opacity,

\[
\tau_{\rm IGM}(\lambda_o,z)=\int_0^z
\rho_d(z')\kappa_{\rm ext}[\lambda_o/(1+z'),z']
\frac{c_{\rm light}\,dz'}{(1+z')H(z')}. \tag{O04}
\]

\(\rho_d\kappa\) has inverse-length units. This term is a possible propagation extension, **not a fitted intergalactic-dust model present in this workspace**. Cosmological transport does not establish its amplitude.

For real supernova catalogues, \(z_{\rm HEL}\), \(z_{\rm CMB}\) and \(z_{\rm HD}\) differ. The spectral/phase engine uses measured heliocentric redshift; the released-distance approximation uses \((1+z_{\rm HEL})D_M(z_{\rm HD})\). O01–O03 were derived with one unperturbed redshift; substituting two frames is a specified peculiar-motion convention, not an exact relativistic treatment of all Doppler/beaming and gravitational perturbations. See the cosmology chapter for that boundary.

## 2. Photons, filters and magnitudes

For a photon-counting detector with collecting area \(A_{\rm tel}\), exposure duration \(\Delta t\), and dimensionless throughput \(T_b(\lambda)\), the expected source photoelectrons are

\[
N_{e,b}=A_{\rm tel}\int_{\rm exposure}dt_o\int d\lambda_o\,
\frac{\lambda_o}{hc_{\rm light}}T_b(\lambda_o,t_o)f_{\lambda_o}(t_o). \tag{O05}
\]

The factor \(\lambda/(hc_{\rm light})\) converts energy to photon number. Throughput must include quantum efficiency if \(N_e\) denotes electrons. Use a consistent wavelength unit for \(hc\) and \(f_\lambda d\lambda\); for Å, \(hc\simeq1.98644586\times10^{-8}\) erg Å. A filter already tabulated with photon weighting must not receive the factor twice. Holding flux/throughput fixed during a short exposure reduces the time integral to \(\Delta t\). This approximation is excellent only when the transient and atmosphere are sufficiently stable during the exposure.

Define \(I_b[f]=\int\lambda T_b(\lambda)f_\lambda d\lambda\). For a primary standard spectrum \(f^P_\lambda\) with assigned band magnitude \(m^P_b\), calibration is

\[
m_b=m^P_b-2.5\log_{10}\frac{I_b[f]}{I_b[f^P]},\qquad
F_{{\rm cal},b}=10^{0.4(Z_b-m_b)}
=10^{0.4(Z_b-m^P_b)}\frac{I_b[f]}{I_b[f^P]}. \tag{O06}
\]

Area, exposure time and \(hc\) cancel in the calibrated ratio when source and standard definitions use the same response. Here the saved DES flux convention is \(Z_b=27.5\); this number sets numerical units, not a physical law or a universal value for every SNANA dataset. The AB reference has constant \(f_\nu\simeq3631\) Jy and thus \(f_\lambda=c_{\rm light} f_\nu/\lambda^2\), with wavelength-unit conversion. A constant \(f_\lambda\) is not an AB standard. The [primary photometric definitions of Hogg et al.](https://arxiv.org/html/astro-ph/0210394v1) specify the band/standard distinction.

The usual magnitude-flux differential and distance modulus are

\[
\delta m=-\frac{2.5}{\ln10}\frac{\delta F}{F},\qquad
\mu=5\log_{10}\frac{D_L}{10\,\mathrm{pc}}
=5\log_{10}\frac{D_L}{\mathrm{Mpc}}+25. \tag{O07}
\]

The differential is a small-error approximation; \(m=-2.5\log_{10}F+Z\) is defined only for positive flux. A noisy background-subtracted flux may be negative without negative emitted radiation. Such measurements belong in a flux likelihood rather than being dropped or logarithmically transformed.

For an observed band X and emitted band Y, write \(m_X=M_Y+\mu+K_{XY}\) in the absence of dust. Eliminating luminosity normalization using O02 and O06 gives an explicit photon-counting cross-band correction:

\[
K_{XY}=-2.5\log_{10}\!\left[
\frac{1}{1+z}
\frac{\int\lambda_o T_X(\lambda_o)L_{\lambda_o/(1+z)}d\lambda_o}
 {\int\lambda_o T_X(\lambda_o)f^0_{X,\lambda_o}d\lambda_o}
\frac{\int\lambda_e T_Y(\lambda_e)f^0_{Y,\lambda_e}d\lambda_e}
 {\int\lambda_e T_Y(\lambda_e)L_{\lambda_e}d\lambda_e}
\right]. \tag{O08}
\]

\(f_X^0,f_Y^0\) are zero-magnitude reference spectra in their respective systems. The emitted spectrum in both integrals is evaluated at the same corresponding rest-frame phase. This derivation permits different bands and standards. Same-band constant \(L_\nu\) with an AB standard gives \(K=-2.5\log_{10}(1+z)\), a useful sign test. Dust must be consistently included either in the spectrum or as separately defined band attenuation. SALT forward integration already performs the redshift/filter mapping; adding a second K correction would count it twice.

**Implementation verdict.** `scripts/phase2/independent_flux/engine.py:Engine.prepare/flux` uses \(\lambda T\), the primary-SED denominator, \((1+z)^{-1}\), observer-frame MW extinction and observer time divided by \(1+z\). Those factors are consistent with O02–O06. It uses the local SALT normalization \(10^{-12}\) and `MAG_OFFSET=.27`; these are calibrated model conventions, not constants derived from an explosion. `sncosmo_flux` provides a separate integrator using shared assets. Agreement checks numerical implementation, not the correctness of the training data or the underlying population.

**Correction made in this audit.** The direct integrator previously discarded nonzero filter transmission outside the rest-wavelength model grid and returned a partial-band flux as though it were a full-band prediction. It now rejects an evaluated band lacking full spectral support, while allowing unused filters to be prepared. Zero-throughput endpoints are harmless. Tests include a partially covered g band at \(z=1\), an entirely unsupported g band at \(z=2\), and all 528 saved reference epochs. Extrapolating an unknown spectrum is a scientific choice requiring a model, not a numerical default.

## 3. What actually lands on a detector

For exposure e and pixel p, a minimal linear image model is

\[
\Lambda_{ep}=F_e P_{ep}(\boldsymbol x_{\rm SN})+
 (P_e*G)_p+B_{ep}+D_{ep},\qquad \sum_pP_{ep}=1, \tag{O09}
\]

where quantities are expected electrons, \(P_e\) is the normalized point-spread function including pixel integration, \(G\) is the constant host scene, B sky and D dark current. Exposure/photometric scaling is implicit in the scene term here and must be explicit when exposures have different scales. The raw charge obeys the ideal photo-detection approximation

\[
n_{ep}\sim\mathrm{Poisson}(\Lambda_{ep}),\quad
 d_{ep}=n_{ep}+\epsilon_{{\rm read},ep},\quad
\mathrm{Var}(d_{ep})=\Lambda_{ep}+\sigma_{{\rm read},ep}^{2}. \tag{O10}
\]

This assumes independent photon arrivals and ideal detector response. Gain converts electrons to ADU; variance in ADU divides by gain squared. Nonlinearity, flat-field uncertainty, correlated electronics, imperfect PSF/host fits, and resampling require additional terms. Calibrated flux estimates are not raw Poisson counts. Gaussian errors become an approximation after extraction. DES's [scene-modelling photometry paper](https://arxiv.org/abs/1811.02378) identifies the fitted transient-plus-host structure; the [five-year release](https://arxiv.org/abs/2406.05046) distinguishes SMP from difference-image products.

A linear resampling or flux estimator \(\widehat{\boldsymbol f}=W\boldsymbol d\) propagates noise by

\[
C_f=W C_d W^T. \tag{O11}
\]

If difference images reuse a reference image n, then \(d_e=s_e-a_e n\) has \(\mathrm{Cov}(d_e,d_{e'})=a_ea_{e'}C_n\) for otherwise independent exposures. This shows why repeated epochs or a shared template are not automatically independent measurements. Simultaneously fitting a host also couples estimated epochs. An empirical excess scatter or error floor does not specify which of these mechanisms caused it.

`pixel_noise.py` and `noise_tests.py` inspect calibrated data and simple apertures; neither implements a complete image-level likelihood. The annular median-background uncertainty currently uses a Gaussian approximation with a \(\pi/2\) factor. It is not a general physical detector law. The repository therefore cannot certify O09–O11 end to end from its flux fits alone.

## 4. Where SALT enters and where physics stops determining the formula

The implemented rest-frame spectral model is

\[
S_\lambda(p)=x_0\,[M_0(p,\lambda)+x_1M_1(p,\lambda)]
10^{-0.4c_{\rm SALT}CL(\lambda)},\qquad
p=(t_o-t_{0,o})/(1+z_{\rm HEL}). \tag{O12}
\]

The matrices \(M_0,M_1\) and colour law CL are learned surfaces; \(x_1\) is dimensionless and \(c_{\rm SALT}\) a fitted colour coordinate. This phase origin is fitted B-band maximum, not the explosion time used in radioactive-decay and ejecta equations. The [SALT3 implementation specification](https://sncosmo.readthedocs.io/en/stable/api/sncosmo.SALT3Source.html) establishes the executable convention. A fitted \(x_0\) absorbs distance and luminosity normalization. The engine does not independently insert physical \(L\), \(D_L\), nickel mass and an ejecta diffusion time into O12.

In the local colour-law convention, CL(4302.57 Å)=0 and CL(5428.55 Å)=−1. Thus positive colour changes V relative to the B anchor at fixed \(x_0\). It is not equivalent to applying a positive dust optical depth to every wavelength. The alternative `family='f99'` subtracts the F99 value at the B anchor and permits signed colour. It is an **empirical colour-shape replacement with free amplitude**, not an absolute passive dust screen. The `phase_colour` alternative adds \(d\tanh(p/20\,\mathrm{day})\) to the colour parameter; that is a diagnostic basis, not a transport solution.

The saved SNANA convention is

\[
m_B=-2.5\log_{10}x_0+10.635,\qquad
\delta m_B=-\frac{2.5}{\ln10}\delta\ln x_0. \tag{O13}
\]

The additive 10.635 depends on the pinned model/calibration normalization. Substituting an unrelated sncosmo amplitude zero point would be an error. For an absolutely calibrated rest-frame passband, a physical magnitude instead follows O06–O08 applied to a spectrum at 10 pc. A SALT coordinate bearing the name mB must not be assumed to measure bolometric output, nickel mass, or a universal physical B-band luminosity without that calibration link.

The distance estimator is typically written, with an explicitly defined correction sign,

\[
\widehat\mu=m_B+\alpha x_1-\beta c_{\rm SALT}-M
+\Delta_{\rm host}-\Delta_{\rm bias},\quad
\Delta_{\rm bias}=\mathbb E[\widehat\mu_{\rm pre-bias}-\mu_{\rm truth}\mid\text{simulated selected sample}]. \tag{O14}
\]

Published columns may use a different named/sign convention; reconcile them algebraically before applying anything. \(\alpha,\beta,M\), a host step and a bias calibration are empirical relations/estimators. They do not follow from radioactive decay, Maxwell's equations or Einstein's equations. A bias derived under one population/selection model is conditional on that generator, and may require recalculation when the population changes.

A useful *illustrative* intrinsic-plus-dust magnitude model gives

\[
\begin{split}
m_B&=\mu+M_0-\alpha_{\rm true}x_1+\beta_{\rm int}c_{\rm int}+A_B+f(T,Z,\mathcal B)+\epsilon,\\
c_{\rm SALT}&\simeq c_{\rm int}+E,\qquad A_B\simeq R_BE,\\
\widehat\mu-\mu&\simeq(M_0-M)+(\alpha-\alpha_{\rm true})x_1
 +(\beta_{\rm int}-\beta)c_{\rm int}+(R_B-\beta)E
 +f+\Delta_{\rm host}-\Delta_{\rm bias}+\epsilon. \tag{O15}
\end{split}
\]

Here T is progenitor formation-to-explosion delay, Z progenitor birth metallicity and \(\mathcal B\) binary/channel properties; f is an additional intrinsic magnitude response. The approximations are essential: SALT colour is not an exact physical B−V excess and broad-band extinction varies with spectrum/phase. After conditioning out stretch and other brightness terms (or assuming they and the luminosity residual are uncorrelated with total colour), two independent centred latent colours give a regression coefficient

\[
\beta_{\rm eff}=\frac{\beta_{\rm int}\mathrm{Var}(c_{\rm int})+
R_B\mathrm{Var}(E)}{\mathrm{Var}(c_{\rm int})+\mathrm{Var}(E)}. \tag{O16}
\]

Correlation adds cross terms, measurement error changes regression, and selection changes moments. Therefore changing the colour mixture changes fitted \(\beta\) at fixed microscopic grains. The local hierarchy explicitly uses this effective-coordinate interpretation; [Mandel et al.](https://arxiv.org/abs/1609.04470) is a primary example of the distinct intrinsic and dust components.

## 5. Propagation perturbations and what uncertainties mean

A weak achromatic magnification produces

\[
\Delta m_{\rm lens}=-2.5\log_{10}\mathcal A,\qquad
\mathcal A=\frac{1}{|(1-\kappa)^2-|\gamma|^2|}\simeq1+2\kappa,
\qquad \Delta m\simeq-\frac{5}{\ln10}\kappa. \tag{O17}
\]

Here \(\kappa\) is lensing convergence and \(\gamma\) shear; this is the geometric-optics thin-lens Jacobian, not an implemented line-of-sight matter reconstruction. The propagation derivation and averaging subtleties are described by [Bartelmann & Schneider](https://arxiv.org/abs/astro-ph/9912508). Even if a specified complete averaging measure conserves mean flux, it need not conserve mean logarithmic magnitude; selection further changes that average. A Gaussian redshift-dependent lensing scatter is a phenomenological summary, not a demonstration that every mean is zero.

At low redshift, \(D_L\simeq c_{\rm light}z/H_0\), so an unmodelled radial velocity dispersion has the familiar leading-order effect

\[
\sigma_\mu\simeq\frac{5}{\ln10}\frac{\sigma_v}{c_{\rm light}z},\qquad
\sigma_{\mu,z}^2=\nabla_{\boldsymbol z}\mu^T C_{\boldsymbol z}\nabla_{\boldsymbol z}\mu. \tag{O18}
\]

The second expression is the general small-error propagation for the specific catalogue's correlated redshift variables. The first omits higher-redshift Doppler/geometry terms and coherent flows. It is not valid to add a second copy to a released covariance already containing it. The phase-2 pre-BBC population data contract explicitly lists redshift, peculiar velocity, lensing and shared calibration as unfinished additions; the released-distance phase-1 covariance is a different data level.

For shared physical nuisance parameters \(\eta\), linearization gives

\[
\delta\boldsymbol f=G\delta\boldsymbol\eta,\quad
C_{f,\eta}=G S_\eta G^T,\quad
C_{\mu,\eta}=D S_\eta D^T. \tag{O19}
\]

A common calibration shift or dust-map normalization produces cross-object covariance. Replacing a shared mode by independent diagonal errors changes the assumed physical uncertainty and lets it average down incorrectly. A covariance also does not correct an unknown nonzero mean shift.

The local response calculation is a weighted linear least-squares identity:

\[
R=(J^TC^{-1}J)^+J^TC^{-1}G,\quad
G_\perp=G-JR,\quad
D=(1,\alpha,-\beta,0)R, \tag{O20}
\]

where \(J=\partial f/\partial(m_B,x_1,c,t_0)\), and the implemented whitened SVD provides the appropriate pseudoinverse. A negligible residual \(G_\perp\) means the perturbation can be hidden in fitted light-curve parameters, not that its standardized distance shift D is negligible. These local derivatives hold the covariance, masks, trained surfaces and \(\alpha,\beta\) fixed. They do not include retraining, changed selection or BBC.

For Gaussian flux errors the normalized likelihood is

\[
-2\log p(\boldsymbol y\mid\theta)=
r^TC(\theta)^{-1}r+\log\det C(\theta)+n\log(2\pi). \tag{O21}
\]

The determinant can be dropped in a parameter comparison only when it is constant. Recomputing flux-dependent noise while minimizing only the quadratic is a different objective, not an equivalent likelihood. `independent_flux/noise.py` intentionally compares both objectives as sensitivity arms; the existence of the unnormalized arm is not itself an undisclosed physics error. SALT training variances, error floors, and intrinsic-scatter prescriptions remain empirical components of C.

Selection changes the physical population seen by the instrument. If \(\epsilon(y,z)\) is a valid measurement-level selection probability, then

\[
p(y\mid z,S=1,\theta)=\frac{p(y\mid z,\theta)\epsilon(y,z)}
 {\int p(y'\mid z,\theta)\epsilon(y',z)dy'}. \tag{O22}
\]

When trigger, host redshift, quality or classification depend on additional latent/observed variables, those variables must be included in the integral. The numerator form is not universally reducible to a magnitude cut. The hierarchy's probit selector is labelled synthetic; it does not establish DES completeness. Population priors and classifier probabilities are not fundamental physical laws.

Finally, arbitrary luminosity drift is observationally degenerate with the inferred distance shape:

\[
m_{\rm std}(z)=M+\mu(z)+e(z),\qquad
\mu_2=\mu_1+g(z),\quad e_2=e_1-g(z). \tag{O23}
\]

No amount of algebraic correctness removes this freedom without information on e. The emitted-energy equations constrain possible explosion models, and multiband/host/independent probes can constrain population changes, but SALT alone does not derive \(e(z)=0\). This is why a correct distance calculation is insufficient to certify the physical assumptions behind a cosmology result.

## 6. Verification boundary

The new `scripts/physics_audit/observation_check.py` checks conservation of integrated spectral energy, wavelength/frequency photon-integral agreement, the flat-\(L_\nu\) K-correction sign, inverse-square scaling, calibrated gray dimming, retained reference predictions and rejection of incomplete band support. Its saved output is `runs/physics_audit/observation-check.json`. These are independent identities and regression checks; they are not a new survey calibration or a rerun of the cosmological posterior.

This chapter has inspected the local integration, fitter, population interface, response and simple detector diagnostics. It has not established the correctness of every external code path, the original detector calibrations, every uncertainty in the released products, or a first-principles explosion prediction for any individual observed SN.
