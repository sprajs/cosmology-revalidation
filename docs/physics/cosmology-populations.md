# Physical equations: cosmology, stellar populations and auxiliary probes

Implementation references in this chapter identify the preserved source snapshot; see [physical foundations](README.md) for the current execution boundary.

This chapter establishes the physical equations used by the local cosmology, age-mapping and directional calculations. Equation identifiers `Cxx` and `Pxx` distinguish this chapter's equations from the original papers' numbering. Equations are written independently, with explicit units, domains and assumptions; agreement with an empirical fit is not treated as a derivation from fundamental physics.

The main result is conditional correctness: the implemented flat, radiation-free FLRW distances, CPL evolution, cosmic clock, BAO distance ratios and normalized cosmic-SFH–DTD convolution have the expected equations and signs. Their astrophysical assumptions are not thereby established. Two actual software faults were found in the piecewise-deceleration implementation: unrequested extrapolation silently saturated the comoving integral, and the public `efunc(..., model='qbins')` dispatched to the wrong model. The audit correction supplies a consistent piecewise expansion function and rejects evaluation beyond the specified redshift domain. These corrections do not change the sampled in-domain physical model.

Scope: the source traced here is [core.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/cosmology/core.py), [reference_chains.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/cosmology/reference_chains.py), [csfh_dtd.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/mapping/csfh_dtd.py), [controlled.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/mapping/controlled.py), [audit.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/directional/audit.py), [core.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/phase2/hierarchy/core.py), and the age-semantics comparison in [pantheon_age_semantics.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/assumption_audit/pantheon_age_semantics.py). The chapter also identifies foundational physics omitted from these compressed models. It does not claim that the full third-party CMB, structure-growth, stellar-population-synthesis, or SNANA implementations have been independently verified.

## C01–C04. Geometry, gravity and energy conservation

Use cosmic proper time \(t\), a dimensionless scale factor \(a(t)\) with \(a(t_0)=1\), speed of light \(c\), and curvature constant \(K\) with units length\(^{-2}\). Homogeneity and isotropy give the FLRW metric

\[
\boxed{ds^2=-c^2dt^2+a(t)^2\left[\frac{dr^2}{1-Kr^2}+r^2d\Omega^2\right].}
\tag{C01}
\]

This is a symmetry assumption about the large-scale geometry. It is not inferred independently for every supernova. Define \(H=\dot a/a\) in inverse time, \(q=-\ddot a/(aH^2)\), and \(j=\dddot a/(aH^3)\). Expansion means \(H>0\); accelerated expansion means \(q<0\). Neither definition requires an assumption about dark energy.

For a perfect fluid, let \(\epsilon\) be **energy density**, in J m\(^{-3}\), and \(p\) pressure, in the same units. With four-velocity \(u^\mu u_\mu=-c^2\),

\[
T_{\mu\nu}=\frac{\epsilon+p}{c^2}u_\mu u_\nu+p g_{\mu\nu},
\qquad G_{\mu\nu}=\frac{8\pi G}{c^4}T_{\mu\nu}.
\tag{C02}
\]

Here a cosmological constant, if present, is included as a fluid with \(\epsilon_\Lambda=\Lambda c^4/(8\pi G)\) and \(p_\Lambda=-\epsilon_\Lambda\); do not add it again on the left. Evaluating the temporal and spatial Einstein equations in C01 yields

\[
\boxed{H^2=\frac{8\pi G}{3c^2}\epsilon-\frac{Kc^2}{a^2},
\qquad \frac{\ddot a}{a}=-\frac{4\pi G}{3c^2}(\epsilon+3p).}
\tag{C03}
\]

Both sides have units s\(^{-2}\). If instead \(\rho=\epsilon/c^2\) denotes mass-equivalent density, the acceleration combination is \(\rho+3p/c^2\), not \(\rho+3p\). Thus the shorthand in `docs/first-principles.md` has been corrected to make the density and \(c\) convention explicit. Curvature occurs in the first Friedmann equation but has no separate term in the acceleration numerator.

Local stress-energy conservation follows from \(\nabla_\mu T^{\mu\nu}=0\). Equivalently, energy in a comoving volume obeys \(d(\epsilon a^3)=-p\,d(a^3)\):

\[
\boxed{\dot\epsilon+3H(\epsilon+p)=0,
\quad \epsilon(a)=\epsilon_0\exp\!\left[-3\int_1^a[1+w(a')]d\ln a'\right],
\quad w=p/\epsilon.}
\tag{C04}
\]

The second relation holds for each **separately conserved** component. Interacting components obey \(\dot\epsilon_i+3H(\epsilon_i+p_i)=Q_i\), with \(\sum_iQ_i=0\), and need different evolution laws. Nonrelativistic matter gives \(\epsilon_m\propto a^{-3}\); freely propagating radiation gives \(\epsilon_r\propto a^{-4}\); a cosmological constant is constant. These results combine GR and the matter model, not kinematics alone. See [Carroll's cosmological-constant review](https://ned.ipac.caltech.edu/level5/Carroll/Carroll2.html) for the standard dynamical framework.

## C05–C09. Redshift, dark energy and deceleration

For comoving emission and reception in FLRW, successive wave crests along the same null path satisfy \(dt_o/a_o=dt_e/a_e\). Consequently

\[
\boxed{1+z=\frac{\lambda_o}{\lambda_e}=\frac{a_o}{a_e},
\quad dt_o=(1+z)dt_e,
\quad \frac{dt}{dz}=-\frac1{(1+z)H(z)}.}
\tag{C05}
\]

The source time is therefore shorter than the observed light-curve interval by \(1+z\). This redshift is the background cosmological redshift; measured redshift also includes local motions and gravitational perturbations.

With \(E=H/H_0\), critical energy density \(\epsilon_{c0}=3c^2H_0^2/(8\pi G)\), and \(\Omega_K=-Kc^2/H_0^2\), C03–C04 give

\[
\boxed{E^2(z)=\Omega_m(1+z)^3+\Omega_r(1+z)^4+\Omega_K(1+z)^2
 +\Omega_{\rm DE}f_{\rm DE}(z),}
\quad
f_{\rm DE}=\exp\!\left[3\int_0^z\frac{1+w(u)}{1+u}du\right],
\tag{C06}
\]

with \(\Omega_m+\Omega_r+\Omega_K+\Omega_{\rm DE}=1\). Massive neutrinos require their evolving energy density rather than a single matter or radiation power across the relativistic transition. The local distance calculator explicitly sets \(\Omega_K=\Omega_r=0\) and \(\Omega_{\rm DE}=1-\Omega_m\). It cannot be extended to the sound horizon or recombination by merely integrating to higher redshift.

The CPL model is an empirical two-parameter choice,

\[
\boxed{w(a)=w_0+w_a(1-a)=w_0+w_a\frac{z}{1+z},
\qquad f_{\rm DE}(z)=(1+z)^{3(1+w_0+w_a)}
 e^{-3w_a z/(1+z)}.}
\tag{C07}
\]

To derive the last expression, write \(w(z)=(w_0+w_a)-w_a/(1+z)\) inside C06; the integral is \(3(1+w_0+w_a)\ln(1+z)-3w_a z/(1+z)\). The exponential sign in `efunc` and the equivalent expression in `csfh_dtd.clock` are correct. CPL is a parameterization, not a microphysical field equation. The cut \(w_0+w_a<0\) makes DE subdominant relative to matter asymptotically in the radiation-free model, but does not validate all early-universe physics. Its future divergence as \(a\to\infty\) prevents treating distant-future extrapolation as an observation. [Linder's original parameterization](https://arxiv.org/abs/astro-ph/0208512) supplies the model context.

Differentiating \(H\), or using C03, gives two independently checkable expressions:

\[
\boxed{q(z)=-1+(1+z)\frac{E'(z)}{E(z)}
=\frac{\Omega_m(1+z)^3+2\Omega_r(1+z)^4+
[1+3w(z)]\Omega_{\rm DE}f_{\rm DE}(z)}{2E(z)^2}.}
\tag{C08}
\]

For flat matter plus DE without radiation,

\[
\boxed{q_0=\frac12+\frac32w_0(1-\Omega_m).}
\tag{C09}
\]

`qvalue` and the low-redshift derived quantities in `reference_chains.py` use precisely these approximations. The public CMB chains were generated with more complete early-time physics; recomputing their \(q\) with C09 is an explicitly radiation-neglecting derived summary, not the exact original Boltzmann-model acceleration at every epoch. Rejecting \((w_0,w_a)=(-1,0)\) does not imply \(q_0>0\). The acceleration criterion is total \(\epsilon+3p<0\), not merely negative \(w\) of one component.

## C10–C15. Distances, flux, curvature and velocity conventions

Radial null propagation gives \(\int c\,dt/a\). Combining it with C05 defines the dimensionless radial integral \(X(z)\) and transverse distance

\[
X(z)=\int_0^z\frac{du}{E(u)},\qquad
D_M=\frac c{H_0}\mathcal S_{\Omega_K}(X),
\quad \mathcal S_{\Omega_K}(X)=
\begin{cases}
\sinh(\sqrt{\Omega_K}X)/\sqrt{\Omega_K},&\Omega_K>0,\\
X,&\Omega_K=0,\\
\sin(\sqrt{-\Omega_K}X)/\sqrt{-\Omega_K},&\Omega_K<0.
\end{cases}
\tag{C10}
\]

An isotropic source emits \(L\,dt_e\). Each photon's energy is reduced by \(1+z\), arrivals are spread over an interval longer by \(1+z\), and the present wavefront area is \(4\pi D_M^2\). Photon-conserving propagation therefore gives

\[
\boxed{F_{\rm bol}=\frac{L_{\rm bol}}{4\pi D_M^2(1+z)^2}
=\frac{L_{\rm bol}}{4\pi d_L^2},
\quad d_L=(1+z)D_M=(1+z)^2D_A.}
\tag{C11}
\]

The inverse-square equation defines luminosity distance. The last equality is distance duality, assuming metric light propagation and photon conservation. Dust, unresolved lensing and source anisotropy must be handled as distinct effects rather than silently redefining background \(H(z)\). The observer's band flux also requires the redshifted source spectrum and band response; bolometric luminosity cannot be substituted directly for a fitted SALT \(B\)-band amplitude. [Hogg's distance derivation](https://ned.ipac.caltech.edu/level5/Hogg/Hogg_contents.html) provides the distance conventions.

The physical distance modulus and the local shape-only quantity are

\[
\boxed{\mu=5\log_{10}(d_L/10\,{\rm pc})
=25+5\log_{10}(d_L/{\rm Mpc}),}
\qquad
\widetilde\mu=5\log_{10}\!\left[(1+z_{\rm HEL})X(z_{\rm HD})\right].
\tag{C12}
\]

`core.mu` returns \(\widetilde\mu\), with \(25+5\log_{10}[(c/H_0)/{\rm Mpc}]\) removed into the fitted intercept. This is intentional; treating it as an absolute physical modulus without restoring that constant would be wrong. The phase-two functions do include the constant, with \(c=299792.458\) km s\(^{-1}\), \(H_0\) in km s\(^{-1}\) Mpc\(^{-1}\), and distances in Mpc. For \(z=0\), \(d_L=0\) and a point-source distance modulus is undefined; distance-shape inference uses \(z>0\).

The exact measured redshift is the ratio of photon energies in emitter and observer frames:

\[
\boxed{1+z_{\rm obs}=\frac{(-u_\mu k^\mu)_e}{(-u_\mu k^\mu)_o}.}
\tag{C13}
\]

For separated cosmological and radial Doppler factors, \(1+z_{\rm obs}=(1+z_{\rm cos})(1+z_{\rm pec})\), where a receding purely radial speed \(\beta c\) gives \(1+z_{\rm pec}=\sqrt{(1+\beta)/(1-\beta)}\). For small speeds this yields \(v_{\rm pec}\simeq c(z_{\rm obs}-z_{\rm cos})/(1+z_{\rm cos})\). Simply subtracting \(cz-H_0D\) is a low-redshift approximation. [Davis and Scrimgeour](https://arxiv.org/abs/1405.0105) analyze that distinction.

The release convention used here is \(d_L^{\rm release}\propto(1+z_{\rm HEL})X(z_{\rm HD})\), with the release's velocity corrections and covariance. It agrees with the intended released-product calculation. It is not an exact solution for arbitrary moving sources/observers, which also requires aberration and flux transformation. The directional audit deliberately varies `zHEL`, `zCMB`, `zHD` and its reproduced Local-Group convention; these alternatives do not become independent distance measurements. In `load_c1`, a projected speed is inserted into the one-dimensional relativistic Doppler expression. That is a reproduction of the chosen convention and agrees to first order in speed, not a complete three-dimensional Lorentz boost at second order. No historical reproduction input is silently replaced by this audit.

Unmodeled gravitational lensing changes apparent flux by a magnification \(\mathcal A\), not extinction:

\[
\boxed{\mathcal A=\frac1{\left|(1-\kappa)^2-|\gamma|^2\right|},
\quad \Delta m=-2.5\log_{10}\mathcal A
\simeq-\frac5{\ln10}\kappa\quad(|\kappa|,|\gamma|\ll1).}
\tag{C14}
\]

Here \(\kappa\) is convergence and \(\gamma\) shear, using the image-plane Jacobian of the lens mapping in geometric optics. The absolute determinant gives positive image flux magnification; its original sign describes image parity. The linear magnitude expansion is restricted to the weak, positive-parity branch. Unresolved strong-lensing images require summing positive magnifications and accounting for time delays. See [Bartelmann and Schneider](https://arxiv.org/abs/astro-ph/9912508). This chapter does not independently implement the simulations' lensing distribution or its covariance. A redshift-dependent Gaussian magnitude scatter, where imported, is a statistical approximation to line-of-sight structure, not a physical lensing solution. Likewise a low-\(z\) velocity-scatter estimate \(\sigma_\mu\simeq(5/\ln10)\sigma_v/(cz)\) is a linear propagation approximation, not a complete correlated velocity-field likelihood. Both terms are already partially compressed into release covariances.

Curvature matters when differentiating distances. Define \(D=H_0D_M/c\). On the expanding branch before a closed-universe distance turning point,

\[
E=\frac{\sqrt{1+\Omega_KD^2}}{D'},\qquad
\boxed{q=-1+(1+z)\left[\frac{\Omega_KDD'}{1+\Omega_KD^2}-\frac{D''}{D'}\right].}
\tag{C15}
\]

This follows from differentiating C10, using \((d\mathcal S/dX)^2=1+\Omega_K\mathcal S^2\), then C08. Flatness removes the first term. Thus inferring acceleration from distances needs two derivatives and a geometry assumption. A fit over a finite bin is not an unconstrained point measurement of \(q(0)\).

## C16–C19. Kinematic models and directional cosmography

Integrating the first equality of C08 gives the exact kinematic relation

\[
\boxed{E(z)=\exp\!\left[\int_0^z\frac{1+q(u)}{1+u}\,du\right].}
\tag{C16}
\]

For `model='kinematic'`, \(q(z)=q_0+q_1z/(1+z)\), hence

\[
E(z)=(1+z)^{1+q_0+q_1}e^{-q_1z/(1+z)}.
\tag{C17}
\]

No matter/DE split or Friedmann acceleration equation is needed for this identity, but FLRW geometry remains. For piecewise constant \(q_i\) on \([z_i,z_{i+1}]\), continuity fixes \(E_i=E(z_i)\) recursively and gives

\[
\boxed{E(z)=E_i\left(\frac{1+z}{1+z_i}\right)^{1+q_i},
\quad E_{i+1}=E_i\left(\frac{1+z_{i+1}}{1+z_i}\right)^{1+q_i},}
\]
\[
\boxed{\int_{z_i}^{z}\frac{du}{E(u)}
=\frac{1+z_i}{E_i}
\begin{cases}
\dfrac{1-[(1+z)/(1+z_i)]^{-q_i}}{q_i},&q_i\ne0,\\
\ln[(1+z)/(1+z_i)],&q_i=0.
\end{cases}}
\tag{C18}
\]

`expm1` or an analytic series supplies the removable \(q_i=0\) singularity. The integrand is positive when \(H>0\); \(X\) must increase. Before correction the local implementations stopped adding path length after the last bin (2.5 for NumPy, 1.3 for JAX), although a finite \(q\) implied continuing expansion. This violates C18 outside the declared interval. Corrected behavior rejects the NumPy call and returns an invalid value under JAX tracing, with sampler input validation; it does not invent an extension of the physical model. The newly explicit NumPy piecewise \(E\) now agrees with C16 and C18.

Expanding \(a(t)\) about today and eliminating lookback time yields

\[
\boxed{\frac{H_0d_L}{c}=z+\frac{1-q_0}{2}z^2
-\frac{1-q_0-3q_0^2+j_0-\Omega_K}{6}z^3+O(z^4).}
\tag{C19}
\]

The code's cubic coefficients are correct for \(\Omega_K=0\). They are a truncated series, not an exact luminosity-distance law over arbitrary \(z\). The equality's domain is controlled by the discarded terms, not by positivity of the cubic alone. The existing directional tests preserve approximation errors and parameter biases up to \(z=0.8\). See [Visser's cosmographic derivation](https://arxiv.org/abs/gr-qc/0309109).

The directional ansatz replaces \(q_0\) in this polynomial by \(q_m+q_d(\hat n\cdot\hat n_d)e^{-z/S}\). Its dimensionless parameters \(q_m,q_d,S\) define a phenomenological angular/redshift fit. Substituting a varying \(q\) into a fixed-\(q_0,j_0\) Taylor series is not a derivation of a self-consistent anisotropic spacetime, and is not equivalent to integrating C16 along each direction. In particular, \(dq/dz\) contributes to \(j=q(2q+1)+(1+z)dq/dz\) for an isotropic FLRW history. An angular acceleration claim needs a specified metric, null propagation and velocity treatment; the reproduced fit alone supplies none of those missing foundations.

## C20–C22. BAO and the boundary of the CMB calculation

The distance combinations in the local BAO likelihood are

\[
\boxed{D_H=\frac c{H(z)},\qquad
D_V=[zD_M^2D_H]^{1/3},\qquad
\frac{D_M}{r_d}=\frac c{H_0r_d}X(z),\quad
\frac{D_H}{r_d}=\frac c{H_0r_dE(z)}.}
\tag{C20}
\]

\(D_V\) is the conventional isotropic volume-distance combination, not a third independent radial null path. All ratios are dimensionless; the two displayed expressions using \(X(z)\) assume the local flat geometry, while a curved model must use the transverse distance in C10. `hrd` means \(H_0r_d\), in km s\(^{-1}\), **not** the dimensionless reduced Hubble parameter \(h\) times \(r_d\). Keeping this product free avoids fixing the ruler's length from a CMB fit. It does not remove assumptions in the measured compressed BAO statistic.

Before decoupling, photon pressure supplies the restoring force and photons plus baryons supply inertia. In the tightly coupled, adiabatic limit \(\delta\rho_b/\rho_b=3\delta\rho_\gamma/(4\rho_\gamma)\), and \(p_\gamma=\epsilon_\gamma/3\), yielding

\[
\boxed{c_s^2=\frac{c^2}{3(1+R_b)},\quad
R_b=\frac{3\epsilon_b}{4\epsilon_\gamma},\quad
r_d=\int_0^{t_d}\frac{c_s(t)}{a(t)}dt
=\int_{z_d}^{\infty}\frac{c_s(z)}{H(z)}dz.}
\tag{C21}
\]

The drag epoch \(z_d\) depends on baryon-photon momentum coupling and the ionization history; it is not freely interchangeable with the last-scattering redshift \(z_*\). Radiation, baryons, neutrinos and recombination affect C21. The local `BAO.prediction` does not calculate it. [DESI DR2 equations 1–6](https://arxiv.org/html/2503.14738v3) specify the measured distance combinations and the ruler convention.

For context, the CMB acoustic angle and, for a flat scalar adiabatic calculation, linear-theory angular spectra have the schematic physical form

\[
\theta_* = \frac{r_s(z_*)}{D_M(z_*)},\qquad
C_\ell^{XY}=4\pi\int d\ln k\,\mathcal P_{\mathcal R}(k)
\Delta_\ell^X(k)\Delta_\ell^Y(k).
\tag{C22}
\]

Here temperature/polarization anisotropies are normalized to the mean CMB temperature, so the displayed \(C_\ell\) and transfer functions are dimensionless; spectra in kelvin squared require the corresponding temperature factor. The transfer functions \(\Delta_\ell\) require coupled Einstein–Boltzmann perturbation evolution, recombination, reionization and the specified neutrino model; lensing remaps the primary fields. Foreground and instrumental likelihoods introduce further ingredients. The local work reads posterior chains produced by that machinery and computes derived late-time quantities. It does not solve or verify all those equations, and a one-number prior on \(r_d\) or \(\theta_*\) is not a general substitute for their full likelihood. [CAMB's original implementation paper](https://arxiv.org/abs/astro-ph/9911177) is a primary reference for the numerical transfer calculation.

## P01–P03. Cosmic time and what a host age means

For an expanding Big-Bang model, elapsed cosmic time is

\[
\boxed{t(z)=\int_z^\infty\frac{du}{(1+u)H(u)}
=\frac1{H_0}\int_{-\infty}^{-\ln(1+z)}\frac{d\ln a}{E(a)}.}
\tag{P01}
\]

`clock` evaluates the second integral from \(\ln a=-16\), adding the matter-dominated boundary \(H_0t(a_{\min})\simeq2a_{\min}^{3/2}/(3\sqrt{\Omega_m})\). That boundary is consistent with the three declared radiation-free cosmologies and negligible at the adopted start. It would fail for arbitrary early-dark-energy domination. Physical conversion is \(H_0^{-1}=977.79222168/H_0\) Gyr when the numerical \(H_0\) is in km s\(^{-1}\) Mpc\(^{-1}\). Cosmic age, lookback time \(t_0-t(z)\), a star's age and progenitor delay are distinct.

For a galaxy forming mass at rate \(\Psi_G(t_f,Z)\) per metallicity interval, a star born at \(t_f\) has age \(\tau=t-t_f\). The three common age moments are

\[
\boxed{\langle\tau\rangle_W=
\frac{\int_0^t d\tau\int dZ\,\tau\,\Psi_G(t-\tau,Z)W(\tau,Z)}
{\int_0^t d\tau\int dZ\,\Psi_G(t-\tau,Z)W(\tau,Z)}.}
\tag{P02}
\]

Formed-mass weighting has \(W=1\); surviving-mass weighting uses the fraction \(f_{\rm surv}(\tau,Z)\) remaining in the specified stellar/remnant inventory; luminosity weighting uses band luminosity per formed mass \(\ell_b(\tau,Z)\). Those weights are not interchangeable. The observed host flux constrains these ages only through population-synthesis spectra, IMF, metallicity, attenuation and SFH assumptions. \(0\le\langle\tau\rangle_W\le t(z)\) for nonnegative histories and weights; any violation is a useful diagnostic but agreement alone is insufficient validation.

The mass-loss approximation actually used in `csfh_dtd.calculate` is

\[
\boxed{f_{\rm surv}(\tau)=1-0.046\ln\left(1+\frac{\tau}{0.000276\,{\rm Gyr}}\right).}
\tag{P03}
\]

The conversion \(0.276\) Myr \(=0.000276\) Gyr is correct. This is an empirical stellar-evolution/IMF approximation, not a universal conservation law, and should not be extrapolated to ages where it becomes negative. The calculation uses it for an illustrative cosmic formed-mass distribution, not a measured SN-host distribution. [Childress et al., appendix A](https://arxiv.org/html/1409.2951v1) specifies the Chabrier-IMF mass-loss prescription. The MD14 CSFH normalization uses a different reference IMF; a constant mass rescaling cancels in a normalized age moment, but an age/metallicity-dependent mismatch does not.

## P04–P07. Star formation and delay-time distributions

The two cosmic star-formation-rate densities implemented here are

\[
\boxed{\psi_{\rm B13}(z)=
\frac{0.180}{10^{-0.997(z-1.243)}+10^{0.241(z-1.243)}}}
\quad [M_\odot\,{\rm yr}^{-1}\,{\rm Mpc}^{-3}],
\tag{P04}
\]
\[
\boxed{\psi_{\rm MD14}(z)=
0.015\frac{(1+z)^{2.7}}{1+[(1+z)/2.9]^{5.6}}}
\quad [M_\odot\,{\rm yr}^{-1}\,{\rm Mpc}^{-3}].
\tag{P05}
\]

These are observationally calibrated fits, not consequences of Friedmann equations or nuclear burning. Their numerical coefficients, reference IMF and redshift extrapolation require empirical assessment. The stable `logaddexp` implementation in the code is algebraically equal to P04. References are [Behroozi et al., equation F1 and table 6](https://arxiv.org/abs/1207.6105) and [Madau and Dickinson, equation 15](https://arxiv.org/abs/1403.0007); the exact acquired text is retained under `runs/mapping/sources-2026-09-20/`.

A DTD \(\Phi(\tau,Z,B)\) is the expected explosion rate per **formed** stellar mass at delay \(\tau\), optionally conditional on metallicity and binary parameters. Its units are events mass\(^{-1}\) time\(^{-1}\). The smooth shape and truncated power-law shapes used by this workspace are

\[
\boxed{\Phi(\tau)=A\frac{(\tau/t_p)^\alpha}{1+(\tau/t_p)^{\alpha-s}},}
\qquad(t_p,s,\alpha)=(0.3\,{\rm Gyr},-1,20),
\tag{P06}
\]
\[
\boxed{\Phi(\tau)=A(\tau/\tau_*)^s\,\mathbf1_{\tau\ge\tau_{\min}},}
\quad(\tau_{\min},s)=(0.3\,{\rm Gyr},-1)
\ \text{or}\ (0.04\,{\rm Gyr},-1.13).
\tag{P07}
\]

Here \(A\) carries the rate units and \(\tau_*\) is a fixed reference time. The local code omits common factors such as \(A\) and \(\tau_*^{-s}\) because it normalizes the final delay distribution. This is correct for age moments; the resulting unnormalized array cannot be treated as an absolute SN event rate. The smooth function tends to \(\tau^\alpha\) at short delays and \(\tau^s\) at long delays. Its \(t_p\) is a transition scale, not an exact hard cutoff or exactly the maximum: \(\tau_{\rm peak}=t_p(-\alpha/s)^{1/(\alpha-s)}\) for \(\alpha>0,s<0\), approximately \(0.346\) Gyr for the adopted case. Small nonzero support below a physically possible white-dwarf formation time is a smooth approximation, not a claim of instantaneous explosions.

The shape in P06 is correctly transcribed from the Childress DTD family; the chosen sharpness 20 is a declared reconstruction choice. Neither P06, P07 nor their minimum delays follow from fundamental physics without a stellar/binary population model. A single \(\tau^{-1}\) law over \(0<\tau<\infty\) is not normalizable; finite cosmic time and a nonzero lower cutoff matter.

## P08–P11. How physical progenitor delays enter an event sample

One motivation for a late-time \(\tau^{-1}\) DTD comes from gravitational-radiation-driven circular binaries. With orbital separation \(b\), constant masses \(m_1,m_2\) and total \(M=m_1+m_2\), Newtonian orbital binding energy and leading quadrupole loss give

\[
E_{\rm orb}=-\frac{Gm_1m_2}{2b},\qquad
\dot E_{\rm orb}=-\frac{32G^4m_1^2m_2^2M}{5c^5b^5},
\quad \dot b=-\frac{64G^3m_1m_2M}{5c^5b^3},
\]
\[
\boxed{t_{\rm merge}=\frac{5c^5b_0^4}{256G^3m_1m_2M}.}
\tag{P08}
\]

Integration of \(\dot b\) supplies the factor 256. A population with \(dN/db_0\propto b_0^{-1}\) then gives \(dN/dt_{\rm merge}\propto t_{\rm merge}^{-1}\) by change of variables. This depends on an assumed separation distribution, circular point masses, negligible other torques, and merger-driven explosions; total progenitor delay also includes stellar evolution before a double-white-dwarf binary exists. It is therefore a physical motivation, not a derivation of a universal SN Ia DTD or the particular 40/300 Myr cutoffs. [Peters's original radiation-reaction calculation](https://journals.aps.org/pr/abstract/10.1103/PhysRev.136.B1224) supplies the orbital result. P08 is explanatory context, not an implemented binary-evolution solver in this repository.

For star formation \(\psi\), independent progenitor cohorts add linearly. A mass \(\psi(t-\tau)d\tau\) formed at the appropriate time contributes its DTD-weighted event rate now:

\[
\boxed{R_{\rm Ia}(t)=\int_0^t\psi(t-\tau)\Phi(\tau)d\tau,
\qquad p(\tau\mid t,\mathrm{SN})=
\frac{\psi(t-\tau)\Phi(\tau)}{R_{\rm Ia}(t)}.}
\tag{P09}
\]

With consistent time units the rate has units events time\(^{-1}\) volume\(^{-1}\), and the PDF has inverse-time units. The local script integrates over delay in Gyr while CSFH amplitudes are printed per year; its common factor \(10^9\) cancels in the normalized moments. It must be restored for absolute rates. The code evaluates formation redshift at cosmic time \(t(z)-\tau\); using \(t(z)+\tau\) would form progenitors after their explosion and would be a sign error. No extra \(dt/dz\) factor belongs inside this integral because its integration variable is already time. Converting the integration variable to formation redshift would require that Jacobian.

For an observed sample, the cosmic PDF is insufficient. Write \(n(G|z)\) for the galaxy abundance distribution and \(S\) for discovery, classification, light-curve and host/redshift selection. Then

\[
\boxed{p(\tau|z,S=1)\propto
\int dG\,dZ\,dB\,dI\,dD\,
n(G|z)\Psi_G[t(z)-\tau,Z]\Phi(\tau,Z,B)
p(I,D|\tau,Z,B,G)P(S=1|I,D,G,z).}
\tag{P10}
\]

Normalize over \(\tau\). \(I\) represents intrinsic emission and \(D\) dust/geometry. The cosmic-SFH convolution follows only when DTD and selection assumptions permit the integrals to collapse. Galaxy number weighting, formed-mass weighting and SN-event weighting give different distributions. For a list already sampled as SN events, multiplying each host again by its SN rate double-counts event weighting. The local convolution has no survey selector; its output is a cosmic-volume sensitivity, not a recovered selected-sample progenitor distribution.

If absolute event counts were modeled, observer-time rates would also contain

\[
\frac{dN}{dz\,d\Omega\,dt_o}=
\frac{R_{\rm Ia}[t(z)]}{1+z}\,
\frac{cD_M(z)^2}{H(z)}\,\varepsilon_{\rm sel}(z),
\tag{P11}
\]

for a comoving rate density and a consistently population-averaged efficiency. The \(1/(1+z)\) is time dilation, and \(cD_M^2/H\) is comoving volume per redshift/solid angle. A power-law event-rate input in a simulator is an empirical population model, not a substitute for deriving P09 or validating the selection function. Conditioning a distance analysis on observed redshift can omit the event-count likelihood, but it does not remove flux-dependent selection.

## P12–P15. Age-to-luminosity transport and an actual host-SFH mismatch

Suppose the **remaining standardized absolute-magnitude** response is \(M_{\rm rem}(\tau)=M_0+s_\tau\tau\), with \(s_\tau\) in mag Gyr\(^{-1}\). The predicted mean population shift is

\[
\boxed{e(z)=s_\tau[\langle\tau\rangle_{z,S}-\langle\tau\rangle_{0,S}].}
\tag{P12}
\]

For \(s_\tau=-0.030\) mag Gyr\(^{-1}\) and younger high-redshift progenitors, \(e(z)>0\): their standardized absolute magnitudes are fainter, and subtracting \(e(z)\) from apparent magnitude reduces the inferred distance modulus. `correction_subtracted_mag_*` stores \(+0.030[\tau(0)-\tau(z)]\); the sign is consistent with this convention. **The value and universality of \(s_\tau\) are empirical hypotheses, not consequences of energy conservation or white-dwarf burning.** The response refers to what remains after the release's existing corrections; applying a pre-correction population effect again can double count it.

Linearity selects the mean exactly, since \(\mathbb E[s_\tau\tau]=s_\tau\mathbb E[\tau]\). Replacing it with the median is a different assumption. For example, constant SFH and \(\Phi\propto\tau^{-1}\) on \([a,b]\) give

\[
\boxed{\langle\tau\rangle=\frac{b-a}{\ln(b/a)},\qquad
\tau_{1/2}=\sqrt{ab}.}
\tag{P13}
\]

At \(a=0.04,b=13\) Gyr these are 2.24073 and 0.72111 Gyr. The script intentionally saves both quantities as separate sensitivity templates. A median template is not the physical mean prediction of a linear magnitude law unless an additional justified relationship links their evolution. For a nonlinear emission law, integrate the whole response: \(e(z)=\int M_{\rm rem}(\tau)p_z(\tau)d\tau-\int M_{\rm rem}(\tau)p_0(\tau)d\tau\); neither one representative age nor one slope generally suffices.

An observed host-age slope transfers to delay only under additional conditions. If a stable conditional map has \(\mathbb E(\tau|H,z,S)=a+bH\) and residuals have zero conditional mean, then

\[
\boxed{b_{YH}=s_\tau\frac{\operatorname{Cov}(H,\tau)}{\operatorname{Var}(H)}=s_\tau b,
\qquad b_{YH}\Delta\langle H\rangle=s_\tau\Delta\langle\tau\rangle.}
\tag{P14}
\]

This is a conditional statistical identity, not a physical law of stellar populations. Evolving map intercept/slope, age error, confounding by dust/metallicity, or selection can invalidate its application while leaving every equation of stellar dynamics intact. Also, \(t\propto H_0^{-1}\) at fixed \(E\), but a fixed physical DTD cutoff changes its position relative to that clock. A magnitude correction expressed in Gyr is therefore not automatically independent of \(H_0\) merely because the SN intercept is marginalized.

The separately preserved MC-Age/FSPS audit identifies a real age-definition mismatch, not a new gravitational law. Let \(u\) be time since star formation began, \(U\) its current duration, \(T\) a transition time and \(\tau_s\) the early SFH scale. In a common normalization \(K=T e^{-T/\tau_s}\), the fitted historical FSPS history is

\[
\boxed{\Psi(u)\propto
\begin{cases}
u e^{-u/\tau_s},&0\le u\le T,\\
\max\{0,K[1+s(u-T)]\},&T<u\le U,
\end{cases}\quad
H_{\rm formed}=U-\frac{\int_0^Uu\Psi(u)du}{\int_0^U\Psi(u)du}.}
\tag{P15}
\]

\(s\) is a fractional slope per Gyr. The preserved MC-Age age postprocessor instead uses \(\max[0,K+s(u-T)]\): it omits multiplication of the late slope by \(K\), so it computes an age for a different SFH than the SED model. The independent semantics script correctly contrasts these two definitions, including the zero-SFR stopping point, using analytic moments and quadrature. This audit does not alter upstream historical chains or promote a postprocessing fix to a photometric refit. An updated interpretation of published ages still requires their actual chain provenance, fitted spectra and consistent mass-return conventions. See `docs/assumption-audit/pantheon.md` and [age-semantics.json](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/assumption_audit/pantheon/age-semantics.json) for the source-level evidence and scope.

## Verification record and unresolved assumptions

Reproduce the new targeted checks with

The original execution commands are preserved with the linked historical calculation records. Use the [workflow guide](../workflows.md) for current commands.

The script uses the separate pinned hierarchy Python environment for JAX checks. Results and exact source hashes are in [`cosmology-check.json`](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/physics_audit/cosmology-check.json). On the audited source:

| Independent check | Result |
|---|---:|
| CPL distance versus Astropy, three cosmologies, \(0.001\le z\le2.5\) | Maximum \(1.93\times10^{-15}\) mag difference |
| \(q\) versus finite-difference \(d\ln H/dz\) | Maximum \(1.84\times10^{-10}\) |
| CPL energy-conservation identity | Maximum \(4.24\times10^{-10}\) |
| Piecewise \(E\) versus directly accumulated segment expression | Maximum \(1.34\times10^{-15}\) |
| Constant-\(q\) analytic comoving integral | Maximum \(5.64\times10^{-13}\) dimensionless |
| JAX constant-\(q\) distance | Maximum \(9.03\times10^{-13}\) mag |
| JAX gradient exactly at \(q=0\) | Finite; invalid redshifts rejected as NaN under `jit` |
| BAO ratios versus independent Astropy distances | Maximum \(2.14\times10^{-14}\) |
| Cosmic ages versus Astropy over the three declared clocks | Maximum \(9.96\times10^{-8}\) Gyr |
| B13/MD14 formula transcription | Maximum \(1.39\times10^{-17}\) in printed SFR units |
| C14×B13 delay integrals versus adaptive quadrature at \(z=0,1\), shared declared clock | Mean difference \(<10^{-10}\) Gyr; median difference \(<7.22\times10^{-7}\) Gyr |

All actual checked Pantheon points lie below \(z_{\rm HD}=2.26137<2.5\); the three prepared hierarchy datasets lie below 1.12132<1.3, and current synthetic/coverage files below 0.95. Thus the domain repair does not invalidate prior fits through a changed in-domain distance law. Previously archived generated results retain their original source provenance; rerunning inference is not needed merely to relabel a new boundary guard.

The JSON also records a transparent radiation-omission example with \(H_0=70\), \(\Omega_m=0.3\), \(T_{\rm CMB}=2.7255\) K and massless neutrinos, compared with the radiation-free clock. Including radiation changes the modulus by \(-0.000420\) mag at \(z=2.5\), and the cosmic age by approximately \(-0.00529\) Gyr at \(z=0.001\). These approximation differences are separate from the much smaller integration errors above. This example is not a fitted radiation correction or a license to omit radiation in C21.

The major unresolved physical questions remain the actual standardized luminosity response to progenitor properties, the transport from selected host ages to progenitor delays, dust/metallicity/SFH degeneracies in host spectra, representative DTD/CSFH weights, flux-dependent selection, and a self-consistent recomputation of any age correction when cosmology changes. A fixed-template correction freezes those choices. Numerically accurate FLRW integration cannot establish them, and no universal sign of present acceleration follows without specifying them.
