# Physical luminosity, energy release, and the empirical SALT interface

Implementation references in this chapter identify the preserved source snapshot; see [physical foundations](README.md) for the current execution boundary.

Audit date: 2026-09-21. This chapter establishes the emission equations needed by this repository's Type Ia supernova analysis. Equation labels `L01`–`L35` are stable references. It distinguishes conservation laws, constitutive approximations, measured nuclear data, and empirical light-curve relations. The derivations below are independent algebra; linked primary papers supply historical formulations, laboratory input, or physical counterexamples.

**Finding:** the active project predicts calibrated flux from empirical SALT spectral surfaces and models the distribution of SALT parameters. It does not calculate an explosion, radioactive energy deposition, or radiation transport from ejecta properties. Its existing SALT documentation largely acknowledges this correctly. The absence of a physical emission calculation is an untested physical bridge, not evidence that every fitted light curve is wrong. This chapter supplies that bridge as equations with explicit assumptions; it does not silently substitute a simplified explosion model into cosmology fits.

## 1. Quantities that must not be interchanged

Use source rest-frame time \(t\) measured from explosion for heating and expansion. SALT phase \(p\) is measured from fitted maximum light, so \(t=p+t_{\rm rise}\), with a generally uncertain rise time. The SALT parameter called `t0` is an observer-date maximum, not the gamma escape time used below. Use \(t_\gamma\) for that escape scale and \(\tau_m\) for the diffusion scale.

| Symbol | Meaning and units (cgs unless stated otherwise) |
|---|---|
| \(M_{\rm ej},M_{\rm Ni,0}\) | Total ejected mass and initial radioactive nickel mass, g |
| \(N_j,\lambda_j=1/\tau_j\) | Number of nuclei and decay rate, s\(^{-1}\); \(\tau_j\) is a mean life |
| \(Q_j\) | Energy per decay, erg; channel and neutrino convention must be stated |
| \(\dot Q_{\rm rad},\dot Q_{\rm dep}\) | Generated non-neutrino radioactive power and deposited thermalizing power, erg s\(^{-1}\) |
| \(E\), \(e\), \(u\) | Stored internal energy (erg), specific internal energy (erg g\(^{-1}\)), energy density (erg cm\(^{-3}\)) |
| \(L_{\rm bol},L_\lambda,L_\nu\) | Emitted luminosity: erg s\(^{-1}\), per cm or per Å, and per Hz; choose one wavelength unit consistently |
| \(F\), \(I_\nu\) | Local radiative energy flux, erg s\(^{-1}\) cm\(^{-2}\), and specific intensity, erg s\(^{-1}\) cm\(^{-2}\) Hz\(^{-1}\) sr\(^{-1}\) |
| \(\kappa\), \(\chi=\rho\kappa\) | Mass opacity, cm\(^2\) g\(^{-1}\), and extinction coefficient, cm\(^{-1}\) |
| \(x_0,x_1,c_{\rm SALT}\) | Empirical amplitude, shape, color coordinates; none is automatically \(M_{\rm Ni}\), diffusion time, or physical reddening |

Here \(c_{\rm light}\) denotes the speed of light whenever SALT color might cause confusion. Energy and photon flux differ by photon energy, \(h\nu=hc_{\rm light}/\lambda\).

The bolometric quantity in the diffusion equation is the reprocessed thermal radiation emerging from the ejecta, ideally integrated over its complete spectrum. Unthermalized escaping decay gamma rays belong to a separate output channel. A UV/optical/IR integral with missing far-UV or IR is a *quasi-bolometric* measurement until the missing flux is accounted for.

\[
L_{\rm bol}(t)=\int_0^\infty L_\nu(t)\,d\nu
             =\int_0^\infty L_\lambda(t)\,d\lambda,
\qquad
L_\nu=L_\lambda\left|\frac{d\lambda}{d\nu}\right|
       =L_\lambda\frac{\lambda^2}{c_{\rm light}}.
\tag{L01}
\]

This is a spectral definition, not a blackbody assumption. The last equality uses \(L_\lambda\) per unit length in the same units as \(\lambda\) and \(c_{\rm light}\). Per-Å versus per-cm conversion introduces \(10^8\).

## 2. What would be required to predict the explosion

The physical system before free expansion obeys mass, momentum, and energy conservation. A nonrelativistic, inviscid, unmagnetized, diffusion-limit formulation is

\[
\frac{D\rho}{Dt}=-\rho\nabla\!\cdot\mathbf v,\qquad
\rho\frac{D\mathbf v}{Dt}=-\nabla(P_{\rm gas}+P_{\rm rad})-\rho\nabla\Phi,
\qquad \nabla^2\Phi=4\pi G\rho.
\tag{L02}
\]

The radiation-pressure form assumes a nearly isotropic comoving radiation field. Outside that limit, its force requires the radiation momentum equation or a transport calculation, and cannot be replaced by an arbitrary pressure. Nuclear abundances and specific thermal energy obey

\[
\frac{DX_i}{Dt}=\dot X_i(\rho,T,\{X_j\}),\qquad
\frac{De}{Dt}+P\frac{D(1/\rho)}{Dt}
 =\dot q_{\rm heat}-\frac{1}{\rho}\nabla\!\cdot\mathbf F_{\rm rad}.
\tag{L03}
\]

For thermal energy excluding nuclear rest masses, a reaction network's nuclear heating follows mass-energy conservation:

\[
Q_r=\left(\sum_{\rm reactants}m_i-\sum_{\rm products}m_i\right)c_{\rm light}^2,
\qquad
\dot q_{\rm nuc}=
\frac{1}{\rho}\sum_r Q_r\,\mathcal R_r-\dot q_{\nu,\rm escape}.
\tag{L04}
\]

\(\mathcal R_r\) counts reaction events per volume per time. Atomic/nuclear masses, electrons, positrons, and neutrinos must be treated with the same convention on both sides; otherwise annihilation or rest-mass energy can be counted twice. Reaction rates, electron degeneracy, ionization, gravity, and flame/detonation evolution close these equations. Conservation alone does not specify those functions, ignition conditions, or the nickel yield. The light-curve work here starts after these unknowns have already determined the ejecta.

Once pressure/gravity acceleration becomes small, ballistic homologous expansion gives

\[
\mathbf r=\mathbf v t,\quad \nabla\!\cdot\mathbf v=3/t,
\quad \rho(v,t)=t^{-3}f(v),
\quad
M_{\rm ej}=4\pi\!\int f(v)v^2dv,
\quad E_{\rm kin}=2\pi\!\int f(v)v^4dv.
\tag{L05}
\]

These last integrals assume spherical symmetry. They follow by substituting \(r=vt\), \(dr=t\,dv\) into mass and kinetic-energy integrals. Homology does not require uniform density. For a uniform-density sphere with outer speed \(v_{\max}\),

\[
R=v_{\max}t,\qquad \rho=\frac{3M_{\rm ej}}{4\pi R^3},
\qquad E_{\rm kin}=\frac{3}{10}M_{\rm ej}v_{\max}^2.
\tag{L06}
\]

Consequently the coefficient is not \(1/2\) if the quoted speed is the outer homologous speed; \(E_{\rm kin}=M v_{\rm rms}^2/2\) instead uses a mass-weighted rms speed. An observed absorption-line velocity need not equal either one. No universal nickel mass, ejecta mass, or progenitor-age luminosity slope follows from L02–L06.

## 3. Radioactive energy: solve the chain before assigning luminosity

For initially synthesized \(^{56}\mathrm{Ni}\), neglecting initial \(^{56}\mathrm{Co}\), the chain is

\[
\frac{dN_{\rm Ni}}{dt}=-\lambda_{\rm Ni}N_{\rm Ni},\qquad
\frac{dN_{\rm Co}}{dt}=\lambda_{\rm Ni}N_{\rm Ni}-\lambda_{\rm Co}N_{\rm Co},\qquad
\frac{dN_{\rm Fe}}{dt}=\lambda_{\rm Co}N_{\rm Co}.
\tag{L07}
\]

Thus, with \(N_0\simeq M_{\rm Ni,0}/(56m_u)\),

\[
\begin{aligned}
N_{\rm Ni}(t)&=N_0e^{-t/\tau_{\rm Ni}},\\
N_{\rm Co}(t)&=N_0\frac{\tau_{\rm Co}}{\tau_{\rm Co}-\tau_{\rm Ni}}
   \left(e^{-t/\tau_{\rm Co}}-e^{-t/\tau_{\rm Ni}}\right),\\
N_{\rm Fe}(t)&=N_0-N_{\rm Ni}(t)-N_{\rm Co}(t).
\end{aligned}
\tag{L08}
\]

Integrating factors give the second line: multiply the Co equation by \(e^{t/\tau_{\rm Co}}\), integrate the Ni source, and impose \(N_{\rm Co}(0)=0\). An initially present Co population adds \(N_{\rm Co,0}e^{-t/\tau_{\rm Co}}\). Mean life and half-life are related by \(\tau=t_{1/2}/\ln2\), and must not be interchanged.

The generated *non-neutrino* power is the energy per decay times the number of decays per time:

\[
\begin{aligned}
\dot Q_{\rm rad}(t)
 &=Q_{\rm Ni}\lambda_{\rm Ni}N_{\rm Ni}
  +Q_{\rm Co}\lambda_{\rm Co}N_{\rm Co}\\
 &=\frac{M_{\rm Ni,0}}{56m_u}\left[
 \frac{Q_{\rm Ni}}{\tau_{\rm Ni}}e^{-t/\tau_{\rm Ni}}
 +\frac{Q_{\rm Co}}{\tau_{\rm Co}-\tau_{\rm Ni}}
  \left(e^{-t/\tau_{\rm Co}}-e^{-t/\tau_{\rm Ni}}\right)\right].
\end{aligned}
\tag{L09}
\]

**A consequential correction to a common shortcut:** \(Q_{\rm Co}(-dN_{\rm Co}/dt)\) is not the Co heating rate. Co is being produced as well as destroyed, so that shortcut even gives negative Co heating early. The correct rate is \(Q_{\rm Co}\lambda_{\rm Co}N_{\rm Co}\). Similarly, writing \(\epsilon_{\rm Co}=Q_{\rm Co}/(56m_u\tau_{\rm Co})\) and multiplying only by the exponential difference omits the buildup factor \(\tau_{\rm Co}/(\tau_{\rm Co}-\tau_{\rm Ni})\). With the lifetimes below, that error underestimates the Co contribution by 7.91%. These are cautions for a future physical implementation; the current project does not contain such a decay implementation to repair.

For a specified conventional laboratory compilation, \(\tau_{\rm Ni}=8.80\) d, \(\tau_{\rm Co}=111.3\) d, \(Q_{\rm Ni,\gamma}=1.75\) MeV, \(Q_{\rm Co,\gamma}=3.61\) MeV, and the branch-averaged Co positron kinetic energy is \(0.12\) MeV. The Co gamma allocation includes annihilation photons, so another \(2m_ec^2\) per positron must not be added. Neutrino energy is excluded. Rounded historical coefficients give

\[
\dot Q_{\rm rad}\simeq\frac{M_{\rm Ni,0}}{M_\odot}
 \left(6.45\times10^{43}e^{-t/(8.80\,{\rm d})}
      +1.45\times10^{43}e^{-t/(111.3\,{\rm d})}\right)
 {\rm erg\,s^{-1}}.
\tag{L10}
\]

These are measured nuclear inputs rather than exact universal decimal constants; the source is [Nadyozhin (1994), equations 6–19 and decay-energy tables](https://articles.adsabs.harvard.edu/pdf/1994ApJS...92..527N). Laboratory electron-capture rates are an approximation for the relevant charge states; unusually ionized conditions require their own check.

An independent integral check follows because each initial nucleus eventually undergoes each transition once:

\[
\int_0^\infty\dot Q_{\rm rad}(t)dt
 =\frac{M_{\rm Ni,0}}{56m_u}(Q_{\rm Ni}+Q_{\rm Co}).
\tag{L11}
\]

This is an energy budget, not the integrated optical radiated energy: escaping gamma rays and work done on expanding ejecta reduce the latter.

## 4. Deposition: radioactive power is not already optical light

In general, deposition must be found by transporting gamma rays, positrons, and their secondary electrons through the ejecta. A useful bookkeeping equation is

\[
\dot Q_{\rm dep}=\sum_{j=\rm Ni,Co}\lambda_jN_j
 \left[Q_{j,\gamma}f_{j,\gamma}(t)+Q_{j,+}f_{j,+}(t)+\cdots\right].
\tag{L12}
\]

Every \(f\) is a deposited *energy* fraction between zero and one, averaged over the corresponding source distribution and spectrum. It is not just an interaction probability if an interaction transfers only part of the photon energy. Prompt thermalization is another assumption; if fast particles or ionization energy store energy for long periods, their stored-energy equations must be added. Positron trapping is sensitive to propagation and magnetic fields, so setting \(f_+=1\) is a declared approximation.

For a straight ray in grey pure absorption,

\[
\tau_\gamma=\int\kappa_\gamma\rho\,ds,
\qquad f_{\gamma,\rm ray}=1-e^{-\tau_\gamma},
\qquad f_\gamma=1-\langle e^{-\tau_\gamma}\rangle_{\rm decay,angle}.
\tag{L13}
\]

The average of the exponential is not the exponential of the average. A scalar grey opacity here represents an effective energy-loss prescription for a problem that includes Compton scattering and energy redistribution. It is not the same opacity used for optical diffusion.

Under L05 with fixed source distribution in velocity coordinates and fixed gamma opacity, \(\rho\propto t^{-3}\) while path length scales as \(t\); hence every homologous column scales as \(t^{-2}\). This motivates a *one-scale approximation*

\[
\tau_{\gamma,\rm eff}=\left(\frac{t_\gamma}{t}\right)^2,
\qquad f_\gamma\simeq1-e^{-(t_\gamma/t)^2},
\qquad
f_\gamma\to\begin{cases}1&t\ll t_\gamma,\\(t_\gamma/t)^2&t\gg t_\gamma.\end{cases}
\tag{L14}
\]

The power \(t^{-2}\) follows from homology; a single exponential interpolation and its normalization do not. For a central source in a uniform sphere, \(\tau_\gamma=3\kappa_\gamma M/(4\pi v_{\max}^2t^2)\); a distributed source or different density profile changes the weighting. Ejecta columns can therefore affect light-curve shape even at fixed nickel mass. A primary treatment of the escape-timescale dependence on ejecta geometry is [Levanon & Soker (2019)](https://academic.oup.com/mnras/article/486/4/5528/5487090).

At late times, if photon and thermal equilibration times are short and stored-energy changes are negligible, \(L_{\rm bol}\simeq\dot Q_{\rm dep}\). Gamma-dominated Co emission in the optically thin limit then behaves approximately as \(e^{-t/\tau_{\rm Co}}t^{-2}\), rather than a pure Co exponential. Positron deposition, other isotopes, or delayed recombination can change that limit.

## 5. Stored energy, expansion work, and diffusion

Integrating L03 over a comoving volume gives

\[
\frac{dE}{dt}=\dot Q_{\rm dep}-L_{\rm bol}
                -\int P\nabla\!\cdot\mathbf v\,dV.
\tag{L15}
\]

This expression assumes no mass exchange across the chosen comoving boundary and includes all material+radiation internal reservoirs in \(E\) and \(P\). External interaction or a central engine contributes additional heating; they cannot be erased by calling the input radioactive power.

For radiation-dominated internal energy, \(u_{\rm rad}=a_{\rm R}T^4\), \(P_{\rm rad}=u_{\rm rad}/3\); with homologous expansion, the integral is \(E/t\), yielding

\[
\frac{dE}{dt}=\dot Q_{\rm dep}-L_{\rm bol}-\frac{E}{t}.
\tag{L16}
\]

Thus \(E\propto t^{-1}\) for an adiabatically expanding radiation field with heating and emission switched off. Gas-dominated monatomic internal energy instead loses \(2E/t\); ionization energy has another effective relation. Confusing these regimes changes the light curve while still leaving superficially plausible dimensions.

Multiplication by \(t\) supplies a useful model-independent-within-L16 check:

\[
tE(t)-t_iE(t_i)
 =\int_{t_i}^{t}t'\dot Q_{\rm dep}(t')dt'
 -\int_{t_i}^{t}t'L_{\rm bol}(t')dt'.
\tag{L17}
\]

The explicit finite start \(t_i\) prevents a singular early homologous approximation from silently removing shock-stored energy. If both boundary terms are negligible, the two time-weighted integrals become equal. This conservation relation does not require grey optical opacity or a prescribed diffusion time. Its observational use does require bolometric coverage and a deposition model, as established by [Katz, Kushnir & Dong (2013)](https://arxiv.org/abs/1301.6766).

In an optically thick medium with a nearly isotropic radiation field, diffusion is the first angular-moment closure:

\[
\mathbf F_{\rm rad}=-\frac{c_{\rm light}}{3\kappa_R\rho}\nabla u_{\rm rad},
\qquad
\frac{1}{\kappa_R}
 =\frac{\int\kappa_{\nu,\rm tr}^{-1}(\partial B_\nu/\partial T)d\nu}
        {\int(\partial B_\nu/\partial T)d\nu}.
\tag{L18}
\]

The Rosseland expression requires an approximately thermal local spectrum and an appropriate transport opacity. In rapidly expanding, line-rich, non-LTE ejecta it is a closure to be checked, not an exact material constant. Diffusion fails near optical depth of order unity and in optically thin regions; unconstrained diffusion can even predict a flux exceeding \(c_{\rm light}u\). Detailed transfer or a justified flux-limited closure is then required. The role of line opacity and time dependence is developed in [Pinto & Eastman, *Opacity and Diffusion*](https://arxiv.org/abs/astro-ph/9611195).

A random-walk estimate gives \(t_{\rm diff}\sim\tau_{\rm opt}R/c_{\rm light}\) up to a structural factor: the mean free path is \(\ell=(\kappa\rho)^{-1}\), about \((R/\ell)^2\) steps are required, and each step takes \(\ell/c_{\rm light}\). At fixed grey opacity and homologous geometry this decreases as \(t^{-1}\). Define the one-zone convention

\[
t_{\rm diff}(t)=\frac{\tau_m^2}{2t},\qquad
L_{\rm bol}=\frac{E}{t_{\rm diff}(t)}=\frac{2tE}{\tau_m^2},
\qquad
\tau_m^2=\frac{2\kappa M_{\rm ej}}{\beta_d c_{\rm light}v_{\rm sc}}.
\tag{L19}
\]

\(\beta_d\) represents the adopted spatial profile/eigenfunction and boundary closure; a commonly used uniform-density grey prescription gives about 13.8. It is unrelated to the Tripp color coefficient \(\beta\). Values of \(\tau_m\) in another convention may differ by \(\sqrt2\), so the stated \(L(E,t)\) closure must accompany the timescale. The random walk establishes scaling, not the numerical prefactor. With \(v_{\rm sc}\propto(E_{\rm kin}/M_{\rm ej})^{1/2}\),

\[
\tau_m\propto\kappa^{1/2}M_{\rm ej}^{3/4}E_{\rm kin}^{-1/4}.
\tag{L20}
\]

A width alone therefore cannot distinguish opacity, mass, and kinetic energy, even within this simplified model.

## 6. Derive Arnett's relation and its limits explicitly

Substitute \(E=\tau_m^2L/(2t)\) into L16 with constant \(\tau_m\). The terms proportional to \(L/t\) cancel, leaving

\[
\frac{dL}{dt}=\frac{2t}{\tau_m^2}(\dot Q_{\rm dep}-L).
\tag{L21}
\]

Multiplying by \(e^{t^2/\tau_m^2}\) and integrating gives

\[
L(t)=L(t_i)e^{-(t^2-t_i^2)/\tau_m^2}
  +\frac{2}{\tau_m^2}\int_{t_i}^{t}
     t'e^{-(t^2-t'^2)/\tau_m^2}\dot Q_{\rm dep}(t')dt'.
\tag{L22}
\]

Taking \(t_i\to0\), with negligible initial stored-energy contribution, gives the familiar radioactive diffusion convolution. It is a history integral: setting \(L(t)=\dot Q_{\rm dep}(t)\) at every time discards storage and the rise entirely. At a differentiable maximum with \(t_{\rm peak}>0\), L21 implies

\[
L(t_{\rm peak})=\dot Q_{\rm dep}(t_{\rm peak})
\quad\text{within the constant-}\tau_m\text{ closure.}
\tag{L23}
\]

The equality is with *deposited* power, not necessarily L09's generated power. It is bolometric, not \(B\)-band peak luminosity, and its time is not automatically the date of \(B\) maximum. The original analytic construction is [Arnett (1982)](https://articles.adsabs.harvard.edu/pdf/1982ApJ...253..785A).

The limit of this derivation can itself be exposed algebraically. If the closure becomes \(t_{\rm diff}=C(t)/t\), then \(E=C(t)L/t\) and

\[
\frac{dL}{dt}=\frac{t}{C(t)}(\dot Q_{\rm dep}-L)
                       -\frac{\dot C(t)}{C(t)}L,
\qquad
L_{\rm peak}=\frac{\dot Q_{\rm dep}(t_{\rm peak})}
                  {1+\dot C(t_{\rm peak})/t_{\rm peak}}.
\tag{L24}
\]

Here \(C\) has units time squared; its derivative divided by \(t\) is dimensionless. This is still a one-zone result, not a general correction formula. If the spatial distribution of heating differs from the assumed radiation profile, there may be no adequate single \(C(t)\). A multizone transfer calculation need not obey L23 even with constant material opacity. Such departures are directly investigated by [Khatami & Kasen (2019)](https://arxiv.org/abs/1812.06522).

The assumptions required before inferring nickel mass from peak light are therefore: a deposited-heating prescription; approximately homologous, radiation-dominated ejecta; the adopted optical diffusion/opacity/spatial closure; stated initial energy; sufficient bolometric coverage; and a peak time relative to explosion. Agreement with an analytic formula does not test all these assumptions independently.

## 7. Why a bolometric model does not determine color or a SALT light curve

The emergent spectrum requires frequency-dependent transport, including Doppler-shifted lines, scattering, fluorescence, and level populations. A flat-space inertial-frame form for radiation propagation is

\[
\frac{1}{c_{\rm light}}\frac{\partial I_\nu}{\partial t}
 +\mathbf n\cdot\nabla I_\nu
 =\eta_\nu-\chi_\nu I_\nu+\mathcal S_\nu[I].
\tag{L25}
\]

\(\eta_\nu\) is true emissivity per steradian; \(\chi_\nu\) includes removal from the ray; \(\mathcal S\) restores scattering into the ray, including redistribution in frequency and direction. Coefficients must be transformed from the moving material frame; a comoving-frame equation has additional Doppler/velocity-gradient terms. Applying static, comoving opacities to L25 without those transformations is not a complete expanding-ejecta calculation. The output luminosity is a surface integral,

\[
L_\nu(t)=\int_{\partial V}\mathbf F_\nu\cdot d\mathbf A,
\qquad \mathbf F_\nu=\int I_\nu\mathbf n\,d\Omega.
\tag{L26}
\]

For physical interpretation of opacities/emissivities, species level populations satisfy rate balance, schematically

\[
\frac{Dn_i}{Dt}+n_i\nabla\!\cdot\mathbf v
 =\sum_{j\ne i}(n_jP_{ji}-n_iP_{ij})
  +\text{ionization/recombination and nuclear terms}.
\tag{L27}
\]

\(P_{ij}\) are transition probabilities per time with collisional and radiative contributions. They depend on electron density, temperature, and the radiation field; abundance conservation and charge neutrality provide additional constraints. Setting the left side to zero is statistical equilibrium; assigning Saha–Boltzmann populations additionally invokes local thermodynamic equilibrium. Neither is automatically valid throughout Type Ia ejecta. Tabulated atomic data and a closure for all these rates are needed before claiming a first-principles spectrum.

For a true opaque thermal surface in LTE,

\[
B_\nu(T)=\frac{2h\nu^3}{c_{\rm light}^2}
       \frac{1}{e^{h\nu/(k_BT)}-1},\qquad
F_\nu=\pi B_\nu,\qquad
L_\nu=4\pi^2R^2B_\nu.
\tag{L28}
\]

The \(\pi\) comes from integrating the projected isotropic intensity over an outward hemisphere, not \(4\pi\). Integration over frequency yields

\[
L_{\rm bol}=4\pi R^2\sigma_{\rm SB}T^4,
\qquad \sigma_{\rm SB}=\frac{2\pi^5k_B^4}{15h^3c_{\rm light}^2},
\qquad a_{\rm R}=\frac{4\sigma_{\rm SB}}{c_{\rm light}}.
\tag{L29}
\]

For a non-blackbody spectrum, \(T_{\rm eff}=[L/(4\pi R_{\rm ph}^2\sigma_{\rm SB})]^{1/4}\) can still be a *definition*, but it does not establish \(L_\nu=4\pi^2R_{\rm ph}^2B_\nu(T_{\rm eff})\). Color temperature, gas temperature, and effective temperature need not coincide. A simple grey photospheric location \(\int_{R_{\rm ph}}^\infty\kappa\rho dr\sim2/3\) is itself a photospheric approximation; line-rich ejecta do not have one sharply defined radius at every wavelength.

An apparent early \(t^2\) rise follows only when an expanding photospheric area \(R_{\rm ph}^2\propto t^2\) is combined with nearly constant relevant surface brightness. Even in a blackbody approximation, \(L_{\rm bol}\propto t^2T^4\), and a fixed optical band need not share that temperature dependence. Thus a \(t^2\) rise is not a universal explosion clock.

The standard empirical decline coordinate is

\[
\Delta m_{15}(B)=m_B(t_{B,\max}+15\,{\rm rest\ d})-m_B(t_{B,\max}).
\tag{L30}
\]

It is a photometric definition. More nickel can alter temperatures, ionization, line blanketing, redistribution into redder wavelengths, and opacity; it is insufficient to infer the \(B\)-band width solely from L20. [Kasen & Woosley (2007)](https://arxiv.org/abs/astro-ph/0609540) demonstrate in their modeled ejecta that ionization/color evolution can dominate the \(B\)-band width–luminosity relation. This supplies a physical mechanism, not a universal fitted slope or proof of no redshift evolution.

## 8. What SALT and standardization actually implement here

The active spectral surrogate is

\[
S_{\rm SALT}(p,\lambda)=x_0[M_0(p,\lambda)+x_1M_1(p,\lambda)]
                  10^{-0.4c_{\rm SALT}CL_m(\lambda)},
\qquad p=\frac{t_{\rm obs}-t_{0,\rm obs}}{1+z}.
\tag{L31}
\]

The tables and color law are learned from observations; they are not solutions of L02–L29. In the notation \(\exp[c_{\rm SALT}CL_e]\), \(CL_e=-0.4\ln10\,CL_m\). In the pinned SNANA convention \(C=-CL_m\), so the factor is \(\exp[(\ln10/2.5)c_{\rm SALT}C]\). The reference anchors are \(C(4302.57\,\mathrm{\mathring A})=0\) and \(C(5428.55\,\mathrm{\mathring A})=1\). Positive color at fixed \(x_0\) leaves the first anchor unchanged and increases red flux. A physical screen instead removes energy along the beam; fitting \(x_0\) absorbs its grey/blue dimming. The empirical model and its training assumptions originate in [Guy et al. (2007)](https://arxiv.org/abs/astro-ph/0701828) and [Kenworthy et al. (2021)](https://arxiv.org/abs/2104.07795).

The executed implementation is [`Engine.flux`](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/phase2/independent_flux/engine.py), using the original DES-SN5YR SALT3 release. [`Engine.prepare`](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/phase2/independent_flux/engine.py) applies observer-frame Milky Way attenuation, wavelength/time redshift transformations, and photon-weighted passband/calibration integration. These operations are dealt with in the extinction/propagation chapters. This emission audit confirms that they integrate SALT \(M_0+x_1M_1\), not a computed radioactive SED.

The amplitude-to-magnitude convention recorded in the local portable fixture is

\[
m_B=-2.5\log_{10}x_0+10.635,
\qquad
\delta m_B=-\frac{2.5}{\ln10}\,\delta\ln x_0.
\tag{L32}
\]

See [`export_portable_pilot.py`](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/phase2/official/export_portable_pilot.py) and the bundled [`SALT2mBcalc`](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/official/portable_pilot/assets/genmag_SALT2.c). The additive constant and `MAG_OFFSET=0.27` are model/calibration conventions; neither converts \(x_0\) directly into a nickel mass. The latter offset is already included in the flux calibration and must not be added twice.

Standardization uses a fitted relation of the form

\[
\widehat\mu=m_B-M+\alpha x_1-\beta c_{\rm SALT}
                 +\Delta_{\rm host}+\Delta_{\rm other}.
\tag{L33}
\]

Correction signs must be defined through their action; a code defining a positive magnitude *bias* often subtracts it instead of using the plus sign above. This relation is empirical. Its linearity, coefficient constancy, environmental terms, residual scatter, and transport to a different selected population are assumptions to test, not consequences of radioactive decay. \(x_1\) is not a measured diffusion time, \(\beta\) is not a nuclear or extinction constant, and \(M\) is not the unprocessed bolometric nickel luminosity.

The hierarchical implementation makes the assumption particularly explicit:

\[
\begin{aligned}
m_B&=\mu+M-\alpha x_1+\beta_{\rm int}c_{\rm int}
                  +R_{B,\rm eff}E+\gamma H+\Delta M(z)+\epsilon_M,\\
c_{\rm SALT}&=c_{\rm int}+E.
\end{aligned}
\tag{L34}
\]

See [`effective_population` and `population_logpdf`](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/phase2/hierarchy/core.py). The explicit luminosity term \(+\gamma H\) in L34 would require \(\Delta_{\rm host}=-\gamma H\) in a distance estimator using the convention of L33. This is an effective SALT-coordinate decomposition. Its \(E\) and \(R_{B,\rm eff}\) do not become physical sightline reddening and an extinction ratio merely because the symbols resemble them. A defensible physical mapping is: compute an intrinsic time series with varied ejecta/abundances; apply a wavelength-dependent dust model; propagate and sample it with the survey's calibration/cadence/noise/selection; refit with the actual SALT model; then measure the mapping to \(x_0,x_1,c\). The [existing SALT assumption report](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/docs/salt-dust-audit/salt-model-assumptions.md) and [hierarchical specification](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/docs/phase2/generative-model.md) already retain this boundary correctly.

For a fixed passband, phase, and magnitude convention, define \(\mathcal L_B\) using the *same detector response as the magnitude*. For a photon-counting system, \(\mathcal L_B=\int [\lambda/(hc_{\rm light})]T_B(\lambda)L_\lambda\,d\lambda\) is the response-weighted photon luminosity (photons s\(^{-1}\)); a common multiplicative normalization cancels below. It is not the energy-weighted integral \(\int T_B L_\lambda\,d\lambda\), whose ratio can differ when the SED changes. Then

\[
\Delta M_B=-2.5\log_{10}\frac{\mathcal L_{B,2}}{\mathcal L_{B,1}},
\qquad \frac{\mathcal L_{B,2}}{\mathcal L_{B,1}}=10^{-0.4\Delta M_B}.
\tag{L35}
\]

This makes a proposed magnitude correction into a required response-weighted luminosity ratio; it does not demonstrate that plausible nickel yield/opacity/spectral evolution produces it. For example, a brighter intrinsic population has negative \(\Delta M_B\). Inferring its distance with a fixed dimmer \(M\) underestimates the distance modulus by that same negative magnitude offset. A drift in bolometric output can differ from a drift in \(B\) output if the SED changes.

## 9. Verification and findings that can be acted on

The reproducible check is [luminosity_check.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/physics_audit/luminosity_check.py); its result and SHA-256 input hashes are [luminosity-check.json](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/physics_audit/luminosity-check.json). The original execution settings are retained in that historical record.

The original execution commands are preserved with the linked historical calculation records. Use the [workflow guide](../workflows.md) for current commands.

| Check | Executed result | What this establishes |
|---|---|---|
| L08 analytic abundances versus independent ODE integration | Maximum normalized abundance difference \(1.62\times10^{-13}\) | Chain solution, signs, buildup factor, and conserved total number |
| L09 time integral versus L11 | Relative discrepancy \(2.3\times10^{-16}\) | Correct channel-energy budget and time-unit conversion |
| L22 quadrature versus L21 numerical ODE | Maximum relative difference \(2.45\times10^{-10}\) at the recorded times | Integrating factor and diffusion-convolution normalization |
| L17 time-weighted conservation | Fractional residual below \(4.11\times10^{-12}\) | Expansion loss and stored-energy bookkeeping within the tested model |
| Planck integral versus Stefan–Boltzmann coefficient | Agreement to floating-point precision | The frequency integral and angular \(\pi\) normalization |
| SALT native-grid nonnegativity, \(-15\le p\le45\) d, \(3500\le\lambda\le8000\) Å, \(-3\le x_1\le3\) | 703 of 27,511 nodes allow \(M_0+x_1M_1<0\) for some allowed \(x_1\) | An empirical linear-model support limitation, not a physical negative luminosity or a measured fitted-distance bias |

The numerical example fixes \(M_{\rm Ni,0}=0.6M_\odot\), \(\tau_m=13\) d, \(t_\gamma=35\) d, and trapped positron kinetic energy. It peaks at 15.2871 d with \(L=1.42760\times10^{43}\) erg s\(^{-1}\); \(L/\dot Q_{\rm dep}=1\), whereas \(L/\dot Q_{\rm rad}=0.994781\). These are synthetic values demonstrating the distinction between deposition and decay generation, not fitted properties of any observed supernova. Using the explicit modern constants in the script and the rounded historical decay energies yields L10 coefficients \(6.44259\) and \(1.44299\), each times \(10^{43}\), and total energy \(1.87747\times10^{50}\) erg per initial solar mass. Small differences from the historical rounded coefficients are convention/rounding differences, not an inconsistent chain.

The SALT support check independently reproduces the prior grid-audit finding. It tests endpoints in \(x_1\), which are sufficient because the surface is linear in \(x_1\). It does not test interpolated passbands at actual fitted epochs or imply 703 affected observations. A positive color multiplier cannot repair a negative monochromatic surface. Before interpreting these regions physically, the actual observing weights and fitted parameter support must be traced; arbitrary clipping would change the model and requires its own validation.

| Project statement or possible assumption | Audit status and correction |
|---|---|
| SALT is an empirical spectral representation | Verified in source; retain this description. |
| The current code derives luminosity from foundational explosion physics | False if asserted; no reaction network, heating/deposition calculation, or ejecta transfer solver is present in the active analysis scripts. |
| The same \(M_0,M_1,CL\) fit light curves, therefore luminosity evolution has been excluded | Does not follow. Free amplitude, limited spectral coordinates, population selection and standardized absolute magnitude remain separate. |
| Decay generation equals observed luminosity | Replace by L12, L15, and a justified transport solution. Equality holds only in specified limits. |
| \(L_{\rm peak}=\dot Q_{\rm rad}\) is a conservation law | Replace by conditional L23 with deposited power, or use L17's broader conservation identity. |
| A blackbody temperature or bolometric diffusion time predicts \(B\)-band standardization exactly | Unsupported; requires frequency-dependent line/ionization transport and photometric mapping. |
| An environmental or age slope uniquely fixes nickel yield evolution | Unsupported; no mapping from age/metallicity/binary history to ejecta yield/structure and then to standardized \(B\) luminosity is implemented. |
| Passing these equation checks validates SALT training or rules out an astrophysical cosmology bias | False; the checks validate algebra and explicitly stated limits only. |

The project should retain its empirically useful standardization while marking its physical origin and population transport as testable hypotheses. No new numerical dust, age, or luminosity correction to the data is justified solely by the equations or synthetic examples in this chapter.
