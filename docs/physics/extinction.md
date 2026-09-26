# Extinction, attenuation and detected light: physical derivation and implementation audit

Implementation references in this chapter identify the preserved source snapshot; see [physical foundations](README.md) for the current execution boundary.

This chapter establishes the physical chain from grain interactions to measured photon counts, then separates that chain from the empirical extinction curves and population assumptions used in this workspace. Equations D01–D28 use LaTeX. They are not a claim that the repository solves microscopic dust physics or the full radiation-transfer equation.

**Finding:** the inspected flux code uses the correct exponential attenuation, observer/rest-frame wavelength distinction and photon-counting wavelength factor. The historical SNANA option called F99 is an approximation to that empirical curve. Both that approximation and the newer F99 spline produce *negative extinction* when extrapolated into part of the released low-R_V population. That result fails the passive-screen interpretation. Replacing the approximation by the spline alone does not repair that population support. A small independent slab calculation also failed at exactly zero optical depth; that local numerical defect is corrected in this audit.

The executable evidence is [extinction_check.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/physics_audit/extinction_check.py), with newly generated results in [extinction-check.json](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/physics_audit/extinction-check.json). Existing audit reports were used to locate relevant code, then the equations, source branches and raw simulation headers were checked again. Only the named routines and bounded tests below were verified; this is not certification of the complete SNANA package, the released light-curve pipeline or its cosmological consequences.

## Quantities and dimensions

| Symbol | Meaning and convention | Dimensions |
|---|---|---|
| \(I_\lambda\) | Specific intensity in a specified ray direction | energy time\(^{-1}\) area\(^{-1}\) solid-angle\(^{-1}\) wavelength\(^{-1}\) |
| \(L_\lambda\) | Emitted spectral luminosity, integrated over emission direction | energy time\(^{-1}\) wavelength\(^{-1}\) |
| \(f_\lambda\) | Spectral energy flux at the telescope | energy time\(^{-1}\) area\(^{-1}\) wavelength\(^{-1}\) |
| \(C_{\rm abs},C_{\rm sca},C_{\rm ext}\) | Single-grain cross sections | area |
| \(\alpha_\lambda\) | Extinction per path length; not SALT stretch coefficient \(\alpha\) | length\(^{-1}\) |
| \(\kappa_\lambda\) | Extinction per dust mass | area mass\(^{-1}\) |
| \(\tau_\lambda\) | Optical depth | dimensionless |
| \(A_\lambda,A_b,E(B-V)\) | Monochromatic extinction, band extinction, color excess | logarithmic magnitudes |
| \(R_V,R_B\) | Ratios of magnitude differences | dimensionless |
| \(T_b(\lambda)\) | Photon detection throughput in passband \(b\) | dimensionless |
| \(\lambda_o,\lambda_r\) | Observer and source-rest wavelength | Å in the inspected code |

The code uses \(10^4/\lambda[\text{Å}]\) for inverse microns. When \(f_\lambda\) is per Å, the photon conversion must use \(hc\) in energy·Å, or convert **both** the wavelength and spectral density consistently. Multiplicative attenuation is dimensionless. Magnitudes are logarithms of flux ratios; they are not energies or cross sections.

## D01–D04: the physical origin of a screen extinction law

For stationary material, elastic scattering and scalar, unpolarized radiation, the local transfer equation is

\[
\boxed{\frac{1}{c_{\rm light}}\frac{\partial I_\lambda}{\partial t}
 +\mathbf n\!\cdot\!\nabla I_\lambda
 =-\alpha_{\rm ext,\lambda} I_\lambda
 +j_{\rm em,\lambda}
 +\alpha_{\rm sca,\lambda}
   \int p_\lambda(\mathbf n',\mathbf n)I_\lambda(\mathbf n')\,d\Omega'.}
\tag{D01}
\]

The scattering phase function is normalized over outgoing directions to one and has units sr\(^{-1}\). The emissivity \(j\) has units \(I\)/length. This is photon/energy bookkeeping along rays; its loss term counts both absorption and scattering **out of** the ray. Its gain term includes scattering **into** the ray. Frequency redistribution, polarization, relativistic material motion and cosmological transport need a more general transfer equation. Those effects are not implicitly included in D01.

For grains of composition \(j\) and radius \(a\), with \(n_j(a,s)\,da\) the number per physical volume in a size interval,

\[
\boxed{\begin{aligned}
C_{\rm ext}&=C_{\rm abs}+C_{\rm sca},\qquad C_i=\pi a^2 Q_i,\\
\alpha_{\rm ext,\lambda}(s)
 &=\sum_j\int n_j(a,s)C_{{\rm ext},j}(a,\lambda)\,da
 =\rho_d(s)\kappa_{\rm ext,\lambda}(s),\\
\tau_\lambda&=\int_{\rm path}\alpha_{\rm ext,\lambda}(s)\,ds.
\end{aligned}}\tag{D02}
\]

Each factor has a physical role: grain dielectric response and shape determine cross sections; the grain size/composition distribution determines the extinction curve; the column density determines its amplitude. Maxwell's equations plus material constitutive relations, scattering boundary conditions and the Poynting energy balance underlie the cross sections. They do **not** determine the astrophysical size distribution or material abundances without observations and additional formation/destruction physics. For passive grains \(C_{\rm abs},C_{\rm sca}\geq0\), hence \(\tau_\lambda\geq0\). The underlying optical constants are themselves measured or modeled inputs. Draine's [optical-property calculations](https://www.astro.princeton.edu/~draine/dust/dust.diel.html) exemplify this physical route; the local F99 routines do not perform it.

Dropping emission and in-scattering **in the measured beam**, D01 reduces to \(dI_\lambda/ds=-\alpha_{\rm ext,\lambda}I_\lambda\). Separating variables gives

\[
\boxed{I_\lambda^{\rm out}=I_\lambda^{\rm in}e^{-\tau_\lambda},\qquad
 f_\lambda^{\rm out}=f_\lambda^{\rm in}e^{-\tau_\lambda},\qquad
 0<e^{-\tau_\lambda}\leq1.}\tag{D03}
\]

This foreground-screen result applies when unresolved scattered light and dust emission make negligible contributions in the measured aperture and time interval. It is not the general solution for mixed emitters and dust, clumpy unresolved regions, light echoes or circumstellar scattering.

The astronomical magnitude definition then supplies, without a fitted physics coefficient,

\[
\boxed{A_\lambda\equiv-2.5\log_{10}
 \frac{f_\lambda^{\rm out}}{f_\lambda^{\rm in}}
 =\frac{2.5}{\ln10}\tau_\lambda
 =1.0857362048\,\tau_\lambda,
 \quad e^{-\tau_\lambda}=10^{-0.4A_\lambda}.}\tag{D04}
\]

Thus the `10**(-.4*A)` implementations are correctly normalized. For a passive screen, \(A_\lambda<0\) is impossible. A negative *fitted correction*, noise realization, differential color warp, or an echo contribution relative to the instantaneous direct light is a different quantity; it is not proof of negative absorption. The SNANA comments allowing negative `AV` for spectral warping do not make negative attenuation physically valid for an actual dust population.

## D05–D08: redshift, dust frames and photons in a filter

For isotropic-equivalent source spectral luminosity and luminosity distance \(d_L\), the relation consistent with \(f_{\rm bol}=L_{\rm bol}/(4\pi d_L^2)\) is

\[
\boxed{f^0_{\lambda_o}(t_o)
 =\frac{L_{\lambda_r}(t_r)}{4\pi d_L^2(1+z)},\qquad
 \lambda_r=\frac{\lambda_o}{1+z},\qquad
 t_r=\frac{t_o-t_{0,o}}{1+z}.}\tag{D05}
\]

The Jacobian \(d\lambda_o=(1+z)d\lambda_r\) cancels the explicit spectral factor when integrated over wavelength. In frequency units the corresponding factor is reversed: \(f_{\nu_o}=(1+z)L_{(1+z)\nu_o}/(4\pi d_L^2)\). Mixing the frequency and wavelength formulae is a common incorrect extra factor of \((1+z)^2\). In the local SALT engine the fitted amplitude absorbs the distance scale; its explicit \((1+z)^{-1}\) is the spectral Jacobian, not an independent computation of inverse-square luminosity dilution.

With a source-host screen and a Milky Way screen,

\[
\boxed{f_{\lambda_o}(t_o)=f^0_{\lambda_o}(t_o)
 10^{-0.4[A_{\rm host}(\lambda_o/(1+z))+A_{\rm MW}(\lambda_o)]}.}\tag{D06}
\]

An intervening screen at \(z_d\) uses \(A_d(\lambda_o/(1+z_d))\). Optical depths, and therefore monochromatic extinctions evaluated along the same observed ray, add. The relevant redshift maps the received spectrum to the dust/source frame; it is not automatically the Hubble-diagram redshift after a peculiar-velocity correction.

**Implementation:** [flux_response.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/salt_dust_audit/flux_response.py), `build_model`, creates separate `effect_frames=['obs','rest']`. [engine.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/phase2/independent_flux/engine.py), `prepare`, evaluates MW F99 at `wave` and the optional host-like shape at `rest=wave/(1+z)`. Pinned `sncosmo/models.py`, `_flux`, transforms phase and wavelength and dispatches each effect in its stated frame. Historical `genmag_SEDtools.c`, lines 2368–2375 and 2466–2472, calls `GALextinct` with `LAMOBS` for MW and `LAMREST` for host. These inspected transformations are correct. At observed 8000 Å, \(z=0.5,E=0.1,R_V=3.1\), the independently evaluated host term is 0.31673744 mag and the MW term is 0.17153438 mag: using one in place of the other is materially different.

For a detector counting absorbed photons, each incident photon carries \(hc_{\rm light}/\lambda\), so the expected source count rate is

\[
\boxed{\dot N_b(t_o)=\mathcal A\int
       f_{\lambda_o}(t_o)T_b(\lambda_o)
       \frac{\lambda_o}{hc_{\rm light}}\,d\lambda_o.}\tag{D07}
\]

\(\mathcal A\) is collecting area, optionally absorbed into the response. Multiply by exposure duration only if the light curve is effectively constant over the exposure; otherwise integrate D07 over time. A bolometer measures an energy-weighted integral instead. This chapter uses the ordinary photon-throughput definition of \(T_b\); one must not insert another \(\lambda\) if a published response already incorporates photon weighting.

For a calibration primary with spectrum \(f_{\lambda,\rm ref}\) and assigned magnitude \(m_{{\rm ref},b}\), the matched photon integrals give

\[
\boxed{\begin{aligned}
m_b&=m_{{\rm ref},b}-2.5\log_{10}
 \frac{\int\lambda T_b f_\lambda\,d\lambda}
      {\int\lambda T_b f_{\lambda,\rm ref}\,d\lambda},\\
F_{{\rm cal},b}&=10^{0.4(ZP-m_b)}.
\end{aligned}}\tag{D08}
\]

Collecting area and \(hc\) cancel. A flat \(f_\nu\) AB primary has \(f_\lambda\propto\lambda^{-2}\), with unit conversion included. SNANA's `FLUXCAL` convention here has \(ZP=27.5\); these calibrated values are not raw photoelectrons. `MAG_OFFSET=.27` is an empirical SALT normalization convention, not an extinction coefficient.

**Implementation:** `engine.py:prepare` uses `wave*trans` in both source and reference integrals and `1/(1+z)` for the source spectrum. Pinned `sncosmo/models.py:_bandflux_single` uses `sum(wave*trans*f)*dwave/HC_ERG_AA`. Historical `genmag_SALT2.c:2921` uses `LAMSED*TRANS` within its own integration/normalization convention. A physically complete band integral must cover **all nonzero throughput** within the model's spectral support; dropping unsupported nonzero filter wings while keeping a full reference denominator is not equivalent to D08. The companion implementation audit corrects that local `engine.prepare` boundary behavior rather than claiming that a mean-wavelength cut proves full support.

## D09–D12: broadband extinction and the meaning of R

Define \(W_b(\lambda,t)=\lambda T_b(\lambda)f^0_\lambda(t)\geq0\). The extinction of a *band measurement* is

\[
\boxed{A_b(t;E,R)= -2.5\log_{10}
 \frac{\int W_b(\lambda,t)10^{-0.4E k(\lambda,R)}\,d\lambda}
      {\int W_b(\lambda,t)\,d\lambda},\qquad
 k(\lambda,R)\equiv A_\lambda/E.}\tag{D09}
\]

This is an integral followed by a logarithm. In general it is neither \(A(\lambda_{\rm eff})\) nor an unweighted mean of \(A_\lambda\). It depends on the source spectrum, phase, band shape, redshift and extinction amplitude. For a nonnegative spectrum and throughput, if every contributing \(A_\lambda\geq0\), D09 proves \(A_b\geq0\). If a trained basis combination predicts negative spectral flux, this probability-weight argument fails; that is another model-domain issue, not a property of physical dust.

Let \(K=0.4\ln10\) and let \(\langle\cdot\rangle_E\) denote the band-weighted average with the already-reddened weight \(W_be^{-KEk}\). Differentiating D09 gives

\[
\boxed{\frac{dA_b}{dE}=\langle k\rangle_E,
 \qquad\frac{d^2A_b}{dE^2}=-K\,\operatorname{Var}_E(k)\leq0.}\tag{D10}
\]

The second derivative is not zero for a broad filter with a varying curve. Equivalently,

\[
\boxed{A_b(E)=E\langle k\rangle_0
 -\frac{K E^2}{2}\operatorname{Var}_0(k)+O(E^3).}\tag{D11}
\]

This derives the limitation of a constant broadband coefficient directly, rather than merely asserting a “bandpass effect.” In a constructed positive \(f_\lambda\propto\lambda^{-2}\) spectrum through a 4000–8000 Å top-hat, the fresh check finds \(A_b(0.1)=0.29594593\) mag versus the linear prediction 0.29895538 mag. The derivative and curvature independently agree with numerical flux perturbations. This illustrative filter is not the DES survey selection and its difference is not a cosmological bias.

For one explicitly specified source spectrum, epoch, pair of passbands and dust amount,

\[
\boxed{E_{B-V}\equiv A_B-A_V,\qquad
 R_V^{\rm band}\equiv\frac{A_V}{E_{B-V}},\qquad
 R_B^{\rm band}\equiv\frac{A_B}{E_{B-V}}
 =R_V^{\rm band}+1.}\tag{D12}
\]

This last identity is exact **by definition**, if the denominator and bands are the same and \(E_{B-V}\ne0\). It is not a derivation of a universal physical value for either ratio. Nor does it establish that the input shape parameter called `RV` in a monochromatic curve, the SALT color coefficient, and a broad-band ratio for every supernova spectrum coincide. At zero reddening use an appropriate limiting derivative if defining a ratio. Grey extinction gives \(E_{B-V}=0\) with nonzero \(A_V\), so this parameterization is singular even though the dust may be physical.

F99 constructs a monochromatic curve informed by synthetic broadband stellar colors. The paper explicitly treats filter bandwidth effects; its [appendix discussion](https://ned.ipac.caltech.edu/level5/Fitzpatrick/Fitz_append.html) explains why a curve evaluated at one filter effective wavelength is not the broadband answer. The local hierarchy already labels `R_B` an **effective SALT-coordinate coefficient** in [generative-model.md](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/docs/phase2/generative-model.md); that qualification is necessary. A fit to that parameter does not measure \(R_V\) merely by subtracting one.

## D13–D17: exactly which empirical extinction curves are implemented

CCM89, O'Donnell94 and F99 are empirical families constrained by sightline measurements. They are not consequences of Maxwell's equations with no additional assumptions. The original [CCM paper](https://doi.org/10.1086/167900) and [O'Donnell study](https://ntrs.nasa.gov/citations/19950037261) describe fitting mean extinction behavior. The equations below record the active families relevant to this workspace and their actual local coefficients; unused newer SNANA alternatives are outside this scoped verification.

Set \(x=10^4/\lambda[\text{Å}]\), numerically in \(\mu\mathrm m^{-1}\). CCM/O'Donnell have the form

\[
\boxed{A_\lambda=A_V\,[a(x)+b(x)/R_V]
       =E\,[R_Va(x)+b(x)].}\tag{D13}
\]

The coefficients below use dimensionless numerical \(x\) expressed in the stated units. For \(0.3\leq x<1.1\),

\[
a=.574x^{1.61},\qquad b=-.527x^{1.61}.
\]

For \(1.1\leq x<3.3\), write \(y=x-1.82\), \(a=\sum_i a_i y^i\), \(b=\sum_i b_i y^i\). The full coefficient arrays in ascending power are

\[
\begin{aligned}
\mathbf a_{\rm CCM}&=(1,.17699,-.50447,-.02427,.72085,.01979,-.77530,.32999),\\
\mathbf b_{\rm CCM}&=(0,1.41338,2.28305,1.07233,-5.38434,-.62251,5.30260,-2.09002),\\
\mathbf a_{\rm O94}&=(1,.104,-.609,.701,1.137,-1.718,-.827,1.647,-.505),\\
\mathbf b_{\rm O94}&=(0,1.952,2.908,-3.989,-7.985,11.102,5.491,-10.805,3.347).
\end{aligned}\tag{D14}
\]

O'Donnell changes this optical/near-UV segment; it is not interchangeable with CCM or F99. For \(3.3\leq x<8\), both use

\[
\begin{aligned}
a&=1.752-.316x-\frac{.104}{(x-4.67)^2+.341}+F_a,\\
b&=-3.090+1.825x+\frac{1.206}{(x-4.62)^2+.263}+F_b,\\
(F_a,F_b)&=(0,0)\quad(x<5.9),\\
F_a&=-.04473u^2-.009779u^3,\quad
F_b=.21300u^2+.120700u^3\quad(u=x-5.9\geq0).
\end{aligned}
\]

For \(8\leq x\leq10\), with \(v=x-8\),

\[
\begin{aligned}
a&=-1.073-.628v+.137v^2-.070v^3,\\
b&=13.670+4.257v-.420v^2+.374v^3.
\end{aligned}\tag{D15}
\]

Outside this numerical domain the inspected legacy branch assigns \(a=b=0\). That software fallback must not be interpreted as physical transparency at all other wavelengths. These branches are at historical `MWgaldust.c:356–420` and current `MWgaldust.c:577–641`.

The implemented F99/FM_UNRED family instead writes \(A_\lambda=E[R_V+k_F(x)]\). Its UV function is

\[
\boxed{\begin{aligned}
k_F(x)&=c_1+c_2x+c_3\frac{x^2}{(x^2-x_0^2)^2+x^2\gamma^2}+c_4F(x),\\
F(x)&=\begin{cases}0,&x<5.9,\\ .5392(x-5.9)^2+.05644(x-5.9)^3,&x\geq5.9,\end{cases}\\
c_2&=-.824+4.717/R_V,\quad c_1=2.03-3.007c_2,\\
(c_3,c_4,x_0,\gamma)&=(3.23,.41,4.596,.99).
\end{aligned}}\tag{D16}
\]

For \(\lambda\leq2700\) Å the UV expression is used directly. For longer wavelengths the local implementation uses a **natural cubic spline** in \(x\) through the following knots of \(k_F=A/E-R_V\), with zero second derivative at the endpoints:

| \(x_i\), in inverse microns | \(k_{F,i}\) |
|---|---|
| 0 | \(-R_V\) |
| \(1/2.65\) | \(-.914616129R_V\) |
| \(1/1.22\) | \(-.7325R_V\) |
| \(1/.600\) | \(-.422809+.00270R_V+.000213572R_V^2\) |
| \(1/.547\) | \(-.051354+.00216R_V-.0000735778R_V^2\) |
| \(1/.467\) | \(.700127+.00184R_V-.0000332598R_V^2\) |
| \(1/.411\) | \(1.19456+.01707R_V-.00546959R_V^2+.000797809R_V^3-.0000445636R_V^4\) |
| \(1/.270\), \(1/.260\) | D16 evaluated at each knot |

\[
\boxed{A_\lambda=E\,[R_V+\mathcal S_{\rm natural}(x;\{x_i,k_{F,i}\})]
 \quad(\lambda>2700\text{ Å}).}\tag{D17}
\]

The updated optical anchors in `extinction`/FM_UNRED are not a literal transcription of the original paper's printed table. The [package documentation](https://extinction.readthedocs.io/en/latest/api/extinction.fitzpatrick99.html) makes that distinction explicit and separates computational wavelength coverage from claimed validity. The current SNANA header permits 912–35000 Å and documents a recommended \(2\leq R_V\leq6\) that it **does not enforce**. The local Python effect exposes about 909–60000 Å. Neither exposed range is proof that the empirical curve is observationally reliable everywhere within it, especially in extrapolated low-\(R_V\) populations. Low \(R_V\) is not itself forbidden by physics; applying this particular empirical family outside its validated support needs separate justification.

**Normalization correction:** current `MWgaldust.c:699` and the historical routine comments call `AV` extinction “at 5495 Angstroms.” The inspected F99 function does not enforce that identity. Direct evaluation at \(R_V=3.1\) gives

\[
 A(5495\text{ Å})=0.9790469793\,A_{V,\rm input}.
\]

The input is the conventional F99 normalization/shape parameter, not an exact pointwise \(A_{5495}\) and not a guaranteed band integral for any SN SED. This chapter corrects that interpretation. Vendor source and calibration assets remain unchanged. Renormalizing the curve to force equality would define a different curve and would require explicit normalization/provenance updates.

## D18–D19: the historical approximation and the physical support failure

The historical DES-linked SNANA `2fe0f56` branch does **not** evaluate D16–D17 for option 99. It calculates

\[
\boxed{A^{\rm old99}_\lambda=A^{\rm O94}_\lambda
 \sum_{i=0}^{10}p_i\left(\frac{\lambda}{1000\text{ Å}}\right)^i.}\tag{D18}
\]

The coefficients, in ascending power, are

\[
\begin{aligned}
\mathbf p=(&8.55929205\!\times\!10^{-2},\ 1.91547833,\ -1.65101945,\ .750611119,\\
&-.200041118,\ .0330155576,\ -.00346344458,\ .000230741420,\\
&-9.43018242\!\times\!10^{-6},\ 2.14917977\!\times\!10^{-7},\ -2.08276810\!\times\!10^{-9}).
\end{aligned}
\]

This fixed wavelength-only ratio cannot in general reproduce the full \(R_V\)-dependent F99 family. It is an additional numerical/empirical approximation, not an alternative derivation of extinction physics. Current SNANA option 99 selects D16–D17; current option -99 restores D18. The released simulation input `PIP_D5YR_SIM_V2_DATADESSIM_4D_P21.input:17` requests option 99, whose meaning therefore depends on the executable version. The fresh compiled check establishes that current -99 and historical 99 agree exactly on its grid.

At positive \(E=0.1\) mag and \(R_V=0.4\), the spline gives

\[
\boxed{A_{6500}=-0.01847554\text{ mag},\quad
 A_{8000}=-0.01842972\text{ mag},\quad
 10^{-0.4A_{8000}}=1.01711928>1.}\tag{D19}
\]

At 8000 Å the positive-\(E\) zero crossings occur at \(R_V=0.65323402\) for D18 and \(R_V=0.66338574\) for D17. These numbers follow from the inspected coefficients and an independent root solve. They are not fitted real-SN extinction measurements.

The fresh raw-header census contains 71,946 written simulated objects in 25 files, including 40,369 with `SIM_RV<2`. At the 8000 Å sign boundary, 4,601 have negative historical extinction and 4,791 would have negative extinction under the spline applied to the same `SIM_AV,SIM_RV` values. An example is simulated `SNID=431914`, `RV=.5475049615`, `AV=.004549657926`, giving historical \(A_{8000}=-.0006214315\) mag. The JSON records its full source path and file hashes. These counts refer to written simulation truth, not independent empirical dust observations or to objects necessarily passing every downstream fit/BBC cut.

The source path passes this quantity directly into the host transmission: historical `genmag_SEDtools.c:2469–2472` computes `pow(10,-.4*XT_MAG)`, and `genmag_SALT2.c:2921` multiplies the spectral contribution by `HOSTXT_FRAC`. There is no nonnegativity correction between these inspected points. The direct-beam dust component is therefore allowed to amplify light. This contradicts D02–D04 for passive-screen dust even when its final synthetic summaries happen to fit observations.

**Required physical correction before interpreting this as dust:** choose and document a nonnegative curve over the entire nonzero passband support and the intended population domain; refit the population, then regenerate and refit matched simulations, selection and bias corrections. A hard \(R_V\) truncation or `max(A_lambda,0)` changes the model and is not a derived unique solution. Merely upgrading option 99 changes the curve but retains the demonstrated support failure. This audit preserves the historical assets and does not invent a cosmological correction from these monochromatic values.

## D20–D22: attenuation geometry and the corrected slab calculation

For unresolved emitters distributed through an absorption-only homogeneous slab, assume the same intrinsic emissivity profile at both wavelengths, no scattering, and let \(u\in[0,1]\) be the fraction of the total optical depth between an emitter and the observer. Uniform emissivity gives

\[
\boxed{\frac{F_\lambda}{F^0_\lambda}
 =\int_0^1e^{-\tau_\lambda u}\,du
 =\begin{cases}(1-e^{-\tau_\lambda})/\tau_\lambda,&\tau_\lambda>0,\\1,&\tau_\lambda=0.\end{cases}}
\tag{D20}
\]

Therefore

\[
\boxed{A_{\lambda,\rm slab}
 =-2.5\log_{10}\!\left(\frac{1-e^{-\tau_\lambda}}{\tau_\lambda}\right),\qquad
 A_{\lambda,\rm slab}=\frac{2.5}{\ln10}
 \left(\frac{\tau_\lambda}{2}-\frac{\tau_\lambda^2}{24}
       +\frac{\tau_\lambda^4}{2880}+O(\tau_\lambda^6)\right).}
\tag{D21}
\]

The thin limit averages half the total path depth. For thick slabs \(F/F^0\sim1/\tau\); their attenuation is greyer than a screen with \(F/F^0=e^{-\tau}\). A galaxy-integrated attenuation ratio is consequently not automatically the extinction ratio of a particular SN sightline. This distinction follows even before including real scattering, geometry or stellar-population gradients.

[dust_identifiability.py](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/scripts/standardization/dust_identifiability.py), `slab`, correctly encoded D21 at positive optical depth but formerly returned `NaN` at \(\tau=0\). It now returns exactly zero attenuation there, rejects negative/nonfinite depth, and uses the displayed thin-limit expansion below \(10^{-4}\) to avoid cancellation. The existing positive experiment grid \((.01,.1,.3,1,3,10)\) is bitwise unchanged. Independent numerical integration of D20 agrees to \(1.2\times10^{-16}\) in transmission on the checked grid. Historical output files were not rewritten.

The same toy script assumes a point SN uniformly distributed in \(u\), \(\tau_B=(1+1/R_V)\tau_V\), and a sole detection requirement \(A_{B,\rm SN}<A_{\rm cut}\). Then

\[
\boxed{\begin{aligned}
A_{B,\rm SN}(u)&=\frac{2.5}{\ln10}\tau_Bu,\\
f_{\rm selected}&=\min\!\left[1,\frac{A_{\rm cut}}{(2.5/\ln10)\tau_B}\right],\\
\langle E_{B-V,\rm SN}\rangle_{\rm selected}
 &=\frac{2.5}{\ln10}\frac{\tau_V}{R_V}\frac{f_{\rm selected}}2.
\end{aligned}}\tag{D22}
\]

At zero depth use \(f_{\rm selected}=1\) for a positive cut and zero mean reddening. The script's published `result()` grid uses strictly positive depths; this audit corrected the reusable `slab` primitive, not the whole toy-selection API. D22 is algebraically correct under its declared narrow-band/two-opacity and toy-selection assumptions. It does not describe a calibrated magnitude-limited survey with varying intrinsic luminosity, redshift, cadence or signal-to-noise.

## D23–D24: what a dust map and a scattering alternative add

Diffuse far-infrared dust emission, in an optically thin approximation, has the physical form

\[
\boxed{I_\nu^{\rm dust}\simeq\int\rho_d(s)\kappa_{\rm abs,\nu}(s)
 B_\nu[T_d(s)]\,ds.}\tag{D23}
\]

This thermal emission is not already \(E(B-V)\). Inferring reddening requires dust temperature, opacity/emissivity, calibration against reddened sources, sky/background separation and assumptions about dust variation. The inspected SFD/CSFD code consumes calibrated maps; it does not derive reddening from D23. Its `0.86` factor is an empirical normalization correction associated with the [Schlafly–Finkbeiner calibration](https://arxiv.org/abs/1012.4804), not a physical constant. Historical `MWgaldust.c:279` applies it under the designated map option. `map_dust.py:54` applies it once to the difference of two maps retained on the same original SFD scale. `flux_response.py` uses already-scaled header values. Every new input needs its map normalization identified; applying the factor twice would be an error.

If scattered circumstellar light matters, a useful schematic for a time-variable source is

\[
\boxed{F_\lambda(t)=F^0_\lambda(t)e^{-\tau_{\rm ext,\lambda}}
 +\int_0^\infty\Psi_\lambda(\Delta t)
                  F^0_\lambda(t-\Delta t)\,d\Delta t
 +F_{\lambda,\rm thermal}(t).}\tag{D24}
\]

The transfer kernel \(\Psi\geq0\) contains geometry, phase function, path attenuation and aperture acceptance, and has units of inverse time for this normalization. It cannot be inferred from one \(R_V\) value. Global energy bookkeeping must account for energy absorbed and reradiated, all angles and delayed times. The [Goobar circumstellar calculation](https://arxiv.org/abs/0809.1094) illustrates how geometry and scattering can produce low effective \(R_V\); it does not justify letting a passive-screen coefficient become negative. D24 is an explicit missing model if this mechanism is to be invoked, rather than something silently included in F99.

## D25–D28: microscopic electromagnetic and thermal foundations

To make the physical origin of D02 explicit, use **SI units in this section**: distance in metres, power in watts, \(\epsilon_0\) and \(\mu_0\) the vacuum permittivity and permeability. Maxwell's equations for the microscopic fields and total charge/current are

\[
\boxed{\nabla\!\cdot\!\mathbf E=\rho_q/\epsilon_0,\qquad
 \nabla\!\cdot\!\mathbf B=0,\qquad
 \nabla\!\times\!\mathbf E=-\partial_t\mathbf B,\qquad
 \nabla\!\times\!\mathbf B=\mu_0\mathbf J+\mu_0\epsilon_0\partial_t\mathbf E.}
\tag{D25}
\]

For continuum grain optics, microscopic currents are represented by a material response. In a linear, local, isotropic, nonmagnetic grain, the frequency-domain constitutive relations are \(\mathbf D=\epsilon_0\epsilon_r(\omega)\mathbf E\), \(\mathbf B=\mu_0\mathbf H\). Tangential E and H match at a surface without imposed surface current; normal D and B match without imposed free surface charge. Incident-wave plus outgoing-wave boundary conditions complete the scattering problem. The complex dielectric function \(\epsilon_r\), grain shape and orientation are physical inputs, not universal constants or quantities specified by Maxwell alone. Magnetic or anisotropic materials require corresponding tensor/magnetic response.

Dot the two curl equations with E and B and use the divergence identity for a cross product. This gives Poynting's energy balance and, for harmonic convention \(e^{-i\omega t}\), the passive dielectric heating rate:

\[
\boxed{\begin{aligned}
\partial_t\!\left(\frac{\epsilon_0E^2}{2}+\frac{B^2}{2\mu_0}\right)
 +\nabla\!\cdot\!\left(\frac{\mathbf E\times\mathbf B}{\mu_0}\right)
 &=-\mathbf J\!\cdot\!\mathbf E,\\
\langle P_{\rm abs}\rangle
 &=\frac{\omega\epsilon_0}{2}\int_{\rm grain}
       \operatorname{Im}\epsilon_r(\omega)\,|\mathbf E(\mathbf r)|^2d^3r\geq0.
\end{aligned}}\tag{D26}
\]

The first line accounts for reversible energy exchange as well as dissipation; its field energy is the microscopic vacuum-field energy, not a proposed \(\epsilon(\omega)E^2/2\) storage formula for a dispersive medium. The second line assumes passive response, \(\operatorname{Im}\epsilon_r\geq0\), and that conductivity is already included in \(\epsilon_r\); adding the same conductivity loss again would double count it. Integrating the surface flux separates absorbed power from radiation scattered away. This supplies the positivity in D02–D04.

For a plane wave, \(F_{\rm inc}=\epsilon_0c_{\rm light}|E_0|^2/2\). Define \(C_{\rm abs}=P_{\rm abs}/F_{\rm inc}\) and \(C_{\rm sca}=P_{\rm sca}/F_{\rm inc}\); each has units m\(^2\). A useful derived example is a homogeneous sphere of radius \(a\) in vacuum, small enough that \(ka\ll1\) and \(|\sqrt{\epsilon_r}|ka\ll1\), with \(k=2\pi/\lambda\). Write its electric dipole as \(\mathbf p=\epsilon_0\alpha\mathbf E_0\), so **\(\alpha\) has units m\(^3\)**, not the alternative SI polarizability units that include \(\epsilon_0\). The electrostatic sphere boundary problem and dipole radiation give

\[
\boxed{\begin{aligned}
\alpha_0&=4\pi a^3\frac{\epsilon_r-1}{\epsilon_r+2},\qquad
\alpha\simeq\frac{\alpha_0}{1-i k^3\alpha_0/(6\pi)},\\
C_{\rm ext}&=k\operatorname{Im}\alpha,\qquad
C_{\rm sca}=\frac{k^4}{6\pi}|\alpha|^2,\qquad
C_{\rm abs}=C_{\rm ext}-C_{\rm sca},\\
C_{\rm abs}&\simeq4\pi k a^3\operatorname{Im}\frac{\epsilon_r-1}{\epsilon_r+2},\qquad
C_{\rm sca}\simeq\frac{8\pi}{3}k^4a^6
       \left|\frac{\epsilon_r-1}{\epsilon_r+2}\right|^2.
\end{aligned}}\tag{D27}
\]

The last line gives leading Rayleigh terms away from regimes requiring higher-order corrections. Its familiar \(a^6/\lambda^4\) scattering scaling does not imply an exact universal \(\lambda^{-4}\) dust curve: dielectric dispersion and the size distribution also vary. The radiation-reaction denominator is necessary for consistent dipole energy accounting: a lossless grain scatters and extinguishes, but has zero absorption. Subtracting the scattering term from a purely real *uncorrected* electrostatic polarizability would spuriously produce negative absorption. Finite-size dynamic corrections and higher multipoles are omitted here. [Draine's radiative-reaction derivation](https://adsabs.harvard.edu/pdf/1988ApJ...333..848D), Sections II–III, provides the relevant primary treatment; its Gaussian-unit polarizability differs from the explicitly declared SI-volume convention here by \(4\pi\).

Finally, Kirchhoff's absorption/emission relation and energy conservation establish the thermal basis of D23. For an orientation-independent cross section in radiation with mean intensity \(J_\nu=(4\pi)^{-1}\int I_\nu d\Omega\), an isothermal grain with internal energy \(U_d\), neglecting collisional heating and nonthermal energy losses, obeys

\[
\boxed{\begin{aligned}
\frac{dU_d}{dt}&=4\pi\int_0^\infty C_{\rm abs}(\nu)
           [J_\nu-B_\nu(T_d)]\,d\nu,\\
\text{radiative equilibrium:}\qquad
\int C_{\rm abs}(\nu)J_\nu\,d\nu
 &=\int C_{\rm abs}(\nu)B_\nu(T_d)\,d\nu,\\
B_\nu(T)&=\frac{2h\nu^3}{c_{\rm light}^2}
                 [e^{h\nu/(k_BT)}-1]^{-1}.
\end{aligned}}\tag{D28}
\]

Both sides of the first line are powers; the \(4\pi\) integrates radiance over solid angle. **Absorption**, not total extinction, heats the grain. Summing the emitted \(C_{\rm abs}B_\nu\) over a dilute grain population yields the emissivity used in D23; integrating it along an optically thin sightline gives that equation. Very small grains can undergo temperature fluctuations after individual photons and require an energy-distribution treatment rather than one equilibrium temperature. These regimes are discussed in [Draine's dust-emission review](https://vo.ned.ipac.caltech.edu/level5/Sept19/Draine2/paper.pdf), Section 9.

D25–D28 establish the microscopic laws and a controlled limiting example. This repository does not supply a dielectric-function calculation, a Mie/discrete-dipole solution, a grain-formation/size-distribution theory or a stochastic-heating calculation. Consequently its empirical F99 coefficients and dust populations cannot be described as having been uniquely built from those foundations.

## Verification verdict and remaining physical choices

| Component | Verdict from this audit | What remains assumed or unverified |
|---|---|---|
| Exponential screen factor and magnitude conversion | Correct in the inspected code; derived from D01–D04 | Negligible received scattered/thermal light; appropriate aperture and time regime |
| Observer/rest dust frames and spectral Jacobian | Correct in named SNANA, sncosmo and local paths | Correct redshift metadata for each path; no audit of every vendor call site |
| Photon-weighted calibrated band integral | Correct factors in inspected implementation | Throughput convention, full spectral support, trained-SED validity and calibration assets |
| F99/CCM/O'Donnell numerical equations | Independently reconstructed and compared against freshly compiled routines | Empirical applicability to each SN host; mean law does not eliminate sightline diversity |
| Historical F99 option name | Version-dependent approximation, now explicitly documented | Original executable/configuration must remain part of provenance |
| Negative red-optical extinction in low-\(R_V\) support | Verified physical inconsistency for a passive-screen interpretation | No measured downstream distance or cosmology bias established |
| \(R_B=R_V+1\) | Exact only for consistently defined band quantities | SALT/effective coefficients cannot inherit it without a mapping calculation |
| `AV` as exactly \(A_{5495}\) | Incorrect comment for inspected F99 implementation | A deliberate monochromatic convention would need renormalization and regenerated products |
| Uniform slab | Analytically correct under stated geometry; zero-limit bug fixed | Galaxy geometry, scattering, source distribution and survey selection are illustrative |
| Map normalization | Empirical calibration, applied once in inspected paths | Map-specific spatial errors, contamination and covariance need data-level validation |

The fresh numerical grid comprises 1,803 wavelengths from 1000–10000 Å and eight \(R_V\) values, including deliberately extrapolated cases. Maximum absolute discrepancies in mag per unit \(E\) are \(2.85\times10^{-14}\) for F99 versus fresh C, \(5.81\times10^{-11}\) versus `extinction` 0.4.9, \(3.1\times10^{-14}\) for CCM/O'Donnell versus C and \(3.63\times10^{-12}\) for D18 versus historical C. These are implementation checks on common empirical equations, not independent confirmation of dust physics. A separate positive-support grid across 1000–35000 Å and \(2\leq R_V\leq6\) found no negative values; a finite grid is not proof for every possible input.

To rerun the bounded evidence without changing historical products:

The original execution commands are preserved with the linked historical calculation records. Use the [workflow guide](../workflows.md) for current commands.

The JSON records input hashes, installed numerical versions, actual simulation source paths and all stated numerical diagnostics. Resolving physical population support and cosmological impact requires the full regenerated/refitted analysis described above; a correct transmission formula by itself cannot establish those inferences.
