# From measurements to acceleration

## Metric and kinematic assumptions

In a homogeneous and isotropic FLRW geometry, define H=ȧ/a, 1+z=a₀/a(t), and q=−aä/ȧ². Expansion is H>0; acceleration is q<0. Differentiating redshift gives dt/dz=−1/[(1+z)H], hence

\[
q(z)=-1+(1+z)\frac{H'(z)}{H(z)},\qquad
E(z)=\frac{H(z)}{H_0}=\exp\left[\int_0^z\frac{1+q(u)}{1+u}\,du\right].
\]

Let χ=∫₀ᶻdu/E(u). With Ωk=−kc²/(a₀²H₀²), the transverse comoving distance is (c/H₀)Sₖ(χ), where Sₖ=χ for Ωk=0, sinh(√Ωkχ)/√Ωk for positive Ωk, and sin(√−Ωkχ)/√−Ωk for negative Ωk. Transparent photon propagation and distance duality give dL=(1+z)DM. In the released peculiar-velocity-corrected SN convention, replace the prefactor by (1+zHEL) and use the cosmological/Hubble-diagram redshift in the integral. These assumptions must not be silently carried between redshift frames.

For dimensionless D=H₀DM/c, E=√(1+ΩkD²)/D′. Consequently

\[
q=-1+(1+z)\left[\frac{\Omega_k D D'}{1+\Omega_kD^2}-\frac{D''}{D'}\right].
\]

At zero curvature this becomes −1−(1+z)D″/D′. Reconstructing acceleration therefore requires two distance derivatives and a curvature assumption (or independent curvature information). Noisy, finite-redshift SN measurements cannot supply an arbitrarily resolved q(z), especially a point value at z=0, without regularity or parametric assumptions. A q-bin fit reports averages under its parameterization; it does not directly observe q(0).

## Dynamics adds assumptions

GR plus separately conserved nonrelativistic matter and dark energy gives, neglecting radiation at the fitted SN redshifts,

\[
E^2=\Omega_m(1+z)^3+\Omega_k(1+z)^2+
\Omega_{DE}\exp\left[3\int_0^z\frac{1+w(u)}{1+u}du\right].
\]

For flat CPL w(z)=w₀+wₐz/(1+z), the dark-energy factor is (1+z)^{3(1+w₀+wₐ)}exp[−3wₐz/(1+z)]. It follows that q₀=½+3w₀(1−Ωm)/2, and q(z)=½[Ωm(z)+(1+3w(z))ΩDE(z)]. Radiation would add Ωr(z); curvature contributes zero to the acceleration numerator. The low-redshift calculator deliberately neglects radiation, whereas the public CMB chains include early-universe physics. This numerical approximation is separate from the much stronger physical assumption of constant SN luminosity after standardization.

The acceleration condition ρ+3p<0 follows from the GR acceleration equation. It is not needed to define or reconstruct kinematic q. Rejecting the fixed point (w₀,wₐ)=(−1,0) is different from rejecting q₀<0: an evolving dark energy model can accelerate today, and a currently decelerating model can have accelerated previously. CPL extrapolations into the future are model extrapolations, not observations of future dynamics.

## What supernova magnitudes identify

At the released distance-product level, write m_std(z)=𝓜+5log₁₀[(1+zHEL)Sₖ(χ)]+ΔM(z)+ε. The intercept 𝓜 includes the unknown standardized absolute magnitude and H₀. With covariance C, integrating a flat prior over 𝓜 yields the quadratic form rᵀAr, where

\[
A=C^{-1}-\frac{C^{-1}11^TC^{-1}}{1^TC^{-1}1}.
\]

The omitted normalization is parameter-independent when C is fixed. That is why an uncalibrated SN-only distance-shape fit does not determine H₀. It does not imply that an astrophysical correction expressed in Gyr is independent of H₀: the cosmic age scale still contains H₀⁻¹, and a DTD with a fixed physical delay adds another clock.

For any two expansion histories with distances d₁ and d₂, the replacement

\[
\Delta M_2(z)=\Delta M_1(z)+5\log_{10}[d_1(z)/d_2(z)]
\]

leaves every predicted apparent magnitude identical. This is an exact observational non-identifiability if arbitrary luminosity evolution is admitted. It is not evidence that a physical SN population realizes the required function. Age, dust, selection and calibration measurements constrain that freedom only through their own measurement/selection/transport models. A low-redshift host-age correlation alone does not identify its mean high-redshift continuation.

A universal slope error σs in a correction s f(z) contributes σs² f fᵀ, a correlated, rank-one covariance, or equivalently a single shared latent slope. Adding σs²fᵢ² only to each diagonal spuriously makes one physical uncertainty average down with SN count. Our sensitivity fits use a common amplitude where slope uncertainty is enabled. Neither option represents the entire uncertainty in SFH, DTD, dust or the transfer from hosts to progenitors.

## Additional probes

BAO measures combinations DM/rd, DH/rd=c/(Hrd), and DV/rd=[z DM²DH]^{1/3}/rd in a specified reconstruction/likelihood. A low-redshift BAO fit can leave H₀rd free without a CMB ruler calibration; its q₀ still extrapolates below the first BAO measurement. Calibrating rd with the CMB imports an early-universe and gravity model, recombination and neutrino assumptions, and correlated likelihood choices. A distance-modulus plot constructed by fixing H₀rd to a preferred BAO+CMB model is therefore not an independent, model-free agreement test.

Constant-w fits to CMB, BAO and SNe are not localized measurements of w at each probe's characteristic redshift. Their distance integrals and parameter degeneracies overlap; ordering their fitted w values with probe redshift is not itself a measured w(z) history. Likewise, a reduced one-dimensional KL divergence after adding data does not by itself validate a correction or establish cross-probe concordance. Joint likelihood contributions and posterior predictive discrepancies provide more direct checks.

## Numerical record

The independent implementation is `scripts/cosmology/core.py`; `scripts/cosmology/run.py validate` checks Astropy distances, analytic constant-q limits, q from a numerical H derivative, offset invariance, covariance handling, and compressed versus full likelihoods over broad parameter draws. The validation manifests record input and code hashes. The declared reconstruction assumptions are part of the result, not removable caveats.

Primary context: [Son 2025](https://doi.org/10.1093/mnras/staf1685), [Pantheon+ cosmology](https://arxiv.org/abs/2202.04077), [DESI DR2](https://arxiv.org/abs/2503.14738), and [Sah 2026](https://arxiv.org/abs/2606.09650). Equations above are derived directly rather than inferred from agreement between those papers.
