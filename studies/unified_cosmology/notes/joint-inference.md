# A shared cosmological measurement

The supernova, baryon-acoustic-oscillation (BAO) and cosmic-microwave-background (CMB) observations constrain different parts of the same expansion model. Their combination requires one background cosmology and one sound horizon. It does not require treating a published sound-horizon estimate as an additional independent observation.

## Observations and sample identity

The main supernova input is the current Dovekie release: **1,820 distance rows**, comprising 1,623 DES, 117 Foundation and 80 other low-redshift objects. The packed matrix is a **total precision matrix**. We reconstruct its symmetric matrix, invert it to obtain covariance, and preserve the exact distance-row ordering. Subsetting a precision matrix would condition on the omitted observations; the appropriate marginal subset comes from the covariance instead. The source and independent calculation checks are in [the survey note](survey-selection.md).

The 13 DESI DR2 measurements retain their full released covariance and observable identities. Pantheon+ and the historical DES3YR distance likelihood are alternatives, not extra independent supernova samples. Their shared events and calibration prevent multiplying them with Dovekie without a cross-covariance. DES3YR's 20 released bins represent 329 supernovae and contain 18 informative bins; they are not 20 individual supernovae.

The original CMB target is full Planck 2018 TT/TE/EE, low-multipole temperature and polarization, and lensing. A separate contemporary comparison adds the explicitly specified ACT and SPT likelihoods with their stated Planck cuts. These targets have different observations and nuisance assumptions. The [external-probe note](external-probes.md) distinguishes the paper description, released-chain configuration and directly evaluated public implementation.

## Equations and units

For a spatially flat background, a radial comoving distance is

\[
D_M(z)=c\int_0^z\frac{dz'}{H(z')},\qquad D_A(z)=D_M(z)/(1+z).
\]

The released supernova convention requires the cosmological-frame redshift in the distance integral and the heliocentric redshift in the photon redshift factor:

\[
D_L=D_A(z_{\rm HD})(1+z_{\rm HD})(1+z_{\rm HEL}),\qquad
\mu=5\log_{10}[D_L/{\rm Mpc}]+25.
\]

For BAO, the predictions are $D_M/r_d$, $D_H/r_d=c/[H(z)r_d]$, or $D_V/r_d=[zD_HD_M^2]^{1/3}/r_d$, as indicated by each released row. The late-time calculation leaves $H_0r_d$ free, in **km/s**. The CMB calculation derives $r_d$ in Mpc from the same physical baryon and matter densities that produce its temperature, polarization and lensing spectra. It does not multiply in a second Planck-derived ruler prior.

The late-time calculation uses matter plus dark energy over $z\leq2.33$ and omits radiation. The CMB calculation retains radiation, neutrinos and recombination through CAMB. These are separately labelled approximations rather than silently identical backgrounds.

## Absolute brightness and unmeasured evolution

The standardized supernova prediction contains an unknown common brightness offset $M$ and an optional residual mean luminosity term $B(z)$:

\[
m_i^{\rm pred}=M+\mu(z_i)+B(z_i).
\]

Let $r=\mu+B-m^{\rm observed}$, $P=C^{-1}$ and $u=P\mathbf1$. Integrating the unknown offset with a flat prior gives

\[
\chi^2_{\rm SN}=r^TPr-\frac{(r^Tu)^2}{\mathbf1^TP\mathbf1}.
\]

Subtracting a constant from $r$ before this calculation improves numerical stability without changing its value. The likelihood retains the covariance determinant and the offset-integral normalization. Consequently, a software field named `chi2__released_sn`, which includes those constants, is not the conventional quadratic statistic. The separately exported `sn_chi2` is that statistic. The unspecified normalization of the improper offset prior precludes absolute evidence claims.

The baseline sets $B=0$ after the released corrections. Sensitivity calculations allow either

\[
B(z)=\epsilon\frac{\log(1+z)}{\log2},\qquad
\epsilon\sim U(-0.5,0.5)\ {m mag},
\]

or four natural-cubic shape coefficients at $z=0.1,0.4,0.8,1.3$, with zero value at $z=0$ and independent Gaussian widths of 0.1 or 0.3 mag. For the latter model, Gaussian integration is exact: replace $C$ by $C+\sigma_B^2FF^T$, where $F$ contains the basis functions. These widths are **chosen sensitivity priors**, not host-age measurements. Positive $B$ means a dimmer standardized supernova at fixed distance.

Arbitrary grey evolution is exactly degenerate with the supernova distance curve: $\mu\rightarrow\mu+f(z)$ and $B\rightarrow B-f(z)$ preserve the prediction. BAO and CMB add independent physical constraints on the expansion history, but they do not by themselves identify age as the cause of a fitted brightness term. The [host observations](host-likelihood.md) and selection-conditioned population likelihood remain necessary for that causal interpretation.

## Expansion and the rate of change of acceleration

We report

\[
q=-\frac{\ddot a}{aH^2}=(1+z)\frac{H'}H-1,
\qquad
j=\frac{\dddot a}{aH^3}=q(2q+1)+(1+z)\frac{dq}{dz}.
\]

Negative $q$ means accelerating expansion. If $q<0$, negative $j$ means that the positive scale-factor acceleration is decreasing with time. A decreasing acceleration is not yet deceleration. The sign of the time derivative of $q$ is a third, distinct diagnostic because the normalization $aH^2$ also evolves.

For the late-time CPL model, $w(z)=w_0+w_a z/(1+z)$ and

\[
E^2=\Omega_m(1+z)^3+(1-\Omega_m)(1+z)^{3(1+w_0+w_a)}
\exp[-3w_a z/(1+z)].
\]

Writing $f_{\rm DE}$ for the dark-energy fraction at that redshift gives

\[
q=\tfrac12+\tfrac32 f_{\rm DE}w,
\quad
j=1+\tfrac92 f_{\rm DE}w(1+w)+\tfrac32 f_{\rm DE}\frac{w_a}{1+z}.
\]

The full CMB analysis instead differentiates its actual CAMB background. Independent analytic and numerical checks test the sign, units and derivative formulas. In matter-plus-Λ cosmology $j=1$ is imposed by the model; its zero posterior uncertainty is not a separate observational measurement.

## Priors, support and numerical reliability

The exact parameter ranges and decisions precede the corresponding posterior interpretation in [the original design](../code/inference/design.json) and [the contemporary extension](../code/inference/modern-inference-design.json). Spatial flatness, general relativity, fixed neutrino mass and primordial-spectrum assumptions are part of the measurement. It is not a simultaneous test of every cosmological theory.

CAMB 1.6.6's ordinary CPL interface rejects $w_0+w_a>0$. This is an effective support restriction even if no separate prior names it. The late-time calculation has no such cut. A separate table-based background experiment can evaluate some of that region, but fixed-parameter failures to fit the CMB are not a profiled exclusion of the entire region.

Four independently initialized chains, rank-normalized split diagnostics, effective sample sizes, prior boundaries and numerical-domain failures accompany posterior intervals. Interacting ensemble walkers are not mislabelled independent chains. An independently seeded nested-sampling calculation checks difficult late-time tails. Optimized points alone do not establish posterior uncertainties or a significance against ΛCDM.

Any interpolated CMB spectrum is only a computational proposal. Exact-provider checks, untouched validation spectra and direct CAMB evaluation of selected posterior points precede qualified intervals. Untrimmed density-ratio weights correct the proposal to the stated numerical target; weight concentration and tail diagnostics test that correction. Finite checks cannot prove that an unvisited distant mode does not exist.

The remaining physical boundary is explicit: a released-distance joint fit conditions on the supernova reduction, standardization, contaminants, selection and covariance. Successful numerical inference does not turn those upstream assumptions into observations, nor does it establish a measured age-dependent replacement correction.
