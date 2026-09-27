# Nebular emission lines and the residual age association

Measured hydrogen emission lines do not settle whether an additional age correction is required. In the 138 hosts with both Balmer lines detected above formal signal-to-noise 3, the conditional age coefficient changes from **−0.0265 ± 0.0212 to −0.0129 ± 0.0233 mag/Gyr** when the measured line ratio is included. Both estimates are uncertain. The direction and size of the change depend on line selection, and the added age term does not improve heldout prediction after the nebular measurement is included.

This exploratory extension was specified after the earlier age and spectral results, before fitting these comparisons. It is not an independently preregistered test. The exact individual ZTF light-curve amplitude blinding transformation remains unknown, so these coefficients cannot be certified as an unblinded cosmological measurement.

## What was measured

We used the existing 165-object ZTF/TITAN spectroscopic brightness cohort and joined on exact supernova, physical-host and 64-bit spectrum identifiers. All 165 joined uniquely. Two objects, ZTF19aamhqow and ZTF20acyvzbr, have invalid negative Hα flux errors; they were excluded from line-based comparisons. The other **163** retain all signed measurements, including **five** nonpositive Balmer measurements. There are 138 hosts with Hα and Hβ S/N>3, 122 with both S/N>5, and 75 with all four diagnostic lines above S/N=3 and a star-forming BPT classification.

MPA-JHU line fluxes are measured after stellar-continuum subtraction and already corrected for Milky Way extinction. Their extraction therefore retains stellar-model dependence, although the observables differ from broad-band SED ages. The central 3-arcsec aperture does not encompass 99 of the 138 primary supernova positions. We applied no second foreground correction. Nebular gas, stellar-continuum attenuation and supernova line-of-sight dust are distinct quantities. See the [SDSS measurement description](https://www.sdss4.org/dr17/spectro/galaxy_mpajhu/).

For detected lines the regressor is

$$D=\ln\!\left[\frac{F_{\mathrm{H}\alpha}}{2.86F_{\mathrm{H}\beta}}\right].$$

The constant 2.86 sets the regression zero point; changing it is absorbed by the intercept. We do not impose nonnegative attenuation or convert the ratio into a supernova extinction correction. First-order error propagation uses

$$\operatorname{Var}(D)=
\frac{s_\alpha^2}{F_\alpha^2}+
\frac{s_\beta^2}{F_\beta^2}-
\frac{2\rho s_\alpha s_\beta}{F_\alpha F_\beta},$$

with primary $\rho=0$ and sensitivities −0.5 and +0.5 because inter-line covariance is unavailable.

The broad 163-host check uses the finite signed contrast

$$Q=\frac{F_\alpha-2.86F_\beta}
{\sqrt{F_\alpha^2+(2.86F_\beta)^2+s_\alpha^2+(2.86s_\beta)^2}}.$$

This retains weak and negative measurements without logarithms. It is a sensitivity, **not a physical attenuation estimator**. First-order uncertainty propagation is imperfect: 10,000 Gaussian perturbations around measured fluxes give Monte Carlo/delta-method variance ratios with 5th/50th/95th percentiles 0.79/1.05/2.01 among 25 weak or nondetected objects. Maximum local contrast bias is 0.21. These perturbations diagnose the approximation; they are not a latent-flux posterior or a coverage proof.

## Common-cohort results

Every before/after comparison below uses identical objects. The baseline fits width, supernova colour, linear redshift, a probabilistic host-mass step, local host colour, and quadratic width/colour terms. Its variance retains the released light-curve covariance, redshift uncertainty and fitted residual scatter. Nebular measurement variance is propagated through its coefficient. The age parameter is the public TITAN global mass-weighted age median, conditioned upon rather than re-inferred.

| Cohort | Hosts | Age slope before nebular term | Age slope after nebular term | After nebular term and separate SED controls |
|---|---:|---:|---:|---:|
| Hα and Hβ S/N>3 | 138 | −0.0265 ± 0.0212 | −0.0129 ± 0.0233 | +0.0045 ± 0.0241 |
| Hα and Hβ S/N>5 | 122 | −0.0257 ± 0.0234 | −0.0227 ± 0.0251 | −0.0019 ± 0.0260 |
| BPT star forming | 75 | −0.0276 ± 0.0319 | −0.0254 ± 0.0319 | −0.0215 ± 0.0319 |
| All valid signed lines, contrast Q | 163 | −0.0221 ± 0.0192 | −0.0241 ± 0.0198 | −0.0005 ± 0.0206 |

Slopes and conditional standard errors are in mag/Gyr. The final column includes continuous SED mass, stellar Av and metallicity as distinct controls. On the primary common cohort, the corresponding SED-controlled age coefficient **before** nebular adjustment is −0.0042 ± 0.0233 mag/Gyr.

For the primary 138 hosts, adding age after nebular control changes the five-fold heldout log predictive score by **−1.34**, with conditional host-bootstrap interval **[−4.15, +1.09]**. Adding nebular information to a model already containing age gives **+1.48 [−2.23, +5.60]**. Neither demonstrates predictive improvement. With separate SED controls, adding age after nebular control gives −1.93 [−3.71, −0.26]. Physical hosts stay together in each fold. Bootstrap intervals condition on the fitted folds and omit uncertainty from retraining the entire procedure.

The primary age coefficient's 95% normal interval after nebular adjustment is **[−0.0586, +0.0327] mag/Gyr**. Conditional normal-approximation power for an absolute 0.03-mag/Gyr coefficient is only **25%**, falling to **16%** in the BPT sample. After baseline and nebular adjustment, the primary age dispersion is 0.72 Gyr, below the median marginal age uncertainty of 0.94 Gyr. Absence of predictive improvement therefore does not establish absence of an age effect.

The before/after coefficients are correlated. We have not estimated their paired difference uncertainty; describing the point-estimate change as a measured fraction “absorbed by dust” would be unjustified. Stronger-line and BPT samples show much smaller changes, while the signed-contrast sample changes in the opposite direction.

## Error and selection sensitivities

The MPA-JHU authors warn that formal errors can understate repeat-observation scatter. Historical recommended multipliers are 2.473 for Hα, 1.882 for Hβ, 1.566 for [O III] and 2.039 for [N II]. We treated these as a sensitivity, not verified recalibration for every CAS measurement. On the same 138 hosts, the nebular-adjusted age coefficient becomes −0.0094 ± 0.0226 mag/Gyr, or +0.0081 ± 0.0237 with SED controls. Selecting using scaled-error S/N gives 116 Balmer and 67 BPT hosts; neither yields positive heldout age evidence. See the [author error documentation](https://wwwmpa.mpa-garching.mpg.de/SDSS/DR7/raw_data.html).

Changing the line-error correlation from −0.5 to +0.5 changes the primary nebular-adjusted age point estimate from −0.0120 to −0.0139 mag/Gyr. Independently propagating the marginal age error gives −0.0131 ± 0.0236. That sensitivity omits unavailable joint age–dust–mass–metallicity covariance and is not a complete errors-in-variables model.

[Groves, Brinchmann & Walcher (2011)](https://arxiv.org/abs/1109.2597) identified a roughly 0.35 Å Hβ equivalent-width bias in the historical DR7 catalogue. A conditional correction of this size, applied only to the 136 primary objects with emission-sign equivalent widths and refitting the same 136-object comparison, gives −0.0176 ± 0.0242 after nebular adjustment. We do not assume that a population-average historical correction exactly applies to this release or to each object.

## What this establishes

Nebular measurements add a separate dust-related comparison, but this sample does not discriminate “already corrected” from “additional correction required.” It is small, selected for spectroscopy and line strength, and often samples gas away from the supernova. Conditioning on colour, host properties and nebular emission may condition on physical mediators. Joint line/continuum uncertainties, a progenitor-to-environment model and selection modelling remain necessary. A physically interpretable luminosity correction also requires verified unblinded light-curve amplitudes.

All 88 full-data fits passed an independent SLSQP objective check; the largest −2-log-likelihood difference was 2.35×10⁻⁸. All 48 heldout model evaluations completed without rank or optimizer failures. Independent checks verify exact joins, physical-host folds, signed-flux retention, proxy gradients, flux-unit invariance and predictive-score arithmetic.

The [frozen design](../code/nebular_dust/design.json), [compact results](../results/nebular_dust/summary.json), [numerical validation](../results/nebular_dust/validation.json) and [reproduction instructions](../code/nebular_dust/README.md) retain assumptions, source hashes and sensitivities. Downloaded inputs and generated row-level tables remain outside Git.
