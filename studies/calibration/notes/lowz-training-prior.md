# A low-redshift training prior is not exactly independent of expansion shape

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The [existing BayeSN source review](../../dust/notes/standardization-dust.md)
identifies an important alternative-model assumption. G26 uses an informative
cosmological distance prior below z=.08 and a flat distance prior above it.
The rationale is that the low-redshift distances are approximately sensitive
only to H0. [G26 methodology, section 2.3](https://arxiv.org/html/2606.19429v1#S2.SS3).

We tested that approximation analytically, without fitting supernova data.
In flat geometry the dimensionless luminosity distance is

`D_L(z)=(1+z) integral_0^z dz'/E(z')`.

For constant deceleration q, `E(z)=(1+z)^(1+q)`, so
`D_L=(1+z)[1-(1+z)^(-q)]/q`; the q=0 limit is
`(1+z) ln(1+z)`. H0 is an overall multiplicative scale and drops out of the
distance ratios. A single absolute-luminosity/intercept adjustment cannot
remove their remaining dependence on redshift.

The following changes compare two declared, nonaccelerating kinematic histories
to flat LCDM with Omega_m=.3. The second column uses the same H0; the third
subtracts the difference at z=.02 to remove one common gray intercept.

| Alternative at z=.079 | Delta distance modulus (mag) | Delta after anchoring at z=.02 (mag) |
|---|---:|---:|
| q=0 | -.043521 | -.031821 |
| q=+.5 | -.084667 | -.062225 |

The archived BayeSN code defaults to Omega_m=.28. Using that reference gives
anchored changes -.033739 and -.064144 mag. The .28-to-.30 reference change
itself is only .001919 mag in this anchored contrast. Thus the approximation
can be adequate within a narrow familiar cosmological family while failing
as an exact argument across substantially different expansion histories.

For scale, the code's 150 km/s peculiar-velocity term alone corresponds to
.013743 mag at z=.079. **This is not the full prior width and these ratios are
not significance levels.** Spectroscopic redshift errors and latent intrinsic
scatter matter; the inspected archived training code combines `sigma0` with
the redshift term in its `Ds` prior. Its default source does not identify the
exact G26 training execution or hybrid-prior implementation.

This calculation establishes a finite prior-shape dependence, not an observed
G26 calibration error. No actual training-sample weighting, SED retraining,
host/dust inference, held-out prediction or corrected cosmology was performed
here. Before treating G26 as an expansion-independent discriminator between
accelerating and nonaccelerating models, a useful sensitivity is to repeat
training with alternative low-redshift distance shapes or to infer a flexible
low-redshift distance relation jointly, then compare on unchanged held-out
photometry with the same calibration and selection treatment.

[Executable calculation](../../light_curve_fitting/code/lowz_distance_prior.py),
[numeric results](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/lowz_distance_prior/result.json)
and [full fixed grid](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/lowz_distance_prior/distance-shape.csv)
are preserved with hashes and exact archived code lines. Independent Astropy
distances agree with direct quadrature to relative error 2.1e-13; the q=.5
closed form is independently checked. These are kinematic examples, not
estimates of the universe's actual expansion.
