# Conditional Gaussian check of nonlinear nuisance fitting

Could nonlinear fitting alone create a residual along the frozen observer-band
direction? We tested that restricted mechanism at the same six previously fixed
redshift ranks used for the calibration-response check. This is separate from
the native P21/G10 simulation and selection test.

For each object, generate 128 light curves from its independent sncosmo SALT3 mean plus
Gaussian noise with its exact exported SNANA covariance. Hold the accepted
epochs and covariance fixed. Refit brightness, width, colour and time; compare
the observer score at the new fitted tangent with the generating tangent's
linear projection, using identical random noise. The observer coefficients
remain the original 43-object values. No raw observed flux is used as simulated
truth, and no coefficient is learned from the simulated residuals.

The generating and refitted means use the **same sncosmo implementation**;
they are not instrumented SNANA mean evaluations. Their largest fractional
mean difference from the saved SNANA reference is .00059684. Internal mean
consistency is appropriate for this isolated curvature check, but is not an
independent verification of the SNANA optimizer. The original frozen protocol's
word “native” referred loosely to the existing SALT implementation; this
paragraph supplies the precise implementation distinction without modifying
the preserved numerical artifacts.

All **768** fits succeeded, with no boundary hits. All **48** prescribed
second-start checks succeeded and agree in objective within 8.22e-11.
Pooled nonlinear chi-square per degree of freedom is .99697, compared with
.99743 for the paired linear control.

An [independent Astra calculation](../../runs/research_2026_09_26/astra_design/gaussian-control-independent.json)
reproduces the means and Monte Carlo errors from all saved draws within 1.8e-16.

| Six-object summed quantity | Nonlinear minus linear mean | Monte Carlo SE |
|---|---:|---:|
| Matched product along the frozen vector | -.023499 | .004507 |
| Squared prediction norm (information) | +.002639 | .004331 |
| Fixed-vector log-score gain | -.024818 | .004955 |

The generating information is 2.84390. Dividing the paired matched-product
change by that information gives a descriptive amplitude shift of about -.0083.
This small effect points opposite the positive real-data pattern in this
particular control. The unpaired nonlinear matched mean is .1060±.1516 Monte
Carlo SE, consistent with zero. The paired comparison is more precise because
the two calculations share each noise realization.

These Monte Carlo errors describe the finite experiment, not uncertainty in a
survey bias. This does not validate the real covariance, a physical intrinsic
scatter law, native SNANA priors/optimizer, model-dependent covariance,
clipping, detection or classification. The six fixed cadences are not a
weighted representation of all 1,020 validation objects. The broader native
simulation control remains necessary.

The [script](../../scripts/research_2026_09_26/gaussian_fit_control.py) runs with
`phase2/env-official/bin/python`, which supplies the already installed sncosmo
dependencies. An initial launch in `.venv` failed at import before any numerical
output and is [retained](../../runs/research_2026_09_26/gaussian-fit-wrong-environment.log).
The successful [result](../../runs/research_2026_09_26/gaussian_fit_control/result.json),
[all draws](../../runs/research_2026_09_26/gaussian_fit_control/draws.csv), protocol,
executed source and hashes are preserved. The output directory cannot be
overwritten by the script.
