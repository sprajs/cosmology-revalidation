# Constructed alternatives under the existing shared calibration model

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The twelve-mode predictive covariance retains some information distinguishing
the two constructed mean changes. Under the same discovery-conditioned null,
the full-amplitude observer probe scores 6.18 nats above the spectral probe,
but **both lower the predictive score relative to adding neither**. This is
a conditional probe of fixed offsets, not evidence for a physical calibration
correction or against general spectral evolution.

The [protocol](sed-shared-protocol.md) was frozen before these scores. It uses
all 1,020 validation objects and the saved nonlinear fits to a **noiseless**
target. No observed-target fitted means are used as candidate predictions.
The existing four-parameter nuisance complement projects these response vectors
into the same coordinate system as the twelve released shared modes. The shared
posterior mean and covariance are the calibration-only values learned on the
43 discovery objects. Nothing is optimized on validation.

For predictive mean m and covariance V, define the fixed probe gain
`G=v^T V^-1(y-m) - (v^T V^-1 v)/2`. The inherited mode scales, nominal per-object
covariance and shared discovery-to-validation correlation are unchanged.

| Fixed probe | Information I | Centered matched product M | Predictive gain G, nats |
|---|---:|---:|---:|
| Nonlinear noiseless observer response | 47.192095 | 22.292612 | -1.303435 |
| Nonlinear noiseless spectral response | 50.779827 | 17.910006 | -7.479908 |
| Original infinitesimal observer control | 46.602541 | 22.076688 | -1.224582 |

Before shared uncertainty, the squared norm of the difference between the two
nonlinear responses is 6.298645. Under V inverse it is 5.633516. Each response's
own squared norm falls much more: observer 332.814212 to 47.192095, spectral
315.442203 to 50.779827. The shared calibration family covers much of their
common response while covering less of their difference. Thus the earlier
1.90% unshared squared mismatch alone does not describe identification after
shared uncertainty is admitted.

The extra nominal projection removes squared norms .076418 and .070207 from
the noiseless nonlinear responses, and .003052 from their difference. These
small but retained terms distinguish this calculation from a full nonlinear
shared likelihood. It is not the earlier observed-flux maximized-score
comparison, nor a replacement for the twelve-mode integrated comparison.

The original observer score is reproduced within 1e-9. A separate thin-QR
covariance solve verifies the Woodbury calculation within 1e-10; the shared
covariance can only decrease the information metric, also checked numerically.
Every object's saved input and prediction hash, epoch order, covariance and
convergence flag is checked. The design vectors and their manifest were written
before the residual was read by the score phase.

A separate Astra calculation uses an SVD covariance solve, checks 4,088 hashes,
and reproduces scores and Gram matrices within 1.1e-12. It finds conditional
vector cosine .94313 and difference information equal to 11.94% of the observer
information. Its [independent record](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/astra_design/sed-shared-independent.json)
retains the inputs and method.

The fixed full-amplitude offsets are deliberately not rescaled to improve these
scores. A proper comparison of physical models would recondition calibration
coordinates jointly on discovery under each specified alternative, supply an
independently justified spectral/population prior, and propagate covariance,
training and selection. The result leaves those questions open. It demonstrates
that accounting for already included corrections affects interpretation; it
does not identify a unique correction from the flux residual.

Artifacts: [design and geometry](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_shared_probe/design.json),
[conditional scores](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/sed_shared_probe/score.json),
[source](../code/sed_shared_probe.py).
