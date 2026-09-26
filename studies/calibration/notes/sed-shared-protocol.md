# Shared-systematic discrimination of constructed mean changes

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

Frozen before this comparison's conditional scores are read. Earlier validation
and nonlinear results are known, so this is an exploratory identification check.

Use all 1,020 validation objects, exact existing rows, nominal local nuisance
projection and the twelve-mode shared calibration-only posterior learned on
the original 43 discovery objects. Retain its original mode scales and all
shared discovery-to-validation covariance. Do not refit any spectral coefficient,
amplitude, nuisance prior, or calibration posterior on validation.

First construct two data-independent vectors from the resolved full nonlinear
fits to the **native noiseless target**: observer-filter change and the fixed
discovery-designed .05-mag-RMS spectral alternative. Project their saved fitted
mean changes into the existing nominal four-parameter complement, aligning
original epochs exactly. This retains the distinction between a nonlinear
noiseless response and a new fit to observed data. Record the norm lost under
this extra nominal projection. Also retain the original infinitesimal observer
direction as an exact reproduction control. Do not substitute observed-target
fitted response vectors into this design.

The governing predictive mean and covariance are `m=A mu_D` and
`V=I+A C_D A^T`, where the twelve-mode posterior comes only from discovery.
For each fixed vector v report I=v^T V^-1 v, M=v^T V^-1(y-m), G=M-I/2.
First report the outcome-free separation `(v_observer-v_SED)^T V^-1
(v_observer-v_SED)` and its counterpart with V=I. Compare shared and unshared
geometry; shared uncertainty cannot increase separation in this metric.
Verify Woodbury arithmetic by an independent thin-QR covariance solve.

These are probes obtained by adding each fixed vector to the **same** predictive
distribution. They are not jointly refitted physical models, integrated evidence
for spectral evolution, or posterior probabilities of calibration failure.
Both directions were chosen through known discovery geometry. Conditioning
the shared posterior anew under each physical alternative, covariance feedback,
retraining, selection and BBC remain separate. No cosmology is fitted here.

A design command writes and hashes the projected vectors and covariance inputs
without reading validation residuals or observed flux. A subsequent score
command first verifies those hashes, then reads the saved residual and posterior
mean. Preserve both phases and an exact original-observer score reconstruction.
