# Smooth rest-frame SED identification geometry

**The small rest-frame comparator leaves a serious identification gap.** A
prespecified smooth spectral/phase family with a smooth redshift-dependent
mean can closely imitate the observer-band response on both cohorts. A
constructed 0.05-mag RMS SED change chosen using discovery **design geometry
only** reproduces 98.12% of the frozen observer vector's validation information
after finite broadband SED reintegration. This is an existence example inside
an empirical model, not evidence for actual evolution, an inferred correction,
or a calibrated physical population prior. The finite SED calculation still
uses a **linear SALT nuisance projection**; nonlinear nuisance refitting,
retraining, classification, selection and BBC have not been performed.

## Frozen design, before evaluating the new geometry

This bounded calculation asks whether shared rest-frame spectral/phase mean
changes can imitate observer-band offsets after fitting each object's
amplitude, stretch, colour and peak time. It reads no observed flux or residual
array. All three observer contrast directions are primary; the already frozen
discovery observer vector is a secondary diagnostic. This is a local response
and information calculation, not another residual significance test.

Write the multiplicative perturbation as
`F_new(p,lambda,z) = F_nominal(p,lambda,z) exp[-kappa Delta_m(p,lambda,z)]`,
where `kappa=0.4 ln(10)`, wavelength is rest Angstrom, phase is rest days, and
`Delta_m` and its coefficients have units mag. Positive Delta_m dims the SED.
Define `x=2 ln(lambda/2000)/ln(11000/2000)-1`, `u=tanh(p/20)`, and
`v=ln((1+z)/1.5)/ln(1.5)`. The latter is a smooth selected-population drift
coordinate, not a physical evolution law or an inferred parent population.

The **small basis**, fixed before the calculation, consists of `P1(x)`,
`P2(x)`, `P3(x)`, and `P1(x)u`. The **broad sensitivity** consists of
`Pn(x)` for n=1..8, `Pn(x)u` for n=0..8, and `Pn(x)P2(u)` for n=0..8:
26 columns. The small basis is nested inside the broad basis. For each, compare
a shared mean and mean plus its product with v: 4, 8, 26 and 52 coefficients.
No knots or orders are selected using the frozen observer vector or residuals.
The broad basis tests smooth representation flexibility; it is not a calibrated
prior over physical supernova populations. Pure phase-independent gray terms
are excluded because each object's fitted amplitude absorbs them exactly.

Derivatives must use the actual released SALT3 DES5YR spectral surface,
redshift, phase, Milky Way attenuation and DES photon-weighted passbands.
Use the native integration grid of at most 5 observer Angstrom, with a
2.5-Angstrom convergence check on every fiftieth object. Reproduce the archived
native integrated mean. The native fractional spectral response is mapped onto
the official mean, which enforces the exact gray-amplitude null. This is an
explicit independent spectral derivative approximation at official fitted
coordinates, using the already archived native x1/c/t0 tangent and exact
SNANA frozen covariance; it is not a native SNANA SED derivative export.

For each object compute `S=Q^T L^-1 G_SED` and
`O=Q^T L^-1 G_observer`, where `LL^T=C` and Q is orthogonal to `L^-1 J`.
Stack separately over the 43 discovery and 1,020 validation objects. Report
three canonical correlations of the column spaces, singular values/rank and
the remaining observer information eigenvalues `1-rho^2`. These describe
best possible linear imitation within the specified family, not actual
imitation in the data. Fit the response-to-response mapping on discovery and
transport it unchanged to validation as an additional geometry check.

Amplitude accounting uses the RMS of Delta_m on a fixed product quadrature:
uniform x over [-1,1], uniform phase over [-15,45] days, and uniform v over
[-1,1] for drift families. Observer contrasts use an orthonormal, zero-sum
griz basis; a coefficient vector of norm .02 therefore has griz Euclidean
norm .02 mag. Report SED RMS and maximum magnitude for each least-squares
canonical match. Fixed sensitivity penalties are `(SED_RMS/sigma)^2`, with
sigma=.01,.02,.05,.10 mag; they are explicit regularization scales, not
measured uncertainty. Report deterministic mismatch norms, not p-values.

The basis is physically defined as a change in the emitted SED and can
represent shared training error or a smooth change in selected-population mean.
It is not identified with dust, metallicity, host mass, age or explosion
physics. Empirical SALT dimensionality and spectral population dependence
motivate admitting such alternatives, but do not fix our basis or amplitudes:
[SALT3](https://arxiv.org/abs/2104.07795),
[SALT3 host spectral model](https://arxiv.org/abs/2209.05584), and
[SALT3+ dimensionality study](https://arxiv.org/abs/2502.09713).

Artifacts and results are appended after execution. No observer residual
coefficient will be re-estimated, and no fit, covariance, clipping, classifier,
selection or BBC product will be changed.

## Executed geometry on both cohorts

All 43 discovery objects / 1,649 epochs and 1,020 validation objects / 39,606
epochs are retained, with their original repeated-row multiplicities. The
residual-space dimensions are 1,477 and 35,526. Exact frozen SNANA covariance
and the archived nuisance Jacobian are used; only its amplitude column is
set to the official mean's exact derivative. Numerical ranks are respectively
4, 8, 26 and 52 in both cohorts.

The three canonical correlations below compare the **whole observer
contrast space** with each rest-frame space. They are not correlations with
observed residuals. A value near one means little separation along that
particular pair of response directions, without restricting SED amplitude.

| Rest-frame family | Columns | Discovery correlations | Validation correlations |
|---|---:|---|---|
| Original two dust laws + phase-colour comparator | 3 | .97046, .45470, .32065 | .97276, .47759, .29931 |
| Small shared mean | 4 | .99719, .76633, .34242 | .99721, .76219, .36356 |
| Small mean + redshift drift | 8 | .99872, .98838, .78060 | .99857, .98746, .76932 |
| Broad shared mean | 26 | .99889, .93158, .47755 | .99876, .91195, .51197 |
| Broad mean + redshift drift | 52 | .99967, .99827, .96252 | .99950, .99723, .95267 |

The validation remaining-information eigenvalues are .00558/.41906/.86783
for the small shared family and .00100/.00554/.09241 for broad plus drift.
These are fractions relative to each corresponding observer direction's
information, not common fractions of the total dataset. Broader shared mean
alone leaves .00247/.16835/.73789. Permitting population drift changes the
identification problem much more than increasing a shared spectral mean's
dimension alone.

This has a simple analytic explanation. A smooth observer spectral response
`h(log lambda_obs)` becomes `h(log lambda_rest + log(1+z))`. A quadratic h
contains a rest-wavelength quadratic, a wavelength-times-redshift term and
a redshift-only gray term; the last is absorbed by each object's amplitude.
Thus smooth rest-frame drift and smooth observer response can be equivalent
even before introducing arbitrary rest-frame flexibility. Actual griz band
offsets are integrated step-like responses, so the numerical broadband test
is still necessary; this argument does not replace the integrations.

The unbounded broad family has weakly measured spectral directions. For
example, its unregularized validation best match to the frozen vector has
77.4-mag full-grid RMS without drift and 21.8 mag with drift. Those solutions
are outside a credible local physical regime. Canonical overlap by itself
therefore cannot establish a plausible finite explanation.

## Fixed-amplitude frontier and finite SED check

The reporting amendment converts the four already specified amplitude scales
into hard full-grid RMS bounds. It also adds supported-wavelength RMS and
finite integrations; it changes no basis, cohort or target. The original
output and source are preserved. No observed-flux, residual, chi-square, or
residual-score array was accessed; the manifest records all archive keys read.

For the secondary frozen observer vector, the table reports the fraction of
its validation information that remains unmatched. These are deterministic
responses at official fitted coordinates, not observed residual scores.

| Full-grid SED RMS limit, mag | Small + drift, 8 columns | Broad + drift, 52 columns |
|---:|---:|---:|
| .01 | .8161 | .3705 |
| .02 | .6791 | .1642 |
| .05 | .4666 | .01729 |
| .10 | .2946 | .00365 |

This table uses validation **geometry** to obtain the nearest admissible
response. A more stringent transfer calculation learns the response mapping
from discovery geometry and transports its coefficients unchanged. At the
.05-mag RMS bound, the latter leaves .01897 of validation information
unmatched in the linear calculation, versus .01884 after the exact finite
SED multiplier is integrated. Its supported-wavelength throughput RMS on
validation is .04263 mag, full-grid maximum is .1921 mag, and maximum on the
actual evaluated wavelength/phase support is .14762 mag. The finite response
differs from its linear prediction by .37566 in whitened norm, or 2.129% of
the linear response norm. The unmatched finite response norm is about 2.49;
no p-value is assigned to this constructed composite alternative.

For clarity, the 98.12% value is `1-||O c-S theta||^2/||O c||^2` (using the
finite spectral response for S theta), and is not a coefficient of
determination for actual observations. The spectral coefficients minimize
response-to-response discrepancy, not an outcome likelihood. The full
observer space is primary because the previously frozen c itself encodes
discovery outcomes. The 52-column family and .05-mag budget remain arbitrary
declared sensitivity choices; no independent observation supplies their
population probability or prior mass.

“Supported-wavelength RMS” gives each accepted epoch equal weight and uses
normalized positive `lambda*T(lambda)` weights inside its actual passband.
It is neither detected-photon nor signal-to-noise weighting. It prevents
assessing amplitude solely at unobserved rest wavelengths, while preserving
the full-grid RMS/max diagnostics. Neither norm justifies the physical
population model. The passband integration itself retains **signed** native
SALT spectral flux throughout; it never clips negative spectral cells.

## Constructed local distance ambiguity: a further gate is required

Using the same verified 255 lowest / 255 highest zHEL membership as the
existing twelve-mode distance analysis, define compensation for a hypothetical
mean change by `delta_p=-(L^-1 J)^+ L^-1 delta_f`, and standardize with
`(1,.16087,-3.1178,0)`. The frozen observer vector gives a high-minus-low
fixed-reference pre-BBC contrast of **+.060683 mag**. The discovery-designed
.05-mag broad-drift SED example gives **-.186324 mag** in the linear response,
or **-.179120 mag** with finite SED integration and still-linear parameter
compensation. The latter differs from the observer example by -.239803 mag.

These are constructed tangent-order responses, not distance estimates or
measured biases. Similar projected residual patterns can absorb very different
colour/amplitude shifts. The approximately .24-mag separation is large enough
that full nonlinear nuisance refits are a necessary next gate before treating
it as a quantitative distance ambiguity beyond tangent order. Finite SED
integration does not satisfy that gate. Refitting alpha/beta, spectral training,
selection, classifier/BBC and external population support are further distinct
requirements. Arbitrary redshift-dependent gray luminosity evolution remains
unconstrained by the entire projected test.

## Verification, model support and reproducibility

Independent QR canonical correlations agree with the SVD calculation to
1.26e-13; the constrained amplitude solutions satisfy their primal KKT
equations to relative 7.12e-15. Native integrated means reproduce the archive
to relative 1.10e-13. Observer response matrices reproduce their archived
values to 8.00e-15, and the gray projected norm is below 5.53e-14. The 5 to
2.5-Angstrom integration check on every fiftieth object changes the projected
broad response by at most 2.10e-6 relatively. Independently integrating the
linear spectral multiplier reproduces the matrix response to 3.72e-15.

Native SALT has negative spectral cells carrying more than 1e-6 of absolute
band photon weight in 16 discovery epochs and 782 validation epochs. Ten
validation **native** integrated means are negative; nine corresponding
official means are exactly zero and one is +.000053. They are early g/r/i
epochs at rest phase -14.81 to -12.83 days. Their maximum absolute native
mean divided by quoted error is .632; their exact ledger is saved. Negative
spectral cells can occur even when the band-integrated mean is positive.
These are empirical-model support warnings, not physical negative photons.
No rows or spectral values were removed or clamped. A separate native-J /
row-removal audit was assigned to the simulation reviewer after this result.

Using the raw native spectral derivatives instead of mapping their fractional
response to the official mean changes the broad matrix norm by .0262% on
discovery and .0989% on validation. Validation broad-drift canonical
correlations become .99952/.99705/.95338, versus .99950/.99723/.95267.
This sensitivity does not remove the geometric ambiguity; it does not certify
the archived native nuisance tangent at the problematic early epochs.

The covariance here is the existing per-object exact frozen covariance. The
calculation does not integrate the twelve shared calibration modes again,
reproduce the 8.68/8.30 log-score increment, validate the full selection
likelihood, or rank mechanisms. A fresh SED-versus-observer predictive
comparison needs independent spectral/population restrictions and a frozen
selection-aware protocol. It must preserve the broad-family sensitivity
rather than treating the old three-column comparator as all rest-frame physics.

Main [script](../../scripts/research_2026_09_26/sed_identification.py),
[numerical result](../../runs/research_2026_09_26/sed_identification/result.json),
[design matrices](../../runs/research_2026_09_26/sed_identification/design-matrices.npz),
[finite responses](../../runs/research_2026_09_26/sed_identification/finite-sed-responses.npz),
[distance response ledger](../../runs/research_2026_09_26/sed_identification/constructed-distance-responses.csv),
[nonpositive model epochs](../../runs/research_2026_09_26/sed_identification/nonpositive-model-epochs.csv),
[independent check](../../runs/research_2026_09_26/sed_identification/independent-geometry-check.json),
and [archive-key/input manifest](../../runs/research_2026_09_26/sed_identification/manifest.json)
retain the calculation. The initial rank-reporting field accidentally counted
rows rather than columns; its preserved output is explicitly labelled, and
the corrected report uses 4/8/26/52. No numerical SVD result changed.

Reproduce with the existing environment and preserved output directory:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python scripts/research_2026_09_26/sed_identification.py --compute --amplitude-amendment
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/env-official/bin/python runs/research_2026_09_26/sed_identification/independent_geometry_check.py
```
