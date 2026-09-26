# Smooth rest-frame SED identification geometry

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
