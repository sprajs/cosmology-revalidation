# Independent review of the broadband response pilot

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../../README.md) for current execution status.

Reviewed 2026-09-21 against the script and saved outputs hashed in `flux-response-review.json`. The reviewer did not modify `scripts/salt_dust_audit/flux_response.py`. The independent check is reproducible with:

```bash
phase2/env-official/bin/python scripts/salt_dust_audit/theory_flux_response_review.py
```

**The response calculation is internally correct for its stated synthetic, local, frozen-covariance problem.** It is not an independent inference of actual dust or an official SNANA/BBC rerun. The review covers all 64 saved objects and both weighting choices; independent wavelength integration additionally covers the lowest-, middle-, and highest-ranked selected redshifts.

## Mathematical and numerical checks

The forward derivative J is with respect to additive amplitude magnitude, x1, c, and observer-frame t0. The nuisance derivative G is with respect to host or MW E(B−V), or an additive photometric magnitude. For `C=L L^T`, weighted least squares gives `R=(L^-1 J)^+ L^-1 G`. The SVD implementation, parameter covariance `V diag(s^-2) V^T`, and residual projection in the reviewed script implement this correctly.

An independently evaluated pivoted-QR least-squares solution reproduces all saved response entries to a maximum absolute difference of `2.884e-14`. The inverse-covariance identity residual is `2.443e-13`; the largest difference in standardized-magnitude response is `1.199e-14`. Every four-parameter tangent has full numerical rank. The whitened residual is orthogonal to the fitted tangent space to relative norm `6.702e-16`. Direct evaluation of sncosmo's interpolated quadratic surface variance **before its negative-value floor** finds zero negative values at all 2,581 pilot epochs. This check is specific to sncosmo and does not instrument SNANA's path.

Adding the same photometric magnitude to every filter is exactly degenerate with changing the amplitude magnitude. Summed calibration columns recover `(dm,dx1,dc,dt0)=(1,0,0,0)` to `1.415e-9`. No fitted physical information should be inferred from this gray mode. There are no nonpositive integrated model-flux values among the saved pilot epochs, although the separate source-grid audit finds some negative monochromatic mean-surface values elsewhere.

The photon-integral identity

\[
\frac{1}{f_X}\frac{\partial f_X}{\partial E}
=-0.4\ln(10)\,
\frac{\int d\lambda\,\lambda T_X(\lambda) F_\lambda k(\lambda)}
     {\int d\lambda\,\lambda T_X(\lambda) F_\lambda}
\]

was evaluated independently with dense trapezoidal wavelength integration. Host k is evaluated at `lambda_obs/(1+z)`, and MW k at `lambda_obs`. For five extinction columns across the three checked objects, it agrees with the saved finite differences to a maximum relative difference of `1.745e-6`. This independently verifies the sign, wavelength frame, and E versus A_V unit conversion; it still uses the same `extinction.fitzpatrick99` implementation for the law itself.

The finite-step checks show nonlinear-minus-linear standardized shifts ranging from `−8.78e-7` to `−4.95e-9` mag at E=0.001, `−8.77e-5` to `−4.92e-7` mag at E=0.01, and `−0.002177` to `−0.0000114` mag at E=0.05. These are selected deterministic injection checks, not uncertainty estimates or a statement about the full sample's finite corrections. They demonstrate why a derivative per unit E must not be extrapolated automatically to an E=1 physical change.

## Sign and normalization interpretation

- Positive host/MW nuisance here means **extra attenuation in the synthetic target**, with the fit's dust configuration held fixed. Increasing assumed MW attenuation in the model while keeping observed flux fixed has the opposite local response.
- The calibration columns mean a positive additive **photometric magnitude**, so flux decreases. A positive zero point in the equation `FLUXCAL=counts*10**(0.4*(27.5−ZP))` also decreases FLUXCAL at fixed counts, but conventions that directly increase calibrated flux have the opposite sign. The explicit dimming definition is authoritative; a bare “zero-point error” label is insufficient.
- `dm_amplitude_mag` changes `x0` by `exp[−0.4 ln(10) dm]`. Its sign is therefore the usual magnitude sign. The standardized response `[1,0.16087,−3.11780,0] R` uses the published fixed alpha/beta, and is the change in **pre-BBC standardized magnitude**. It does not include a response of alpha, beta, host correction, population fitting, or bias correction.
- Applying `MAG_OFFSET=0.27` and each KCOR `Primary Mag` before weighting is appropriate: SNANA `genmag_SEDtools.c` adds the primary magnitude to its model zero point, which produces the corresponding `10**(−0.4*PrimaryMag)` FLUXCAL factor. Constants cancel from simple relative derivatives but matter for absolute model-flux covariance relative to reported flux errors.
- The data/model median flux-ratio diagnostic is only a sanity check. Across objects it ranges from approximately 0.872 to 1.080, with median 0.994; this does not establish point-by-point SNANA agreement. Existing differences in interpolation, error floors, and covariance treatment remain explicit.

## Conditions that must accompany scientific interpretation

The source color model, reference coordinates, accepted epochs, spectroscopic redshift, and foreground map value are fixed. This is a tangent to **synthetic model data centered on published coordinates**, not the derivative of the actual observed-data optimum. Nonzero observed residuals, priors, iterative rejection, and parameter-dependent covariance can change the latter. Frozen-covariance weighting also omits derivatives of both the quadratic form and log determinant in a full varying-covariance likelihood.

`measurement_only` uses a diagonal measurement covariance. Any unprovided correlated photometric/background uncertainties are consequently outside it. `measurement_plus_model` includes sncosmo's SALT3 per-object model covariance, evaluated at the reference point; it does not include full shared training/calibration covariance, a complete intrinsic population model, or SNANA's separate MW covariance. No cross-SN covariance should be inferred from a block stored for one object.

The absorption fraction compares a nuisance flux direction to the span of **all** fitted light-curve directions, including freely fitted gray amplitude. Large absorption establishes local difficulty distinguishing a perturbation from changes in SALT coordinates. It does not by itself measure intrinsic-versus-dust identifiability, prove there is dust, or bound the cosmological bias. Inspect the surviving residual modes, induced standardized magnitudes, and redshift-dependent population distribution together.

The four host-E columns at different R_V values are alternative local dust-law hypotheses. Treating them as four independent physical nuisance parameters would require an explicit physical model. They are not four dust measurements. A nonzero R_V derivative around an extincted SN would require a nonzero host E reference and its own perturbation; R_V has no first-order effect at E=0.

Finally, `R S R^T` is an uncertainty covariance only when S has a justified definition, amplitude, and correlation structure. These response matrices intentionally do not supply that missing physical prior. Joint shared modes must remain shared across SNe. Given these restrictions, the matrices are useful conditional diagnostics and inputs to a later explicitly specified sensitivity analysis.
