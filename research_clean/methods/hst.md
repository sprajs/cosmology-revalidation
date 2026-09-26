# HST signed pixels and quoted variance

## Scope

These workflows use eight current WFC3/IR F160W calibrated FLTs for one search visit and one template visit, plus the official pixel-area map and the header-named flat references. They study signed repeat photometry and variance accounting. They do not reproduce the historical RAISIN image reduction or measure a cosmological correction.

Exposure order is fixed by `EXPSTART` and filename. Ordinals 1 and 3 define source masks; 2 and 4 provide held-out aperture brightness sums. DQ and finite-value/error eligibility enter support, but held-out brightness amplitudes do not select positions.

## Geometry and extraction

`hst-geometry` rebuilds a 6-arcsecond tangent-plane lattice, excludes the SN neighbourhood, constructs design-exposure source masks, and integrates fractional-pixel aperture/annulus weights. The aperture radius is 0.4 arcsec and the annulus spans 1.2–2.0 arcsec. WCS curvature, doubled-resolution closure, constant-image cancellation and support checks are retained. Pixel coordinates use the original one-based convention consistently.

The strict branch requires every contributing pixel to be valid in every exposure. Its support is zero: no strict estimate is produced. The separately specified masked branch requires at least 90% PAM-weighted aperture coverage and 75% native annulus coverage in all eight exposures. It retains 170 positions across 16 spatial tiles from a 385-position geometry grid. Rebuilt operator arrays must equal the frozen arrays exactly before that replay is declared successful.

The original 120-second geometry benchmark cap is retained, together with the bounded working-array design. A slow machine can stop before completion; partial geometry is never accepted as complete. The bundled frozen operators remain available for the repeat calculation. This is a resource boundary, not a reason to change the sample.

For retained pixels, define aperture weights `a`, annulus weights `b`, pixel-area map `P`, validity `v` and relative PHOTFLAM scale `k`. The signed linear weights are

`w = k * [v*a*P - sum(v*a*P)/sum(v*b) * v*b]`.

The extracted flux is `w^T SCI` and the diagonal reference variance is `sum(w² ERR²)`. No positive-flux selection or full-aperture renormalization is applied. `hst-repeat` recomputes these sums from the FITS pixels, preserving the frozen support.

## Repeat statistics

Within a visit, form the ordinal-2 minus ordinal-4 flux difference and the sum of their quoted variances. Report signed standardized means, the uncentered second-moment ratio and the weighted-intercept-centred variance ratio. Spatial-tile resampling is descriptive and requires at least 12 occupied tiles with at least three positions each; 2,000 replicates are used by default.

Expect uncentered `sum(difference²)/sum(quoted variance)` of approximately `0.7082533` in the search visit and `0.7287246` in the template visit. These values do not identify which part of the quoted error, detector covariance, reference processing or aperture model explains the difference. No error array is rescaled.

## Flat-reference variance

`hst-flat` propagates the exact header-named PFLT and DFLT reference ERR planes on the same support. Unit LTM and the five-pixel LTV difference determine the reference-to-science coordinate mapping. The diagonal flat contribution is

`sum[(w*SCI)² * ((ERR_P/P)² + (ERR_D/D)²)]`.

It also calculates covariance from the same flat-reference detector pixels shared by two exposures. Sky-position coincidence alone does not imply equal detector pixels. Under independent reference pixels, expect flat marginal fractions of roughly `0.00005258` and `0.00005141` of quoted pair variance. The corresponding fractions removed by same-pixel shared-flat covariance are about `0.00000340` and `0.00000274`.

Those particular supplied terms are much smaller than the observed repeat-variance difference. Unknown spatial/temporal detector covariance and other reference terms remain unresolved. The separate [dark controls](hst-dark.md) now reproduce raw paired moments and a frozen native dark slope comparison, but those total moments mix several detector processes. Neither the dark controls nor a source-formula comparison under independent reads justifies changing ERR without identifying the actual covariance relevant to science photometry.

```bash
uv run --frozen python research.py run hst-geometry --name hst-geometry
uv run --frozen python research.py run hst-repeat --name hst-repeat
uv run --frozen python research.py run hst-flat --name hst-flat
```
