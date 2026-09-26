# CSP source lineage and filter response

## Literal photometry lineage

`csp-lineage` reads the original CSP DR3 archive and its SNooPy tables directly. Decimal strings preserve the published magnitude, uncertainty and time literals. Multisets retain repeated rows; joins use object identity, band and time rather than brightness.

The raw physical labels `Jdw` and `Hdw` become `J` and `H` in the intermediate tables. Expect 5,491 raw NIR rows and exactly the same 5,491 intermediate rows after those relabellings, with no unmatched archive values. This establishes a label transformation, not an empirical passband correction.

The release comparison then applies the stated converter equations:

`flux=10^[-0.4*(mag-27.5)]`

`error=flux*(10^[0.4*magerr]-1)`.

Values are compared at the exported five-decimal scientific formatting. The workflow retains all LISTed objects and all mismatches. Expect 73 of 76 objects to match that complete converter identity. Of 520 raw WIRC-J rows, 303 released rows have an unambiguous WIRC physical label. A missing/ambiguous physical label is not filled by selecting the closest magnitude.

## Conditional photon response

`csp-passbands` uses the archived KCOR SN SED, its BD17 primary spectrum and the pinned CSP physical throughput tables. It calculates photon counts as `integral lambda*T(lambda)*f_lambda(lambda/(1+z))/(1+z) d_lambda`.

The natural-reference comparison sets equal differential BD17 magnitude, with no fitted offset. Phase and redshift contrasts cancel static zero-point conventions. Signed tabulated transmission tails are retained. The grid contains seven phases, five redshifts and four filter pairs, producing 140 rows by default.

Piecewise-linear spectrum/throughput on their union knots makes the photon integrand cubic on each interval. Two-point Gaussian quadrature is compared against four-point quadrature and an independent 0.5-Angstrom trapezoid. Required maximum magnitude gaps are `1e-10` and `1e-4` respectively. Out-of-support wavelength or phase requests stop the calculation.

```bash
uv run --frozen python research.py run csp-lineage --name csp-lineage
uv run --frozen python research.py run csp-passbands --name csp-passbands
```

The output quantifies response conditional on the chosen spectral model, reference spectrum and filters. It does not establish the exact historical calibration execution, observed SN spectra, trained-model independence or a corrected distance. The native-fit passband experiments and their branch ambiguity are not silently replaced by this spectral integral.
