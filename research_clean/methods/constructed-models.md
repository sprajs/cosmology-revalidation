# Explicit dust and selection models

## Dust geometry

`dust` distinguishes absorption along a point-source sightline from spatially averaged galaxy attenuation. In the homogeneous absorption-only slab, transmission is `(1-exp(-tau))/tau` and attenuation is `-2.5 log10(transmission)`. The code implements the transparent limit and a small-tau expansion, checked against numerical integration.

Supernovae are placed uniformly in depth; the illustrative detection rule is `A_B < 0.5 mag`. The output evaluates three microscopic R_V values and six optical depths. It demonstrates that integrated-light attenuation R_V can differ from individual SN extinction R_V. Scattering, actual host geometry, luminosity selection, cadence and redshift dependence are absent by construction.

It also evaluates the exact reparameterization `E' = E+k*age`, `c_int' = c_int-k*age`, `b' = b-(R_B-beta_int)*k`. Observed colour and brightness remain unchanged although the explicit age coefficient changes. Expected numerical identity errors are below `1e-12`. This identifies a model degeneracy; it does not show that any particular physical decomposition describes the observations.

## Correlated sign selection

`lib/selection.py` provides Gaussian likelihoods for `C=diag(d)+v v^T`. The latent scalar factor reduces orthant probabilities to one-dimensional integration. The code distinguishes:

- all signed measurements;
- values conditional on all retained measurements being positive;
- observed positive values plus known negative signs for the other eligible epochs.

The censored likelihood includes the probability of the missing signs. It must not receive an additional truncation normalizer. Its full eligible schedule and omitted signs are required inputs, not inferred automatically from an archive of positive detections.

For multirow correlated orthants the validated numerical domain requires `sum(v²/d) <= 0.5` and latent log-density mode within `[-4,4]`. Scalar and independent cases are analytic. Unsupported factor strength or tails raise errors; agreement between two quadrature orders alone is not accepted as proof.

`sign-selection` checks the kernel against adaptive quadrature and a dense multivariate Gaussian density, then runs a clean seeded recovery experiment: 128 draws, 24 fixed epochs, amplitude 0.4, unit independent variance and loading 0.06. It compares full signed, naive positive-only and correctly censored fits on the same simulated draws. This compact generator is fully specified in the code; it is not presented as an exact replay of all historical simulation designs. Bias, Monte Carlo error and optimizer boundaries are reported.

```bash
uv run --frozen python research.py run dust --name dust
uv run --frozen python research.py run sign-selection --name sign-selection
```

The expected result is numerical likelihood closure and a demonstration that sign-dependent omission can affect amplitude estimation under the specified generator. No real survey noise, censoring law or cosmological correction is estimated here.
