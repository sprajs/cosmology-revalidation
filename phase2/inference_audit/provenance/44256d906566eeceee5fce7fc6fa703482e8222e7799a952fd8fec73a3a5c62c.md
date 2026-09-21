# Independent engine: validation and deliberate failures

These tests use no DES observations and no SNANA likelihood. They establish numerical and statistical behavior under stated synthetic laws; they do not validate those laws as descriptions of real supernovae.

`scripts/phase2/hierarchy/validate.py` compares the integrated Gaussian-plus-exponential likelihood with independent adaptive quadrature in one, two and three dimensions. The maximum log-density difference across 36 cases is 1.2e-14. Independent selection integration agrees to 2.8e-16. Constant-q distances agree with analytic de Sitter, coasting and Einstein–de Sitter limits to 5.1e-13 magnitude. These comparisons test algebra and numerical implementation, not astrophysical correctness.

The selected 1,000-object Bayesian pilot generates latent stretch, intrinsic effective colour, exponential reddening, luminosity scatter and correlated measurement errors, then applies a known probability of detection to each draw. Redshift, host class and measurement covariance are conditioned on; this is not a volumetric survey simulation. The selection probability is deliberately strong at high redshift. The fitted model integrates over its latent population and divides by the probability of selection. Two independent chains, each with 600 warmup and 1,000 retained draws, have no divergences and R-hat below 1.002. Changing the tiny-matrix implementation and sampler mass adaptation gives consistent results. The faster run takes 81.5 seconds. The first realization's dust scale differs from truth by about 2.5 posterior standard deviations; that result is preserved, not tuned away.

The prespecified repeated-recovery experiment uses 30 new seeds, 1,000 selected objects per seed and an independently implemented NumPy generator with SciPy distance quadrature. Maximum-likelihood estimates and observed-Hessian errors are used here: this is a frequentist recovery diagnostic, **not** Bayesian simulation-based calibration. All 30 fits converge from two starts to the same objective within 3e-8; their Hessians are positive and no parameter reaches a declared bound. Across the eleven parameters, nominal 95% intervals include truth in 27–30 of 30 repetitions. For the dust scale they include truth in 30/30, with a binomial 95% coverage interval of 0.884–1.000. Thirty realizations cannot tightly certify 95% coverage. The retained pilot discrepancy is not accompanied by evidence of a repeated dust-scale bias in this experiment.

Four further seeds per case let the expansion parameters vary, with generating q0=-0.5 and q1=1.2 in q(z)=q0+q1*z/(1+z). These small samples are counterexamples and debugging diagnostics, not precision estimates of bias:

| Analysis/generating law | Recovered q0, four seeds | Typical reported standard error |
|---|---|---|
| Correct likelihood and selection | -0.311, -0.497, -0.321, -0.667 | 0.11–0.12 |
| Incorrectly omit selection normalization | -0.995, -1.111, -0.955, -1.308 | 0.11 |
| Same-variance t4 luminosity scatter, fit Gaussian scatter | -0.591, -0.455, -0.460, -0.694 | 0.13–0.15 |
| Equal-mean two-exponential dust, fit one exponential | -0.396, -0.276, -0.527, -0.516 | 0.11 |
| Unmodelled +0.30*z/(1+z) magnitude evolution | -0.700, -0.766, -0.646, -0.901 | 0.11 |

The rows with different generating laws consume random numbers differently, so they are seed-matched experiments, not identical latent realizations. The correct-versus-omitted-selection rows use the same saved data. All of these fits can converge well while estimating the wrong acceleration under a misspecified law. The sign and magnitude of the failures are specific to these synthetic scenarios; they are not measurements of a bias in DES. Heavy tails alone do not create an equally large shift in this limited experiment.

Full configurations, input realizations, per-parameter estimates/errors, Hessians' minimum eigenvalues, optimizer messages, objective gaps and hashes are preserved under `phase2/hierarchy/coverage/`. The exact loaded source is snapshotted under `phase2/hierarchy/provenance/`. The experiment is reproducible with:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 uv run --project phase2/hierarchy --frozen python scripts/phase2/hierarchy/validate.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 uv run --project phase2/hierarchy --frozen python scripts/phase2/hierarchy/coverage.py
```

The next gate is empirical: a DES selector must reproduce the actual trigger, host-redshift and fit/classification quality decisions or quantify its approximation. The toy selector is not silently applied to DES. Actual parameter covariance must also come from a coherent refit, and model comparison must retain the same held-out objects and a common flexible distance relation. Source extraction, SALT training and shared calibration assumptions remain outside the scope of this summary-observable engine and require the separate observation-level audit.
