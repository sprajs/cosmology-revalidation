# Improving precision of the frozen infrared predictive integral

A direct bridge between the optical posterior and an auxiliary optical-plus-infrared posterior is a defensible next integration method. It estimates exactly the same predictive density as the original calculation. It has not been run on BayeSN, and it is not yet known whether it can meet the unchanged 0.05 Monte Carlo error limit economically. The [failed direct score](../results/bayesn-synthetic-dense99-case00.json) remains unchanged; no additional physical chains or observed-data fits were launched for this investigation.

For either proper-distance arm, let θ be the **same 47 latent coordinates**, π(θ) its normalized prior, O the optical vector, and N the complete six-measurement J/H vector. Define

$$q_0(\theta)=\pi(\theta)L_O(\theta),\qquad
q_1(\theta)=q_0(\theta)L_N(\theta),\qquad
p_i=q_i/Z_i.$$

The desired conditional density is

$$p(N\mid O)=r=Z_1/Z_0=E_{p_0}[L_N].$$

For any suitable common bridge function h, integration of the same product gives

$$r=\frac{E_{p_0}[q_1h]}{E_{p_1}[q_0h]}.$$

Using the asymptotically optimal independent-draw bridge, with fixed sample fractions s₀ and s₁, reduces this to the scalar equation

$$r=\frac{\operatorname{mean}_{p_0}[L_N/(s_0r+s_1L_N)]}
{\operatorname{mean}_{p_1}[1/(s_0r+s_1L_N)]}.$$

This is the common-support normalizing-constant identity and mixture bridge of [Meng and Wong (1996), §§3–4](https://www3.stat.sinica.edu.tw/statistica/j6n4/j6n43/j6n43.htm), specialized here to an added likelihood factor. Their optimality statement assumes independent draws; correlated MCMC requires separate error accounting. The stable implementation solves for log r with bounded logistic terms. A root solver converging does not prove adequate overlap.

No Gaussian approximation to the actual 47-dimensional posterior is needed. The common prior, optical likelihood and any common coordinate-transform Jacobian cancel inside q₁/q₀. The added factor still includes every Gaussian normalization constant and sums all six residual terms **within the same θ**. Multiplying six separately averaged predictions would estimate a different quantity.

Fitting the auxiliary distribution does not retune the frozen predictive model: it supplies integration draws in regions where Lₙ contributes strongly. Its posterior mean or fitted best light curve must never replace the predictive integral. The original optical posterior, proper priors, M20 surfaces, residual-SED covariance, calibration, passbands, epochs, uncertainty model and score definition all stay fixed. A model change, outcome-driven noise inflation or selection of favorable cases would invalidate this equivalence.

The existing `Kernel` supports arbitrary supplied row vectors and one shared 47-coordinate input. Case 0 has 15 optical and six infrared rows. A separate executor could concatenate the vectors and call the unchanged `numpyro_model` with a combined-row kernel, preserving the single D, Aᵥ, Rᵥ, shape, timing and 42 residual-SED coefficients. Its observation-site name is cosmetic; the new executor would explicitly label the distribution auxiliary. It must first check that its log density minus the optical log density equals the normalized six-point NIR likelihood at independent physical and transformed points. The existing optical worker and its file-isolation guard should remain untouched. Only a separate, sealed **synthetic** payload would be available to an auxiliary worker after optical qualification.

The [conjugate test](../code/bayesn_bridge_gaussian.py) uses a 47-dimensional Gaussian optical posterior with six informed linear directions and a correlated six-dimensional observation. Its other 41 directions integrate out exactly. All auxiliary samples are independent analytic Gaussian draws, so this is not a BayeSN runtime benchmark. Across 64 fixed-seed repetitions per scenario, log-density RMSE was:

| Prespecified Gaussian case | Direct 4,000 draws | Bridge 4,000 + 4,000 draws | Direct 64,000 draws |
|---|---:|---:|---:|
| Moderate concentration | 0.0220 | 0.0162 | 0.00643 |
| Narrower likelihood | 0.0626 | 0.0353 | 0.0189 |
| Shifted likelihood | 0.1167 | 0.0405 | 0.0283 |

The [test record](../results/bayesn-bridge-gaussian.json) retains every replicate. Constant-likelihood, reciprocal-ratio and large additive log-density controls agree exactly at recorded precision; independent scalar conjugate quadrature agrees within 2.3 × 10⁻¹⁶ in log density. The bridge improved precision per total draw in the narrower and shifted examples; direct 64,000 draws were more precise. Analytically cheap auxiliary draws in this test do not establish a physical speedup.

For a future physical bridge, write a=Lₙ/(s₀r+s₁Lₙ) and b=r/(s₀r+s₁Lₙ). The derivative magnitude of its estimating equation is

$$K=s_0E_{p_0}[ab]+s_1E_{p_1}[ab].$$

The leading variance of log r is the sum of the two **Monte Carlo variances of the sample means**, divided by K². Those variances must preserve within-chain dependence; the independent-draw standard error in the Gaussian helper is not an acceptable physical-MCMC error estimate. The future design should retain the larger between-chain and contiguous-block estimates, inspect block-size sensitivity, and propagate the independently seeded two-arm errors into the final score difference. The original error threshold stays **<0.05**, alongside every original auxiliary-chain convergence diagnostic. Poor overlap, an unstable K, remaining optical-sample error or an unconverged auxiliary posterior must remain a failure. The frozen optical sample imposes a precision contribution that arbitrarily many auxiliary draws cannot remove for a fixed bridge.

The measured direct difference MCSE is 0.1549. Sixteen times as many equally effective optical draws would reduce it to about 0.0387 under a stationary square-root scaling assumption; that is a forecast, not a guarantee. Repeating four chains with 1,000 warmup and 16,000 retained draws would require 8.5 times the present iteration count. Scaling the observed 168.1 combined worker-minutes gives roughly **24 worker-hours**, with considerable uncertainty from adaptation, autocorrelation and CPU contention.

An initial auxiliary campaign would instead add four chains of 1,000 warmup and 1,000 retained draws per arm, subject to the unchanged gates. It would process 21 rather than 15 epochs and may have tighter or more difficult geometry. A linear epoch-count forecast gives approximately four combined worker-hours, but it is too crude to set a launch budget or promise success; the broad arm alone could exceed the previous two-hour cap. A short, separately declared gradient/runtime probe is therefore needed before allocating auxiliary chains. Neither that probe nor the chains have been launched.

The [preparatory design](../code/bayesn-bridge-integration-design.json) fixes the conjugate controls and states the prerequisites for a later executor. The source paper was retrieved from its [university-hosted PDF](https://dept.stat.lsa.umich.edu/~xuanlong/courses/stat700-f09/meng-wong-96.pdf), SHA-256 `d3ed5b55fac1437a5de7ffa7674d7bd55de0784536449e44606aaa67a1682d9e`; it is an ignored reference download. This remains method development for one synthetic case, with no age-correction or observed cosmological conclusion.
