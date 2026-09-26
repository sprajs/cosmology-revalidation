# Kernel v2 review and a frozen minimal recovery design

The original v2 snapshot passes this independent bounded review. Of 103
prespecified weak-factor, mode-boundary and native-covariance synthetic
cases, 63 are accepted and 40 correctly refused by the actual computed
precision/mode conditions. The accepted 64/128/256-node log probabilities
differ from independent adaptive integration by at most **5.69e−14**.
Three synthetic censored-density checks also agree within 5.69e−14.
The adaptive integrator's largest reported relative error estimate is
2.00e−13. No observed flux or fitted mean enters these tests.

## Equations and numerical scope

For `a_i=s_i μ_i/sqrt(d_i)` and `b_i=s_i v_i/sqrt(d_i)`, the latent log
integrand is `ell(z)=−z²/2 + sum log Φ(a_i+b_i z)` plus a constant.
Its derivative is `−z + sum b_i λ(a_i+b_i z)`, with
`λ(t)=φ(t)/Φ(t)`. Because `0<λ(t)[t+λ(t)]<1`,

`−[1+sum b_i²] <= ell''(z) <= −1`.

The unique mode is therefore inside [−4,4] precisely when the endpoint
scores have the stated signs. The factor-precision bound 0.5 limits curvature
to [−1.5,−1]. These are valid mathematical guards. They do **not** by themselves
prove a universal 1e−10 Hermite error bound for every finite input.

The independent reference uses a stable `erfcx` Mills ratio, brackets the
mode, and adaptively integrates `exp[ell(mode+u)−ell(mode)]` over u∈[−12,12].
Strong log concavity bounds omitted scaled tails by Gaussian tails beyond12;
this is negligible at the reported precision. It avoids underflow for tiny
orthant probabilities and does not reuse Hermite nodes. The design spans
dimensions2/4/16/74, nominal total factor precisions0.01/0.25/0.5,
standardized means−3/0/+3, all-positive/alternating signs, target modes±3.99
and±4.01, and 15 additional native-C-only synthetic sign/mean patterns.
Some nominal precision0.5 cases round just above the exact code threshold
and are conservatively refused; no tolerance was relaxed.

The v2 review is against unchanged snapshot hash
`ed4c52369914f8cdbb291605de414c9aeae63af6bca4ff4506cda0f438ea26a2`.
Its extreme-tail stress cases are refused, but the old Mills calculation
emits overflow warnings. The parent subsequently hardened the active kernel
with `sqrt(2/pi)/erfcx(−t/sqrt(2))` for negative t. That formula is exact;
setting the positive-tail value to0 at t≥38 is negligible at the declared
floating-point accuracy, not an exact mathematical equality.

This review separately checks that active hardening against the direct
moderate-tail formula (max absolute discrepancy5.16e−14) and requires the
−1e12/−1e150 mode refusals with warnings treated as errors. Those pass.
The reviewed active kernel is copied to `recovery-kernel.py`, SHA256
`52ac4f277e96d61b6b3a32f10417fb013d4416b09a5432a38d2e513b7f4668dc`.
No original v2 artifacts are replaced or relabeled as v2.1 results.

## Smallest useful paired recovery screen — design only

The frozen `recovery-protocol.json` and `recovery-design.npz` specify one
74-coordinate synthetic amplitude problem, **not a supernova template**:
`h_i=0.5 sqrt(C_ii)`, true amplitude a=1. C is the saved B-anchor covariance,
used only as a fixed mathematical noise design. Its native model-error/MW
components and fitted-state provenance do not establish extraction truth.
The synthetic signal is not taken from observed flux or a fitted light curve.

Generate128 complete paired draws `y=h+v Z+sqrt(d)*epsilon`, one shared
latent plus74 independent errors. Use each exact same full vector for four
estimators:

1. Full signed Gaussian on all74 coordinates.
2. Naive Gaussian on the newly positive coordinates A, with marginal C_AA.
3. Joint censored known-cadence likelihood: `p(y_A, Y_N<0)`.
4. Full-pattern conditional likelihood:
   `p(y_A,Y_N<0)/P(Y_A>0,Y_N<0)`.

The fourth denominator is over the **whole** realized sign pattern. It is
not the kernel's `log_positive_conditioned`, which conditions only the named
retained A rows. Correlation makes these different. The joint censor arm
does not receive an additional normalizer. For comparison, a single separate
old-noise seed creates an archived-positive mask; each new signed vector on
that fixed mask receives an ordinary marginal Gaussian fit. This cheap
analytic control distinguishes fresh truncation from fixed-cadence loss.
No draw is regenerated because its number of positives is inconvenient.

The amplitude domain is explicitly [0,4], including the zero boundary as a
limiting signal. The three Gaussian cases use the closed-form amplitude and
report unconstrained and boundary-constrained results. The two proper
likelihoods use nested129/257-point grids, every sampled local bracket and
endpoint, followed by continuous polishing. All-negative full-pattern
conditional data contain no retained values and have identically flat
conditional likelihood; flag nonidentification, never select an arbitrary
amplitude or silently drop the draw. Preserve every zero, upper boundary,
tie, quadrature refusal and optimization failure.

The benchmark uses the first8 predeclared realizations; they remain in the
128, not a tuning set. The entire screen is capped at10 minutes on one
worker. It entails256 one-dimensional likelihood fits and384 analytic
amplitude calculations, with roughly a few hundred million scalar CDF terms
before cache/polish savings. No SNANA or HMC runs are needed. If the measured
runtime forecast exceeds the cap, retain the partial benchmark and report
the resource limit; do not alter truth, masks or tolerances.

At true amplitude and all candidate optima, require 64/128/256-node
log-likelihood agreement within1e−8 and independent adaptive checks.
Coarse/fine best log likelihood must agree within1e−7 and identified
amplitude within1e−5. Central score steps1e−4 and5e−5 must agree within1e−6.
Report amplitude bias/RMSE, paired changes and Monte Carlo standard errors;
also report the score at truth. Correct likelihoods have zero expected
score, but their finite-sample MLEs need not have zero bias. A score departing
from zero by more than4 Monte Carlo standard errors is a diagnostic to
investigate, not permission to change the frozen truth or generator.

Report `deltaD=−2.5 log10(a)` separately. At a=0 it is undefined/infinite;
report that mass and the exact finite subset rather than computing a
silently conditioned average. The nonlinear logarithm has Jensen bias even
when unconstrained amplitude is unbiased. No dust law, shape fitting,
population evolution, survey event selection or cosmology is inferred.

This is a kernel/estimator recovery screen before a native-template,
multiple-shape-branch experiment. The latter still requires its own frozen
truths, cadence/duplicate identity gates and covariance assumptions.

Recovery protocol SHA256:
`39e07714c9c9de2b61d1a3ec3378335e14e9d3c6822202da1c8543d9fffb4bf0`.
No recovery draws or estimator fits have been launched by this task.
