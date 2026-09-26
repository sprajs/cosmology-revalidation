# Using known omitted signs without discarding correlated errors

The first RAISIN native fixed covariance has the exact numerical form
`C = diag(d) + v vᵀ`, with positive diagonal remainder. The reconstruction
residual is `2.76e-17` relative to the matrix norm; its largest pairwise
correlation is 0.002695. This permits a one-dimensional latent-Gaussian
integral for sign probabilities. It does not prove that the real extraction
noise follows this covariance or a Gaussian distribution.

The [kernel](../../scripts/research_2026_09_26/onefactor_sign_likelihood.py)
implements the conditional model

    Y_i = m_i + v_i U + sqrt(d_i) epsilon_i,
    U, epsilon_i independently standard normal.

For known signs `s_i` (positive=+1, negative=−1),

    P(s_i Y_i > 0 for all i)
      = integral phi(u) product_i Phi[s_i (m_i + v_i u)/sqrt(d_i)] du.

The product is inside the shared-latent integral. Multiplying independent
marginal sign probabilities would generally lose the correlation.

For measured values on rows A and only signs on rows N, the joint censored
density is `p(y_A) P(sign_N | y_A)`. Conditioning on the measured values gives

    Var(U | y_A) = 1 / [1 + sum_A v_i²/d_i],
    E(U | y_A) = Var(U | y_A) sum_A v_i(y_i − m_i)/d_i.

The remaining sign integral uses this updated Gaussian latent distribution.
This differs from a likelihood for retained positive values conditioned on
those retained rows being positive, which divides their joint Gaussian
density by their joint positive-sign probability. Conditioning on the
**entire** known sign pattern requires a denominator involving both retained
and omitted indices. The code keeps these estimands separate; no
normalization is applied twice.

The [original frozen validation](../../runs/research_2026_09_26/onefactor_sign_validation/protocol.json)
passes its synthetic checks: independent-error limit, all 16 sign patterns
in a four-dimensional example, Gaussian density against an independent
implementation, and a single censored row against its analytic
conditional-normal formula. Errors are below `9e-16`. Gauss–Hermite orders
64, 128 and 256 agree to approximately `2e-15` for the tested log probabilities.
With 100,000 synthetic draws, a four-dimensional sign frequency is 0.082680
against 0.082850 predicted (Monte Carlo SE 0.000872). Using native covariance
but an entirely synthetic mean gives 0.629370 against 0.630759 predicted
(SE 0.001526). These are recovery checks, not observational significances.

The [independent review](../../runs/research_2026_09_26/raisin_profile_solver_review/onefactor-review/README.md)
verified the formulas and native weak-factor examples, but found a generic
quadrature failure: a very strong shared factor can make all three quadrature
orders agree on an inaccurate answer. The original results/source are
preserved. The [v2 amendment](../../runs/research_2026_09_26/onefactor_sign_validation/v2/amendment.json)
therefore makes scalar and independent cases analytic, and explicitly
restricts multirow quadrature to total factor precision `sum(v²/d) <= 0.5`
and a latent log-density mode within [−4,+4]. The mode check uses the
derivative of the log-concave integrand. Unsupported strong/tail cases raise
an explicit error; quadrature-order agreement alone is not treated as proof.
Sparse/two-row covariance factor reconstruction also refuses its unsupported
case explicitly rather than returning an invented factor.

The [v2 checks](../../runs/research_2026_09_26/onefactor_sign_validation/v2/result.json)
repeat the original recovery, reproduce the strong scalar analytic result
and verify all three refusal paths. The native matrix's total factor
precision is 0.01934. These checks cover a stated numerical domain, not every
possible Gaussian parameter combination. No observed flux likelihood has
been evaluated or fitted.

The independent v2 review compared 103 prespecified cases to separate adaptive
quadrature: 63 accepted cases agree within `5.69e-14` in log probability;
40 outside the supported domain were correctly refused. A subsequent
[tail-hardening check](../../runs/research_2026_09_26/onefactor_sign_validation/v2_1/tail-guard.json)
uses the scaled complementary error function for the negative-tail inverse
Mills ratio, preventing warning/overflow in the remote-mode refusal path.
The original and amended sources remain separate. The
[frozen paired recovery design](../../runs/research_2026_09_26/raisin_profile_solver_review/onefactor-review/v2/recovery-protocol.json)
tests estimators on a synthetic signal. The [128-draw execution](onefactor-recovery.md)
is complete. A separate root verifier imports neither the kernel nor the
execution code: it directly integrates the latent Gaussian times observed
normal densities and missing-row CDFs. All 2,205 candidate likelihoods agree
within `1.14e-13`, all 384 matrix amplitude solutions within `1.34e-15`,
and 16 independent bounded searches reproduce the optima within `1.36e-7`.
Those searches include all three zero-boundary cases; their log distances
remain undefined. Summary arithmetic agrees within `8.33e-17`.
The [review record](../../runs/research_2026_09_26/onefactor_recovery/root-independent-review/result.json)
preserves its explicitly retrospective selection and scope.

Only fixed Gaussian covariance and a known observation process are covered.
Jointly estimated errors, unknown rejection causes, non-Gaussian extraction,
model-dependent covariance, event/follow-up selection and population
normalization remain separate. Available signed measurements should use
their full signed density; deliberately censoring them discards information.
The kernel supplies a controlled comparison for information lost through
known sign omission before a physical bias correction is considered.
