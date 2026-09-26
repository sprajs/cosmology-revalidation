# One-factor sign/censor kernel: independent mathematical review

The formulas in the reviewed kernel are correct for the explicitly fixed
Gaussian covariance `C = diag(d) + v vᵀ`, with all `d_i > 0`. The reviewed
source SHA256 is
`9fe224b5d79d93d423a2a8824d4453fd607a002b97c56c982be5b95f0b60a308`.
No observed flux or fitted native mean was used by this review. The native
covariance alone enters synthetic checks.

Write `Y_i = μ_i + v_i Z + sqrt(d_i) ε_i`, with independent standard-normal
`Z, ε_i`. For observed index set A, the posterior latent distribution is

`V_Z = [1 + Σ_A v_i²/d_i]⁻¹`,
`m_Z = V_Z Σ_A v_i (y_i−μ_i)/d_i`.

The joint density of the observed values and known omitted-row signs is

`p(y_A, signs_N) = p(y_A) E_{Z|y_A}[ ∏_N Φ(s_i(μ_i+v_i Z)/sqrt(d_i)) ]`.

This exactly matches `log_censored`. Its outer `log_gaussian` correctly uses
the marginal `C_AA`, including its determinant. The conditional latent
loading becomes `v_N sqrt(V_Z)` and the conditional mean shifts by `v_N m_Z`.

The normalizer depends on the sampling statement:

- With known full eligible cadence and retained positive values plus omitted
  negative signs, use that joint censored density, restricted to `y_A>0` by
  the sign-selection caller. Do not divide by another sign normalizer.
- If explicitly conditioning on the entire realized sign pattern, divide
  the joint censored density by `P(Y_A>0, Y_N<0)`. This changes the experiment
  by discarding the pattern's parameter information.
- `log_positive_conditioned` instead gives `p(y_A | Y_A>0)` for the named
  retained rows. Its denominator is only `P(Y_A>0)`, with their marginal
  covariance. It does **not** condition on or model the omitted-row pattern.
  That is a different valid conditional sampling statement, not a substitute
  for the full known-cadence censor likelihood.
- A frozen archived-positive mask with newly drawn signed noise has no fresh
  sign truncation. It uses the ordinary marginal Gaussian on that fixed mask.

The generic censor kernel permits arbitrary signs and observed values;
positive-selection callers must enforce the positive support for retained
values themselves. Estimated errors, native flags, actual measurement noise,
event selection and a moving parameter-dependent C are separate modeling
requirements. The rank-one algebra does not establish them.

## Independent synthetic checks

`check.py` adds direct one-dimensional numerical integration of a bivariate
censor likelihood against the exact conditional-normal expression. It verifies
that integrating its retained positive value gives the sign-pattern
probability, while dividing by that pattern probability gives unit integral.
The retained-only positive conditional density also integrates to one.
All-observed and none-observed limits, a full Gaussian comparison, common
loading-sign reversal and permutation invariance pass. Maximum core error
is 4.45e−16.

Using only the saved native covariance, 15 synthetic combinations of
standardized mean −10, −3, 0, 3, 10 and all-positive, alternating, and first-nine
negative signs have maximum 64/128/256-node log-probability spread
4.55e−13. These are synthetic standardized means, not observed truth or
estimated light curves.

## Changes or restrictions needed before generic reuse

Agreement among the three quadrature orders is insufficient as a universal
accuracy gate. A deliberately strong-factor scalar check with μ=−0.2, d=1,
v=1000 returns approximately −0.69314718 for *all* three orders. The exact
`log Φ(μ/sqrt(d+v²))` is −0.693306770125, so every order misses by 1.596e−4.
The narrow CDF transition lies between quadrature nodes. This does not
invalidate the weak-factor native covariance results, but generic use needs
an analytic scalar case and an explicit strong-factor restriction or a
separately validated adaptive integration strategy. Merely increasing and
comparing even Hermite orders is not a proof of accuracy.

`decompose_covariance` also raises an uninformative `StopIteration` for a valid
correlated 2×2 covariance and for a 3×3 covariance with only two active
loadings. Their diagonal-plus-one-factor representation is not uniquely
recoverable from three nonzero off-diagonal entries. An explicit unsupported
domain error, or a documented robust construction, is appropriate. The
current full native covariance does not hit this edge case. Factorize the
supported full covariance once and restrict its `(d,v)` arrays for row subsets
where possible; do not refactor a two-row submatrix unnecessarily.

The parent was notified of both limitations before kernel recovery use.
`protocol.json`, `result.json` and `check.py` retain all synthetic values,
exact comparisons and source hashes. No source repair, observed score,
observed fit or new HMC was performed here.
