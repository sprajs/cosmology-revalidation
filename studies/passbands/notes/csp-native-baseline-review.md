# CSP native baseline: reproduction and convergence are separate

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The original 42-object, three-iteration baseline fails the frozen cohort gates.
Those failures remain preserved. A separately specified twelve-iteration,
current-input comparison now passes numerical gates for all 42 objects and
gives a mean filter-label response of **+0.000483432 mag**; see the
[complete response and independent reviews](csp-stable-filter-response.md).
It does not retrospectively repair the original historical reproduction.

The independently reviewed pilot changes only six known WIRC J labels for
2004ef to the paper-supported RC2 approximation. Its raw fitted distance changes
by **+0.0009077025 mag**; unchanged control 2005hc is exact. Actual post-initialization
starts separated by 0.4 mag converge within 4.30e-8 mag. All 42 exported pilot
callback quadratics close independently within 2.28e-13, and the expected
multiplicative brightness law closes within 2.06e-9 relative error. The six
default callback physical-row multisets agree after only the declared relabel.
This is a conditional processing response, with fixed template, shape,
extinction and timing; it is not a cosmological correction.

## All 42 nominal objects

Nominal and exact-copy FITRES and LCPLOT outputs are byte identical. All 42
objects return native success flags, but:

- 2007A's replay chi-square exceeds the archived value by 0.32136, failing the
  original tolerance. Its distance differs by -0.00028 mag, within tolerance.
  The released peak date differs from the archived initializer by three
  float32 units, or 0.01171875 day. See the
  [source and header audit](csp-peak-header-provenance.md).
- Six objects move by more than 0.001 mag between iterations two and three:
  2007ca, 2008hu, 2004ey, 2006kf, 2008ia and 2009ab. The largest change is
  +0.0101198 mag for 2008ia. All final two row masks agree.
- The largest final covariance change is 0.199023 in relative Frobenius norm
  for 2009ab. The corresponding **inverse-covariance** change is 0.202827.
  Earlier secondary prose conflated these; its original bytes and correction
  history are preserved.
- Nine objects use 24 first-pass rows later than the stored template grid's
  70-day endpoint, reaching 93.915 days. All second- and third-pass rows fall
  within the grid. Full nonzero filter transmission, native wavelength
  limits, shape and positive model means pass for every exported row.
  The source-defined late extrapolation is a separate operation from
  empirically supported interpolation.

Independent root parsing reproduces all 126 callback quadratics within
5.46e-12. The final fixed-covariance amplitude optimum is within 1.87e-5 mag
for every object. Thus the individual frozen-covariance minimizations work;
that does not show that updating the covariance has converged.

## What the next numerical experiment tests

Source inspection finds that the native pipeline stops at its configured
iteration count. Its MINUIT covariance-status check concerns parameter errors,
not convergence of the observation covariance. A separately frozen, nominal-only
9/12-iteration experiment retained all 42 objects, every original header,
the existing objective and the failed historical gates. It tests continuation
stability, genuinely different starts, covariance changes and accepted rows.
Any ensuing filter response must be labelled a response of the **current
released inputs**, not exact reproduction of the discrepant archived input.

The continuation completed and passes independent root review of 2,016
callback states. The first three and first nine states reproduce their shorter
runs exactly. All 42 distances at iteration nine equal those at iteration
twelve at exported precision; different starting distances agree within
2.53e-8 mag. The largest last-update covariance change in the whitened operator
norm is 1.41e-7. The mean all42 change from iteration three to twelve is
**+0.0000413886 mag**, with range **[-0.0000908122, +0.001231436] mag**.
This resolves the tested numerical stopping issue without resolving the
archived-header discrepancy or validating the first-pass extrapolation.
The [independent continuation result](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/csp_native_filter_response/root-state-review/continuation/result.json)
and frozen inputs support a new, separately declared current-input filter test.

The source covariance also depends on model brightness. A fixed point of
iteratively frozen weighted least squares generally differs from minimizing
either the continuously changing quadratic or a normalized Gaussian likelihood.
This mathematical distinction alone does **not** establish estimator bias or
make the Gaussian alternative physically correct. Self-consistent iterative
methods can remove multiplicative-error biases under specified assumptions,
as demonstrated by [Ball et al.](https://arxiv.org/abs/0912.2276).
Supernova comparisons of chi-square and complete-likelihood inference are
also established research, for example
[Lago et al.](https://arxiv.org/abs/1104.2874).
The native estimator and any proposed alternative need recovery tests under an
explicit measurement and scatter model before a difference can become a
correction. [Kessler et al.](https://arxiv.org/abs/1209.2482) provide relevant
primary evidence that the assumed wavelength dependence of intrinsic scatter
can affect recovered supernova distances.

## Artifacts

- [Native cohort baseline report](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/csp_native_filter_response/cohort-execution/baseline-report.md)
- [Root changed-pilot review](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/csp_native_filter_response/root-state-review/changed-pilot-v3/result.json)
- [Root all42 state review](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/csp_native_filter_response/root-state-review/full42-nominal/result.json)
- [Root covariance versus precision changes](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/csp_native_filter_response/root-state-review/full42-nominal/C-versus-W.json)
- [Rowwise source-grid and throughput review](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/csp_native_filter_response/root-state-review/full42-support/result.json)
- [Nominal continuation protocol](../specifications/experiments/filter_response/convergence-diagnostic/protocol.json)

The first two root changed-pilot reviewer attempts failed on path resolution
and an incorrect float32 join for double-precision exported MJDs. Both are
preserved beside the passing third review; neither changed any fit output.
