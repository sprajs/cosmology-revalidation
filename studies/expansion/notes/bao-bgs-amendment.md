# Exploratory BGS extension, frozen after the first cone results

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The registered anisotropic-only cone result was nonrejecting. A subsequent
continuous, nonaccelerating fit predicted BGS DV/rd=8.17963 at z=.295, compared
with 7.94168 +/- .07609; this is a plug-in residual and excludes fit uncertainty.
This motivated the following **exploratory** measured-observable test. It does
not replace the registered result. The six AP contrasts were already inspected;
their maximum deficit was 2.67247 standard errors.

Under the same flat, constant-ruler, positive-expansion and fixed compressed
Gaussian likelihood assumptions, define t=ln(1+z), g=(1+z)DH/rd, d=DM/rd.
No acceleration makes g nonincreasing. At BGS redshift b, d_b/t_b >= g_b.
For every later anisotropic point j, g_b>=g_j. Since

    (DV_b/rd)^3 = [z_b/(1+z_b)] d_b^2 g_b,

the measured quantity

    s_b = (DV_b/rd) / ([z_b/(1+z_b)] t_b^2)^(1/3)
        = ((d_b/t_b)^2 g_b)^(1/3)

satisfies s_b>=g_j. This supplies six *linear* necessary inequalities in
released observables even though DV itself is a nonlinear function of geometry.
It does not require assigning BGS a separate DM or DH measurement.

Freeze a joint maximum-standardized-violation diagnostic across these six,
the six earlier AP contrasts, the five adjacent radial monotonic contrasts,
and the eleven anisotropic interval-integral bounds (28 rows in total).
Keep all correlations induced by shared observables/BGS. Every contrast is
nonnegative under the null and is simultaneously zero for coasting; therefore
the zero-contrast multivariate Gaussian is a least-favourable null for this
specific maximum score. Some rows are redundant: that is allowed and their
correlation is explicitly retained. Do not multiply 28 independent p-values.

Calibrate by 1,000,000 Gaussian realizations using the full original 13-row
covariance; report the binomial interval and all contrasts. Also show the
six BGS family diagnostic separately. The new rule was selected after prior
diagnostics and is exploratory, even though the finite contrast family has
an exact least-favourable calibration. It does not adjust for every possible
research choice, nor independently validate DESI extraction/systematics.

A rejection supports some acceleration somewhere before the tested endpoint,
under these assumptions. It does not determine q(0), dark-energy microphysics,
or the validity of a particular supernova correction. A changed ruler,
curvature or compressed-likelihood failure could change the interpretation.
