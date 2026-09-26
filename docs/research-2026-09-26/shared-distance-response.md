# What the residual test leaves uncertain about distance

The fitted residuals constrain only part of a calibration change. Another part
is absorbed by each supernova's fitted brightness, colour, width and time. We
calculated both parts for all 43 discovery and 1,020 validation objects using
the exact native covariance and saved local derivatives. The
[protocol](distance-response-protocol.md) fixed the target before these results:
the mean pre-BBC standardized response of the highest 255 validation redshifts
minus the lowest 255. The ranges are .61002–1.04763 and .07411–.34336.

The numbers below are **conditional response calculations, not measured biases**.
They propagate the finite Gaussian mode approximation and its discovery-only
posterior. They are neither total distance uncertainty nor a bound on all
calibration or luminosity-evolution errors. Alpha=.16087 and beta=3.1178 remain
fixed; selection and BBC are not rerun.

| Mean-response family | Prior contrast SD (mag) | Discovery-conditioned contrast mean (mag) | Discovery-conditioned contrast SD (mag) | Prior contrast variance remaining |
|---|---:|---:|---:|---:|
| Nine coupled calibration/SALT variations | .013659 | +.011103 | .012683 | 86.21% |
| Twelve variations, adding CALSPEC/MW scale/MW law | .013793 | +.011342 | .012793 | 86.03% |
| Nine variations plus empirical observer-band term | .087327 | +.058297 | .025709 | 8.67% |
| Twelve variations plus empirical observer-band term | .087348 | +.058318 | .025717 | 8.67% |

For the released-variation approximation, discovery residuals leave about 86%
of the prior variance in this redshift contrast. Good residual agreement is
therefore weak evidence that the corresponding distance response is small.
The broader empirical observer family permits a larger response; its positive
mean is not a justified correction, especially because the held-out residual
test does not support applying the original pilot amplitude in full.

An independent follow-up propagates the two already declared prior/weight
sensitivities through these same saved response arrays. The equal-total-power
isotropic observer prior changes the twelve-mode contrast to
**+.062910 ± .026041 mag**, versus +.058318 ± .025717 above, despite similar
flux predictive scores. The alternative paper colour-law weight gives
+.011536 ± .012872 mag for systematics alone and +.058318 ± .025717 for the
inherited observer extension. These remain conditional contributions, not total
uncertainties or detected biases. The
[separate calculation](../../runs/research_2026_09_26/shared_mode_review/distance-prior-sensitivity.json)
agrees with the original arithmetic and retains the unchanged target.

The 86% figure uses **43 discovery objects**, as required by the predictive
test; it is not an information limit for all available light curves. A separate
[design-only calculation](../../runs/research_2026_09_26/shared_distance_information/result.json)
uses the information matrices of all 1,063 objects without reading residual
values or fitting a new mean. Within the twelve-mode Gaussian approximation,
the contrast SD would be .00942 mag, leaving **46.6%** of its prior variance.
With the empirical observer extension it would be .01207 mag (.01236 with the
isotropic prior). These still exclude unknown modes and population/selection
uncertainty.

The exact CALSPEC-versus-observer attribution degeneracy has effectively zero
projection on this particular contrast (squared target projection below 6e-31):
its gray component cancels in high minus low. That does not constrain arbitrary
redshift-dependent gray luminosity evolution, which is outside these finite
mode families. Reproduce this information check with
[shared_distance_information.py](../../scripts/research_2026_09_26/shared_distance_information.py).

The calculation retains shared covariance: if `h` is the vector of mode
responses for the high-minus-low contrast, its variance is `h.T V_D h`.
Treating each supernova's calibration error as independent would lose these
cross-object terms and incorrectly suggest that averaging 255 objects removes
the systematic uncertainty.

For a mode's nominal-unit flux change `delta F`, local refitting at fixed data
gives `delta p=-pinv(L^-1 J) L^-1 delta F`, where
`p=(mB,x1,c,t0)` and `C=L L.T`. The standardized response is
`delta mB + alpha*delta x1 - beta*delta c`. CALSPEC is first mapped back to
nominal flux units; coupled SALT variants already contain their released
training response. This does not reproduce a full refit of alpha/beta or
preserve an invariant physical interpretation across retrained coordinates.

As an exact tangent check, uniform one-magnitude model dimming has zero
projected residual and is compensated by `delta mB=-1`, with the other three
coordinates unchanged. The largest gray projected norm is 5.53e-14. The unit
amplitude is a derivative check, not a proposed one-magnitude finite correction.
Weighted normal equations independently reproduce the SVD parameter responses
within 5.2e-15; tangent-plus-orthogonal decomposition agrees within 1.1e-15.
An independent Astra review confirmed the sign, shared propagation and limited
interpretation.

Reproduce using [shared_distance_response.py](../../scripts/research_2026_09_26/shared_distance_response.py)
with `--modes 9` or `--modes 12` into new preserved output directories.
[Nine-mode results](../../runs/research_2026_09_26/shared_distance_response_9/result.json)
and [twelve-mode results](../../runs/research_2026_09_26/shared_distance_response_12/result.json)
include checks and exact contrast definitions. Each run also saves every
object/mode response, cohort membership, arrays and hashes of all consumed
numerical inputs and outputs. The original results are immutable by the script.
