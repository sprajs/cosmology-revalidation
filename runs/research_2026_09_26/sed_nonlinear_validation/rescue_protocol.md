# Numerical refinement after preserved primary run

The frozen full run reports two-start numerical gradient failures. Preserve
all primary outputs and failed gates. For each failing object, rerun the same
targets/models from the identical original two starts, bounds and covariance,
with all finite-difference Jacobian steps reduced by a factor of ten. The
amplitude derivative remains analytic. Do not change the perturbation or
loosen any convergence threshold. Save refined results separately.

In addition to the unchanged gates, at every returned start solution measure
the Fisher-coordinate gradient using the refined, half-refined and
quarter-refined Jacobians. Require all three norms below the original1e-4
threshold. Independently differentiate the actual half-chi-square objective
along the four local Fisher-normalized parameter directions using central
steps1e-3,3e-4,1e-4,3e-5. Report all values and require the vector norms at
the last two steps below1e-4. Thus a small derivative reported by one finite
step alone cannot make a failed fit pass.

Only if every refined fit passes both the original and extra checks may its
object checkpoint replace the failed numerical checkpoint in a separately
saved resolved view. Retain all1,020 objects. Compare all refined object
responses and objectives to primary values. Preserve primary and refined
source/input hashes, the initial failure count, resolved count and remaining
failure count. This is numerical closure of an unchanged objective, not an
outcome-driven change to the scientific alternative.
