# Official-objective implementation discriminator

Freeze before official43 residual scores. Keep the same43 published
PROB_SNNV19>.999 objects and original published accepted mask. Use existing
SNANA-audit-v3 output-only executable with the conditioned-mask configuration,
no new clipping/phase membership, and each object's published coordinates as
initialization. This is a modern pinned implementation, not a claim to recover
the unknown historical runtime. Match all1649 epochs one-to-one by band, time
and flux, retaining repeated multiplicities. Fail rather than drop objects.

Export actual accepted-epoch predictions, frozen flux covariance/inverse,
fitted coordinates and covariance/prior setup. Verify rTC^-1r+prior closes to
the dumped objective. Check unchanged results against the ordinary executable
on the identical small cohort if needed; the existing patch is output-only.

At the same exported parameters/epochs compare native independent mean flux
with official mean flux under the same exact covariance. Recompute local
native SALT nuisance Jacobian; use the exact amplitude derivative for the
official mean when projecting it. Observer-band templates use the respective
mean flux; rest/phase templates are the previously declared nuisance family.
Record this remaining approximate-tangent boundary explicitly.

Repeat frozen family/scale/object-field prediction and block sign checks on
official-mean/exactC, native-mean/exactC and native-mean/nativeC. The first two
isolate mean implementation; the latter isolates covariance at fixed native
mean. Their comparison is a secondary implementation diagnostic after the
observed64/43 signal, not new pristine discovery significance. Retain previous
measurement-only results and do not select the most favorable covariance.

Before interpreting a field/band residual as photometric calibration, quantify
the projected official-minus-native residual subspace and its ability to
explain the apparent band coefficient. If historical F99 remains a plausible
source, run a separately named matched OPT_MWCOLORLAW=-99 arm after verifying
the built source's option mapping. No data or cosmology correction follows
from these fixed-mask linearized tests alone.
