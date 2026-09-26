# Public SALT training-input overlap

The transferred residual is not confined to the identifiable public training
inputs. Of 1,020 validation supernovae, 130 are in the public K21 DES input list
and 890 are absent. Their original fixed-pattern amplitudes are almost equal:
.582329 and .582468. The absence group retains positive matched products in
all ten DES fields. This is a descriptive overlap sensitivity, not a fully
independent training validation.

The [DES release documentation](https://des-sn-dr.readthedocs.io/en/latest/2_LCFIT_MODEL.html)
identifies the nominal SALT3.DES5YR training sample as the K21 sample, with
Fragilistic recalibration and observer-U removal. The same statement is in the
pinned [model README](../../sources/repos/des-science__DES-SN5YR/2_LCFIT_MODEL/README.md).
The exact accepted DES training execution remains unavailable. Consequently,
public input-list membership is the stated scope of this check.

We acquired `docs/_static/SALT3TRAIN_K21_PUBLIC.tgz` from the already pinned
SALTShaker commit `5183c2dcd8f0a6fba61e5a2137bcf186120ef43f`. Its 43,947,844 bytes
match the pinned Git tree blob SHA-1, and a separate SHA-256 is recorded.
Reading the thirteen lists named in `Train_SALT3_public.conf` yields exactly
1,083 source files and 1,083 unique string SNIDs. Source headers, survey names,
sky coordinates and hashes are retained in the [roster](../../sources/updates/2026-09-26-training-roster/roster.csv).
No downloaded training code was executed.

All 206 listed DES IDs match the pinned released DES HEAD by ID **and exact
sky coordinates**. Six of the 43 discovery objects and 130 of the 1,020
validation objects have these matches. The [protocol](training-overlap-protocol.md)
and membership were frozen before computing any membership-stratified score.
The original discovery vector is unchanged; discovery objects were not removed
or refitted after learning their membership.

| Validation membership | Objects | Matched product M | Information I | Fixed gain G, nats | Amplitude M/I |
|---|---:|---:|---:|---:|---:|
| Present in public K21 DES input | 130 | 37.717742 | 64.770523 | 5.332481 | .582329 |
| Absent from that input | 890 | 154.683541 | 265.565617 | 21.900733 | .582468 |
| All original validation objects | 1,020 | 192.401284 | 330.336139 | 27.233214 | .582441 |

All original objects appear exactly once and the two strata sum to the original
M, I and G. The present group has positive matched products in nine of ten
fields; the absent group in ten of ten. No significance test is assigned to
these post-result descriptive partitions.

All observations still share trained surfaces and calibration anchors. Other
survey aliases, training-time cuts and the exact executed calibration/training
chain are not fully linked. “Absent” therefore means absent from this precise
public DES input list, not proven absent from every training or calibration
dependency. This result neither establishes overfitting nor supplies a physical
correction.

Artifacts: [archive acquisition and blob verification](../../sources/updates/2026-09-26-training-roster/acquisition.json),
[roster and source hashes](../../sources/updates/2026-09-26-training-roster/roster-manifest.json),
[membership and sky closure](../../runs/research_2026_09_26/training_overlap/design.json),
[scores and field contributions](../../runs/research_2026_09_26/training_overlap/score.json).
