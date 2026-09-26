# Final static release review

**PASS for root release of Stage A only.** Independent byte hashing verifies all23 frozen input/source records and sizes. No executor was imported or executed, and no scientific pixel arrays or aperture outcomes were evaluated by this reviewer.

- Input freeze SHA-256: `4bd87011fe9b519090e65cef6f020d62502b87b185ebd2abc7725266003fdb6f`.
- Executor SHA-256: `c5cfac0efa5a53fd6895f0251c519cafac84c3a218838ac08e80a063c34dbb01`.
- Geometry SHA-256: `1d7e0c9b4990d551da563b483d89d45d374a57ea19175c4a04f9066f18648bc4`.

The source uses only design exposures1&3 for source detection. Stage A accesses measurement SCI solely for finite-value validity, while ERR/DQ determine the declared accepted-pixel support; no held-out aperture or pair sum occurs. Source-mask selection, data-dependent DQ and ERR remain explicit conditioning, not a claim of statistically independent preprocessing.

The native operator has correct PAM/PHOTFLAM units, complete aperture/background weights and retained-area masked normalization. Stage B carries shared-template covariance through `H D Hᵀ`, retains signed results, and includes the separately frozen raw-moment diagnostics. The algebraic checks do not prove a Gaussian measurement likelihood or unbiased quoted variances.

Prefreeze issues are resolved: fewer-than30-position branches stop before sums; insufficient tile support suppresses resampling without deleting positions; numerical geometry/synthetic failures stop rather than prune candidates; DQ decoding uses the actual16stored bits; finite positive variance is required; partial ledgers/operators are retained; and a complete primary checkpoint precedes dependent sensitivities. Radius/reversal diagnostics retain explicit dependence and cannot replace the primary.

Resource safeguards use a50MB retained-operator-field ceiling, eight-pixel quadrature chunks at up to128² samples, and one concatenated field at a time. These bound the algorithm's new-array allocation design; this is not a measured process-RSS certification. Stage A has the120-second initial benchmark cap. Stage B requires its own root release bound to the completed Stage A hash and one cumulative cap covering its sensitivities.

`final-static-review.json` preserves all23 verification records and the detailed checklist. `final-run.py`/`final-geometry.py` are immutable reviewed snapshots; `verify_final_freeze.py` reproduces the independent hash audit without importing the executor. The earlier draft review and source snapshots remain present.

This recommendation authorizes no action by itself. Root controls release. The experiment concerns current2026 public products, not historical RAISIN reproduction, CRNL correction, total-source flux recovery under masking or cosmology.
