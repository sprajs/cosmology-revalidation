# Initializer versus actual optimizer start

The frozen gate correctly stopped the branch. All eight `joint_minus` NML initializers and first `PROSP_PEAK_BOUNDS` values are exactly **−2 days** from nominal. Native `FITINI_ADJUST` moves all eight by **+2 days** to the same grid winner, MJD57707.80078125, before `CSP_ENTRY` and `MNPARM`. The actual optimizer-entry separation is therefore **zero for all eight**, rather than the prescribed −2days. The parameter was not ignored by the reader, and the physical peak intervention was not tested yet: this run perturbs the initialization grid, not the optimizer entry.

`per-object.csv` and `result.json` give independently parsed numerical evidence. Both branches have22656 saved CSP records. They are not bitwise identical. The largest final distance difference is4.5793555614e−9mag and peak difference1.8982245820e−7day. First-iteration prior centers are exactly equal; later centers differ by at most0.0002444760612day because they follow the slightly different preceding fits. These small final differences are evidence of stability to this initialization-grid change, not evidence that the required optimizer-entry multistart was performed. No paired NIR timing effect has been computed.

The native default grid is `[-4,-2,0,+2,+4]` days around the NML initializer. The nominal candidate MJDs are57703.80078125 through57711.80078125 in2day steps; the minus grid spans57701.80078125 through57709.80078125. They share four candidate coordinates, including the winning57707.80078125. The logs print `Adjust PKMJD=57705.801 ==> 57707.801` for each minus object and `57707.801 ==> 57707.801` for each nominal object. Full precision comes from the bound and entry records, not these rounded display lines.

Source mapping in the frozen restricted build:

- `snlc_fit.car:8360`: applies the user NML peak initializer.
- `snlc_fit.car:8424`: logs/asserts the domain and initializer before adjustment.
- `snlc_fit.car:7038–7040`: defines the default peak grid.
- `snlc_fit.car:9308–9368`: saves the initializer, evaluates the grid, and replaces INIVAL by the grid winner.
- `snlc_fit.car:9516`: even its fallback retains the selected peak grid value.
- `snlc_fit.car:2130–2139`: logs CSP_ENTRY after adjustment.
- `snana.car:9228`: passes INIVAL to MNFIT_DRIVER; `snana.car:36586` passes its value to MNPARM.

A literal post-grid `INIVAL(IPAR_PEAKMJD)+=±2` would have an additional scientific meaning. `FCNCHI2_PRIOR` at `snlc_fit.car:5055` subtracts this same common INIVAL from the trial peak, and `FITINI_PKMJDPRIOR:11078,11159` establishes the five-day prior's squared offset. Such a change moves both the optimizer start and the first-iteration prior center. It should not be described as a start-only check of the same objective. Later covariance iterations recenter INIVAL on the preceding fit, so the trajectory can converge to the same state, but that does not retroactively make the first objectives identical or establish a unique global Gaussian likelihood.

The recommended additive test is a private opt-in hook that changes only a local MNPARM entry value for the peak in the first covariance iteration. Leave common INIVAL, the prior center, all other initial parameters, bounds, covariance initialization, observations and subsequent native recentering unchanged. Log the native source value, requested offset, actual MNPARM value, prior center and fixed bounds separately. Check the finite shifted value against the existing fixed domain before MNPARM. Keep the original±2 entry-separation tolerance and every support/convergence gate; do not widen tolerances or choose an offset from outcomes.

Preserve the current `joint_minus` as an initializer-grid sensitivity, not a displaced-optimizer run. For a new private binary, absent/zero mode must exactly replay the already completed nominal `joint12` and `joint9` saved model/W/objective/FITRES states. Those existing scientific baselines can be reused after exact disabled-path identity is established. Any new ±2 start branches require a separately frozen patch, source/standalone review and root release; no new photons, reseeds or membership changes are needed. All existing11.420924458seconds of the new-estimator native budget remain carried forward.

This diagnosis ran no native process and changed no frozen inputs. The first diagnostic draft incorrectly assumed later recentered prior values would also be bitwise equal; its assertion failed before results were written. It is preserved as `initial_overstrict_prior_identity.py`; the final diagnostic explicitly distinguishes first and later prior centers without changing any native gate.
