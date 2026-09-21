# Saved forward-comparison checkpoint

Primary work complete: nine simulated selection stages, paired generated-attempt and LIBID-cluster uncertainty; actual measured-fit conditional forecasts on817train/202test DES objects; all9 physical/noise arms;300 joint bootstrap replicates; actual Hessian observables; preregistration and six numerical falsifiers. Report: `docs/phase2/forward-discrimination.md`. No cosmology/physical population winner inferred.

Exact commands, repository root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/phase2/forward_discrimination/test_numerics.py
.venv/bin/python scripts/phase2/forward_discrimination/selection.py --require-fits
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/phase2/forward_discrimination/compare.py --bootstrap 300 --sensitivities
```

The last command reruns the deterministic primary and runs all seven registered sensitivity arms. Each primary run took well under a minute on this host; the complete sensitivity suite was not launched before the user-requested pause. It intentionally uses300bootstrap replicates throughout; do not silently present a smaller pilot as the registered result.

Pending before strong interpretation:

1. Run registered sensitivities and compare on identical original test CIDs within each common-support cohort; report changed support counts, especially strict cuts and narrower bandwidth. Across-arm pairs within a run already share exact rows. Across-sensitivity score deltas require their common intersection.
   Final audit identified a field-provenance discrepancy: real field is frozen fold table first-PHOT-epoch field; simulated field is FITRES. Real compound FITRES labels can therefore enter under their simple first-epoch field, contrary to the preregistration assumption. Preserve primary; add a separately named consistency sensitivity using first-PHOT-epoch fields on both sides or excluding real multi-field objects. This was identified after outcomes and must be labelled accordingly.
2. Check the classifier reconstruction agent's validation. Only if real probabilities reproduce can a newly preregistered matched-SNN simulation selection sensitivity close part of that mismatch. HEAD and fitted simulated CIDs are unique; all-generated CIDs are reused, so retain generated_attempt_index.
3. Ask the independent audit to inspect energy-score implementation, weighted-empirical V-statistic finite-simulation effects, bandwidth approximation, nuisance alignment and bootstrap. The numerical tests verify formula/invariance properties but not physical fidelity or exact conditional calibration. The unadjusted many-diagnostic intervals do not establish a physical winner.
4. Full BBC/systematics membership, nonIa mixture, shared calibration/model uncertainty, independent-seed coverage and W22 SN_age are still absent. These counts do not normalize the exact real cohort or arbitrary continuous population alternatives.
5. The paired selection contrast is a clear actionable empirical result: +0.10mag reduces quality count1314→1193; true noise×1.2 changes written2678→2682 but quality1314→1227. Preserve source hashes and stage distinctions.

Files owned by this workstream: `scripts/phase2/forward_discrimination/`, `phase2/forward_discrimination/`, `docs/phase2/forward-discrimination.md`. No phase1 files or other workstream outputs changed. Primary process completed normally; no long-running forward-comparison process remains.
