# Saved joint12 final-iteration retry diagnosis

The `bounds iteration coverage` stop is caused by a checker assumption that is narrower than the unchanged native estimator: the native engine revisits its final covariance iteration once when MINUIT reports a parameter covariance status below 3. CID22 and CID64 each execute `1,...,12,12`; the other62 execute `1,...,12`. These are real calls with repeated preparation, covariance construction and minimization, not duplicate output lines.

A separate numerical warning remains material: CID22 and CID64 return MINUIT code4 on every iteration4–12 and on the retry. MINUIT documents code4 as abnormal termination, including nonconvergence. Exact repeated states do not establish peak stationarity or a global optimum. The original frozen checker/failure and all scientific thresholds remain untouched; this review is not a declaration that the joint stage passed.

## Source mechanism

All citations below refer to the frozen private start-only source in `restricted-peak-engineering/start-only-hook/build/src/`; exact excerpts and hashes are in `source-evidence.txt` and `input-hashes.json`.

- `snlc_fit.car:2324–2335`: `OPT_MNSTAT_COV>0`, final iteration, `MNSTAT_COV<3`, and not already repeating cause `ERRFLAG_REPEAT_ITER=-1`. The default option is1 (`7173`). The NML does not override it.
- `snana.car:9365–9369`: error−1 decrements the iteration number and sets `LREPEAT_ITER`. The main loop increments it again and repeats preparation/minimization. The `not already repeating` condition prevents a second covariance-status retry.
- `snlc_fit.car:8087–8131`: a repeated ITER12 gets starting values, step sizes and the recentered peak-prior center from ITER11. `FITINI_COV` also receives `FITVAL(:,ITER−1)` (`2066`). This retry is not an additional iteration13 or an independent displaced start.
- `snana.car:36684–36688`: the printed “MIGRAD returns” value is actually the return from a `MINIMIZE` command. `minuit.F:3696–3704,3782–3799` defines code4 and the MIGRAD/SIMPLEX fallback.
- `snana.car:36745–36746`: `MNSTAT_COV` is MINUIT's parameter-uncertainty status. It is distinct from the observation covariance C used in rᵀC⁻¹r. Status1 is an approximate, inaccurate matrix.
- `snana.car:36770`: the subsequent `EXIT` command overwrites the driver's IERR. Its normal return11 (`minuit.F:3978`) is not a successful-minimization code. All770 callbacks have EXIT11, including the20 MINIMIZE4 callbacks.

Final repeated `MNSTAT_COV` itself is not in this FITRES schema and is not separately printed: after the one permitted retry, that repeat condition is disabled. It is therefore not valid to infer an accurate final parameter Hessian from `ERRFLAG_FIT=0`, a positive observation C, or the absence of a second “Bad COV” line. Native DLMAGERR/PKMJDERR accuracy is not certified here. The planned Monte Carlo summaries use the empirical64-draw variation, not those parameter errors.

## Independent saved-state checks

`review.py` streams the existing log, preserves every occurrence, reconstructs all770 C matrices/objectives independently of the frozen stage checker, and writes a callback ledger. It performs no native execution.

| Quantity | Result |
|---|---:|
| Physical objects / callbacks |64 /770|
| Bound assertions / callbacks |770 /770|
| Native MINIMIZE return0 / return4 |750 /20|
| Physical mean calls, no support violations |4,263,658|
| Maximum absolute objective reconstruction error |4.690e−13|
| Minimum observation-C eigenvalue |0.4858107|
| Saved native process time |17.949420339 s|

The20 abnormal returns all belong to CID22/64, iterations4–12 and repeated12. Both initial final12 attempts explicitly print `MNSTAT_COV=1` before the retry. Their repeated final12 means, all row tokens, full W, full C, objective components, numerical entry values and bound tokens are exactly identical. Only the repeat flag changes. Their repeated final12 D/peak/C also equal ITER11 exactly.

CID22 final D=42.002880601977495 and peak=57707.804007839368; CID64 final D=42.005624834644372 and peak=57707.879802205491. These numbers are recorded only to identify the repeated native state, not as timing-effect results. Both first and repeated occurrences remain in the saved arrays/ledger. There has been no9-iteration comparison, displaced-MINUIT comparison, NIR intervention, paired effect estimate, refit, new photon or reseed in this review.

## Proposed additive checker correction, not executed

The bookkeeping correction is justified, provided it preserves occurrences instead of dropping duplicates. A future version can key records by `(CID, nominal_iteration, occurrence)` and require the exact source-supported pattern `1..N` with at most one extra finalN visit explicitly linked to a preceding covariance-status retry and `LREPEAT_ITER=T`. It must reject arbitrary repetitions, missing visits, unknown retry causes, extra objects and unmatched bound/objective/start records.

The bound-assertion count should equal actual preparation occurrences; nominal iteration coverage must still be exactly1..N. Apply every existing row, finite/SPD, objective, phase/domain and fixed-nuisance check to every occurrence. Preserve comparisons of the final retained state to the preceding distinct iteration and separately to its retry predecessor using the existing thresholds. The MNPARM/MNPOUT checker also needs occurrence keys, with requested starts applied only at ITER1 and exactly zero thereafter. Do not silently overwrite repeated keys. Full-null comparisons must retain exact multiplicity and order.

Record the MINIMIZE and EXIT returns separately for every occurrence and report the final parameter-covariance uncertainty explicitly. Such an additive correction would allow the already-prescribed9/12 and actual-start diagnostics to be evaluated under a reviewed release without repeating joint12 or changing the estimator. It would not erase the MINUIT warning or itself prove convergence. If independent point-estimator convergence remains unresolved after those diagnostics, withhold the64 paired summary; do not select another branch, discard CID22/64, refill or weaken thresholds. No new native execution is authorized by this note.
