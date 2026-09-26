# Independent design review: native peak restricted to tabulated support

The proposed restriction is scientifically defensible as a **new conditional estimator**, not a repair of the failed original result. The original twelve-iteration joint branch remains failed: CID4 made 304 unsupported mean calls and CID7 made147; final supported states cannot erase those evaluations. No fits, builds, photometric scores or distance outcomes were produced or inspected in this review. `review_metadata.py` is independent of the experiment executor and reads only epoch/band/redshift/input-peak metadata plus structural support records and source files.

## Definition and independent metadata check

For fixed observation time t_i and fixed heliocentric redshift z, the native rest phase is (t_i−τ)/(1+z). Requiring phase a_i≤phase≤b_i implies

    τ ∈ [t_i−(1+z)b_i, t_i−(1+z)a_i].

Intersect these intervals over **every fixed epoch**, including every band involved in the observer/rest transformation. Do not recompute the domain from a flux-selected or parameter-dependent accepted subset. With the active SNooPy grid [−20,+70] days and KCOR [−20,+85], the common phase range is [−20,+70]. Shape, redshift, model and filters stay fixed under this bounded experiment. This establishes support of the tabulated numerical operator; it does not establish empirical accuracy of the SN model throughout that range.

All eight HEAD/PHOT objects have the same117 epoch/band multiset as the frozen SIMLIB. Their actual FITS1D times span57686.007 to57773.053; fixed zHEL is0.453000009059906. Consequently the common absolute peak interval is

    [57671.34299936581, 57715.067000181196] MJD,

with width43.724000815389445 days. The lower endpoint is set by the last z-band epoch and the upper endpoint by the first g-band epoch. Relative to the declared engineering center it spans[−36.457781884,+7.266218931] days; these offsets describe the interval but do not define it. All three declared starts and their ±4-day initialization grids remain inside. The original final |peak−truth|≤4-day engineering gate remains a separate unchanged gate, as explicitly directed by root; it is neither the fitting bound nor proof of support. A failure retains the object and stops the branch rather than selecting a survivor cohort.

Native precision matters. FCNSNLC uses R8EP_MJD, populated directly from SNLC8_MJD, and evaluates `(MJD−MJDOFF)−peak_offset` in REAL*8. The displayed/R4 epoch array is not its input. LOAD_EPALL separately computes a REAL*4 time difference, redshift factor and phase for masks. `metadata-result.json` records both routes at raw endpoints and one-REAL*8-step inward endpoints; all are supported for these inputs. An implementation must store the absolute domain, convert consistently to native offset coordinates, and assert its constancy at every iteration. If endpoint rounding requires an inward move, use the smallest documented representable move that passes actual arithmetic, not an arbitrary scientific margin. No timestamp replacement or clipping is required here.

## Required source semantics

The inspected native source initializes peak INIBND to initial±100 days in FITINI_PARVAL, passes INIBND directly to MINUIT MNPARM, and normally rejects outside-bound parameters through FCNCHI2_PRIOR before the mean loop. SIGMA_ONLY bypasses that prior; diagnostic calls can also bypass the large-prior return. Therefore an INIBND change alone does not prove the all-call property.

A separate default-off implementation can replace the first-iteration peak bounds with this metadata-derived interval before initialization-grid or covariance mean calls. The absolute bounds must not follow the previous fitted peak, its error, the moving Gaussian-prior center, or subsequent phase masks. On later iterations assert and retain the same bounds. Keep the existing peak prior and covariance updates unchanged so this is an isolated domain change. FITINI_EPVAR itself initializes arrays without invoking model means; FITINI_COV and FITINI_ADJUST do invoke downstream evaluation routes and need the initialized domain in place.

Add a premodel refusal on every USRFUN path, before the rest-frame template or KCOR calls. An unsafe call must leave an explicit failure diagnostic and stop; do not silently clip the phase, clamp a parameter, replace the model value, skip its epoch, or count a rejected call as a supported one. The guard must cover direct, covariance and SIGMA_ONLY routes as well as MINUIT. Existing all-call finite/model-phase/shape/redshift tracing remains enabled. The initialization grid is already supported for this cohort, so no grid pruning is needed. Plots remain disabled.

The disabled hook must reproduce the original binary's science outputs exactly on the declared identity control; active noiseless runs must recover the same supported input signal under existing gates. Hash the patch and private source/build; retain the failed branch. The physical interval is broad, so passing support does not certify numerical uniqueness. Keep the existing12-versus9 iteration, actual-start, covariance, state, stationarity and fixed-row gates rather than replacing them with ERRFLAG=0.

## Interpretation of a bounded outcome

Report the minimum distance to each bound and any active-bound flag. A boundary optimum obeys a one-sided/KKT condition, not an interior zero derivative; Hessian-based peak uncertainty can be misleading there. Do not discard boundary draws or report a timing-bias average from successful survivors. No new inference gate is proposed here beyond reporting and the existing all-eight conditions.

The pure-peak comparison remains conditional on the same generated photons, same fixed shape/AV/RV, same measurements/errors, same native prior/covariance recipe and same NIR epoch set. Its estimated joint peak remains correlated with the NIR noise used again in the distance fit; replacing that peak with independent Gaussian jitter would change the experiment. Converting the estimated peak into the NIR input's REAL*4 header is a separate declared quantization operation; verify the resulting NIR model calls and preserve it identically in null adapters. The domain restriction itself does not identify a real-population distance bias, establish survey selection closure, or justify a cosmology correction.

`metadata-result.json` hashes13 inputs, records all eight metadata checks and the exact original failure support counts. Source excerpts underlying the derivation are in snlc_fit.car (FCNSNLC/FCNCHI2_PRIOR/USRFUN/FITINI_PARVAL/LOAD_EPALL) and snana.car (MNFIT_DRIVER), under the hashed private support build. Pending implementation review will be added separately; this note does not authorize execution.
