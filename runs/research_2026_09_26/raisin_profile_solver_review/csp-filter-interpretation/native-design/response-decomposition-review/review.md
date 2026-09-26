# Independent review of the mean/covariance bookkeeping

**PASS, all42 objects.** Recalculation from the final native nominal and changed means, measured fluxes and full covariance reproduces every reported term to **9.63e−16 mag**, and all aggregate means to better than1e−15mag. No native fit or model probe was made. The root protocol was fixed before the full42 changed-filter outputs, but after the two-object pilot; this is an explanatory decomposition, not a new discovery test.

| Cohort | Native total ΔD | Changed mean at fixed old covariance | Remaining term |
|---|---:|---:|---:|
| All42 | +0.483432007 mmag | +0.634403783 mmag | −0.150971776 mmag |
| Metadata-affected32 | +0.634504509 mmag | +0.832655303 mmag | −0.198150794 mmag |

The independently checked convention is m1,0=m1 exp[K(D1−D0)], K=ln(10)/2.5. Thus m1,0 is the changed-filter mean brought back to the nominal distance. With u=L0⁻¹m1,0 and v=L0⁻¹y, C0=L0L0ᵀ, the positive best amplitude is a0=(uᵀv)/(uᵀu) and ΔDmean=−ln(a0)/K. The last column is defined by subtraction ΔDnative−ΔDmean. Applying the same amplitude check with the changed final covariance meets the original0.001-mag optimum-gap gate for every object.

This implementation uses this agent's existing native parser and pre-fit raw-row map. It constructs a new accepted-row multiset alignment, reverses occurrence assignment relative to the root's forward queues, and whitens with the covariance Cholesky factor instead of evaluating the root's direct precision-matrix quadratic products. It requires that all measured MJD/flux/error coordinates match and that the resulting permutation is complete. Duplicate multiplicity is retained; mixed-target partial copies would fail, and covariance symmetry under duplicate swaps is required if necessary. There are **no identical (MJD,flux,error,old-band) duplicate groups in these final42 masks**, so no tested duplicate ambiguity is claimed. This does not say that the upstream release is globally duplicate-free.

The ten controls have **exactly zero native total response**. Their two component terms are generally small, equal and opposite: the largest absolute component is9.96842e−9mag. These equal the nominal model's independently evaluated fixed-C analytic optimum gap exactly in this implementation. They arise because the native finite numerical solution need not be the exact analytic amplitude optimum. Calling them arbitrary arithmetic roundoff alone is too strong; native stopping/precision and exported means can contribute. They are negligible relative to the declared gate, but are preserved, not zeroed.

For clarity, an alternative analytic-to-analytic definition would subtract that nominal optimum gap from ΔDmean; all ten controls are exactly zero under that alternative in the independent calculation. That would change the prescribed component definition, so it is **only a notation check** saved in the independent ledger, not a replacement for the primary reported decomposition.

The remaining term includes the path from fixed old covariance to the converged changed native recipe and numerical optimum residue. It is order dependent; reversing the path or choosing another intermediate covariance changes the components. It is not an independent physical bias, dust or calibration attribution. Changing the weighting can cancel much of the mean-only response for an individual object (for example2008hu), without making either component an error bar or a likelihood-ratio result. Source covariance, training and selection assumptions are still inherited.

The root's figure was visually checked here: filled points show total native shifts, open points the fixed-old-covariance term, and diamonds the unchanged controls. Connecting lines are accounting differences, not uncertainty intervals. The axes, mmag conversion and42-object mean agree with the recomputation.

![All42 processing response](/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/csp_native_filter_response/root-state-review/figure/csp-filter-response.png)

Artifacts: this directory's `check.py`, `result.json`, `check.log` and manifests. The original root results and protocol remain under `csp_native_filter_response/root-state-review/response-decomposition/` and `response-decomposition-protocol.json`, SHA256 `ba0dc23d96b63dd3660594b88adc959fb1405d885f870685e1e3f6936fc1f959`. The [main stable response review](../stable-filter-response/review.md) and [identification/priority review](../interpretation-next-steps/review.md) retain the broader scientific limitations.
