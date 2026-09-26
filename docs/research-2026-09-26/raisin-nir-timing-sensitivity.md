# Ten DES16 NIR distances under fixed header-date changes

The exact historical-source baseline reproduced all ten archived author NIR-only `DLMAG` values at printed precision, with all ten `NDOF` values equal. This permits a narrow conditional sensitivity test: change the photometry `PEAKMJD` header while retaining the same released J/H measurements and fixed SNooPy shape/AV settings. It does **not** establish that any released header date is wrong or provide a distance correction.

The [original protocol](../../runs/research_2026_09_26/raisin_nir_timing_sensitivity/protocol.json) was frozen before a fit (SHA-256 `3cb9ed02a920ddda6dd9742c607d2a49f672729c9a8b5e2d36c1d0fffcb14aeb`). The first baseline process aborted before producing FITRES, in `check_file_docana` while reading the absolute VPEC override path. The failed NML, log and execution record remain under `fits/baseline/`. A [prescore technical amendment](../../runs/research_2026_09_26/raisin_nir_timing_sensitivity/protocol-shortvpec-amendment.md) used a work-local, byte-identical copy of that same VPEC file; its [frozen protocol](../../runs/research_2026_09_26/raisin_nir_timing_sensitivity/protocol-shortvpec.json) has SHA-256 `e42ceec5781fa2ae699ba51ac73d1a6cfe2b0c5d72b9e2bbf6ea5cebc7e489dd`. No fit setting, photometry, source, KCOR or model was changed by that repair. SNANA v11_03c commit `06f2ccfdbf99d62c23dec66d0bf7b7e9452604d2` was built in an isolated clone using project-local compiler/libraries. The build-generated HBOOK/ROOT preprocessor diff, build log and binary hashes are recorded. This is source/version matching, not the original historical executable.

Both independent baseline directories produced byte-identical FITRES and LCPLOT files. The [baseline gate](../../runs/research_2026_09_26/raisin_nir_timing_sensitivity/baseline-gate.json) passed: ten `ERRFLAG_FIT=0` records, fixed peak/stretch/AV/RV, exact author DLMAG and NDOF at printed precision, and `NDOF = accepted J/H epochs − 1` for every object. Only then were the five timing variants run. An [independent input check](../../runs/research_2026_09_26/raisin_nir_timing_sensitivity/independent-verification.json) found that each of the 50 shifted object files differs from its released source only in the two matching `PEAKMJD` header lines. All 70 baseline/variant fitted peak values follow the header after single-precision conversion and four-decimal FITRES printing. The two baseline directories have identical photometry bytes.

| Imposed header time relative to released | Mean fitted ΔDLMAG (mag) | Largest absolute object Δ (mag) | Changed accepted J/H masks |
|---|---:|---:|---:|
| −1 day | +0.01305 | 0.03790 | 0/10 |
| −0.5 day | +0.00598 | 0.01886 | 0/10 |
| +0.5 day | −0.00546 | 0.01992 | 0/10 |
| +1 day | −0.01025 | 0.03992 | 0/10 |
| Author `raisin_t0` values (−1.13 to +1.48 days) | +0.00500 | 0.04284 | 0/10 |

The largest absolute response in each condition is DES16S1bno, which has three accepted J/H epochs at baseline. All 50 variant fits return `ERRFLAG_FIT=0`, and no accepted-row multiset or NDOF changes in this specified range. Exact per-object peaks, distances, data `FITCHI2`, errors, and mask comparisons are in the [ledger](../../runs/research_2026_09_26/raisin_nir_timing_sensitivity/per-object-response.csv) and [result](../../runs/research_2026_09_26/raisin_nir_timing_sensitivity/result.json). The fitted-distance changes are conditional finite responses, not a common-likelihood optimization proof across shifted model phases or a physical bias estimate. No amplitude-start variation was performed; the source can reinitialize distance, and this was declared before fits.

The available simulation path provides a separate timing provenance clue. The pinned DES simulation input sets `USE_SIMLIB_PEAKMJD:1` and `GENSIGMA_SEARCH_PEAKMJD:0.01`; v11_03c `snlc_sim.c:1745–1748` parses this into the peak-smear scale, `:7872–7876` chooses the Gaussian branch, `:13129–13153` forms true simulated peak plus a Gaussian draw, and `:20685–20688` stores it as `SEARCH_PEAKMJD`. Its text writer emits that value as `PEAKMJD` (`sntools_dataformat_text.c:116–117`), and the reader maps that header back to `SEARCH_PEAKMJD` (`:1847–1849`). The two newly acquired author simulation NMLs fix peak, shape and AV and fit J/H only; their blob hashes and source checks are in the [timing-source review](../../runs/research_2026_09_26/raisin_nir_timing_sensitivity/timing-source-review.json). This is an available simulated-header recipe, not a demonstrated creator of the **observed** DES16 photometry headers or a certified effective bias-correction execution. The author `raisin_t0` file and the observed header peaks are distinct for all ten objects.

Optical sign restoration could affect the NIR fit through a changed upstream observed header estimate, but the header-producing operation and actual timing changes are not established. Optical/combined fit shape and AV are not imported into this fixed NIR branch. Detection/follow-up and simulation bias selection also require their own execution-linked rerun; none was performed here.

A separate [root verifier](../../scripts/research_2026_09_26/verify_raisin_nir_timing.py)
imports none of the execution/summary code. It verifies 46 input/output hash
entries, independently parses all 70 native records, checks all 70 input
header-only comparisons and unchanged epoch multisets, and reproduces every
condition's summary to `1e-13` or better. Its
[result](../../runs/research_2026_09_26/raisin_nir_timing_sensitivity/root-independent-review/result.json)
and [finite response differences](../../runs/research_2026_09_26/raisin_nir_timing_sensitivity/root-independent-review/timing-response-differences.csv)
are preserved. Those derivatives describe timing sensitivity only; they are
not assigned a timing-error distribution or converted into a new error budget.
