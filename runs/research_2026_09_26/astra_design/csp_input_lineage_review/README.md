# Independent CSP magnitude/flux lineage and prepared-input review

**PASS.** This review independently reparses the original CSP DR3 archive and source RAISIN photometry, then reconstructs the proposed filter-token changes. It imports none of the root lineage or preparation code and performs no native fits or outcome scoring. The result validates the concrete processing intervention, not a physical calibration correction, fitted distance or cosmology.

## Archive to released data

- All **5,491 NIR rows** have identical decimal magnitudes/errors and multiplicities in `SN_photo.dat` and the per-SN SNooPy files after the two declared label mappings, `Jdw→J` and `Hdw→H`.
- Of the 76 LISTed RAISIN CSP photometry objects, **73 occur in the original archive**: 65 have NIR rows and eight have no NIR rows in either corresponding source/release representation. Presence is determined from the SNpy file roster, independently of whether a dictionary of NIR rows was populated.
- The eight without NIR are 2005bl, 2005bo, 2005ir, 2005w, 2006fw, 2006py, 2007mm and 2007ol.
- The three additions absent from this archive—2012fr, 2012ht and 2015F—remain unmatched. They are not falsely treated as zero-row exact reproductions.
- For every archive-present object, the released magnitude/error multiset and the printed flux/error values match. Fluxes were independently recomputed with 50-digit Decimal arithmetic from `F=10^[−0.4(m−27.5)]` and `σ_F=F(10^[0.4σ_m]−1)`, then compared at the converter's printed precision. No imported converter function is used.
- Metadata alone yields **303 released rows uniquely attributable to raw Jdw** and **288 to raw Hdw**. Physical-label candidates use object, mapped band and time; magnitude values never resolve a mixed-instrument match. Duplicates are retained.

The first two root reporting revisions remain intact. Their presence-detection and no-NIR reporting errors do not recur in this independent reconstruction.

## Prepared intervention

All **163 input hashes** recorded by the preparation protocol verify. The full cohort has exactly the same 42 CSP objects as the selected 79-object published RAISIN cohort.

For every object, the review starts from original source bytes and identifies eligible rows directly from the original archive metadata. It replaces only the single `J` byte with `j` in uniquely attributable Jdw observations. The expected whole file is then compared against each prepared arm:

- `nominal` and `nominal_copy` are byte-identical to source for all objects.
- `known_Jdw_to_j` differs by **188 filter-token bytes across the full cohort**, or **six in the two-object pilot**. All other bytes, including headers, fluxes, errors, magnitudes, times and row order, are unchanged.
- All **ten unaffected controls** are byte-identical in the changed arm.
- All **four ambiguous 2006kf observations** remain `J`: source lines 151/152/170/171 at MJD 54037.810/.820. Their source metadata admits both J and Jdw, so no arbitrary assignment is made.
- Within each scope, prepared fit NMLs differ only in the explicit data-arm path. The independently reconstructed changed-row ledger exactly matches the saved ledger.

This does not certify native accepted masks, rest-filter/K-correction selection, optimizer convergence, uncertainty estimates, source-instrument interpretation, or the historical execution. Those remain separate gates. In particular, changing an observer-filter identifier can alter the native model's operations; byte-local input changes alone do not prove a physically isolated output response.

## Reproduce and inspect

Run from the repository root:

```sh
python runs/research_2026_09_26/astra_design/csp_input_lineage_review/verify.py
```

`result.json` contains the verdict and counts. `independent-objects.csv`, `independent-physical-rows.csv`, `independent-changed-tokens.csv`, `independent-file-gates.csv` and `independent-ambiguous-unchanged.csv` retain the joins and exact changes. `input-manifest.json` hashes every directly consumed file. `manifest.json` hashes this review's closed artifacts.
