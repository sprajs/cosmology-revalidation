# Optical V labels follow the physical filter dates

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

The raw `SN_photo.dat` label `V0` is misleading as a universal physical-filter identifier, but the per-SN SNooPy files resolve it correctly by observation date. **All 2,528 rows match the documented date-specific filter choice:** 760 retain V0/LC3014 and 1,768 become V/LC9844. The tagged SNANA converter and released KCOR input preserve these distinctions. This check does not identify an optical analogue of the J WIRC → RC1 ambiguity, and no optical correction or new fit is warranted from this label observation.

A separate completeness gap remains: **21 published LC3009/V1 observations on 12 objects are absent from both formats in the DR3 convenience archive.** This is a source-coverage finding, not a fitted bias estimate.

## Source interpretation

The primary [Krisciunas et al. 2017 paper, §6.1.1](https://pure.au.dk/ws/files/118728791/Krisciunas_2017_AJ_154_211.pdf) describes LC3014 breaking on January 14, 2006, brief use of LC3009, then LC9844 from January 25. LC3014 and LC9844 use the same tertiary-star natural magnitudes because their stellar color terms agree to approximately .002; LC3009 receives a separate reduction. This does **not** equate the two physical transmission functions or require applying a stellar color correction directly to supernova magnitudes. Its Table 9 contains the intermediate-filter photometry. The [official CSP filter page](https://csp.obs.carnegiescience.edu/data/filters) and archive README give the same dates.

The institutional PDF was acquired directly, 4,221,158 bytes, with hash and headers in `paper-pdf-acquisition.json`. The failed arXiv HTML request is preserved separately. The 2017 and 2020 errata were not reinterpreted by this bounded optical audit.

## Metadata-only date test

`date-protocol.json` was saved before the new date comparison. It fixes the civil-UTC boundaries MJD 53749 and 53760 from the documented dates, with a separate literal-integer-JD sensitivity. The parser reads only object, filter and time tokens from the immutable 134-object DR3 archive; it never parses the photometric magnitude/error fields.

| Packaged SNooPy label | Rows | MJD range | Physical source designation |
|---|---:|---|---|
| V0 | 760 | 53249.79–53748.77 | LC3014 |
| V | 1,768 | 53760.64–55305.53 | LC9844 |
| V1 | 0 | none | LC3009 |

Every raw V0 object/time tuple has exactly one V/V0/V1 candidate in the per-SN file with the same multiplicity. There are **zero civil-date mismatches**. The literal rounded-JD reading falsely assigns four January 13 observations to the intermediate filter; those four are retained in the sensitivity ledger rather than hidden by changing the time convention. No observation occupies the intermediate civil-date interval.

This independently verifies the metadata-specific optical lineage refinement; the separate magnitude/error equality result is not used to select this date rule. The original failed universal V0→V hypothesis remains preserved under `runs/research_2026_09_26/csp_optical_magnitude_lineage/`.

## Registry and converter mapping

The official SNooPy source snapshots have the following registrations:

| SNooPy label | Registry throughput | Tagged SNANA label | Released CSP KCOR throughput |
|---|---|---|---|
| V0 | `V_LC3014_tel_ccd_atm_ext_1.2.dat` | m | LC3014 |
| V1 | `V_LC3009_tel_ccd_atm_ext_1.2.dat` | o | LC3009 |
| V | `V_tel_ccd_atm_ext_1.2.dat` | n | LC9844 |

The generic V file was checked, not identified from its name: **its complete numeric array is exactly equal to `V_LC9844_tel_ccd_atm_ext_1.2.dat`** in the locally acquired 2018 CSP throughput collection. This holds both for the pre-2017 snapshot `5576598c4603c48fe195c4901f572671feed785f` and the current pinned snapshot `3b97ef08e48b4d535f587d19a14eb2d0ea3a7f33`. All four corresponding V0/V1 historical/current array comparisons also match exactly. No synthetic SN or stellar photometry was evaluated.

The pre-2017 repository tree stores V0/V1 curves under a different top-level path from its generic V curve. The initial 404 and the exact-tree-path correction are recorded; this is a source-asset comparison, not proof that the old snapshot was a fully installed/executed package. The 2018-01-08 registry change at `c95928e059145b1f7f8ea88985b49a23d8db79b0` changes the V1 Vega reference magnitude from .0096 to .015 while leaving the V/V0 curve registrations intact. No V1 observations occur in this archive, so that fact does not imply a correction to these rows.

Relevant local execution inputs are `sources/repos/RickKessler__SNANA@v11_04k/util/translate_CSPDR3.py:75` and `sources/repos/djones1040__RAISIN_DataRelease@a383c4b/kcor/kcor_CSPDR3_BD17.input:16`. These establish the saved mapping; they do not prove a particular historical preprocessing run beyond the separately verified released-data joins.

## Published intermediate-filter data missing from the archive

The source-only parser extracts the 21 dates and object names from published Table 9. A same-object V/V0/V1 join within .01 day finds **zero matches** in either `SN_photo.dat` or the individual SNpy files. All dates lie in the interval where the archive has no V observations. The ledger records JD and the explicit `MJD=JD−2400000.5` conversion. This check reads no Table 9 magnitudes/errors and does not restore any measurements.

This limits the README claim that the convenience file concatenates all of Tables 7–12. It does not show why those observations were omitted, whether they passed the analysis selection, or their impact on a fitted distance. Any restoration experiment would require a separately frozen source/magnitude/error/quality bridge, the LC3009 calibration convention, and native baseline closure. The current audit stops before that step.

## Evidence and scope

Owned evidence is under `runs/research_2026_09_26/astra_design/csp_optical_filter_identification/`: source PDF/HTML, source-acquisition hashes, `date-label-ledger.csv`, `date-result.json`, six pinned curve files, `registry-comparison.json`, and `table9-missing-metadata.csv` / `table9-result.json`. `manifest.json` covers the closed sources, outputs, this report and the read-only shared source snapshots.

This source audit resolves a label concern while exposing a smaller completeness question. It makes no statement about gray luminosity, dust population, extinction correction, distance bias or cosmology. No native fit, observed-flux score, package change, or alteration of existing files was performed.
