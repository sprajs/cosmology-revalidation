# Observation-level DES data audit

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

This audit examines every released real-data SNANA HEAD and PHOT row in the acquired DES-SN5YR tag 1.3, the acquired Dovekie release, and the public DES difference-imaging release. It operates below fitted distances. **These are calibrated flux measurements, not detector pixels or uncalibrated images.** No source was edited, no cosmological conclusion was fitted here, and no metadata value was silently imputed.

The local [preregistration](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/preregistration.json) was written before computing the row-audit outcomes. It is a prospective local register, not an external preregistration. It acknowledged prior release-level knowledge. Separate follow-up registers explicitly disclose the duplicate findings known before designing their proposed refits.

## Versions and physical coverage

The original release is `des-science/DES-SN5YR` tag 1.3, commit `e3493cb3b9505fc1f3f392d887364b50dae21439`. The acquired Dovekie repository is commit `c9a4fcafc4cbd19bd750dee47fc76194a45c181f`. The public DIFFIMG files retain the previously acquired Zenodo metadata and checksums. [Input hashes](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/inputs_manifest.json), [output hashes](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/outputs_manifest.json), and the [machine-readable summary](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/summary.json) identify the exact bytes.

| Physical photometry table | Objects | PHOT rows including delimiters | Actual observation rows | HEAD columns |
|---|---:|---:|---:|---:|
| DES SMP | 19,706 | 1,798,736 | 1,779,030 | 96 original; 99 current |
| Foundation | 185 | 6,195 | 6,010 | 95 |
| Combined LOWZ | 342 | 36,974 | 36,632 | 246 |
| DES DIFFIMG | 31,636 | 18,106,216 | 18,074,580 | 122 |

Every PHOT table has 20 columns. The distinct files contain 19,896,252 non-delimiter measurement records. This is a count of records in different reductions, **not** independent exposures or independent supernovae. SMP and DIFFIMG use the same observing programme. Current/original releases likewise must not be counted twice.

All three PHOT FITS files are byte-identical between original and current releases; Foundation and LOWZ HEAD files are identical too. Current DES HEAD renames `IAUC` to `NAME_IAUC` and adds `NAME_TRANSIENT`, `LENSDMU`, and `LENSDMU_ERR`. Its only changed shared columns are the four host magnitudes and their errors, involving approximately 16,704 rows. Host mass, redshift, and source pointers are unchanged. A different calibration model or fitting input can change inferred SALT parameters without changing the stored flux table. Conversely, changed host magnitudes are not proof that corresponding stored stellar masses were re-inferred.

## Integrity, missingness, and flux-quality findings

All seven version/dataset combinations have unique SNIDs, valid 1-based inclusive pointers, matching NOBS and pointer lengths, non-overlapping object spans, and the expected MJD=-777 row after every object. No non-delimiter row is unowned. Every non-delimiter measurement has finite positive MJD, finite calibrated flux, finite positive flux uncertainty, and a nonblank band. A full independent readback compared **every physical HDF5 column and row** with its original FITS column; [verification passed](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/verification.json). Repeated PHOT groups are HDF5 hard links to the same verified arrays.

These checks do not establish that every measurement is scientifically usable. DES SMP contains 517,801 negative-flux records, DIFFIMG 8,486,366, Foundation 137, and LOWZ none. Negative flux is retained in a signed-flux likelihood. An entirely positive LOWZ release must not be treated as if it sampled nondetections identically to DES. The maximum full-catalogue SMP SNR is extreme (~2.5 million), so numerical validity alone cannot replace sample/type and quality selection. The [field profiles](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/field_profiles.csv.gz) count blanks, nonfinite values and common sentinels for every actual field; sentinel counts are descriptive and do not classify every negative number as missing.

Among all 19,706 DES headers, 11,398 have the -999 heliocentric-redshift sentinel and 8,643 have the -999 host-mass sentinel. Host colour uncertainty is unavailable for **all** objects: 16,704 entries are -999 and 3,002 are -9999. Host SFR, specific SFR and colour also have substantial sentinel populations. The host-mass-error documentation identifies statistical-only uncertainties; they are not full stellar-population modelling errors. A generative host/colour model cannot manufacture missing joint PDFs from these columns.

DES SNANA epoch cuts use PSF **FWHM in arcseconds**, computed as `PSF_SIG1 * PIXSIZE * 2.355`; directly cutting `PSF_SIG1` in pixels gives the wrong answer. Electron zero-point is `ZEROPT + 2.5 log10(GAIN)`, with the source's treatment of extremely small gain. This was checked in the pinned SNANA source. In SMP, 5,088 rows fail the .5–2.75 arcsec PSF interval and 2,285 fail the electron-zero-point interval. The initial in-progress audit's PSF-unit mistake was corrected, logged, and the complete audit rerun; no initial count is treated as a result.

DIFFIMG includes 4,307,599 rows with one of its documented error/quality bits (8,16,32,64,128,256,512,1024,2048). That descriptive mask is deliberately not imposed on SMP: the published SMP fit accepts many observations with bits 1024 and 2048. Band/field/flag histograms and signed-SNR distributions are in the machine summary. Per-object `27.5 - 2.5 log10(5*FLUXCALERR)` summaries are labelled **error-derived depth proxies**; they are not independently measured detection efficiencies.

## Selection membership and the public fit mask

Literal IDs connect every released Hubble-sample object to photometry: original 1,635 DES + 118 Foundation + 76 LOWZ = 1,829; current 1,623 DES + 117 Foundation + 80 LOWZ = 1,820. The current classification table has 17,733 objects, while original nominal probabilities cover 1,635. The original no-redshift classifier table has 3,547 rows, 1,322 matching the original Hubble table. These are different conditioning sets, not interchangeable classifications or extra independent samples.

The partial DES cut reconstruction checks available redshift, redshift error, MW reddening, epoch count, PSF/zero-point, two-band SNR and approximate rest-phase coverage. It retains 5,144 original or 5,145 current DES objects and includes every corresponding published DES Hubble member. It uses published fitted PKMJD for matched objects and a header peak estimate otherwise. Therefore it is a **diagnostic proxy**, not an exact eligibility classifier or a selection likelihood. Low-z proxy columns use DES thresholds and are explicitly inapplicable to low-z scientific selection. Rest-wavelength cuts, original iterative outlier handling, classifier conditioning, BBC validity, exact unobserved denominators and complete survey efficiencies require separate reconstruction.

The published original LCPLOT includes 62,702 accepted and 15,231 excluded data rows. All match the public calibrated FITS when one accounts for its single-precision MJD serialization: cast original MJD to float32 before the public three-decimal formatting. Matching the original double-precision times with a naive half-milliday tolerance falsely misses many valid rows.

[Flag matching](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/published_lcplot_flags.json) finds no bits 8,16,32,64,128,256 or 512 in unambiguously matched accepted observations, but does find bit 2 (612), bit 1024 (2,526), and bit 2048 (7,982). Ninety-four accepted records have ambiguous source flags from repeated flux records; two allow a bit-32 or a non-bit-32 source candidate. This caveat prevents a stronger universal claim. The official reproduction team used the bit-32 evidence to test an omitted historical mask.

For pilot CIDs 1442085, 1896213, 1289306 and the discrepant 1256426, every public LCPLOT flux agrees with public PHOT to floating-point/formatting precision, and uncertainty differences are below 5e-5 fraction. There is no evidence of a missing per-point flux correction in those examples. Five extra pilot epochs absent from public LCPLOT all carry bit 32. For 1256426 after that mask, the available 27 rows match but iterative clipping chooses six instead of seven accepted points, with seven acceptance switches and very different residual chi-squares. This is a fitting/rejection-baseline discrepancy, not an established photometric correction or a cosmology result. See [targeted checks](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/lcplot_audit.json) and [the mask32 comparison](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/lcplot_raw_1256426_mask32.csv).

## Repeated observations: evidence and limits

The DES SMP file contains 10,026 rows in exact CID/MJD/band repetition groups, affecting 1,566 catalogue objects. Original Hubble members include 118 affected objects and 724 such records; current Hubble members include 115 and 720. DIFFIMG has no exact CID/MJD/band repetitions, although image-ID keys are not universally unique there.

In the public original fit, 94 repeated raw groups on 68 objects contain multiple accepted records: 374 accepted records among 376 raw records in those groups. This reaches the fit input and deserves a controlled sensitivity test. However, a full-column audit prevents treating every four-row group as four copies of one identical image measurement. Most groups have the pattern **A,B,A,B**: A and B differ in positions by tens of pixels, PSF, sky noise and often gain or flags, despite the shared image/time/band key. Some also differ in flux. Collapsing A and B without resolving exposure metadata would overstate what is known.

Within those groups, full equality of **all 20 photometry columns**, checked separately within CID, identifies 4,328 exact duplicate pairs (8,656 rows) on 1,550 catalogue objects. There are 168 exact pairs in groups containing accepted public measurements; source-to-LCPLOT ambiguity means this is not itself an exact accepted-pair count. [Group identities](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/duplicate_identity.json), [exact row mapping](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/exact_duplicate_rows.csv.gz), and [accepted-group mapping](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/duplicate_exposure_published_acceptance.csv.gz) preserve the evidence.

The [targeted refit register](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/duplicate_refit_preregistration.json) and [conservative amendment](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/duplicate_preregistration_amendment1.json) specify nominal-versus-dependence sensitivities after the common baseline gate. The conservative arm removes exact full-row repeats while preserving distinct A/B records. Broader same-key correlation alternatives remain hypothetical. No rows have yet been changed here, and no claim about an acceleration bias follows until fit, selection and cosmological covariance effects are quantified.

## Fitted parameters and feasible hierarchy

The original metadata has x0–x1 covariance rounded to exactly zero in 1,825/1,829 objects and x0–colour covariance zero in 1,828/1,829. Combining its quoted marginal errors and cross-columns yields five non-positive-semidefinite correlation matrices. The current metadata preserves small cross-terms but yields eleven incompatible matrices by the same construction. These are mathematical diagnostics of this attempted Gaussian representation; they are not proof of an invalid published likelihood.

The release enables MINOS. The pinned SNANA `FIX_COVAR_LCFIT` implementation explicitly discusses mixing MINOS marginal errors with parabolic MNEMAT cross-covariances and provides optional repair conventions. Replacing these matrices silently with nearest-PSD matrices would add an unstated modelling choice. The [covariance archive](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/fitted_parameter_covariances.npz) preserves unrepaired x0/x1/c matrices, their Jacobian transform to mB/x1/c and PSD flags. Independent re-fitting is needed for an internally coherent full measurement covariance.

The public data support a direct calibrated-flux likelihood with signed fluxes, passbands, redshift/MW information, explicit population assumptions, quality selection and an independently reconstructed fit covariance. They also support conditional alternative colour-law prediction tests. They do **not** directly provide detector-pixel calibration likelihoods, full host SED/age posteriors for all objects, progenitor ages, an unconstrained separation of intrinsic colour and dust, or the full unseen survey denominator. The [deferred flux-colour holdout design](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/data_audit/flux_colour_holdout_preregistration.json) keeps those identification and data-leakage constraints explicit.

## Data contract and reproduction

[contract.json](../specifications/data_audit/contract.json) defines `photometry_audit.h5`. Groups are `original/{DES,Foundation,LOWZ}`, `current/{DES,Foundation,LOWZ}`, and `diffimg/DES`; each contains every actual `head/` and `phot/` column. `head_row_zero_based` maps observations to headers; -1 identifies unowned delimiters. External `phot_row` and HEAD pointers are one-based. `audit_flags` describes problems without dropping rows; `basic_eligible` retains negative fluxes and has no SNR threshold. objects.csv.gz (`phase2/data_audit/objects.csv.gz`; historical local reference) preserves all header fields plus audit summaries and membership; it is keyed by dataset and HEAD row/SNID.

Run from the repository root:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python scripts/phase2/data_audit/audit.py
.venv/bin/python scripts/phase2/data_audit/finalize.py
.venv/bin/python scripts/phase2/data_audit/verify.py
.venv/bin/python scripts/phase2/data_audit/lcplot_audit.py
.venv/bin/python scripts/phase2/data_audit/lcplot_audit.py --mask32
.venv/bin/python scripts/phase2/data_audit/accepted_flags.py
.venv/bin/python scripts/phase2/data_audit/duplicate_identity.py
.venv/bin/python scripts/phase2/data_audit/exact_duplicates.py
.venv/bin/python scripts/phase2/data_audit/hash_outputs.py
```

The first completed row-audit run encountered a final schema-comparison exception because current DES renamed IAUC. All row products were already written; schema-aware finalization recovered that final comparison without discarding data or hiding the exception. The fixed main script and exact executed code snapshot are retained. Readback verification is independent of that recovery. The HDF5 and decompression cache are large local artifacts and should not be mistaken for small Git-tracked evidence.
