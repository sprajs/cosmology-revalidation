# CSP-II spectra identify passband-shape effects, not a distance correction

The later FIRE spectra add a real observational constraint beyond the fixed Hsiao/SNooPy templates: a same-spectrum ratio measures the effect of changing a physical passband without requiring a supernova distance or luminosity normalization. The frozen test finds a larger early-to-late change for **J WIRC versus RC1 (−19.97 mmag)** than for **J WIRC versus RC2 (+2.70 mmag)**, using the same 13 objects. This supports treating the instrument identity as consequential. It does not establish which transformation was applied to archived photometry, and it is not a corrected distance likelihood.

This report contains the source audit, outcome-free feasibility inventory, and the subsequently authorized, separately frozen spectral experiment. No Lu template was downloaded or evaluated, and no supernova fit or cosmology was run.

## Primary inputs and what is independent

The [official CSP release](https://csp.obs.carnegiescience.edu/data) supplies the [spectra archive](https://csp.obs.carnegiescience.edu/data/CSPII_NIR_Ia_spectra.zip) and a separate template archive. The former was acquired under the explicit 40 MB resource amendment: 31,173,962 bytes, SHA256 `b43532916bb3defb6e041655954dc838e652716d5a6d191274d155e9f3b6d0f5`. Its 1,027 ZIP entries were individually hashed and safely extracted. The 153,724,421-byte template remains deferred; bounded HTTP ranges recovered its member inventory and README without executing its pickle or installing BYOST.

[Lu et al. (2023), §§II–IV](https://arxiv.org/html/2211.05998v2) describe A0V-based flux/telluric calibration and explicitly leave the measured spectra unmatched to photometric colors. Their phases and color stretches come from multiband SNooPy `max_model`. Their color check uses SNooPy interpolation and has about 0.2 mag RMS. The template's K-correction assessment color-matches spectra and randomly withholds **spectra**, rather than whole supernovae; it excludes relative flux-calibration uncertainty. It therefore cannot supply independent broadband or object-held-out validation here. The released 339 spectra include eight host-contaminated spectra excluded from the 331-spectrum template sample.

Our exact-name comparison finds no overlap of the 98 objects with either the 134-object DR3 archive or the 79-object M20 paper Table 1. Three names occur in the 76-object RAISIN CSP photometry list—SN2012fr, SN2012ht and SN2015F—but none in the selected 79-object RAISIN cosmology cohort. This is a normalized official-name join, not a sky/alias-complete training audit. The 2011–2015 observations postdate the named Hsiao07 model; the exact training/execution provenance of every stored Hsiao/SNooPy component remains unclosed. All host-clean spectra are Lu training inputs: these data are not independent validation of that new template.

The minimum future independence requirement is optical-only timing/shape determination from documented light curves, plus object-disjoint spectral calibration/validation. The current supplied phase labels do not satisfy that timing requirement. Absolute spectrophotometric standards, per-region telluric acceptance, calibration covariance and the execution-linked SNooPy training roster are not in this archive.

## Schema and feasibility gates

The table and every TXT/FITS pair were checked before integrals. Results are in `runs/research_2026_09_26/astra_design/csp_spectral_identification/metadata-result.json` and the two CSV ledgers alongside it.

- 339 spectra of 98 objects; 331 host-clean spectra. Redshift 0.0036–0.086, phase −13.3166 to +96.3217 days, color stretch 0.481–1.300.
- Wavelength is **observer-frame micrometres**, while FITS `YUNITS` identifies flux per Å. TXT has wavelength, flux, error; full spectral covariance is absent. Counts below use Å consistently.
- The CSV repeats the `SNRH` header twice; both columns are identical. Positional parsing preserves this instead of silently overwriting a field.
- Each FITS file has one extra final row absent from TXT; its flux is nonfinite. Every retained TXT value exactly matches the wavelength-joined FITS value, including NaNs. The first failed whole-array equality assertion and corrected join are preserved.
- TXT error arrays contain 1,758 nonfinite pixels across 331 spectra, with no finite nonpositive errors. FITS `ERRSCALE` spans 1.307–9.507; these recorded scales were **not applied again**. All pixels contributing to the executed actual-redshift contrasts have valid errors, but this does not establish independent-pixel covariance.
- Recomputed phase from rounded table values differs from the supplied phase by at most 0.00764 day. Supplied phases are retained.

The metadata-only protocol was frozen before support counting: early [−7,+7], late [+10,+20] rest days, host flag No, color stretch 0.75–1.18, quoted peak error ≤0.5 day and band median S/N ≥10. It yields **17 Y / 13 J / 12 H paired objects**; without SN2012fr, 16/12/11. The declared 15-object population-feasibility threshold is met only in Y; J/H support a bounded operator check, not a population conclusion. No threshold was relaxed after counting. A later [+20,+35] metadata sensitivity has 14 paired objects in each band and was not substituted for the primary window.

## Frozen passband experiment

The follow-up protocol is `passband_test/protocol.json`, SHA256 `a7823c8572eefd8f682f77a16f50bc9dae28cc1d7c0686f723a0428043ecfefa`. Before integration it fixed 138 spectra from 53 objects, producing 503 spectrum/passband pairs within phase [−7,+30]. J WIRC–RC1 and WIRC–RC2 are primary; H is secondary and Y a control. The exact 2018 atmospheric-throughput tables and BD17 spectrum are shared with the parent's `csp_passband_contrast.py` inputs, all hashed.

For spectrum `f`, define photon counts `I_b(f)=∫λ T_b(λ) f_λ(λ)dλ` and

`C_W−R(f) = −2.5 log10[(I_W(f)/I_W(BD17))/(I_R(f)/I_R(BD17))]`.

Gray flux normalization and distance cancel exactly. BD17 has zero assigned differential magnitude; actual natural-system zero-point constants are not being inferred. They also cancel from the same-object late-minus-early difference. The primary integrates each delivered observer-frame spectrum directly, without assuming whether a foreground correction has already been applied. It retains signed spectral values and throughput tails. It uses equal spectrum means within each phase window, then equal weight per object.

The secondary redshift grid 0, .01, .03, .05, .08 transforms wavelengths by `(1+z_target)/(1+z_source)` and transforms flux inversely. This moves the **delivered spectral shape**, including any remaining extinction/calibration imprint. It is not a dust-free intrinsic SED, observed high-redshift sample, or population-transport validation.

All 3,018 integrals passed coverage, finite-flux and positive-count gates without extrapolation. Gauss2/Gauss4 agreement is ≤1.61×10⁻¹⁵ mag; independent 0.5 Å quadrature differs by ≤1.43×10⁻⁵ mag. Gray rescaling changes a contrast by ≤1.47×10⁻¹⁵ mag. A separate implementation using direct 0.25 Å trapezoids reproduces all 503 actual-redshift contrasts within 8.23×10⁻⁷ mag, and the 55 object/passband phase differences within 8.69×10⁻⁷ mag. All summary reaggregations agree exactly. The parent also independently reconstructed all 3,018 cases with exact endpoint/midpoint Simpson integration, importing none of this runner: maximum difference 1.08×10⁻¹⁴ mag, with actual-redshift phase summaries agreeing to 10⁻¹⁵. Evidence: `runs/research_2026_09_26/csp_spectral_root_review/result.json` and `scripts/research_2026_09_26/verify_csp_spectral_contrasts.py`. These checks establish computation, not spectral calibration.

## Observed spectral-shape results

These are **descriptive means and sample dispersions**, not errors on a population mean or discovery significances.

| Late-minus-early passband contrast | Objects | Equal-object mean (mmag) | Sample SD (mmag) | Individual range (mmag) |
|---|---:|---:|---:|---:|
| J WIRC − RC1 | 13 | −19.966 | 12.426 | −35.514 to +0.762 |
| J WIRC − RC2 | 13 | +2.699 | 3.004 | −0.513 to +8.043 |
| H WIRC − RetroCam | 12 | +4.314 | 1.785 | +1.595 to +6.985 |
| Y WIRC − RetroCam | 17 | −6.225 | 18.059 | −37.559 to +29.402 |

Removing SN2012fr gives J means −20.383/+2.418 mmag. Shifting the supplied peak by ±its quoted error, reassigning windows within the frozen spectrum cohort, gives J WIRC–RC1 means −19.741/−19.563 mmag and WIRC–RC2 +2.655/+2.585 mmag. This is a limited boundary sensitivity, not an independently measured clock-error distribution.

The prescribed smooth ±0.2-mag color tilt changes any individual J phase contrast by at most 0.540 mmag; conditional F99 foreground removal changes it by at most 0.246 mmag. These specified perturbations are not a prior or exhaustive calibration bound. Fine spectral-response errors can survive a smooth-tilt check.

Telluric treatment is the larger caveat. Removing observer-frame 1.35–1.45 and 1.80–1.95 μm from **both filters and reference integrals** changes some individual-spectrum J WIRC–RC1 contrasts by 26.40 mmag, although the maximum change of a paired-object phase contrast is 3.984 mmag (WIRC–RC2: 1.738 mmag). This common-support result has a different estimand; it does not reconstruct the missing full-band flux. No primary row was removed because of this outcome. Missing per-region quality masks prevent a claim that every full-band contrast is free of telluric systematics.

The same delivered shapes, transported to target z=0 and .08, give mean J WIRC–RC1 phase contrasts −15.109 and −53.758 mmag, versus WIRC–RC2 +4.240 and −2.859 mmag. This exposes passband/redshift dependence without requiring a luminosity distance, while retaining the stated spectral-calibration and transport limitations.

![Spectral passband contrasts](/home/szymon/Documents/ChatGPT/supernova/runs/research_2026_09_26/astra_design/csp_spectral_identification/passband_test/spectral-phase-contrast.png)

## What this can and cannot identify

The experiment independently tests the **shape-sensitive part** of mapping a measured spectrum between physical instruments. It shows why a passband-label ambiguity need not behave like one common magnitude offset. It does not test the gray luminosity law, absolute distance, or cosmology; those directions cancel by construction. It also does not isolate phase-template error from inaccurate phase labels, dust, or wavelength-dependent spectral calibration. The measured phase contrasts cannot be inserted as a global photometric correction.

A minimally defensible next calibration comparison must first close the raw WIRC → published/SNpy → SNANA label and magnitude/error lineage. For affected real epochs, use the verified physical throughput, independently timed spectral observations or an object-held-out spectral model, and a covariance for residual spectral/throughput uncertainty. Refit affected light curves with the appropriate instrument operator and validate on overlapping actual WIRC/RetroCam measurements, preserving source selection and flags. Any trained templates tied to the previous convention and the bias simulation must then be updated consistently. This route could yield a supported distance-likelihood change; renaming a filter or applying the table's phase mean alone cannot.

The immediate high-information control is therefore the parent's magnitude/error lineage test plus a bounded same-object, dual-instrument observation check. The separate archived-simulation timing intervention remains more direct for quantifying the known near-truth-peak assumption. The large Lu template would add model flexibility at this stage without resolving the missing calibration/timing evidence, so its download remains deferred.

## Reproduction and preserved evidence

All owned scripts and ledgers are under `runs/research_2026_09_26/astra_design/csp_spectral_identification/`. Run `inventory_v2.py` with the workspace `.venv/bin/python`; run `passband_test/run.py` and `check_direct.py` with the existing `phase2/env-official/bin/python`. The executed runner refuses to overwrite completed spectrum results. No packages or shared inputs were changed.

`spectra-acquisition.json` records every ZIP member hash. `metadata_protocol.json` is the original support freeze. `inventory-attempt-ledger.json` records the TXT/FITS trailing-row correction. `template-readme-acquisition.json` records the initial local-ZIP-header size failure and the verified central-directory/CRC extraction. The initial scripts and raw HTTP ranges remain intact. `passband_test/manifest.json` preserves original result hashes; a final manifest additionally covers the independent check, closed log, figure and this report.
