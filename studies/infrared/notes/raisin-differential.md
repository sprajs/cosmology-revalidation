# RAISIN: paired optical–NIR distance audit

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

Status: **bounded descriptive comparison complete; a full paired likelihood is not identified by this release.** No cosmology, dust correction or population model was fitted here. The protocol and physical membership were frozen before computing paired redshift outcomes. Pinned inputs were not changed.

For the same 79 supernovae, the released optical-minus-NIR high-redshift minus low-redshift contrast is **+0.00124 mag**. Its descriptive paired scatter error is **0.04497 mag**. Undoing both exported bias and host-mass terms changes the contrast to **+0.07541 mag**. Thus agreement of the final distance branches depends materially on their existing corrections. This does not establish a new correction, a discrepancy significance, or an independent confirmation of the luminosity–distance relation.

## Frozen question and result

The common sample is exactly 42 low-z CSP, 19 high-z PS1 and 18 high-z DES objects in identical nominal FITRES order. There are no branch-specific removals. Low means zHD<0.1 and high means zHD>0.2; no object is between these ranges. Each CID is joined to its released light curve and sky position.

For branch b and object i, define d_i = μ_b,i − μ_NIR,i and

\[
T_b=\frac{1}{37}\sum_{i\in H}d_i-\frac{1}{42}\sum_{i\in L}d_i=a^T d.
\]

The same physical object's distance–redshift relation cancels in d_i. Any branch-constant zero point cancels in T because a^T1=0. This cancellation does not remove branch-dependent selection, standardization, calibration, extinction or population assumptions. The test is unweighted, so the result is not selected by potentially different branch error models.

The immutable [protocol](../specifications/experiments/raisin_differential/protocol.json), SHA256 `caa811c0c0c22e055b8971bfcd8802a9cf8a496b404d54ce9fe4c942672e0567`, precedes the outcomes. The [frozen membership](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/frozen-membership.csv) preserves every ID, position, redshift and host split.

| Existing additive terms retained | Optical−NIR T (mag) | Descriptive paired SE | Combined−NIR T (mag) | Descriptive paired SE |
|---|---:|---:|---:|---:|
| Released bias and mass terms | +0.001236 | 0.044973 | −0.027098 | 0.027341 |
| Bias only; undo mass | +0.030859 | 0.045271 | +0.029898 | 0.028559 |
| Mass only; undo bias | +0.045788 | 0.044204 | +0.017825 | 0.026631 |
| Undo both | +0.075412 | 0.044510 | +0.074821 | 0.027998 |

These arithmetic removals leave the original light-curve fits and sample fixed. “Undo both” remains an empirical SNooPy fitted-distance comparison, not a raw flux result or a validated alternative selection correction. The optical and combined branches are not independent replications.

The differential bias term changes the optical–NIR contrast by −0.044553 mag; the differential mass term changes it by −0.029624 mag. Their sum explains the change from +0.075412 to +0.001236 exactly. In the combined–NIR case the corresponding changes are −0.044923 and −0.056996 mag.

Prespecified descriptive splits of the released optical–NIR result are:

| Split | High / low count | T (mag) | Paired scatter SE |
|---|---:|---:|---:|
| PS1 high / CSP low | 19 / 42 | −0.000670 | 0.059174 |
| DES high / CSP low | 18 / 42 | +0.003248 | 0.062230 |
| Host log10(M/Msun)≥10 | 11 / 35 | +0.031827 | 0.051504 |
| Host log10(M/Msun)<10 | 26 / 7 | −0.050523 | 0.078882 |

These are descriptive; the two survey contrasts reuse all 42 low-z objects. The host threshold was fixed from the documented imbalance, not optimized on these outcomes. The small seven-object low-mass low-z cell prevents a precise host-matched comparison. Full per-object and per-variant values are in [paired-object-ledger.csv](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/paired-object-ledger.csv) and [results.json](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/results.json).

## Signed distance accounting

The audit uses only `distances/w`, not the separate H0/calibrator branch. For each branch, the exported CosmoSIS vector satisfies

\[
m_b^{\rm cosmosis}=\mathrm{DLMAG}-19.36
\]

to its six-decimal printing precision. Both CosmoSIS `zcmb` and `zhel` contain FITRES **zHD**. Its name field is always 0.0, so it is not an object identifier. Membership/order must come from FITRES and the matching numerical columns.

The released arithmetic is, up to a branch-constant normalization,

\[
\mathrm{DLMAG}=\mu_{\rm fit}-\mathrm{DLMAG\_biascor}+\mathrm{MASS\_CORR}+c_b.
\]

In FITOPT004,005,023–026 the light-curve photometry/model fit is unchanged. Across these variants, `δDLMAG + δDLMAG_biascor − δMASS_CORR = 0` closes at approximately 10^−14 mag. This establishes the sign of the released correction increments; mass terms can also be re-estimated after a bias variant.

Independent closure of the nominal columns gives

\[
\mathrm{mures}=\mathrm{DLMAG}-\mathrm{MASS\_CORR}-\mu_{70,0.3}(z_{HD})-k_b,
\]

where the flat ΛCDM reference is used **only to identify the exported column algebra**, not to estimate this contrast. The branch constants k_b are −0.0279324554 (NIR), −0.1055609249 (optical) and −0.0619375466 mag (combined); closure is around 10^−13 mag. Consequently `mures` is not the final mass-corrected distance residual. Subtracting `mures` across branches would silently select a different correction state.

Jones2022 Appendix E says the table's bias is added to its “raw distance.” This conflicts with the released-column sign. Three independent rounded Table7 entries instead satisfy `DLMAG + biascor − MASS_CORR − table_raw = 0.16030…0.16064 mag`, constant within the table's 0.001-mag rounding. We preserve this documentation inconsistency and use the directly verified release algebra. The original pre-correction FITRES normalization is not supplied. See [algebra audit](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/algebra-covariance-audit.json) and [paper provenance/sign check](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/primary-paper-provenance.json). [Jones et al. 2022, v2](https://arxiv.org/html/2201.07801v2)

## Covariance: what closes and what does not

Each `*.covmat` has a zero diagonal. Its companion `lcparams` dmb **already includes that systematic group's diagonal**. The consistent reconstruction is

\[
C_{{\rm sys},g}=C_{{\rm off},g}+\operatorname{diag}(dmb_g^2-dmb_{stat}^2),\qquad
C_{{\rm total},g}=C_{{\rm off},g}+\operatorname{diag}(dmb_g^2).
\]

The `stat.covmat` is exactly zero. Treating a systematics `.covmat` alone as a full covariance gives an indefinite matrix; adding the diagonal a second time also gives the wrong errors. Restored systematic matrices have tiny negative eigenvalues of order 3×10^−7 mag², consistent with six-decimal dmb rounding; no eigenvalues were silently clipped. The scalar variances used below are positive.

No fitted optical–NIR statistical or intrinsic-scatter cross block is supplied. `model/OIR.J22/OIR.INFO` is a seven-band **simulation dispersion model** with empirical band correlations. It is not an optical–NIR distance cross-covariance: using it would additionally require the actual joint fit response, peak-time propagation and selection model.

For the primary optical–NIR contrast:

* The descriptive SE, sqrt(s_H²/37+s_L²/42), is **0.04497 mag**. It measures paired object scatter under independent-object sampling, not a full uncertainty including shared calibration/population terms.
* Published statistical marginal errors with per-object cross-branch correlation ρ=0 give **0.04991 mag**, only an illustration. Allowing each ρ_i∈[−1,1], while retaining independence across different SNe, gives SD bounds **0.00932–0.06997 mag** from Σa_i²(σ_opt,i∓σ_NIR,i)².
* The supplied `all` systematic marginals alone give contrast SDs **0.04212 mag** (optical) and **0.05125 mag** (NIR). Unknown cross covariance gives the Cauchy bounds **0.00913–0.09337 mag** for their difference.
* Using each branch's full statistical-plus-systematic marginal covariance and allowing any compatible cross block gives the broader SD bounds **0.00996–0.11698 mag**. These are conditional covariance bounds, not confidence limits or an estimated error bar.

There is also a real export-reconstruction gap. After projecting out a constant intercept, several documented signed variants reproduce their marginal groups with unit outer-product weights: e.g. NIR low-z calibration and both stretch variants close to relative precision ≲1.3×10^−4. However:

* The NIR `lcfitter` covariance does not match the supplied FITOPT006 distance-response direction (correlation magnitude only 0.014). The subsequent author-code acquisition finds its covariance bytes identical to the author's `tmpl` covariance. It is not licensed to infer a shared cross block from that group.
* The pecvel covariance is not reproduced by either the distance change or its exact reference-Hubble-residual change; even the latter leaves roughly 27–30% relative matrix mismatch after centering.
* The two released mass variants are collinear after intercept removal. Their separate weights cannot be identified from the marginal matrix alone; numerical nonnegative weights are not unique physical priors. More concretely, FITOPT004 is labelled a mass-divide change to 10.44 dex, yet all three branches retain the nominal `MASS_CORR` sign for every object, including the 16 hosts between 10 and 10.44. Only the two-level amplitude changes. The source-linked follow-up below explains the propagation defect and independently establishes both unit weights.
* The category names do not provide a complete disjoint partition of `all`. For example, optical `lowzcal` is zero although the optical CSP calibration FITOPT010–014 responses are nonzero. Summing the apparently disjoint named groups misses `all` by 72%, 32% and 66% in Frobenius norm for NIR, optical and combined respectively. The subsequently acquired source defines `lowzcal` as only CSP Y/J/H and `photcal` as Y/J/H plus the HST/CALSPEC mode, even when processing other branches. This explains part of the category mismatch; it is not proof that the `all` marginal is invalid.

The [systematic reconstruction ledger](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/systematic-reconstruction.json) records those successes and failures. It does not select a covariance that makes the paired result more or less significant.

The physically common signed-input mapping is available for calibration/MW variants and the labelled population variants. Their **partial, conditional** root-sum-square responses in the primary contrast are 0.02762 mag (calibration), 0.01444 (MW), 0.03841 (stretch populations) and 0.00754 (extinction populations), assuming unit variance and independence of the documented variant amplitudes. These are sensitivity scales, not a recovered complete paired uncertainty. All 28 signed changes are preserved in [systematic-contrast-shifts.csv](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/systematic-contrast-shifts.csv).

## Inputs, fitting conventions and signed modes

The pinned commit is `a383c4bd03c9fbfd64bf5bfda525aec38d32b39c`. All **490 files** pass both acquisition SHA256 and Git-blob identity checks. The [provenance ledger](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/provenance-verification.json) and saved primary-paper PDFs preserve the versions used.

The release contains 118 light-curve files representing 117 physical objects: 77 CSP, 21 PS1 and 19 DES before final w-sample selection. Two files for calibrator 2007sr differ in peak metadata; that object is outside this 79-object test. FLUXCAL has zeropoint 27.5. **Correction to this report's initial prose:** the photometry inventory records zero negative fluxes; the later [signed-flux audit](raisin-flux-sign-audit.md) verifies that all 23,007 observations across the 117 physical files are positive and establishes sign-dependent omission in the matched DES ancestors. The earlier sentence claiming negative flux was retained was wrong. No nonpositive quoted flux errors were found. The three nominal KCOR FITS include filter transmissions, primary spectra, SN SED/K-correction tables; nine example NMLs cover the three surveys and three branches. The SNooPy_B18 FITS contains the empirical light-curve grid. The OIR dispersion file is separate. Inventory/checks are in [FITS](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/fits-inventory.json), [filters](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/kcor-filter-inventory.csv) and [photometry](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/photometry-inventory.json). The inventory itself did not have this prose error and remains unchanged.

All nominal NIR rows have stretch=1, AV=0, and zero exported errors on those quantities and peak time. NIR NMLs fix all three; only amplitude is fitted. They use observer JH for high-z HST and YyJjH for CSP, empirical RV=1.518, and rest phases −15 to +45 days. Optical/combined examples fit shape, AV and peak, use phases −7 to +45 days, and set `OPT_COVAR_FLUX=0`.

The paper describes the fixed NIR peak as optically determined; the repository README says optical+NIR. The fixed NIR peak agrees with the supplied photometry header within 0.0117 day, but generally differs from the exported final optical and combined fit peaks (median differences 0.3008 and 0.2930 day). The precise upstream peak-producing fit is therefore not execution-linked here. Either description implies shared information rather than independent NIR timing. See [peak ledger](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/peak-origin-audit.csv).

Nominal NIR FITCHI2/NDOF ranges from 0.043 to 47.48, median 2.76. The supplied DLMAGERR spans 0.115–0.195 mag and is not a naive free-parameter Hessian error; fixing shape/color deliberately leaves real diversity in the light curves, with distance dispersion treated empirically. The large χ² values alone do not justify rejecting these objects or rescaling their errors by sqrt(χ²/NDOF). [Jones et al. 2022, sections III.2–III.5](https://arxiv.org/html/2201.07801v2)

The machine-readable [signed NML ledger](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/nml-signed-variants.json) retains exact commands. Key points:

| FITOPT | Signed meaning and status |
|---|---|
| 001 | MW E(B−V) multiplied by **0.95**. |
| 002 | Shared wavelength-dependent magnitude shift `0.00714 × wavelength_in_microns`; not solely an HST-only offset. |
| 003 | Supplied shifted peculiar-velocity header; keep common redshift change paired. |
| 004 / 005 | Nominal documentation: mass divide 10.44 / mass-step size shift; placeholders in fitting NML, downstream corrections. |
| 006 | NIR BayeSN alternative; SNooPy bias corrections reused in the paper; not a regenerated independent bias pipeline. |
| 007–014 | CSP Y +0.03, J/H +0.02, B/V/g/r/i +0.01 mag in the specified filters. |
| 015–018 | PS1 g/r/i/z +0.003 mag. |
| 019–022 | DES g/r/i/z +0.006 mag. |
| 023 / 025 | Low-z / high-z stretch population shifts. README mistakenly labels 025 low-z; NML explicitly says HIGHZ. |
| 024 / 026 | Global / low-z extinction-population variant. Paper specifies reducing exponential E(B−V) scale 0.13→0.10, approximately −0.05 in mean AV for RV=1.52. |
| 027 / 028 | K-correction / template-flux variants. Released shifted distance vectors exist, but the three referenced `_sys.fits` KCOR files are absent. |

NMLs also name an RV=3.1 FITOPT029, while this release supplies only FITOPT000–028. Historic environment paths and example NMLs are not a complete executable pipeline. The subsequent public-author-repository acquisition below supplies processing code and additional outputs, but its complete covariance reconstruction still fails. No generic replacements were silently introduced.

## Shared observations and training information

At a strict 1-arcsec coordinate threshold, **47/79** RAISIN objects match Pantheon+ (29 CSP, 11 PS1, 7 DES), producing **76 distance rows** because some physical objects have multiple measurements. Those rows are not independent new supernovae. **17/18** RAISIN DES objects match DES5YR photometry, **15** occur in its distance table, and **14** occur in our existing validation1020 cohort. DES16C2cva has no close photometry match; it was not force-joined. DES16X3cry's matched DES record lacks a usable redshift, while DES16S1bno's matched record has a discrepant redshift; both are outside the DES5YR distance subset. See [sky-match ledger](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/overlap-matches.csv) and [nearest DES metadata](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/des-nearest-metadata.csv).

The newly pinned, config-listed public K21 SALT3 training roster contains matches to **43/79** RAISIN objects: 22 CSP-selected objects, 14 PS1 and 7 DES. Of these, 41 have sky separation≤1 arcsec; two additional exact-name matches have offsets between 1 and 3 arcsec, consistent with rounded training coordinates and explicitly flagged. A RAISIN CSP object can enter the training roster through another observing survey. This is an **input-roster overlap proxy**, not proof of final execution acceptance or exact photometric equality after recalibration. The DES documentation identifies the K21 sample, Fragilistic recalibration and exclusion of CFA-U data; an execution-linked accepted training list remains absent. [DES model documentation](https://des-sn-dr.readthedocs.io/en/latest/2_LCFIT_MODEL.html), [training overlap](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/training-input-overlap-summary.json)

RAISIN adds real rest-frame red/NIR photometry with useful wavelength sensitivity. It does not add 79 independent optical-selected objects to DES/Pantheon. SNooPy templates and its empirical dispersion model use low-z CSP information. The 2024 BayeSN reanalysis uses the same sample, a model trained on low-z optical–NIR data including CSP, and external distances computed under fixed flat ΛCDM (H0=73.24, ΩM=0.28), with inherited intrinsic-SED/gray-scatter structure. Its dust-population posterior is conditional on those assumptions and cannot be imported as an independent cosmology-free constraint. [Thorp et al. 2024, v2](https://arxiv.org/html/2402.18624v2)

## Source-linked FITOPT004 follow-up and frozen counterfactual

A bounded follow-up located the public author repository [RAISIN_cosmo at b888214a5cbae38ac0bf886488ce7734f5b77e87](https://github.com/djones1040/RAISIN_cosmo/tree/b888214a5cbae38ac0bf886488ce7734f5b77e87), dated 2022-05-05. Relevant code and three SUBMIT.INFO files were downloaded with exact Git-blob and SHA256 verification. **All 171 released w distance/covariance/lcparams files have exact blob matches in that repository**, including all three baseline/004/005 distance branches. The release's `massstep` covariance is the author's file named `massdivide`. The [cross-repository mapping](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/mass-threshold/cross-repository-all-distance-matches.json) and [acquisition record](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/mass-threshold/code-acquisition.json) preserve this linkage.

This is a concrete source/output propagation defect rather than an unexplained header alone:

1. `SUBMIT.INFO` maps FITOPT004 to `MASS_DIVIDE` and FITOPT005 to `MASS_STEP`, with baseline light-curve fits reused.
2. `nuisance_params` sets `msteploc=10.44` for MASS_DIVIDE and uses it in the Gaussian-CDF mixture membership probabilities that estimate the step amplitude. [cosmo_sys.py, lines 550–566](https://github.com/djones1040/RAISIN_cosmo/blob/b888214a5cbae38ac0bf886488ce7734f5b77e87/raisin_cosmo/cosmo_sys.py#L550-L566)
3. The final DLMAG and MASS_CORR assignment instead uses `HOST_LOGMASS > 10` and `<=10` in both branches of its conditional. The fitted divide is not passed to this application. [lines 619–631](https://github.com/djones1040/RAISIN_cosmo/blob/b888214a5cbae38ac0bf886488ce7734f5b77e87/raisin_cosmo/cosmo_sys.py#L619-L631)
4. The caller bias-corrects each fit option, estimates nuisance parameters, writes `_new.FITRES`, and then constructs covariance. Each selected systematic contributes one unit-weight outer product after subtracting that variant's inverse-DLMAGERR² weighted mean shift. [lines 639–680](https://github.com/djones1040/RAISIN_cosmo/blob/b888214a5cbae38ac0bf886488ce7734f5b77e87/raisin_cosmo/cosmo_sys.py#L639-L680), [caller](https://github.com/djones1040/RAISIN_cosmo/blob/b888214a5cbae38ac0bf886488ce7734f5b77e87/raisin_cosmo/cosmo_sys.py#L819-L831)

This behavior agrees with the exact released rows. It conflicts with interpreting FITOPT004 as applying the step at the stated new divide. It does not by itself establish which historical cosmology chain used which output revision.

The separately frozen [counterfactual protocol](../specifications/experiments/raisin_differential/mass-threshold/counterfactual-protocol.json), SHA256 `48de7b316671d0d7b664cd9b27a147739c271dc1f7450dcfb83cf67e48c6cf41`, retains all fitted amplitudes, errors, nuisance parameters, nominal distances and other variants. If A4 is the published FITOPT004 half-step amplitude, it changes only

\[
c_{4,i}^{cf}=A_4\,[2\mathbf{1}(M_i>10.44)-1],\qquad
\mu_{4,i}^{cf}=\mu_{4,i}^{pub}+c_{4,i}^{cf}-c_{4,i}^{pub}.
\]

The exact `>` / `<=` semantics are retained; no object lies exactly at 10 or 10.44. No amplitude is refitted and no probabilistic host reassignment is invented. The affected 16 objects are:

| CID | Stratum | log10 host mass |
|---|---|---:|
| 2004ef | low CSP | 10.439 |
| 2004gu | low CSP | 10.092 |
| 2005hc | low CSP | 10.307 |
| 2006hx | low CSP | 10.105 |
| 2007ca | low CSP | 10.118 |
| 2008ar | low CSP | 10.002 |
| 2007as | low CSP | 10.220 |
| 2009aa | low CSP | 10.430 |
| 2005iq | low CSP | 10.318 |
| 2007A | low CSP | 10.094 |
| 2009ab | low CSP | 10.360 |
| PScJ440236 | high PS1 | 10.385 |
| PScH540087 | high PS1 | 10.092 |
| PScF520062 | high PS1 | 10.070 |
| PScB480464 | high PS1 | 10.283 |
| DES16X3zd | high DES | 10.053 |

Every one has a positive published FITOPT004 correction, but the declared-threshold application would make it negative. The full fitted step amplitudes (2A4) are 0.02594360, 0.05431113 and 0.07287397 mag for NIR, optical and combined respectively. Thus every affected row changes by the negative of that branch's full step. The exact nominal/actual/expected signs and distance increments are in the [48-row, three-branch ledger](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/mass-threshold/affected16-by-branch.csv).

| Branch | Published 004 high−low response | Conditional 004 high−low response | Change | Published mass-group SD | Conditional mass-group SD |
|---|---:|---:|---:|---:|---:|
| NIR | −0.001934 | +0.001355 | +0.003289 | 0.025295 | 0.025258 |
| Optical | +0.012484 | +0.019369 | +0.006885 | 0.028004 | 0.031678 |
| Combined | +0.029906 | +0.039144 | +0.009238 | 0.035337 | 0.043435 |

All entries are mag. The covariance columns concern each single branch's high−low contrast, not an invented optical–NIR cross covariance. The paired optical−NIR FITOPT004 response changes **+0.014418→+0.018014 mag**, a **+0.003596 mag** propagation sensitivity. Combined−NIR changes **+0.031840→+0.037789 mag**, or **+0.005949 mag**. These alter a systematic variant, not the nominal +0.001236-mag paired result.

For each branch, the unchanged source formula for FITOPT004+005 reproduces the released mass-group covariance within the predeclared literal export-rounding gate: maximum off-diagonal error ≤4.65×10^−9 mag² and diagonal error ≤3.09×10^−7 mag². The counterfactual mass-group matrix replaces only the old004 outer product with its new004 counterpart. After intercept projection, the response rank changes **1→2**. No source matrices were overwritten or projected to positive definiteness.

**The full `all` matrix does not pass that gate.** Unit-weight reconstruction through FITOPT028 has relative matrix errors 2.755 (NIR), 0.00568 (optical) and 0.00265 (combined); including the author's separately acquired FITOPT029 does not fix it. Therefore no full repaired `all` covariance is released here; this closure result does not support a full corrected cosmology calculation. The [result](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/mass-threshold/counterfactual-result.json) and arrays (`runs/research_2026_09_26/raisin_differential/mass-threshold/counterfactual-arrays.npz`; historical local reference) contain reviewable old/new mass-group matrices and the candidate component change, with full-matrix replacement explicitly withheld. A later audit must identify the exact vector/weight revisions that generated `all` before carrying this propagation change into a full likelihood.

The bounded Git-history follow-up does not remove this blocker. The last `cosmo_sys.py` change is [91c7c632a5a1, 2022-01-02](https://github.com/djones1040/RAISIN_cosmo/commit/91c7c632a5a17c76e372c9f9a1182ba17acefa27), byte-identical to May5 HEAD. The preceding Dec28 and Nov12 source versions have **identical parsed all-covariance inner loops**; their differences concern preparation of fits/errors and category lists, not all-matrix weights or exclusions. The three branch-specific CSP SUBMIT.INFO files last changed Dec27 and consistently list 30 options, with 004=MASS_DIVIDE, 005=MASS_STEP and 029=RV. The matching NIR `all` and FITOPT004 outputs last changed together [May2 at 5f12c6228f2a](https://github.com/djones1040/RAISIN_cosmo/commit/5f12c6228f2ad93dc458afc1c1d2895235c81c08), after a preceding redshift-fix commit; there is no intervening source-file revision supplying a new all-covariance recipe. This is a bounded provenance audit, not a search over outcome-fitted weights. [History and exact option mapping](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential/mass-threshold/history-audit-result.json)

An independent [mass-propagation checker](../code/verify_raisin_mass_repair.py), reading raw FITRES0/4/5 and released covariance/lcparams rather than importing the analysis code, reproduces the 16 flips, 11-low/5-high counts, ranks, shifts and old/new mass-group arrays. Maximum array disagreement is **7.39×10^−15**. It retains the protocol's rounded published covariance and swaps only FITOPT004, preserving the original rounding residual. [Independent record](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_mass_independent_verification/result.json)

## Next decisive test and reproduction boundary

The next useful experiment is a **joint optical/NIR flux fit and forward-selection test with one physical latent object per CID**. First recover execution-linked upstream peak fits, accepted epochs, covariance vector/weight revisions and missing systematic KCOR inputs. Reproduce zero-perturbation branch distances and the `all` marginal exports. Perturb shared optical flux/calibration and refit the timing/selection quantities before carrying them into NIR; pair the same noise/population realization across branches to estimate the statistical/intrinsic cross block. Retain signed flux, host-dependent selection, non-detections and template-subtraction errors. Validate any dispersion law against the released residual correlation data rather than treating OIR.J22 as an exact distance covariance.

Only after that closure should an independently constrained luminosity/color/extinction-population alternative be compared through the same selection procedure. A near-zero final paired contrast cannot establish arbitrary gray evolution's absence: a shared gray luminosity change cancels in the optical–NIR difference itself. Nor does the +0.0754-mag correction-removal contrast measure physical drift.

Reproduction uses [raisin_differential.py](../code/raisin_differential.py) with `phase2/env-official/bin/python`, `OPENBLAS_NUM_THREADS=1`, commands `freeze`, `audit`, `contrast`, followed by the bounded follow-up's `freeze-mass` and `mass-counterfactual`. Both freeze commands refuse to overwrite their original protocols. All numerical outcomes are local derived artifacts in [raisin_differential/](https://github.com/sprajs/cosmology-revalidation/tree/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_differential). Primary source files were rehashed after the work; no pinned input was edited.

An independent [checker](../code/verify_raisin_contrasts.py), which does not import this analysis script, reproduced all 40 contrasts, paired empirical SEs and supplied-marginal covariance bounds directly from FITRES/lcparams/cov files: maximum difference **1.64×10^−15**, with the export's tiny negative systematic eigenvalues preserved. See its [verification record](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/raisin_independent_verification/result.json). This validates the arithmetic, not the absent cross-branch covariance or physical correction model.
