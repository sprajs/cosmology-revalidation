# Original DES-SN5YR: executable reconstruction from calibrated fluxes

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

This is a **partial reproduction with an explicit failed exact-match gate**, not a claim that the complete original cosmological pipeline has been regenerated. The initial bounded wave now runs the public SNANA fitter and simulator, refits the calibrated light curves, exports a consistent full measurement covariance, and identifies specific pipeline gaps. The released-distance cosmology result from phase 1 is a different, already completed target.

## What is fixed, acquired and executable

The original target is DES-SN5YR tag 1.3, commit `e3493cb3b9505fc1f3f392d887364b50dae21439`. Current Dovekie products were not substituted. The plan was frozen in `phase2/official/PLAN.md` before fitting outcomes; later evidence-led changes are in `AMENDMENTS.md`.

The primary scientific specification is [Vincenzi et al. 2401.02945v2](https://arxiv.org/html/2401.02945v2), especially sections 2.7, 3.1–3.6 and 5; the released Pippin YAML and SNANA namelists specify actual settings. The [official release documentation](https://des-sn-dr.readthedocs.io/en/latest/) supplies the data hierarchy. Local text versions of the paper were read in addition to the README files.

The public calibrated SMP flux files already contain chromatic/DCR corrections and the anomalous host-surface-brightness flux-error correction. The June 12, 2024 namelist comment explicitly disables their old runtime application because they are now in the release. Reapplying them would double-count them. Calibrated flux is `mag = 27.5 - 2.5 log10(FLUXCAL)`; it is not a detector-image pixel-level remeasurement.

The exact [July 3, 2024 SNDATA archive](https://zenodo.org/records/12655677) was downloaded (1,653,076,974 bytes), hashed and fully extracted under `phase2/official/inputs/SNDATA_ROOT`. It provides:

- `kcor/DES/DES-SN5YR/calib_DES-SN5YR_{DES,LOWZ,Foundation}.fits.gz`, including the exact filter transmissions, primary SEDs and zero-point offsets;
- DES DECam response curves; DES Fragilistic AB offsets g −0.0008, r −0.0116, i −0.0090, z +0.0071 mag;
- the original SALT3 model, SFD Galactic-extinction maps, survey SIMLIB/HOSTLIB/efficiency files, population PDFs and released classifier weights;
- original Pippin, light-curve-fitting and BBC inputs.

The old README's [2020 SNDATA link](https://zenodo.org/records/4015340) predates DES-SN5YR; the acquired 2024 record supplies the actual missing assets. `inputs/sndata-inventory.txt`, `source-manifest.json` and the downloaded Zenodo record preserve provenance.

SNANA current snapshot `886408a4e171896db5eaa97e735a655f50cec2db` (2026-09-18) and historical snapshot `2fe0f564361a873860661ff61080b4db9c607edf` (2023-05-11) both compile locally. The historical pin is supported by the released simulation HEAD version, **not proven to be the real-data light-curve-fit executable**. The full source archive and commit metadata are retained. Both `snlc_fit.exe` and historical `snlc_sim.exe` now run. No global packages were installed: the existing system GCC/CFITSIO plus locally extracted matching Arch GCC Fortran and GSL packages were used. `build-manifest.json` and `snana-output-only.patch` identify the executable and modifications. An isolated Python 3.12 environment is pinned in `inputs/requirements.lock.txt`; sncosmo was installed as a fallback, but the reported fits use actual SNANA.

## Official fitting choices retained

The DES fit uses SALT3.DES5YR; fixed spectroscopic redshift; Galactic reddening from a fresh SFD map lookup multiplied by 0.86 exactly once; Fitzpatrick 1999 extinction with Rv=3.1; filter mean rest wavelengths 3500–8000 Angstrom; fitted rest phases −15 to +45 days; three iterations; MINOS; a first-iteration 0.05 peak-flux error contribution; and observation rejection at delta-chi-square 10. The default covariance includes SALT model and Milky Way extinction uncertainty. PSF cuts are on FWHM in arcseconds, not on the PHOT table's sigma in pixels. `ZPNPE=ZEROPT+2.5 log10(GAIN)`.

The actual SALT amplitude convention is `mx=-2.5 log10(x0)`, with SNANA's tabulated `mB=mx+10.635`. The SALT model's `MAG_OFFSET=0.27` is internal to predicted photometry and must not be applied as a second correction to the released distance.

## Distance algebra closes after resolving a source-precision trap

For all 1829 SNe, the released equation is

`mu = mx + 0.16087*x1 - 3.11780*c + 0.03754*(expit((logmass-10)/0.001)-0.5) - biasCor_mu + 29.95821`.

The metadata CSV rounds host mass to two decimals. Near its very narrow transition this introduces an apparent maximum 0.0171 mag mismatch. Using the unrounded public HEAD host masses reduces the **all-row maximum residual to 0.000266 mag and RMS to 0.000101 mag**, consistent with other printed coefficients/parameters. This tests the published distance algebra; it does not independently fit alpha, beta, gamma, M0 or reproduce the simulation-derived bias correction. The row-by-row ledger is `results/distance-algebra.csv`.

## Calibrated-flux fits and the exact-match gate

A deterministic 12-object pilot was selected by evenly spaced redshift ranks, then all original Hubble-diagram members were attempted. With the literal released base namelist, 1827/1829 fits succeed: 1634 DES plus 193 low-z. `results/all_refit_observables.csv` and `all_refit_covariance.npz` contain the real refits; `all_comparison.csv` and `diagnostics/all_refit_summary.json` retain every discrepancy.

The initial literal-namelist median mB difference is −0.0000096 mag; 90% have |difference| below 0.00289 mag. However, 83 objects differ by more than 0.01 mag and the maximum is 0.71014 mag. There are 271 DES objects with changed degrees of freedom. This cannot be summarized as full reproduction.

Two concrete configuration issues were identified:

1. **Photometry flag rejection:** The DES release lacks a DES `.IGNORE` file and the base namelist does not explicitly reject PHOTFLAG bit 32. Published LCPLOT data exclude that flag: none of 62,608 unambiguously mapped accepted observations has it. Excluding bit 32 exactly restores all pilot observation counts and reduces its maximum mB difference to 0.00223 mag. Other never-observed bits were not automatically declared cuts. This rule is recovered from published accepted observations, not independently documented as the original runtime mask.
2. **Consistent-sample override:** Pippin adds `OPT_SNCID_LIST=1`, whereas the base namelist alone does not. Omitting it loses low-z SN 2007ob at a rest-phase-minimum boundary of approximately five days. The corrected configuration restores the intended override explicitly.

These fixes have separately named outputs and do not overwrite the literal-namelist baseline. The recovered configuration remains subject to the exact-match gate; see its saved diagnostics for the final all-row comparison.

The completed recovered arm yields **1828/1829** fits, including all 194 low-z members. Only DES CID1257112 fails. Its 90th percentile |delta mB| is 0.00150 mag, but 66 objects remain above 0.01 mag. The discrepancy therefore cannot be described solely by its near-zero median. With published `MUERR_FINAL` diagonal weights, the fixed-alpha/beta Tripp differences have weighted RMS **0.0202 mag** and mean +0.00127 mag. A linear flat-LCDM response after profiling one intercept is delta Omega_m **+0.00237** for all matched objects; the high classifier-probability subset gives +0.00884. These are propagation diagnostics using fixed published corrections and weights, not an independently rerun BBC/cosmology fit. The full subset results are in `diagnostics/recovered_weighted_gate.json`.

The largest outlier, CID1256426, has exactly the same 27 available public and reconstructed observations but a different iterative-rejection solution (seven published accepted points versus six). Its published classification probability is zero and its BEAMS distance uncertainty is 185.24 mag: its 0.71 mag parameter mismatch has effectively no published cosmological weight. This does **not** dispose of the general problem: among objects with released PIa>0.999 there are still 32 initial differences above 0.01 mag, up to 0.2102 mag.

Running the historical 2023 fitter on the 85-object outlier/control list does not repair all discrepancies: 60 still have |delta mB|>0.01 mag, with maximum 0.45011 mag. Using published parameters only as an initialization diagnostic also fails to recover all cases. Published per-point fluxes and errors agree with the raw release to float/printing precision, so hidden photometric corrections are not supported as the explanation for the inspected cases. The remaining issue is the exact effective initialization, quality-mask and iterative-rejection implementation/configuration used to generate the released fits. No problem rows were silently discarded.

The historical source archive was checked against the full Git tree, not just its filename: all 223 unchanged source files have the expected Git blob hashes. The sole changed source header has the official build's ROOT/HBOOK preprocessor flags disabled. An anomalous comment mentioning 2026 is already in that pinned upstream blob despite its 2023 commit date; it does not establish a failed short-ref resolution. See `diagnostics/historical-tree-audit.json`.

The additional arm conditioning on the author's published accepted-epoch masks is complete. All **1635 DES objects** fit, and every accepted count agrees with the released LCPLOT: **62,702 observations, zero count discrepancies**. All 194 low-z fits are retained from the recovered arm. Further clipping and dynamic phase reselection are disabled, so this tests fitting conditional on the author mask; it does not independently reproduce clipping. Twelve objects' public LCPLOT counts themselves differ from metadata `NDOF+4` by one to four; `diagnostics/conditioned_mask_counts.csv` retains that release-product discrepancy.

For the conditioned all-survey sample, the 90th percentile |delta mB| drops to 0.000419 mag. Five DES objects initially remain above 0.01 mag. Ten DES objects selected by |delta mB|>0.01 or |delta Tripp|>0.01 were then run from a frozen 3x3 stretch/time starting grid plus a published-reference diagnostic start. The finite-start minimum-chi-square arm is `results/conditioned_multistart_best_refit_*`; it preserves all 1829 objects and each row's source fit label. Four DES brightness outliers remain. This arm is selected by objective, never by closeness to published parameters; it is not proof of a global minimum.

**CID1307748 is a material counterexample to forced agreement.** One start recovers the published-like branch x1=-0.811, t0=56689.307; the other branch has x1=4.526, t0=56694.321. Both use the same 16 accepted observations. The high-stretch branch has lower data chi-square under either common fixed covariance: 3.606 versus 7.316 using its covariance, and 8.512 versus 10.094 using the published-like branch covariance. Prior terms are below 0.000013. Thus the branch ordering is not solely caused by comparing different covariance matrices. The high-stretch branch violates the later BBC |x1|<3 quality requirement even though the object belongs to the published Hubble sample. This exposes a local-solution/sample-selection ambiguity under the reconstructed likelihood; it does not establish that the unknown original executable had an identical objective. No row was removed to conceal it.

## Measurement covariance is now available without CSV precision loss

The published CSV amplitude cross-covariances are rounded to zero in almost every row. Furthermore, SNANA source `FIX_COVAR_LCFIT` documents that its ordinary FITRES scalar errors can be MINOS errors whereas cross-covariances are Hessian values. Squaring those scalar errors does not reconstruct one coherent Gaussian covariance and can produce an indefinite matrix.

The local patch changes output only: it prints float table values with 17 significant digits and dumps the **actual unmodified 4x4 Minuit FITERRMAT** for x0,x1,c,t0 after a final accepted fit. It also exports double-precision fitted parameters. The NPZ preserves the original covariance plus its Jacobian transformation to mx,x1,c,t0. All literal, recovered and conditioned Hessian correlation matrices are positive definite; none was repaired. MINOS error columns remain separately available. This is a coherent numerical export, conditional on the fitter, and not by itself validation of uncertainty accuracy.

The independent flux audit identifies a limitation of the exported Hessian: for CID1896213 its t0 sigma is **0.123924 days**, while SNANA MINOS gives **0.564311** and the release gives **0.5644**. Independent fixed-objective curvature agrees with the wider uncertainty. Three other pilot objects show smaller related t0 discrepancies. A grid/interpolation curvature artifact is a hypothesis under investigation. The independent mB/x1/c block agrees much more closely, but a positive-definite matrix alone must not be presented as validated full-parameter uncertainty.

## Assumptions that numerical reproduction does not validate

`phase2/official/fitter-assumptions.json` records nine implementation assumptions with source file hashes, line references, effective settings and explicit falsifiers. It covers the spectral family, interpolation/extrapolation, calibration/filter response, redshift/extinction, model uncertainty, objective/covariance, priors/initialization, clipping and parameter errors. These are claims to test, not validated facts about supernova physics.

One material detail is that the default fit is iterative generalized least squares: `OPT_COVAR_FLUX=1` holds covariance fixed within each later minimization, while `OPT_CHI2_SIGMA=0` omits the log-determinant term. The covariance includes same-band SALT colour dispersion and a cross-band Milky Way extinction term, as well as diagonal photometry variance. It is not an iid flux likelihood, and it is not a continuously parameter-dependent normalized Gaussian likelihood. An independent implementation must first reproduce this stated objective before comparing a changed likelihood.

The portable 12-object fixture in `phase2/official/portable_pilot/` supplies 1394 original calibrated observations with source-row identifiers, published and refit masks/model curves, true double-precision parameters and Hessian, exact original SALT surfaces, calibration/filter responses, and source hashes. Separate output-only audit builds now additionally export **528 accepted-epoch REAL8 predictions**, exact observation MJDs, frozen flux covariance/inverse, actual per-iteration prior centres and bounds, preceding fit parameters, and the MW/model/data covariance ingredients. Reconstructing r^T C^-1 r plus prior matches the source objective to less than 6e-14, and all three audit builds reproduce unchanged fit parameters/chi-square. The effective t0 prior centre resets to the previous fitted t0 between iterations; it is not always the original flux-search peak.

An initial fixture packaging error included `salt3_lc_model_variance_*` rather than the actual runtime `salt3_lc_variance_*` files. The independent audit caught it; the actual runtime maps are now included and the contract explicitly records the correction. The real SNANA fits always read the full original release directory and were unaffected. Wrong-asset numerical comparisons must not be interpreted as physical evidence. Independent agreement checks implementation; held-out-band/phase residuals, covariance coverage, injection recovery and multi-start stability are needed to challenge assumptions.

## Simulations, selection and what has not yet been regenerated

The first complete released nominal-Ia mock (`...P21-0001`) was actually fitted through the official DES light-curve configuration: 2829 detected objects → 2176 after SNANA preselection → 1936 accepted fits. Applying the published basic BBC quality inequalities to those fits leaves 1384. `results/mock0001_fit_and_quality.csv.gz` retains fit parameters, simulation truth and the mask; `mock0001_all_detected_status.csv` accounts for all 2829 objects.

A newly generated nominal P21 realization (`PH2_pilot02_P21`) also completes the same fitter path: **2678 detected → 2052 preselected → 1825 fitted → 1314 basic-quality**. Its full fit covariance and truth/quality tables are `results/forward_p21_refit_covariance.npz` and `forward_p21_fit_and_quality.csv.gz`. The initial historical simulator run crashed because its 2023 GENPDF loader cannot convert a zero-width analytic alpha distribution to a grid; moving the same constant alpha=0.15 declaration from a derived PDF into the simulation input preserves the distribution and avoids that loader path. The source PDF is unchanged. The forward simulation audit owns the generated-attempt denominator, including reused CIDs, so detected CID alone must not be treated as a unique generation-attempt identifier. Decoded SALT3 model/template/colour/error files in the SNDATA archive and original release were independently checked byte-identical.

This is not the final BBC sample. The remaining cuts include bias-grid coverage, `chi2max=16`, consistent membership across systematic variants and contaminant handling. Detected HEAD files also do not supply the generated-but-undetected selection denominator; the release DUMP/README inputs must be used with the actual simulation generation distribution. The separate selection audit owns that reconstruction.

BBC4D estimates bias as a function of z,x1,c,host mass within survey/field groups. Alpha, beta and the residual host step are jointly fitted with redshift-bin offsets, using the Ia/non-Ia mixture likelihood. Published `PROB_SNNV19` is a classifier probability; `1-PROBCC_BEAMS` is the likelihood-updated Ia probability. They are not interchangeable. Unbinned `MUERR_FINAL` further incorporates the probability-dependent and bin-dependent rescaling from Vincenzi equation 10. Neither renormalization nor contamination is equivalent to adding or removing a single distance offset.

A complete independently regenerated BBC fit, systematic covariance, classifier rerun and cosmological result from these new light-curve fits have **not** been obtained in this wave. The precise gaps are the effective original fit/rejection settings, regenerated bias-correction/CC-prior simulation fit tables and complete sample cut chain; availability of a simulator/model file alone does not close them. Newly inferred physical corrections must stay behind this gate.

## Reproduction commands and scope

From the repository root:

```bash
phase2/env-official/bin/python scripts/phase2/official/prepare_baseline.py
bash scripts/phase2/official/build_current.sh
bash scripts/phase2/official/build_historical.sh
bash scripts/phase2/official/run_snana.sh pilot12_hessian
phase2/env-official/bin/python scripts/phase2/official/audit_fits.py --sample pilot
phase2/env-official/bin/python scripts/phase2/official/audit_fits.py --sample all
phase2/env-official/bin/python scripts/phase2/official/audit_fits.py --sample conditioned_multistart_best
phase2/env-official/bin/python scripts/phase2/official/audit_mock.py
phase2/env-official/bin/python scripts/phase2/official/audit_mock.py --label forward_p21 --source phase2/literature/simulations/outputs/PH2_pilot02_P21
```

Run-specific namelists and source IDs are frozen under `phase2/official/inputs`; the stdout logs contain the complete read-back configuration. The commands expect the acquired archives/extracted project-local build dependencies. Rebuilding does not alter the pinned release or phase-1 outputs.
