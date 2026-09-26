# SNANA dust implementation and simulation audit

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../../README.md) for current execution status.

Audit date: 21 September 2026. These checks establish implementation facts and conditional simulation responses. They do **not** establish bias-free cosmology. No original inputs, other-thread simulations, or released results were changed.

## Findings that affect the analysis

1. The integer `OPT_MWCOLORLAW=99` is version-dependent. It selects the historical polynomial approximation in SNANA `2fe0f564361a873860661ff61080b4db9c607edf`, and the F99 natural cubic spline in current `886408a4e171896db5eaa97e735a655f50cec2db`. Current `-99` restores the old approximation. Although named “MW”, the setting is also used for host extinction, at **rest-frame** wavelength. Changing the executable without freezing its version changes the host-dust generator even if input text is unchanged.
2. The public Dovekie repository's 20 PIPPIN files are **byte-identical** to the original tag-1.3 files. They still specify SALT3.DES5YR, a BBC `x1ERR<1.0` cut, and nine calibration-amplitude weights of 0.3 (sum of squares 0.81). They are not sufficient provenance for executing the final Dovekie pipeline. This does not show that final released Dovekie distances used those stale settings.
3. The released P21 probability maps have appreciable probability outside the F99 implementation's recommended `2<=RV<=6` range. More decisively, actual **released simulated** objects occupy this range: of 71,946 written Ia mocks across 25 original realizations, 40,369 have `RV<2` and 11,965 have `RV<1`; minimum `RV=0.30264568`. The map branch is not restricted by the separately printed default Gaussian `RV` range.
4. At sufficiently small `RV`, both old and exact laws return negative red-optical extinction. Using each released mock's stored `SIM_AV` and `SIM_RV`, historical `A(8000 A)<0` for **4,601/71,946** objects. Substituting exact F99 on the same truth gives **4,791/71,946**. Source tracing shows no clamp before the host multiplier `10^(-0.4 A_lambda)` enters the SED. This is an extrapolated dust component amplifying the source flux, not attenuation by a passive foreground screen. It is not a count of real SNe, a net broadband effect at every redshift, or a measurement of cosmological bias.
5. A controlled actual-SNANA refit of 12 deliberately selected low-RV mocks and 12 comparison mocks found a **small direct response** to flooring negative host extinction: maximum fixed-standardization change **0.000655 mag**, before BBC. Four low-RV objects had affected accepted epochs, eight low-RV objects and all 12 controls were unchanged. This stress result narrows the direct fitting concern for this fixture; it does not test selection/population retuning or the resulting bias grid.
6. BBC mean corrections and error tuning depend on the selected simulated population. They do not marginalize automatically over every plausible dust, intrinsic-color, host-age, calibration, and selection model. Recorded uncertainty scenarios and response matrices must not be mistaken for a proof that these alternatives exhaust the uncertainty.

## Reproducible artifacts

Run from the repository root with the existing pinned official environment:

```bash
phase2/env-official/bin/python scripts/salt_dust_audit/snana_extinction.py
phase2/env-official/bin/python scripts/salt_dust_audit/snana_population.py
phase2/env-official/bin/python scripts/salt_dust_audit/snana_mock_support.py
phase2/env-official/bin/python scripts/salt_dust_audit/snana_config_audit.py
phase2/env-official/bin/python scripts/salt_dust_audit/snana_stress.py prepare
phase2/env-official/bin/python scripts/salt_dust_audit/snana_stress.py run
```

Each output directory below includes SHA-256 manifests. The source snapshots are current SNANA `886408a4...`, historical SNANA `2fe0f564...`, original DES tag 1.3 `e3493cb3b9505fc1f3f392d887364b50dae21439`, and current DES/Dovekie `c9a4fcafc4cbd19bd750dee47fc76194a45c181f`.

| Artifact | Meaning |
|---|---|
| `runs/salt_dust_audit/snana_extinction/wavelength_response.csv` | Exact, old, historical and independent `extinction`-library extinction on 8 RV values and 801 wavelengths |
| `snana_extinction/monochromatic_response_matrix.csv` | `dA/dEBV`, `dA/dRV` at fixed EBV, and exact-minus-old response; **not** distance corrections |
| `runs/salt_dust_audit/snana_population/conditional_moments.csv` | Normalized parent-PDF moments/tail probabilities, by survey/scenario/mass/ZTRUE |
| `snana_population/conditional_extinction_response.csv` | Parent expectation of monochromatic extinction, conditional on independent RV and EBV draws at fixed host properties |
| `snana_population/scenario_response_matrix.csv` | Exact-minus-old plus P21sys1/2/3 and BS21-minus-P21 parent-extinction scenario vectors |
| `snana_population/scenario_cosine_matrix.csv`, `scenario_geometry.json` | Scenario-vector collinearities and singular vectors, with explicit Euclidean grid weighting |
| `runs/salt_dust_audit/snana_mock_support/negative_A8000_simulated_rows.csv` | Actual mock identifiers, stored AV/RV, and evaluated old/exact extinction at 8000 A |
| `runs/salt_dust_audit/snana_config/results.json` | All 20 file comparisons, BBC mask decoding, pinned source excerpts and line numbers |

The independently authored broadband test [`f99_support.py`](../../code/salt_dust_audit/f99_support.py) extends the monochromatic issue using the released SALT3 SED and DES filters; see [`runs/salt_dust_audit/f99_support`](https://github.com/sprajs/cosmology-revalidation/tree/17487bf659fcbdeeea072221492bac14b04a0a85/runs/salt_dust_audit/f99_support). For its explicitly constructed `z=.1`, phase-0, `x1=c=0`, `EBV=.1`, `RV=.4` case, the i-band lies inside the configured mean rest-wavelength fit range and has negative broadband host extinction. It does not provide an occurrence fraction or a fitted-distance bias.

## Version correction: what has been verified

`snana_extinction.py` extracts the contiguous extinction-routine section from actual `MWgaldust.c`; mathematical function bodies are preserved. Only diagnostic/error routines are stubbed. The historical function is renamed to coexist in one shared library. Source/header hashes and generated C are retained. Direct numerical checks give:

- Current approximation (`-99`) versus historical `99`: **exact equality** across the 6,408 grid points.
- Current exact `99` versus independently implemented `extinction.fitzpatrick99` 0.4.9: maximum absolute difference **5.81e-11 mag per E(B-V)** on that grid.

For `EBV=.1`, over rest wavelengths 3500–8000 A, maximum absolute exact-minus-old extinction is 0.00409 mag at `RV=3.1`, 0.01155 mag at `RV=2`, and 0.01592 mag at `RV=1.5`. The last case is outside the recommended domain. These are pointwise extinction changes, not changes in SALT color or standardized distance. The response scales linearly with EBV for fixed RV in these laws. Broadband integration, redshift, SALT refitting and BBC corrections can change the sign and size of a distance response.

The [Dovekie analysis](https://arxiv.org/html/2511.07517v3), Sections 5 and 7 and Appendices A–B, reports replacing the approximation, recalibrating/retraining SALT and rerunning bias corrections; the F99-only intermediate gives `Omega_m=.369 +/- .017`, whereas final Dovekie gives `.330 +/- .015`. It reports correcting calibration weights, changing the stretch-error cut to 1.15, and retaining a modified BS21 systematic whose parameter provenance was lost. Those are author results, not reproduced pipeline products here. The paper's statement about a change in Dust2Dust chi-square does not prove the physical low-RV tail is valid or establish selected-sample coverage for alternatives.

## Map assumptions and collinearities

The nominal DES P21 file's documented latent parameters are intrinsic color mean `-.0736862`, width `.05236237`; low/high-mass RV Gaussian location/width `(3.24740456,.92648388)` and `(1.66218622,.94796795)`; low-redshift EBV scales `(.11078764,.13982172)` and high-redshift scales `(.1224219,.15028344)`; intrinsic beta location/width `(2.07646589,.21728215)`. The simulation uses `alpha=.15`, beta range `.4..3`, and a global `MAG_OFFSET=-.12`. The beta here acts on generated intrinsic SALT color; it is not the BBC fitted effective beta or an independently measured `RV+1` for each SN.

The executable samples **tabulated densities**, not unbounded Gaussians specified only by those header numbers. RV grid spacing is .1, EBV .02, and probabilities are rounded. Normalized interpolated moments therefore differ from Gaussian location parameters. At `ZTRUE=.5`:

| Conditional log mass | Mean RV | P(RV<2) | P(RV<1) | Mean EBV | P(exact A8000<0) |
|---|---:|---:|---:|---:|---:|
|9.0|3.24997|0.08852|0.00682|0.121914|0.00179|
|9.5|2.55777|0.33924|0.08776|0.137289|0.03428|
|10.0/11.0|1.82024|0.60638|0.17400|0.149822|0.06890|

These are conditional **parent-distribution** calculations. The mass variable is mapped to the applicable host-library property; local simulator logs resolve it to `LOGMASS_TRUE`. `ZTRUE` is host-library truth redshift, not automatically an observed fitted redshift. Separate maps condition RV, EBV and stretch on host properties; intrinsic color has its own common PDF. This factorization is an assumption about latent dependence. Selection subsequently induces additional correlations.

The mass PDF grid has 1-dex spacing. It assigns the low-mass shape at mass 9 and the high-mass shape at 10, and SNANA interpolates probabilities between them. Thus the input realizes a mixture across 9–10 dex, not the same sharp boundary as the separately configured BBC mass step. EBV's redshift grid is `[-2.4,.1,2.6]`, so its nominal low/high-redshift header scales must also be interpreted through interpolation rather than as a literal two-bin cutoff. `grid_axes_and_config.json` preserves all coordinates and settings. This is an implementation characterization, not an assertion that the smoothing is unintended.

Integration uses Gauss-Legendre quadrature separately within each piecewise-linear PDF segment, splitting at probability-tail and zero-extinction thresholds. Comparing 4- and 16-point quadrature for the high-mass nominal extinction expectations agrees within `1.22e-11 mag`. Host-color-conditioned P21ur is inventoried but not integrated on the mass grid; no invented color-to-mass mapping is used.

On the explicitly unweighted DES mass/redshift/wavelength grid, the normalized response cosine of P21sys3-minus-P21 with BS21-minus-P21 is **-.93860**, and P21sys1 with P21sys2 is **+.78248**. The five normalized scenario columns have singular values `[1.7341,1.0609,.8943,.2432,.0916]`. These quantify overlapping perturbations in the parent extinction model. They are neither posterior correlations nor independent uncertainty components; adding each column's squared effect as if independent would impose an unsupported covariance model. Actual Hubble-distance propagation requires the measured/released covariance of nuisance parameters, or explicit correlated scenario draws through simulation, fitting, selection and BBC.

## The low-RV path reaches the generator

This is supported by both source and written truth metadata:

1. `snlc_sim.c::gen_RV` delegates to `getRan_genPDF`. The map branch obtains the map's support and samples its interpolated probability, rather than applying `INPUTS.GENGAUSS_RV.RANGE`. The printed default Gaussian range 1–5 is not evidence of a clamp when the map supplies RV.
2. `gen_AV` converts `EBV_HOST * RV` and checks `AV>=0` and `RV>epsilon`. It does not require positive extinction at every wavelength.
3. `genmag_SEDtools.c::fill_TABLE_HOSTXT_SEDMODEL` evaluates the selected law at `LAMOBS/(1+z)` and stores `pow(10,-.4*XT_MAG)` without clamping negative XT_MAG.
4. `genmag_SALT2.c` multiplies `HOSTXT_FRAC` into both photometric and spectral SED integrands. The historical source follows the same path.
5. The release's 25 raw Ia FITS headers provide the numerical AV/RV draws; `snana_mock_support.py` opens those original inputs directly. Per-object threshold checks agree exactly with the evaluated extinction signs.

At `EBV=.1, RV=.4`, exact F99 gives `A6500=-.0184755` and `A8000=-.0184297` mag. The 8000-A zero occurs at `RV=.66338574` (exact) or `.65323402` (old). These roots and all mock counts are reproducible. A passive foreground extinction interpretation fails for these particular extrapolated values, but a phenomenological population model may still reproduce selected optical summaries. Whether that produces a residual cosmological bias is an empirical transport/coverage question; negative extinction alone does not calculate that bias.

## Selection, BBC and errors

The original DES simulation includes DES cadence/depth (`SIMLIB`), Poisson/noise and flux-error corrections, a host library and weight map, detection efficiency plus host-redshift efficiency (`APPLY_SEARCHEFF_OPT=5`), and minimum epoch/SNR/foreground-EBV cuts. Written mocks have already passed generation/write selection. They are not all attempted transients. The 71,946 denominator above must not be reused as a real-survey population or final cosmology-sample denominator.

SALT light-curve fits use fixed MW RV=3.1 and `OPT_MWEBV=3`; code recomputes SFD98 and scales by .86, with .05 fractional reported EBV error. Foreground dust is multiplied at observer wavelength, separate from rest-frame host dust. The fit configuration limits mean rest wavelengths to 3500–8000 A, fit phases to -15..45 days, applies loose x1/color priors and first-pass epoch rejection at delta-chi-square 10. These define the data region in which a response is meaningful. The foreground systematic branch scales EBV by .95; its law branch switches to CCM89/RV=3.0 with amplitude weight .3. This discrete branch is not marginalization over all spatially varying MW laws/maps.

BBC `opt_biascor=4336` decomposes into `16+32+64+128+4096`: multidimensional bias mapping, covariance scale, sample grouping, correction of mu itself, and additive covariance floor. Its table configuration splits survey/deep/shallow samples, has two mass bins with `interp_biascor_logmass=0`, `x1` in [-3,3], color in [-.3,.3], fit-probability/parameter-error cuts and residual chi-square rejection. It simultaneously fits alpha, beta and a residual mass step. Therefore zero residual fitted gamma would not remove mass-dependent dust already present in bias corrections and errors.

The code forms `mu=mB+alpha*x1-beta*c-gammaDM`, subtracts `muBias`, and tunes variance using positive `muCOVadd` or, otherwise, `muCOVscale` while leaving the peculiar-velocity term unscaled. It stores the Monte Carlo `muBiasErr` but sets its local value to zero before variance application, with an explicit comment that separate addition would contradict the covariance-scale correction. This is not evidence that Monte Carlo statistics are ignored everywhere; it does mean that blindly adding the stored column again double counts the intended recipe. Finite bias-grid uncertainty, interpolation support, selection-model uncertainty, and shared simulation errors require their own coherent audit.

[Dust2Dust](https://arxiv.org/html/2112.04456v2) constrains a specified latent dust/color family using several selected-data summaries and forward selection modeling. Successful fits or closure using that same family establish model-conditioned calibration. They do not establish identification against intrinsic spectral/age evolution, alternative dust geometry, selection misspecification, or calibration directions omitted from the model family.

## Controlled nonlinear SNANA refit stress

[`snana_stress.py`](../../code/salt_dust_audit/snana_stress.py) uses original released mock realization 0001. Before examining injected outcomes it chose 12 baseline-quality mocks with RV<.65, equally spaced in redshift rank, and 12 unique nearest-redshift comparison mocks with RV between 2 and 4. The original 24-object selection is retained, including zero-response cases.

For each supported source epoch it integrates the old SNANA law with the released SALT3 SED and filter, then computes the ratio when only negative **host** A_lambda is replaced by zero. This is explicitly an artificial nonnegative floor stress, not a physical dust correction. The measured flux receives `F_true*(ratio-1)`, where `F_true=10**(.4*(27.5-SIM_MAGOBS))` is the original simulator's stored noiseless flux. Thus the original noise residual, cadence, errors, flags and redshifts are preserved, apart from float32 serialization rounding. SNcosmo supplies the broadband integration, calling the extracted actual SNANA extinction routine. SNANA's real `snlc_fit.exe` performs all paired nonlinear fits. No survey selection, parent population, Poisson error generation or BBC grid is rerun.

The fresh baseline matches the existing original-mock refit **exactly** in the compared numerical columns; the independent no-op fixture also agrees exactly. Both baseline and floor variants write all 24 objects and retain their predeclared light-curve quality proxy. The standardization diagnostic holds the release values `alpha=.16087`, `beta=3.11780` fixed: `delta_mu=delta_mB+alpha*delta_x1-beta*delta_c`. It includes no host-step, bias-correction or nuisance-parameter refit.

| Mock CID | zHEL | RV | Affected accepted epochs | Maximum injected accepted magnitude | Fixed-standardization delta_mu |
|---|---:|---:|---:|---:|---:|
|99068|.11588|.40666|11|.0010473|+.0006549|
|169230|.26880|.57772|17|.0000220|+.0000096|
|751585|.36178|.41900|12|.0000717|+.0000379|
|231104|.37702|.41253|11|.0000341|+.0000248|
|Other 8 low-RV mocks|.42443–.98741|.37118–.62774|0|0|0|
|12 comparison mocks|matched redshifts|2–4|0|0|0|

The median response across the low-RV sample is zero. Maximum absolute fitted color and mB changes are .0002481 and .0001030. These are tiny direct responses in this deliberately bounded fixture; no assertion about the survey-wide bias or its uncertainty follows. In particular, wavelengths affected at rest 8000 A often lie outside the observed DES bands at higher redshift. The zero rows demonstrate why the monochromatic low-RV count cannot be translated directly into a distance-bias fraction.

[`fit_response.csv`](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/salt_dust_audit/snana_stress/fit_response.csv) retains every object, fit success and quality outcome. [`epoch_acceptance.csv`](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/salt_dust_audit/snana_stress/epoch_acceptance.csv) maps the actual baseline/floor LCPLOT accepted/rejected flags back to all original source epochs; [`object_epoch_acceptance.csv`](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/salt_dust_audit/snana_stress/object_epoch_acceptance.csv) records the counts, band/phase support and flux/error norms. Epoch references are **zero-based**; original SNANA PTROBS headers are one-based. Rounded LCPLOT MJD values are matched within .0019 day; a repeated same-time/filter exposure pair is distinguished using its observed flux and one-to-one matching. All accepted epochs are inside the injection model's support. No accepted/rejected epoch changes occurred. Native filter means and the fitted baseline phase are recorded for independent wavelength/phase interpretation.

Files outside the source model's phase/wavelength support were deliberately left unchanged and counted; none entered these final fits. Fixtures are reconstructed from original inputs with `prepare`, and their HEAD/PHOT/namelist hashes are included in the manifest. An output-only table-option retry is preserved in the logs. Final fixture repreparation was verified byte-identical to the inputs used for the completed fits. Generated full FITS copies are local, ignored artifacts; a Git-only checkout must reconstruct them before verifying their hashes.

## What would justify correcting the cosmology

The outputs here are response operators and support diagnostics, **not a replacement mu correction vector**. To obtain one, run matched simulated truth/cadence/noise through baseline and alternative valid dust prescriptions; refit/reweight the latent population using independent goodness-of-fit summaries; apply the same detection, host-redshift, fitting/classification and BBC selection; regenerate the bias grid and error tuning; propagate the resulting correlated delta-mu vectors together with calibration/training covariance. Test coverage with held-out injected populations that do not belong to the fitted family. Changes that merely clip negative extinction or impose RV>=2 without retuning the population are useful stress tests, not scientifically justified corrections.

Required unresolved provenance includes exact executed Dovekie job configuration and baseline bias-correction simulations, the relevant Dust2Dust chains/parameter covariance with exact-law/calibration choices, and the pipeline products for the F99-only intermediate. Public post-fit distances and their covariance can support a conditional cosmology analysis, but cannot reconstruct those missing conditionals on their own.
