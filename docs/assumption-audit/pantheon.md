# Pantheon and host-age assumption audit

2026-09-21. This audit preserves all existing results and excludes the separate SALT/SNANA correction work. It inspects the age measurements used alongside Pantheon+, including the original age-fitting dependency. It does not refit host photometry or turn a historical pipeline defect into an unsupported cosmological correction.

## 1. Confirmed: the original R19 code integrates a different SFH from the one fitted to the SED

This is a numerical implementation inconsistency, not just a choice of dust prior. In the preserved [MC-Age revision 92713be](https://github.com/benjaminrose/MC-Age/blob/92713be96a89da991fe53bffcc596a5c0942fc37/calculateAge.py#L99), `runFSPS` passes `sf_slope=tan(phi)` directly into FSPS (`sources/repos/benjaminrose__mc-age/calculateAge.py:108–122`). The posterior conversion again uses `tan(phi)` and passes it unchanged to `integrate_age` (`:763–782`). But the latter integrates a late-time rate `K+s*u` (`:821–839`).

The exact FSPS revision named in [R19](https://arxiv.org/abs/1902.01433), `ae31b2f63d865354ce944e5c22eba6e93e01e67d`, uses **`K*(1+s*u)`**, where `s` is the fractional change per Gyr. In both expressions, `u` is time after the transition and `K=T*exp(-T/tau)` in the common early-phase normalization. Multiplying FSPS's early and late rates by the same `tau` brings the two early rates into identical units and leaves this late-phase discrepancy intact. Neither normalizing total mass nor converting Gyr to years repairs it.

The primary source trace is:

- [Historical `sfhinfo.f90:123–158`](https://github.com/cconroy20/fsps/blob/ae31b2f63d865354ce944e5c22eba6e93e01e67d/src/sfhinfo.f90#L123) computes the relative-slope late rate and mass.
- [Historical `csp_gen.f90:123–147`](https://github.com/cconroy20/fsps/blob/ae31b2f63d865354ce944e5c22eba6e93e01e67d/src/csp_gen.f90#L123) uses that late mass fraction in the actual composite spectrum; this is not an unused metadata formula.
- The same file, `:295–322`, explicitly describes fractional slope, flips the forward/lookback sign, and converts inverse Gyr to inverse years.
- Local immutable copies are `runs/assumption_audit/pantheon/sources/fsps-ae31b2f-{sfhinfo,csp_gen,intsfwght,sfh_weight,sps_vars}.f90`.

The independent script [pantheon_age_semantics.py](../../scripts/assumption_audit/pantheon_age_semantics.py) executes the actual preserved pure MC-Age functions, checks what they pass into a dummy FSPS parameter object, and compares separate analytic integrals and quadrature with all breakpoints supplied. It includes both positive and negative slopes, early and late transitions, a zero-slope control, and a case strictly inside the original prior bounds.

At `z=.05`, R19 Table 5 Model 2 (`tau=.5`, start `1.5`, transition `9`, slope `15`) is inside the priors. MC-Age returns **1.41484 Gyr**, while integrating the fitted FSPS SFH gives **10.68019 Gyr**. The original code's late mass fraction is about 0.998; the actual fitted model's is about 0.00124. Model 7 changes **2.40259 → 6.82501 Gyr**. Zero-slope Model 4 agrees at **4.28826 Gyr** in the independent integrals. These are parameter-level demonstrations, not estimates of catalogue-wide bias. The source's unsplit quadrature differs by `4.3e-5 Gyr` in the zero-slope case; that small numerical issue is distinct from the many-Gyr semantic mismatch.

The authors' [original posterior archive](https://zenodo.org/records/3875482) strengthens this beyond a hypothetical unexecuted-code concern. All eight `circle` validation chains contain 1,020,000 rows each. Recomputing the old age formula reproduces the archived age column with median absolute differences below `5e-7 Gyr` (occasional quadrature outliers are preserved in the output). Thus the released age products actually used the discrepant conversion. Re-expressing each archived SFH under FSPS semantics changes the validation-chain posterior median by up to **+1.29348 Gyr** (Model 7). The large Table 5 parameter-level counterexample and the posterior shift are different quantities and should not be conflated.

The circle archive matches the historical Astropy radiation default, `Tcmb0=2.725 K`. The observed SN5916 chain instead matches the later `Tcmb0=0` default. Both use `H0=70, Omega_m=.27`; accounting for this small version difference is necessary to test the exact age columns. It does not explain the large SFH semantic difference.

**Propagation limit:** [C25 §3](https://academic.oup.com/mnras/article/538/4/3340/8098234) says it follows R19's method and the same metallicity/dust/SFH settings and priors, with an updated FSPS version. Its exact executable age-conversion code and revised SFH posterior draws have not been recovered here. This audit therefore establishes an error in the released R19 dependency and original archived outputs, and a concrete unresolved dependency for C25; it does **not** establish that the revised C25 tables contain the identical error. Both the pro-age and corrected-Pantheon residual regressions using those C25 tables inherit that uncertainty.

The known [July 2026 C25 erratum](https://doi.org/10.1093/mnras/stag1210) fixes a wrong G11 age column plotted in Figure 1 and its wording; the authors state that their quantitative results are unaffected. That plotting correction neither identifies nor resolves the different SFH-to-age inconsistency demonstrated here.

### Original archived posterior support also requires checking

The original global archive has **103 unique hosts**, exactly matching the 103 rows of its deposited age summary and host-photometry file. C25's revised R19 table has 102; original CID15459 is the sole missing object. The reason for that omission is not established here. `r19-sample-crosswalk.json` records the comparison; the full audit keeps the original 103-object denominator.

The first complete stream discovered six hosts with nonfinite reprocessed ages. Inspection showed finite archived numbers but parameter draws outside the stated prior, including **negative tau**. CID18604 alone contains 571 negative-tau draws. The actual preserved `lnprior` returns `-inf` for the inspected example. One such parameter tuple repeats for all 425 retained steps of a walker, while another remains invalid for the first 146 retained steps before moving. This is not a numerical reason to quietly discard a few NaNs.

The archive shape is significant: 102 hosts have `1,020,000 = 2400*(675−250)` rows; CID20048 has 1,019,947, a 53-row shortfall whose origin is not established. The source resets the initial search at `calculateAge.py:378`, creates new unconstrained Gaussian initial positions at `:494–496`, runs 675 steps at `:470–505`, slices off 250 steps at `:547`, and overwrites the trimmed file at `:561` before appending ages at `:804`. These are the intended retained products, rather than the temporary 675-step or initial-search files. The exact execution history of each invalid walker is unavailable, but their inclusion contradicts the support of the stated posterior. There is no final finite-prior filter in that write path.

The authoritative second stream, `pantheon_age_prior_validation.py`, checks **every** draw against all published-code prior bounds, records total/valid/excluded counts and example invalid rows, and reports age changes conditional on that explicitly identified valid subset. It also compares old full-sample and valid-subset medians, so the impact of removing invalid draws is kept separate from the SFH correction. The first stream is preserved under `runs/assumption_audit/pantheon/preliminary/` and is not the authoritative finite summary.

The completed audit covers **105,059,947 archived draws**. It finds **1,288,021 outside the stated prior (1.226%)**, spread across **86 of 103 hosts**. The worst individual host, CID21502, has 13.770% outside support. There are 5,186 negative-tau draws; the largest violation category is the angle bound. Categories overlap and should not be added to obtain the excluded total. No raw parameter or archived-age values were nonfinite; the first-pass NaNs came from attempting the physical analytic formula on those invalid parameter values.

On the 103,771,926 explicitly valid draws, re-expressing the age using the **same fitted FSPS history** gives:

| Catalogue diagnostic | Result |
|---|---:|
| Mean change across host posterior medians | +0.88175 Gyr |
| Median change across host posterior medians | +0.41962 Gyr |
| Range of host median changes | −0.43584 to +5.97984 Gyr |
| Hosts with absolute median change >0.5 / >1 Gyr | 50 / 29 |
| Hosts moving across the 4-Gyr median threshold | 5 (24 below before, 19 after) |
| Largest median change from excluding invalid draws alone | 0.07256 Gyr |

The large SFH effect is therefore not explained by the invalid-draw filtering. CID6962 changes from 3.15998 to 9.13983 Gyr on its valid draws. These are deterministic transformations of original R19 posteriors, not newly validated stellar ages, updated C25 values, or a cosmology result. The [figure](../../runs/assumption_audit/pantheon/r19-age-audit.png) and per-host CSV retain the uncertainty intervals and invalid fractions.

“Valid” here means accepted by the stated prior, not proof of chain convergence or correct posterior weights. Filtering cannot repair the surviving posterior. A [bounded independent review](../../runs/assumption_audit/des/r19-prior-peer-review.json) checked all aggregate summaries, original-source prior boundaries and examples, the raw short CID20048 member, and additional native-FSPS quadrature; it did not reparse the entire 105-million-draw archive.

On valid draws, every host's median absolute reconstruction error for the **old** age formula is below `8.2e-10 Gyr`; the largest host 99th-percentile error is `0.000134 Gyr`. Rare original quadrature outliers reach `0.53236 Gyr`. This precision statement concerns typical archived-column reproduction, not an assertion that every archived numerical age is correct.

The largest outlier is independently reproduced in `r19-quadrature-outlier.json`: CID3674, zero-based row 210652, has archived age **10.04600872948 Gyr**. Executing the preserved source's unsplit quadrature gives **10.04600872958 Gyr**, while both the analytic old-formula integral and separate quadrature with SFH breakpoints give **9.51364893447 Gyr**. This confirms a separate numerical integration error on a valid draw, distinct from the relative-versus-additive slope mismatch. It does not overturn the much smaller typical reconstruction errors or establish a comparable effect on host posterior medians.

The original joint draws also directly exhibit the covariance discarded by an age-only summary: valid-draw posterior dust–metallicity correlations range from −0.893 to +0.528, with absolute correlation above 0.5 in 33 hosts. With the consistent age integral, age–dust correlations range from −0.701 to +0.122, and age–metallicity from −0.769 to +0.417. These are within-host **model posterior** correlations, not between-host population correlations or proof of a physical dust mechanism. They establish why a joint host-property posterior carries information absent from a table of separate errors.

## 2. Confirmed model dependencies: the host-age axis already contains dust and cosmology assumptions

The public source, cited as the C25 methodological dependency, sets:

| Assumption | Exact preserved code | Why it matters |
|---|---|---|
| Young-star dust amplitude fixed at twice the diffuse component; both act on young stars | `calculateAge.py:99–106` | Dust geometry is restricted while age and SFH are inferred. Marginalizing dust amplitude does not marginalize over this geometry or all attenuation laws. |
| `dust2 ~ Normal(.3,.2)` truncated to `[0,.9]` | `calculateAge.py:245–275` | The age posterior is conditional on a particular dust distribution. It is not an independent dust-free measurement. |
| `logZ/Zsun ~ Normal(-.5,.6)` with hard limits | same | Metallicity and age are inferred from the same five-band SED. The R19 paper says sigma `.5`, whereas the preserved implementation uses `.6`; this is a reproducibility discrepancy, not proof of the precise code used for every later table. |
| Fixed `H0=70`, `Omega_m=.27`, cosmological age used in the SED and SFH bounds | `calculateAge.py:46–47,114–122,238–250,663–699` | Ages and their redshift-dependent support already depend on a baseline cosmology. Changing only the Hubble-residual reference cosmology does not test this upstream dependency. |
| Seven parameters fitted to five magnitudes, with bounded SFH and informative dust/metallicity priors | `calculateAge.py:177–212,244–275` | Priors and parameter restrictions supply information needed to break degeneracies. MCMC convergence cannot validate that information. |

The R19 paper itself reports failures to recover mixed and bursty SFHs in its self-consistency experiment (`papers/text/1902.01433v1.txt:454–487`). The conversion error above complicates that validation: some reported input ages are ages of a different SFH from the SED generator. These validation failures cannot simply be interpreted as observational noise or as a validated error model for all host ages.

There is a separate **dust-equation/source discrepancy**. R19 Eq. 3 (`papers/text/1902.01433v1.txt:299–311`) prints young-star attenuation `tau1*(lambda/5500)^(-.7)` and identifies `tau1` with `dust1=2*dust2`. The historical FSPS defaults instead set `dust1_index=-1`, `dust_index=-.7` (`sps_vars.f90:446–448`), and the actual [dust routine](https://github.com/cconroy20/fsps/blob/ae31b2f63d865354ce944e5c22eba6e93e01e67d/src/add_dust.f90#L55) attenuates young stars by **both** components. Its young-star optical depth is `dust1*(lambda/5500)^(-1) + dust2*(lambda/5500)^(-.7)`; old stars see the second term alone (`attn_curve.f90:26–29`). The preserved MC-Age source changes neither exponent. At `dust2=.3`, young-star V-band extinction from the executed expression is 0.977 mag, versus 0.651 mag from the printed expression if its stated parameter mapping is taken literally. This is a difference between the stated and implemented host model, not a measured error in any SN extinction correction. Using the public implementation reproduces the latter model; testing another attenuation law requires refitting the SEDs.

## 3. Foreground extinction propagates to the age axis and should induce shared errors

`sources/repos/benjaminrose__mc-age/util.py:124–178` subtracts Milky Way F99 extinction with fixed `Rv=3.1` from host `ugriz`, using SFD-map reddening at nominal wavelengths `[3543,4770,6231,7625,9134]` Angstrom. It updates the five magnitudes but leaves `err_u` through `err_z` unchanged. `redoGupta.py:203–225` then passes those unchanged errors into a diagonal five-band likelihood (`calculateAge.py:199–212`). No separate foreground-map nuisance or photometric covariance is present in that source path.

The 0.86 map recalibration is **not missing**: [`sfdmap.SFDMap` applies it by default](https://github.com/kbarbary/sfdmap#detailed-usage). The identified issue is unpropagated uncertainty and fixed-law assumptions. Nor does correcting the host and its SN for the same Galactic foreground constitute double correction: they are different flux measurements. It instead creates a possible **age–distance cross-covariance**, because a map/law perturbation changes both the inferred host age and the standardized SN distance. Pantheon distance covariance alone does not contain the missing age derivative.

For a common reddening perturbation `dE`, host magnitudes change as `dm_b=-R_b*dE`. It therefore contributes an outer-product band covariance `R R^T Var(E)`, not independent errors per band. A common map-scale parameter also correlates hosts across the sky and can correlate host-age errors with SN residual errors. A joint propagation needs host SED derivatives/posterior reruns and matching SN systematic-response vectors; the published age medians/errors are insufficient. This is a verified omission in the inspected dependency, with **unknown numerical impact on the C25/Pantheon age slope**.

Using the law at one wavelength per broad filter is another approximation. The [`sfdmap` documentation](https://github.com/kbarbary/sfdmap#how-do-i-get-extinction-at-a-specific-wavelength-or-in-a-specific-filter) explains that true bandpass extinction depends on the source SED and filter integration. Its size in this catalogue has not been computed; it should not be called a detected astrophysical bias from inspection alone.

## 4. Gaussian summaries and uncertainty calibration remain conditional

The C25 paper repeatedly labels the **13.6th–86.4th percentile interval** an average `1 sigma` range (`papers/text/chung2025-published.txt:193,397–398,449`). For an actual normal distribution, those quantiles are **±1.09847 sigma**, not ±1 sigma. R19's public age-summary code uses the conventional 16th/84th quantiles (`calculateAge.py:790–791`; `util.py:22–45`). The unavailable C25 posterior arrays prevent determining whether this is a manuscript typo or the executed calculation. Automatically dividing the released errors by 1.09847 would therefore be unjustified.

Independently, interpreting the 196 matched summary ages as normal posterior distributions assigns more than 5% probability outside `[0, cosmic_age(z)]` to **51 hosts** under either overlap policy. Summed negative-age probability is about `2.70`, and summed probability older than the baseline universe is about `5.11` (G11-first). These sums are expected counts under that **Gaussian approximation**, not observed unphysical hosts or the tail mass of the unavailable real posteriors. The existing LINMIX/full-covariance fits are explicitly approximations and should remain so. Their normal latent-age likelihood is not made physically exact by using a correct SN covariance matrix.

## 5. Existing Pantheon safeguards confirmed, with limits preserved

The existing age investigation already addresses several major traps correctly: joins restricted to SDSS survey 1; 33 G11/R19 overlapping IDs and unique-SN policies; full-vs-diagonal covariance sensitivity; tabulated error/covariance-diagonal disagreement; revised magnitude deltas taken from `m_b_corr`; and distinctions between reversing the whole bias correction and removing a scalar mass step. See [the existing age investigation](../investigations/age-signal.md).

Remaining assumptions are explicit in the current implementation: LINMIX receives no age–HR measurement cross-covariance (`scripts/age_signal/analyse.py:108`; `covariance_linmix.py:19`); the exploratory joint model holds measured mass and SN colour fixed (`fullcov_eiv.py:19,29–35`); and the reversed-bias variants inherit the corrected-distance covariance. These are conditional sensitivity models, not a complete joint age/dust/selection likelihood. A fixed-age GLS or a normal age likelihood cannot determine the missing derivative, SED prior, host-selection normalization, or covariance of the same-sample G11 redshift correction.

## Reproducible evidence

Run `OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/assumption_audit/pantheon_age_semantics.py` for the semantic counterexamples and archived validation checks, then `OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/assumption_audit/pantheon_age_prior_validation.py` for the authoritative full original-global-chain audit. Run `.venv/bin/python scripts/assumption_audit/pantheon_dust_contract.py` for the host-dust contract checks.

The semantic outputs are `runs/assumption_audit/pantheon/age-semantics.json`, `sfh-semantics-cases.csv`, and `age-semantics-manifest.json`. Full valid-subset results are `r19-global-validated-age-audit.csv`, `r19-global-validated-age-summary.json`, `r19-global-prior-violations.json`, and `r19-global-validated-age-manifest.json`. Both original archives are verified against their deposited Zenodo MD5 values in `sources/archive-checksums.json`; `sources/acquisition.json` records URLs and SHA-256. Large original chain bundles remain local and are Git-ignored. No source tables or prior analysis products were replaced.
