# Independent review of the actual-SNANA stress experiment

The written artifacts support the reported **small direct fitting response in this deliberately selected fixture**. The largest fixed-standardization shift is +0.0006548857674 mag before BBC. This review found no discrepancy in the selected IDs, fitting results, injection arithmetic, epoch correspondence, or reported acceptance-mask stability. It does not verify a survey-wide dust correction, regenerated simulation population, BBC bias, or cosmological result.

The review is independent of `snana_stress.py`: it parses the raw FITRES/LCPLOT/FITS products directly and uses a separate matching algorithm. It does not rerun the fitter or independently recompute the broadband dust ratios. [The review JSON](snana-stress-review.json) records the tested script/manifest hashes, numerical results, and the complete independent checking code. No stress source files, fixtures, fit outputs, or manifests were modified.

## Checks against the authoritative artifacts

| Check | Independently verified result |
|---|---|
| Input integrity | Every entry in the stress manifest, intervention contract, and fixture-repreparation check matches its current bytes |
| Selection | Regenerating the baseline-quality selection, 12 evenly spaced low-RV redshift ranks, and unique nearest-redshift controls reproduces all 24 selected CIDs |
| Complete output set | Baseline, no-op, and floor each contain exactly those 24 unique CIDs |
| No-op/reference identity | All **104 shared numeric FITRES columns** are exactly identical between baseline, no-op, and the existing phase2 reference for the selected objects |
| Standardization calculation | Directly recomputing `delta_mB + .16087*delta_x1 − 3.11780*delta_c` agrees with `fit_response.csv` to 6.3×10⁻¹⁷ mag |
| Fit/quality outcomes | All 24 baseline and floor outputs pass the predeclared limited quality proxy |
| Nonzero responses | Only CIDs 99068, 169230, 751585 and 231104 change the reported fixed-standardization quantity |
| Controls | All 12 RV∈[2,4] controls have zero reported fitted response; eight of the 12 low-RV objects also have zero response |
| Actual fitter execution | All three saved `snlc_fit.exe` run records have exit status zero and corresponding fit products |

The fixed coefficients are the original release values. The diagnostic does not refit alpha or beta, change the host step, apply BBC corrections, or estimate a nuisance-parameter posterior. Using differences removes the fixed amplitude normalization offset; it does not turn this expression into a final corrected-distance result.

## Epoch mapping and fit masks

LCPLOT's own column description identifies DATAFLAG=1 as used data and −1 as data rejected in the fit; DATAFLAG=0 denotes a model point and was excluded. For each object and band, this review assigned LCPLOT data rows to original PHOT rows with a global one-to-one linear assignment using both timestamp and measured flux. This differs from the source script's sequential nearest-flux matching. The independent assignment reproduces **every saved baseline/floor epoch flag**.

There are 929 accepted epochs in each fit variant. The intervention affects 58 source epochs in total, of which 51 are accepted, spanning the four responding low-RV objects. No accepted epoch lacks injection-model support. Baseline and floor have exactly the same accepted source-row set: **zero acceptance-mask changes**. The largest match-time difference is 0.001900000003 day, consistent with rounded LCPLOT output. Four time-ambiguous LCPLOT rows across the two variants correspond to the repeated-time/filter case documented by the source agent; measured flux and one-to-one matching resolve the correspondence.

This verifies the frozen fixture's fitter-mask behavior. It does not imply that detection, host-redshift, classification, or BBC selection would remain unchanged in an end-to-end regenerated survey.

## What is preserved, and the numerical limits

Baseline and no-op PHOT table bytes are identical to one another. Their logical columns equal the original source data. All variant HEAD logical columns equal the source as well. FITS rewriting changes some string padding from spaces to null padding, so full table byte identity to the compressed original must not be claimed. This is a serialization distinction, not a change to numerical metadata or interpreted strings.

In the floor fixture, all PHOT columns other than FLUXCAL and SIM_MAGOBS are unchanged, including MJD, BAND, FLUXCALERR, PHOTFLAG, sky/read-noise metadata, zeropoints, gain and PSF properties. Noninjected epochs retain their original FLUXCAL and SIM_MAGOBS. Independently evaluating the intervention formula from original stored SIM_MAGOBS and the recorded band ratio reproduces the new FLUXCAL and SIM_MAGOBS exactly after the script's float32 rounding.

The maximum FLUXCAL rounding error relative to the ideal injected value is 4.58×10⁻⁵ in FLUXCAL units. If the new noiseless flux is reconstructed from the separately rounded, updated SIM_MAGOBS, the largest difference in the inferred noise residual is 9.23×10⁻⁴ FLUXCAL, or 1.18×10⁻⁴ of the preserved flux error. Thus “noise preserving” is accurate for the stated deterministic counterfactual up to small serialization effects; it is not a newly drawn statistically self-consistent observation under a changed signal. Poisson errors, detection probabilities and population selection were deliberately not regenerated.

The band ratio is calculated with the released deterministic SALT3 surface, stored generator parameters, and the historical SNANA host-extinction routine. The intervention floors negative host extinction; it is an artificial stress direction rather than an independently supported physical dust law. The ratio calculation does not reproduce every original simulator detail, including the small C11 scatter realization and late-time extrapolation. Epochs outside its support are left unchanged; none of those epochs enters these final fits. These limitations are declared in the contract and implementation report.

## Interpretation

The monochromatic low-RV finding and this small fitting response are consistent. An extinction anomaly at a specified rest wavelength does not guarantee a comparable change in an observed passband; at higher redshift, the relevant red rest-frame wavelengths can lie outside the observed DES coverage. In this fixture, only four low-RV objects have affected accepted observations. The maximum fixed-standardization shift of about **0.655 mmag** is therefore a result for this intervention and sample, not an upper bound on population evolution, the bias-correction grid, or cosmology.

The source documentation states those boundaries clearly. The next scientific test would retain the paired discipline but regenerate selection and BBC under physically defensible, population-retuned dust alternatives. This review provides no basis for automatically applying the stress shift as a correction or declaring the whole analysis bias-free.
