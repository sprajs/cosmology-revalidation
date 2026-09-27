# Auxiliary background consistency in the native precision check

The full CAMB nonlinear-lensing calculation temporarily enables `WantTransfer` during thermal initialization and restores it afterwards. A plain background calculation from the returned parameters therefore does not necessarily reproduce the same finite-precision drag radius. The [expansion-history investigation](expansion-history.md) established this behavior from CAMB source and four stored full-native controls.

The original 32-point precision worker uses that plain background path for its auxiliary nominal-accuracy comparison. This can fail its sound-horizon consistency gate even when the full-spectrum calculation is intact. Its two full-spectrum likelihood values remain useful for testing accuracy; they must not be replaced merely because an auxiliary reconstruction differs.

The supplemental [review](../code/inference/native_precision_thermal_review.py) preserves the original screen and every native evaluation. It applies the already validated thermal adapter only to new background-only reconstructions at the same nominal and doubled accuracy. It then checks all original nominal rdrag, Ωₘ, q and j tolerances and independently closes the doubled-accuracy H(z) and rdrag against the stored full-native values. Original likelihood differences, priors, coordinates, numerical controls, point selection and spread thresholds remain unchanged. Genuine density failures and large likelihood variation remain failures.

The [independent validation](../results/inference/native-precision-thermal-review-validation.json) uses the four pinned native controls and four separately stored nominal/high native pairs. The original plain-background rdrag error reaches **1.23 × 10⁻⁶ relative**, exceeding the unchanged **10⁻⁷** tolerance for two controls. With the thermal adapter, all four rdrag, Ωₘ, q and j comparisons agree exactly at stored precision. All eight nominal/high H(z) and rdrag references also agree exactly. These checks require **20 background calculations and no new spectra**. They do not test a new cosmological posterior.

Synthetic checks preserve prior and density failures, reject nonfinite quantities before reductions, and reject incomplete or altered records before background work. A 32-record fake consumer verifies parent requalification, input seals, manifest bindings and refusal to overwrite an earlier review. A large likelihood-accuracy variation still requests further investigation after the background repair.

The original screen is never overwritten or silently marked successful. A later supplemental result records the original status and each changed auxiliary closure alongside the re-evaluated finite-screen decision. Passing remains only a bounded 32-point numerical diagnostic, not asymptotic convergence or posterior qualification. No actual 32-point screen has yet been consumed by this supplemental review.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/native_precision_thermal_review_validate.py --output .work/unified-cosmology/native-precision-thermal-validation-repeat.json
.work/unified-cosmology/external-probes/.modern-venv/bin/python studies/unified_cosmology/code/inference/native_precision_thermal_review.py --screen PATH_TO_EXISTING_SCREEN.json --output FRESH_REVIEW.json
```

The actual review first requalifies the original posterior and verifies its sealed 32-point plan, all completed records, logs, spectra, inputs and code. Missing or failed native attempts require their own investigation; a background repair cannot manufacture those evaluations.
