# Joint cosmology and the physical supernova correction

This study combines supernova distances, galaxy BAO and CMB observations in one expansion model, while testing whether the public host and survey observations support a replacement luminosity correction. A shared cosmology is a well-defined conditional measurement. An empirically established age correction requires additional observational and selection tests. Research execution is closed as of 27 September 2026 at the user's request; the [conclusions](../../README.md#9-conclusions) and separate [measurement gaps](../../docs/measurement-gaps.md) distinguish completed evidence from unfinished calculations. No continuation or retry is scheduled.

The joint CMB + BAO + Dovekie flat-ΛCDM fit gives **H₀ = 68.164 ± 0.251 km/s/Mpc**, **Ωₘ = 0.30373 ± 0.00340** and **q₀ = −0.54425 ± 0.00510** at the declared native numerical accuracy. Its sampling and importance-weight checks pass, but the higher-accuracy assessment detects numerical likelihood variation requiring further refinement; these uncertainties are not yet final. Replacing the supernova factor with calibrated Pantheon+SH0ES gives **H₀ = 68.514 ± 0.242**, while withholding its calibrators gives **68.188 ± 0.252**. The withheld calibration has a conditional predictive equal-tail area of **2.56 × 10⁻⁸**; this identifies a model-conditioned discrepancy, not its cause or a frequentist significance. [Joint measurement, calibration alternatives and limits](notes/joint-lcdm-calibration-results.md).

With evolving CPL dark energy and fixed standardized supernova brightness, the same combined probes give **H₀ = 67.460 ± 0.565**, **q₀ = −0.3310 ± 0.0614** and **j₀ = −0.217 ± 0.312**. The q₀ 95% interval is entirely negative, while the jerk interval allows both increasing and decreasing scale-factor acceleration. The dark-energy parameters closely reproduce the published result. Sampling and native-weight checks pass at the declared accuracy, but the higher-accuracy screen detects centered log-likelihood variation of RMS **0.184** across 32 points. These uncertainties are provisional; its parameter-level impact remains unmeasured. [CPL measurement and interpretation](notes/joint-cpl-results.md).

The late-time Dovekie + DESI DR2 fit in flat ΛCDM gives **Ωₘ = 0.30627 ± 0.00768**, **H₀rᵈ = 10086.4 ± 65.1 km/s**, and **q₀ = −0.54059 ± 0.01152**. These are posterior means and standard deviations under the stated model and priors. Its four independent ensembles pass the recorded convergence gates. This row contains no CMB information and does not separately determine H₀. [Numerical record](results/inference/late-lcdm-none-dovekie.json).

The flexible late-time CPL model gives median **q₀ ≈ −0.38**, with central 95% intervals around **[−0.566, −0.203]** across two independent nested calculations. Jerk remains consistent with both signs. A small low-matter-density tail is prior-sensitive and less precisely sampled; the [tail analysis](notes/late-time-tails.md) preserves the separate results and failed initial checks. The CMB calculations are qualified separately. A fitted point, a reproduced published chain, an interpolated spectrum and a converged independently calculated posterior are different forms of evidence.

The released geometric-anchor and Cepheid distance ladder independently reproduces its published matrix result, **H₀ = 73.043 ± 1.007 km/s/Mpc**. Removing its supernova information yields a correlated distance likelihood for **37 host galaxies**. These are conditional calibration results; the H₀ summary is not added as an independent prior to overlapping supernova samples. [Calibration, covariance and remaining interface](notes/distance-ladder.md).

The calibrated Pantheon+SH0ES sample separately gives **H₀ = 73.550 ± 1.017 km/s/Mpc** and **Ωₘ = 0.33245 ± 0.01806** in flat matter-plus-Λ cosmology, retaining the complete calibrator and noncalibrator covariance. An independent matrix and quadrature calculation reproduces the result. [Calibrated cosmology and published-chain comparison](notes/calibration-lcdm.md).

## Scientific reading order

1. [Joint inference](notes/joint-inference.md): observations, equations, nuisance brightness, acceleration and jerk, priors and interpretation.
2. [External probes](notes/external-probes.md): recovered BAO/CMB inputs, likelihood checks, calibration conventions and comparison with the published calculation.
3. [Dependence between probes](notes/probe-dependence.md): measured joint covariance, deliberately omitted correlations and the limits of existing validation.
4. [Survey selection](notes/survey-selection.md): the actual Dovekie sample, covariance semantics, classifier reproduction, host matching and historical DES3YR alternative.
5. [Host observations](notes/host-likelihood.md): aligned host posterior samples, selected optical fluxes, observed spectra and the remaining response/age ambiguity.
6. [Calibrated host spectra](notes/calibrated-hosts.md): 55 matched DESI spectra, signed measurements, overlapping-band covariance and physical limitations.
7. [Physical stellar-age constraints](notes/calibrated-host-physics.md): flexible age, metallicity and dust mixtures, native spectral response, compatible age ranges and the effect of an imposed cosmological clock.
8. [Supernova release comparison](notes/release-comparison.md): separately fitted Dovekie, Pantheon+ and DES3YR alternatives with the same BAO observations.
9. [Absolute-distance calibration](notes/distance-ladder.md): the complete released ladder, independent numerical checks, a Cepheid-only host likelihood and overlap with the cosmology samples.
10. [Calibrated supernova cosmology](notes/calibration-lcdm.md): the separate Pantheon+SH0ES result, its published-chain comparison and the unresolved covariance construction. The [anchored comparison](notes/anchored-target.md) defines its consistent replacement in the shared CMB and BAO target.
11. [Predicting the calibration](notes/calibration-holdout.md): fit CMB, BAO and noncalibrating supernovae, then test their prediction for the withheld calibrator measurements with the full released cross-covariance.

All 17,733 released classification probabilities are reproduced to their printed precision after recovering the measured peak-time input. This resolves a specific classifier-interface discrepancy; it does not reconstruct the complete survey selection or contaminant model. The 1,088 distinct OzDES host spectra add observed high-redshift information, but their lack of relative flux calibration prevents interpreting raw count ratios as calibrated stellar ages. The independent repeat spectra also reject a simple fixed-response, diagonal-noise description for many objects.

Five calibrated spectral bands in 54 of 55 independently matched DESI hosts admit formed-mass stellar ages both below 1 Gyr and above 10 Gyr across the declared flexible mixture family. This limits the age information in those compressed measurements under those assumptions; it does not exhaust the information in the full spectra. It supplies no empirically measured supernova luminosity correction.

## Reproduction

Run commands from the repository root. Source acquisition, ignored execution storage and the dependencies of each lane are specified in the [survey](code/survey_selection/README.md), [host](code/host_likelihood/README.md), [calibrated-host](code/calibrated_hosts/README.md), [external-probe](code/external_probes/README.md) and [inference](code/inference/README.md) instructions. The common interface for distance samples contains `zHD`, `zHEL`, `MU`, `covariance` and `precision`; `MU` denotes the numerical distance-like column and does not imply a Cepheid absolute calibration.

After acquiring and verifying those inputs, the independent late-time calculation can be repeated with:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
.work/unified-cosmology/external-probes/.venv/bin/python \
  studies/unified_cosmology/code/inference/late_geometry.py --validate

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
.work/unified-cosmology/external-probes/.venv/bin/python \
  studies/unified_cosmology/code/inference/late_geometry.py --model lcdm --steps 20000

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
.work/unified-cosmology/external-probes/.venv/bin/python \
  studies/unified_cosmology/code/inference/late_diagnostics.py \
  .work/unified-cosmology/inference/late-geometry/lcdm-none-dovekie
```

Existing chains require explicit `--resume`; they are not silently overwritten. Alternative sample and luminosity assumptions have separate output paths. The CMB target, direct-chain runner, numerical proposal and exact correction live in [the inference code](code/inference), with designs written before the corresponding posterior interpretation. Downloaded observations, third-party implementations, environments, training spectra and posterior chains remain in ignored `.work/`. The repository retains authored calculations, source identities, concise scientific evidence and explicit failed checks.
