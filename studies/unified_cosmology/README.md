# Joint cosmology and the physical supernova correction

This study combines supernova distances, galaxy BAO and CMB observations in one expansion model, while testing whether the public host and survey observations support a replacement luminosity correction. A shared cosmology is a well-defined conditional measurement. An empirically established age correction requires additional observational and selection tests.

The late-time Dovekie + DESI DR2 fit in flat ΛCDM gives **Ωₘ = 0.30627 ± 0.00768**, **H₀rᵈ = 10086.4 ± 65.1 km/s**, and **q₀ = −0.54059 ± 0.01152**. These are posterior means and standard deviations under the stated model and priors. Its four independent ensembles pass the recorded convergence gates. This row contains no CMB information and does not separately determine H₀. [Numerical record](results/inference/late-lcdm-none-dovekie.json).

The flexible late-time CPL model gives median **q₀ ≈ −0.38**, with central 95% intervals around **[−0.566, −0.203]** across two independent nested calculations. Jerk remains consistent with both signs. A small low-matter-density tail is prior-sensitive and less precisely sampled; the [tail analysis](notes/late-time-tails.md) preserves the separate results and failed initial checks. The CMB calculations are qualified separately. A fitted point, a reproduced published chain, an interpolated spectrum and a converged independently calculated posterior are different forms of evidence.

## Scientific reading order

1. [Joint inference](notes/joint-inference.md): observations, equations, nuisance brightness, acceleration and jerk, priors and interpretation.
2. [External probes](notes/external-probes.md): recovered BAO/CMB inputs, likelihood checks, calibration conventions and comparison with the published calculation.
3. [Dependence between probes](notes/probe-dependence.md): measured joint covariance, deliberately omitted correlations and the limits of existing validation.
4. [Survey selection](notes/survey-selection.md): the actual Dovekie sample, covariance semantics, classifier reproduction, host matching and historical DES3YR alternative.
5. [Host observations](notes/host-likelihood.md): aligned host posterior samples, selected optical fluxes, observed spectra and the remaining response/age ambiguity.
6. [Calibrated host spectra](notes/calibrated-hosts.md): 55 matched DESI spectra, signed measurements, overlapping-band covariance and physical limitations.
7. [Supernova release comparison](notes/release-comparison.md): separately fitted Dovekie, Pantheon+ and DES3YR alternatives with the same BAO observations.

All 17,733 released classification probabilities are reproduced to their printed precision after recovering the measured peak-time input. This resolves a specific classifier-interface discrepancy; it does not reconstruct the complete survey selection or contaminant model. The 1,088 distinct OzDES host spectra add observed high-redshift information, but their lack of relative flux calibration prevents interpreting raw count ratios as calibrated stellar ages. The independent repeat spectra also reject a simple fixed-response, diagonal-noise description for many objects.

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
