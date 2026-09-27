# Optical fits predicting held-out near-infrared measurements

The public inputs support a small, conditional prediction test on real supernovae. They do not yet provide a validated near-infrared prediction or an independent host-age correction. The missing J/H forward calculation now passes its numerical checks; converged synthetic recovery, real optical posteriors and the held-out outcome calculation remain unexecuted.

The smallest useful test follows the already frozen pair, DES16E2clk and DES16X3cry. They were selected as the lowest- and highest-redshift members of the ten-object DES16 cohort before predictive outcomes. Their signed optical ancestors, observing metadata, released HST photometry, original passbands and reference spectra, and the pinned M20 BayeSN model are locally available with recorded hashes. Their redshifts are 0.367 and 0.612; the fixed optical masks contain 15 and 14 r/i/z measurements, and the held-out masks contain six and five J/H measurements. J/H refers to the observed filters: the pair samples approximately 0.67–1.25 μm in the rest frame, rather than identical rest-frame infrared passbands. All ten original cohort members have both source photometry files available, covering 175 optical and 36 held-out measurements under the frozen masks.

The [existing design](bayesn-signed-optical-pilot.md) asks whether imposing the paper's external ΛCDM distance prior changes optical dust/SED inference in ways that predict the held-out measurements differently from a proper broad distance prior. Both arms retain the same trained M20 spectral distribution, passive dust support, independent optical-trigger timing prior, fixed masks and calibration. The comparison concerns two model-and-prior predictions for a selected sample. It cannot uniquely distinguish dust, intrinsic colour, timing, calibration or selection failures, or separate grey luminosity evolution from distance. None of the ten names matches the published M20 training-name list; that is a name-based overlap check, not a recovered execution-linked training roster.

The new [source](../code/bayesn_nir_forward_gate.py) and [validation record](../results/bayesn-nir-forward-gate.json) close the first remaining numerical gate recorded in the [previous engineering note](bayesn-signed-optical-engineering.md). They preserve the original scientific inputs and historical incomplete status. The new caller relocates old absolute paths and changes only the wavelength-integration resolution in the author kernel. It does not modify the historical source or any active cosmology inference.

| Numerical check | Result |
|---|---:|
| 1200→2400 integration bins, largest change across five fixed parameter cases | 0.001985 quoted measurement error |
| 2400→4800 bins, largest change; required below 0.1 | 0.00097053 quoted measurement error |
| Dynamic timing versus the official `simulate_light_curve` interface | at most 9.67×10⁻¹⁵ quoted error |
| Eager versus JIT prediction | at most 9.67×10⁻¹⁵ quoted error |
| Changing RV from 1.2 to 6 when AV=0 | exactly zero |
| Expected flux ratio for a +0.15-mag grey shift | relative discrepancy 3.09×10⁻¹⁵ |
| Six phase-knot symmetric steps, ±10⁻⁶ rest-frame day | finite; at most 4.21×10⁻⁷ quoted error |
| Elapsed time / CPU affinity | 193.85 seconds / one CPU |

The five parameter cases comprise the fiducial spectrum, both timing-prior endpoints, zero dust, and a nonzero 42-coordinate residual-SED draw with changed dust and shape. The public static interface shares the underlying author spectral kernel, so agreement checks timing, row layout and parameter delivery; it is not an independent verification of the learned physical model. Finite knot crossings do not establish differentiability at the piecewise-linear Hsiao interpolation knots. These finite tests also do not certify numerical accuracy everywhere in an eventual posterior.

The gate uses epochs and filters from the frozen metadata roster. It reads and hashes the source photometry bytes for identity, but converts only the epoch, filter and quoted-error columns; it never converts or uses the measured flux column. The separately sealed NIR payload is not opened. Quoted errors only set the numerical tolerance here; this check makes no measurement-covariance assumption. All 187 recorded source/input files are rehashed after execution. The result records Python 3.12.14, NumPy 2.2.6 and JAX 0.6.2; execution reports a CPU device. No likelihood, sampler or real-data predictive score is evaluated.

The useful next experiment is therefore an implementation and validation task, rather than another photometry search:

1. Restore the existing optical likelihood in a new execution workspace and perform the frozen synthetic recovery screen: both cadences, three spectral/dust/grey cases and two noise seeds, under both proper distance priors. Preserve failed cases. The previous 80-warmup/20-draw benchmark is a runtime observation, not evidence of convergence; the new forward time likewise does not forecast posterior mixing.
2. Require the existing four-chain diagnostics (Rhat below 1.01, bulk/tail ESS at least 400/200, zero divergences and E-BFMI above 0.3), probability-normalization checks and a stable predictive-score estimator. Verify that swapping held-out values cannot alter the optical inputs, starts, posterior or diagnostics. Propagate the same distance amplitude, timing, dust and residual-SED draw jointly across every held-out epoch and band.
3. Fit only the signed optical rows. Freeze the posterior, input identities and diagnostics before opening the held-out outcomes. Then compare each object's joint NIR-vector predictive log density between the two arms, retaining the declared score MC error below 0.05 and all failures. These released data are held out from this fit, not historically unseen. Extend to the original ten-object denominator only under the same declared rules, rather than selecting favourable objects or bands.

The conditional likelihood has quoted flux errors and the trained residual-SED covariance available. It does **not** have a demonstrated complete measurement covariance across epochs/bands, an execution-linked reproduction of the 2024 Stan dust-population analysis, or a complete population/selection likelihood. The existing design declares diagonal quoted noise as a working model, with separate same-night correlation, baseline-offset and calibration sensitivities. Such assumptions must accompany any future predictive score. Published peak dates are unsuitable independent clocks because optical+NIR fits influenced them; the frozen optical-trigger timing construction avoids that leakage. The unresolved count-rate calibration lineage must not be replaced with an invented universal correction.

Earlier work is not an equivalent held-out test. The ten-object SNooPy signed-versus-positive comparison measures a conditional response to omitted optical epochs with inherited timing. The 79-object optical-minus-infrared distance comparison uses previously fitted products and lacks their full cross-covariance. The old BayeSN execution validated optical forward algebra and short synthetic runtime only. CSP remains useful for separate checks, but 29 of its 42 RAISIN names overlap the published M20 training sample, and its inherited timing and calibration require explicit treatment before calling its NIR flux an external predictive test.

For a repeat of this numerical gate, use a preserved/restored historical workspace root and its isolated BayeSN environment. The root must contain the paths listed in the validation record; `studies/manage.py fetch` can restore individually indexed historical inputs as described in the [study execution guide](../../README.md#execution). The script verifies hashes instead of silently accepting replacement data. It refuses to overwrite either a result or an existing execution design; choose fresh paths for a replay.

```bash
TASK_ARCHIVE=/home/szymon/.local/share/cosmology-revalidation/archive/17487bf659fcbdeeea072221492bac14b04a0a85
"$TASK_ARCHIVE/runs/research_2026_09_26/bayesn_distance_identification/.venv/bin/python" \
  studies/infrared/code/bayesn_nir_forward_gate.py \
  --archive "$TASK_ARCHIVE" \
  --work .work/infrared/nir-forward-replay \
  --result .work/infrared/nir-forward-replay.json
```

The archived environment is not part of a fresh checkout. Reconstructing it and the recorded upstream model/filter assets remains an explicit prerequisite; the current main cosmology environment is not silently substituted. No full synthetic or observed posterior was started by this continuation.
