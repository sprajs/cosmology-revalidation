# Preparing a separately identified native accuracy-2 posterior

This is a **preparation-only deliverable**. It contains a frozen [design](../code/inference/native-accuracy-correction-design.json), tested density and spectrum-conversion contracts, and an immutable numerical-target descriptor. It does not implement or authorize the persistent physical executor or the full observational qualifier. No model, background, spectrum, posterior or new cosmological result was calculated in this preparation.

The proposed target keeps the original physical model, likelihoods, data, priors and all 2,000 selected proposal slots. Only CAMB's `AccuracyBoost`, `lAccuracyBoost` and `lSampleBoost` change from one to two. Each raw importance weight must satisfy both expressions

$$
\log w_2=\log p_2-\log p_{\rm proposal}
=\log w_1+\log p_2-\log p_1.
$$

The old denominator and accuracy-1 records remain immutable. Derived quantities such as the sound horizon, matter density and expansion diagnostics must come from accuracy two. Renaming an accuracy-1 report or passing an accuracy-2 result to a consumer that reconstructs accuracy-one settings would be incorrect.

Before any full calculation, the original parent must freshly qualify. Its fixed 32-point numerical screen must retain genuine likelihood variation after the independent thermal-path review. The same 32 points must then pass the separate accuracy-2 replay and accuracy-2-to-3 screen, including complete records, matching configurations and zero unknown-count attempts. The [read-only prerequisite function](../code/inference/native_accuracy_correction.py) reconstructs these existing checks; it does not turn a 32-point screen into a posterior or authorize further execution.

The 32 existing high-accuracy records can potentially be reused. Their archives contain **Dℓ**, obtained from `get_Cl(ell_factor=True)`, rather than raw Cℓ. Temperature and polarization spectra require a factor $2\pi/[\ell(\ell+1)]$; lensing potential power requires $2\pi/[\ell(\ell+1)]^2$. The latter is potential power, not convergence power. The converter binds the installed Cobaya conversion and unit-handling source, preserves the original archive, and refuses a nonzero temperature/polarization monopole whose units cannot be reconstructed from that archive alone. A failed conversion blocks reuse; it does not silently authorize another native calculation.

A persistent process pool needs another specific validation. The original 32 records came from separate processes, whereas a persistent model can carry numerical or cache state between points. The frozen pilot compares the first and last preselected screen points in both orders and on repeat, using two processes and at most eight log-posterior requests. It must check components, priors, derived quantities, spectra and finalized multipole settings against the existing records and record timing. **That pilot has not run and requires separate authorization.** No unrecorded warmup, cache modification or physical-target change is permitted to force agreement.

If those prerequisites hold, the later executor can retain the 32 verified records and evaluate the remaining 1,968 selected slots using at most four persistent workers, each with one computational thread. Exclusive attempt tickets and claims must prevent retries of interrupted or completed requests. Invocation counts, internal solver counts where available, and unknown interrupted counts must remain distinct. Canonical payload seals and file-hash ledgers must prevent edited records from being resealed as valid.

The eventual typed observational qualifier must recheck all identities and all 2,000 records, then apply the original weight, chain and batch gates unchanged. The current [typed descriptor](../code/inference/native_accuracy_measurement.py) deliberately labels itself `configuration_descriptor_only_not_posterior_qualification`; it cannot substitute for that unfinished qualifier. Downstream calibration, expansion, lensing and luminosity consumers remain unchanged and require separate explicit adaptation later.

## Completed synthetic checks

The [validation report](../results/inference/native-accuracy-contract-validation.json) passes 16 actual configuration-dictionary comparisons, covering ΛCDM/CPL, four luminosity assumptions and CPU/GPU factory routing. It rejects changes to priors or physical assumptions beyond the three boosts. Eight synthetic provider states compare the converter directly with the installed Cobaya method, without constructing a model: maximum relative discrepancy is $2.15\times10^{-16}$ and round-trip error $2.28\times10^{-16}$.

Forty thousand known Gaussian points verify both weight equations, fresh derived values and the unchanged importance-gate implementation. Density accounting closes within $8.89\times10^{-16}$; a deliberately concentrated chain fails the original weight-fraction gates. Mocked evidence-boundary controls reject failed refinement, changed points, bad replay identity, incomplete counts and nongenuine variation. These tests validate the contracts, not a real 2,000-point lineage, process-pool equivalence or numerical convergence across the posterior.

Reproduce the preparation checks with no physical calls:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/inference/native_accuracy_validate.py \
  --output .work/unified-cosmology/native-accuracy-contract-validation-replay.json
```
