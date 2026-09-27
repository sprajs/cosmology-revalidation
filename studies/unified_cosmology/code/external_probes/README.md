# Executable external cosmological measurements

The [scientific note](../../notes/external-probes.md) describes the observations, assumptions, numerical validation and author-configuration discrepancies. The compact records are in [results/external_probes](../../results/external_probes). No posterior inference is established merely by evaluating these likelihoods.

All commands run from the repository root. Downloaded data, third-party code, generated spectra and environments go under ignored `.work/unified-cosmology/external-probes`. Python 3.12.14 was used. Keep OpenMP/BLAS thread counts explicit when running concurrent chains.

## 1. Planck 2018 and DESI DR2

```bash
uv venv --python 3.12.14 .work/unified-cosmology/external-probes/.venv
uv pip install --no-deps --python .work/unified-cosmology/external-probes/.venv/bin/python \
  -r studies/unified_cosmology/code/external_probes/requirements-lock.txt
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.venv/bin/python \
  studies/unified_cosmology/code/external_probes/acquire.py
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.venv/bin/python \
  studies/unified_cosmology/code/external_probes/validate.py
```

`adapter.external_info(mode='full', model='cpl')` provides full Plik TTTEEE, low-T/low-E, native lensing and DESI DR2. `mode='lite'` changes only the high-multipole foreground treatment. `model='lcdm'` fixes w₀=−1,wₐ=0. Callers add their SN likelihood to the same background provider; do not free the BAO sound horizon independently in this CMB combination. The analytic CAMB CPL interface has the recorded w₀+wₐ≤0 domain.

`lensing_precision.py`, `theory_domain.py` and `table_cpl_sensitivity.py` reproduce the native/binary lensing comparison and separately labelled outside-domain diagnostics. They do not alter the baseline target.

## 2. Declared modern Planck–ACT–SPT comparison

Complete step 1 first because Planck/BAO assets are shared. The separate environment keeps these additional likelihood dependencies isolated:

```bash
uv venv --python 3.12.14 .work/unified-cosmology/external-probes/.modern-venv
uv pip install --no-deps --python .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  -r studies/unified_cosmology/code/external_probes/modern-requirements-lock.txt
uv pip install --no-deps --python .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  -r studies/unified_cosmology/code/external_probes/sampling-runtime-lock.txt
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/external_probes/restore_modern.py
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/external_probes/modern_validate.py
OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/external_probes/modern_release_checks.py
```

`modern_adapter.modern_info()` returns the explicitly declared target. It applies each calibration prior once and clears internal candl priors. Shared photons and omitted cross-likelihood covariances remain scientific limitations. `modern_acquire.py` is the initial acquisition recorder; use `restore_modern.py` for repeat acquisition/verification so the frozen pre-training manifest is preserved. Add `--verify-only` to avoid downloads.

The selected SPT datasets use NumPy. The larger installed SPT package declares JAX/GP dependencies for other, unused datasets. Those are deliberately absent; `--no-deps` prevents switching the numerical backend or downloading unpinned extras. The supplemental sampling lock adds MPI and posterior-analysis utilities while preserving the original scientific lock and acquisition identity.

## 3. Proposal matrices and original author chains

After the survey-selection lane has created `normalized/dovekie-total.npz` and acquired its pinned release tree, run `proposal.py`, `sn_proposal.py`, then `modern_proposal.py` in the corresponding environment. They construct proposals, not scientific priors. `sn_proposal.py` uses the observed redshift/covariance design, not measured magnitudes.

Run `author_config.py` to acquire released LCDM/CPL chains, extract exact headers, reproduce their recorded priors and write weighted proposal matrices. Run `author_acquire.py` in the modern environment to recover historical source/likelihood data, then `author_evaluate.py` and `author_density_check.py` to compare five actual retained rows. These are source-informed reconstructions with explicitly recorded runtime gaps, not a second independent cosmological measurement. The shared-calibration diagnostic is separate from the literal public header/source.

## 4. Numerical acceleration checks

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/external_probes/spectral_training.py --workers 4
```

This resumes the frozen 512-training/96-held-out attempted-point design with cached row hashes. Every failure is retained. Its output is exact spectra for a possible emulator, not a posterior or a new prior. Use actual CAMB acoustic-angle coordinates in the row files rather than the requested root-solving coordinates. The root inference lane owns approximation validation and exact posterior correction. `theory_trim_check.py` separately evaluates lower-multipole requests at five fixed points without changing the declared training target.

## Optional exact lensing acceleration

`fast_lensing.py` reassociates the released linear response with its binning matrix. It preserves all response terms and the full joint ACT/Planck covariance. `use_fast_lensing(info)` explicitly opts a caller-owned modern configuration into this wrapper; the default adapter is unchanged. Unsupported calibration/nonlinear variants raise an error.

Run `validate_fast_lensing.py` after the exact training spectra exist. Sixteen stored physical spectra and sixteen algebraic stress curves agree with the unmodified native likelihood to about 3×10⁻¹⁴ in log likelihood. This is matrix reassociation, not a spectral approximation. Its measured native-likelihood median time fell from 0.496 s to 0.00166 s in the recorded one-thread benchmark. It does not accelerate CAMB itself or remove the need for emulator validation and exact posterior correction.

Run `review_fast_lensing.py` in the modern environment for the separately authored 20-case native-response stress check. After all acquisition and validation steps, run `aggregate_validate.py` in that environment. It verifies frozen data/source identities, dependency versions, all 608 attempted-point records and recorded numerical outcomes; it does not rerun cosmological inference or certify an emulator. The complete spectral acquisition contains 509 finite training points, three recorded optical-depth prior exclusions, and 96 finite held-out points.
