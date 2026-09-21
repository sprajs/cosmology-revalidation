# Classification reconstruction: paused handoff

Saved 2026-09-22 approximately 00:02 UTC in response to the user's explicit request to save and pause. **No classifier predictions, probability comparisons, or simulation classification outcomes have been computed.** No background classification or installation processes remain. The independent flux audit is completed and preserved separately.

## Objective and owned paths

Reproduce published SNNV19 probabilities for 1,635 DES objects from public raw calibrated flux, then apply the same classifier and cuts to nine forward Ia simulation variants. Work belongs in `scripts/phase2/classification/`, `phase2/classification/`, and `docs/phase2/classification.md`. Do not modify completed flux results. No global installs, publication, or contact.

## Acquired and verified

* Public nominal classifier directory: `phase2/official/inputs/SNDATA_ROOT/models/classifiers/DES-SN5YR/SNNTRAINV19_z_TRAINDES_V19/`. It contains `model.pt`, `cli_args.json`, `data_norm.json`, configuration, training job, compressed log and simulation predictions. The archive's `hash.txt` is a Pippin task hash, **not a SuperNNova source revision**.
* Source checkout: `phase2/classification/sources/SuperNNova`, https://github.com/supernnova/SuperNNova, detached at `fcf8584b64974ef7a238eac718e01be4ed637a1d`. The primary project README explicitly recommends this commit for reproducing DES analyses; the current README also describes the historical SNANA_DES5yr branch.
* Entirely project-local Python 3.10.21 environment: `phase2/classification/env/`; managed interpreter resides under `phase2/classification/python/`. Packages installed successfully, with `torch==1.13.1+cpu`, `numpy==1.23.5`, `pandas==1.5.3`, `scipy==1.10.1`, `astropy==5.3.4`; all five import successfully. Full freeze is `requirements.lock.txt`; installer log is `install.log`.
* This is **not the exact historical binary environment** (old examples use Python 3.7 / torch 0.4.1). Mathematical inference agreement remains to be tested; no agreement is claimed.
* Inputs and source hashes are recorded in `handoff-manifest.json`.

## Source findings to retain

1. Nominal V19 is a vanilla bidirectional two-layer LSTM, hidden size 32, mean pooling, two classes, dropout 0.05. Full inference sets random length/redshift false, uses eval mode and one sample for vanilla. `supernnova/validation/validate_onthefly.py:classify_lcs` is the simplest official inference entry point and accepts a raw pandas table.
2. Training data job used `--phot_reject PHOTFLAG --phot_reject_list 8 16 32 64 128 256 512 --redshift_label REDSHIFT_FINAL --norm cosmo`. Saved second-stage training cli_args instead has redshift_label `none` and phot_reject null, so blindly using only cli_args for raw preprocessing would omit required transformations. Replace HOSTGAL_SPECZ and HOSTGAL_SPECZ_ERR with REDSHIFT_FINAL and REDSHIFT_FINAL_ERR before inference. Photometry rejection removes observations with any listed flag bit. Verify exact released data settings in Pippin files before asserting reconstruction.
3. No phase/time-window cut is active when `photo_window_files` is null, despite min/max defaults. Classifier uses raw flux, **not SALT fit epochs or SALT fit parameters**.
4. `supernnova/data/make_dataset.py:pivot_dataframe_single_from_df` groups observations into 0.33-day windows starting at each group first epoch, chooses lowest FLUXCALERR for repeated filter within a group, pivots, zero-fills missing filters, computes elapsed time, and casts float64 to float32. Inspect `compute_delta_time` before implementing the wrapper; it has not yet been inspected.
5. Filter order is `[g,i,r,z]`. The 28 all_features are eight flux/error values, delta_time, four host redshift values, then 15 filter-combination indicators. The training feature list is 26 entries (omits host photo-z), but actual code selects indices in **all_features order**, so host spectroscopic z features precede one-hot bands. Do not naively use saved training_features ordering.
6. `norm=cosmo` clips input values to saved minima (flux -2000, errors/time 0), divides all eight flux/error feature entries by the per-light-curve maximum across those entries, and log-standardizes delta_time. Saved delta_time min=0, mean=0.7820349335670471, std=2.481628894805908. Saved flux/error mean/std are not applied by this cosmo branch. Exact normalizer lives in `supernnova/utils/training_utils.py:normalize_arr`.
7. Onthefly get_settings reads cli_args, applies no_dump and points to model directory; check creation behavior before running to avoid historical `/scratch/project2` paths. It applies compatibility fields sntype_var and additional_train_var around ExperimentSettings construction; compatibility is not yet exercised. No model torch.load has yet been attempted.
8. Onthefly packing sorts by sequence length and reverts sorting. It returns SNIDs and predictions shaped n x 1 x 2. Verify class mapping from official code instead of guessing probability column. Batch invariance should be tested on a small subset.
9. Assets documentation describes SNNV19 trained with SALT2 Ia and revised peculiar-Ia/core-collapse populations; parent forward simulations are SALT3 Ia. This is a relevant domain shift, not evidence of misclassification. The release classification README appears to swap prose descriptions of J17 and DESCC columns; resolve against model/source configuration if discussing alternatives.

## Resume sequence (not executed)

1. Read `docs/phase2/plan.md`, current `docs/phase2/selection.md` if present, `3_CLASSIFICATION`, and nominal data/classification Pippin settings. Write and hash a preregistration **before opening any reproduced probability outcomes**. Proposed tests: all-1635 matching count; probability abs/max/RMS and rank agreement; >0.999 and >0.5 membership agreement; probability/logit differences versus z, S/N and sampling; deterministic batch-size invariance. State tolerances before running, and preserve mismatch diagnostics rather than retune against labels.
2. Read official `compute_delta_time`, raw header parsing, model class mapping, and ExperimentSettings path behavior. Build a wrapper using the pinned official onthefly inference. Read HEAD and PHOT through their explicit PTROBS_MIN/MAX bounds; verify unique HEAD CIDs, count, pointer bounds, row boundaries, and no delimiter rows. Do not use early DUMP duplicate CIDs as unique objects.
3. Verify tiny-batch load/inference on CPU with no writes outside owned paths. Use `phase2/classification/env/bin/python` and add `phase2/classification/sources/SuperNNova` to sys.path in the wrapper. Save schema and preprocessing counts. Verify official dataset preprocessing vs onthefly on the same subset if practical.
4. Run 1,635 released DES raw light curves against `sources/repos/des-science__DES-SN5YR@1.3/3_CLASSIFICATION/DES_classification.csv` nominal PROB_SNNV19. Raw files: `sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES/DES-SN5YR_DES_{HEAD,PHOT}.FITS.gz`.
5. If agreement is adequate, classify all detected HEAD objects in each of nine variants P21, BS21, G10, P21_dmplus010, P21_dmminus010, P21_dmz020, P21_rho000, P21_rho090, P21_noisetrue120. Files are `phase2/literature/simulations/outputs/PH2_pilot02_<MODEL>/PH2_pilot02_<MODEL>_{HEAD,PHOT}.FITS`. Export CID, pIa, >.999 flag and preprocessing counts separately by variant. If public pipeline cannot be reproduced, record exact missing/version/statistical barrier and perform only clearly identified feasible checks.
6. Coordinate with `/root/forward_discrimination`: their tables are `phase2/hierarchy/forward/<MODEL>-fitted.csv.gz`, with CID, generated_attempt_index and basic_quality_pass. Join **unique HEAD/fitted CID only**. Their classifier-matched sensitivity should be preregistered before inspecting outcomes. Real 1,063 objects with SNNV19>.999 and perfect-Ia simulated truth are not equivalent cuts.
7. Quantify classifier acceptance and redshift dependence for forward variants. Selected spectroscopic Ia (241+41) and two CC in released Hubble data are not an unbiased probability calibration sample. Same-flux classifier scores are selection information, not a second independent likelihood. No cosmological model or desired residual determines a correction.

## Commands and installation provenance

Already completed:

```bash
git clone https://github.com/supernnova/SuperNNova.git phase2/classification/sources/SuperNNova
git -C phase2/classification/sources/SuperNNova checkout fcf8584b64974ef7a238eac718e01be4ed637a1d
UV_PYTHON_INSTALL_DIR=/home/szymon/Documents/ChatGPT/supernova/phase2/classification/python uv venv --python 3.10 phase2/classification/env
uv pip freeze --python phase2/classification/env/bin/python > phase2/classification/requirements.lock.txt
```

Environment can be recreated with the same local interpreter and `uv pip install --python phase2/classification/env/bin/python -r phase2/classification/requirements.lock.txt --extra-index-url https://download.pytorch.org/whl/cpu`. The CPU wheel index is needed for torch's +cpu version. No outstanding process needs killing; process inspection after install showed only the inspection command itself.

Completed independent flux audit is documented in `docs/phase2/independent-flux.md` with artifacts under `phase2/independent_flux/`. No changes were made there during this classification task.
