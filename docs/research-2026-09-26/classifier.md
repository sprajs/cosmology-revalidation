# Nominal DES SNNV19 reconstruction: failed release reproduction

The released V19 model and available DES raw photometry do **not** reproduce published `PROB_SNNV19` closely enough for forward-simulation selection work. The original, prospectively registered failure is preserved. A later source-motivated diagnostic found an omitted Pippin prediction-stage peak window and narrowed the discrepancy, but still failed the unchanged adequacy gate. No probabilities or classifier acceptance were computed for the nine forward variants; no classifier likelihood was added.

The [preregistration](../../phase2/classification/reconstruction_20260926/preregistration.md) was written before prediction (SHA-256 `925c07d183ae9a4b73cad5337e2210320301133aafa9be8f3145cde1f44591a8`). The [primary output](../../phase2/classification/reconstruction_20260926/des_probabilities.csv) and [metrics](../../phase2/classification/reconstruction_20260926/des_metrics.json) remain unchanged. The attempt used pinned SuperNNova revision `fcf8584b64974ef7a238eac718e01be4ed637a1d` and released model SHA-256 `d02aad4010fd1036bf07b13df0f9a24c0af37097d2802f5a5e4efc166f8af438`. The release FITS copy under `0_DATA` and the Pippin `SNDATA_ROOT/lcmerge` copy are byte identical.

All 1,635 classification CIDs appear uniquely in the 19,706-record DES HEAD. Inclusive one-based PHOT pointers matched `NOBS` for every selected CID. The original protocol retained 163,616 of 165,092 epochs after rejecting 1,476 rows with PHOTFLAG bits 8, 16, 32, 64, 128, 256, or 512. `REDSHIFT_FINAL[_ERR]` replaced `HOSTGAL_SPECZ[_ERR]`. Official onthefly code performed 0.33-day grouping, minimum-error same-filter selection, `[g,i,r,z]` ordering, one-hot features, normalization, packing, and class-0 Ia inference. A 10-object single-versus-batch check agreed within `1.28e-13`.

| Prospective measure | Original attempt | Gate | Windowed diagnostic |
|---|---:|---:|---:|
| RMS probability error | 0.1301 | ≤0.005 | 0.05777 |
| 99th percentile absolute error | 0.8061 | ≤0.02 | 0.32014 |
| Spearman rank agreement | 0.8671 | ≥0.995 | 0.95649 |
| >0.5 membership agreement | 96.76% | ≥99% | 98.84% |
| >0.999 membership agreement | 89.54% | ≥98% | 97.43% |

The primary median absolute error was 0.000152, while the windowed diagnostic median was 0.0000389; the remaining tail is not explained by four-decimal release rounding. Both columns fail every prospective gate measure.

## Documented prediction path and diagnosed omission

The [DES release analysis YAML](../../sources/repos/des-science__DES-SN5YR@1.3/7_PIPPIN_FILES/D5yr_analysis.yml) assigns `DATADES5YRSMP` to `DATADES`, `SNNTESTV19_z_data` to one `vanilla` V19 model with spectroscopic redshift and `cosmo` normalization, and `PROB_SNNV19` to the aggregation alias. Historical Pippin [DataPrep source](https://github.com/dessn/Pippin/blob/bb3ffdd44325f82a86d239c91b386ed0e4d18e44/pippin/dataprep.py) always supplies a `clump_file`, with defaults `OPT_SETPKMJD=16` and `PHOTFLAG_MSKREJ=1016`; its [SuperNNova classifier source](https://github.com/dessn/Pippin/blob/7efa8e9be584bc98d4e4ffab7c59c25f63236918/pippin/classifiers/supernnova.py) passes that file as `--photo_window_files`. SuperNNova's default `photo_window_var` is `PKMJDINI`; the data loader retains epochs satisfying strictly `-30 < MJD-PKMJDINI < 100`, **before** PHOTFLAG rejection. The saved *training* `cli_args.json` has no photo-window file; it does not describe the prediction data stage. The source data job also sets `--redshift_label REDSHIFT_FINAL` and the seven PHOTFLAG bits.

The independently declared [window diagnostic plan](../../phase2/classification/reconstruction_20260926/clump_correction_plan.md) rebuilt the source-defined clump file using local SNANA and applied the window before flags. For the 1,635 CIDs it excluded 61,478 epochs by time, then 752 by flag, retaining 102,862. Its [per-object probabilities](../../phase2/classification/reconstruction_20260926/des_clump_diagnostic.csv) and [metrics](../../phase2/classification/reconstruction_20260926/des_clump_diagnostic_metrics.json) are labelled diagnostic. A local historical SNANA rebuild gave peak times within 0.00005 day of the current-source build and changed zero of these window memberships; [binary, input, and output hashes](../../phase2/classification/reconstruction_20260926/clump_build_provenance.json) are preserved. Neither local binary establishes the exact original-run clump file.

Historical Pippin's [classifier export](https://github.com/dessn/Pippin/blob/7efa8e9be584bc98d4e4ffab7c59c25f63236918/pippin/classifiers/supernnova.py#L393-L422) writes `all_class0` for the vanilla model at four decimal places. Its [aggregation code](https://github.com/dessn/Pippin/blob/2a8f0f85762f34699538e1a7ebddd7ba707a28aa/pippin/aggregator.py#L92-L119) maps task names to probability columns and [merges/renames](https://github.com/dessn/Pippin/blob/2a8f0f85762f34699538e1a7ebddd7ba707a28aa/pippin/aggregator.py#L253-L283) outputs. For `DATADES`, the YAML has one V19 task. This source/config path provides no seed or model ensemble explanation. The `SNNTESTV19_z` second alias is a perfect classifier for BiasCor simulations, not the data probability. Dovekie is not part of this original release YAML path.

## Independent checks and remaining gap

The [first tensor check](../../phase2/classification/reconstruction_20260926/diagnostic_features.json) used five largest original mismatches and five closest agreements declared in the [diagnostic plan](../../phase2/classification/reconstruction_20260926/diagnostic_plan.md). A separate official `make_dataset` pivot and the wrapper produced identical selected tensors. The model's first LSTM layer expects 26 inputs; saved feature indices select flux/error for `g,i,r,z`, `delta_time`, `HOSTGAL_SPECZ[_ERR]`, then 15 filter indicators in saved order. Class index 0 is Ia. Original `HOSTGAL_SPECZ` substitution changes those ten probabilities by roughly 0.0001–0.001, not the large original tail. The [windowed tensor check](../../phase2/classification/reconstruction_20260926/diagnostic_window_features.json) used separately declared five largest corrected mismatches and five nearest agreements: official Astropy FITS loading and the wrapper gave identical retained rows and 26-feature tensors. The excluded photo-z fill differed, but this model does not use photo-z.

An isolated [official SuperNNova test-database build](../../phase2/classification/reconstruction_20260926/test_database_command.json) with the Pippin flags, final redshift, clump, filters, and type map was followed by pinned official [`run.py --validate_rnn`](../../phase2/classification/reconstruction_20260926/official_validate_command.json) using the released weight. Its 17,733-object prediction pickle covers every classification CID. Against windowed onthefly predictions for all 1,635, [comparison metrics](../../phase2/classification/reconstruction_20260926/official_validate_comparison.json) show maximum absolute difference `4.17e-7`, RMS `2.44e-8`, and zero differences above `1e-6`. The official validator's RMS difference from published values is still `0.05777`. This independently validates the local wrapper against the pinned official full dataset path, not against the unavailable original Pippin run.

The test database has its own normalization, but official validation with `--model_files` reloads the *model's saved* `data_norm.json`. A declared [counterfactual test-norm run](../../phase2/classification/reconstruction_20260926/des_test_norm_diagnostic_metrics.json) worsened RMS to `0.06183`; it is not the effective release setting. The archived model's `predictions.csv.gz` contains simulation IDs without their executable raw curves and cannot serve as a same-input DES fixture. The original prediction job's exact source/runtime, data-stage HDF5, clump file, and unrounded DATADES prediction pickle are absent from acquired release assets. The residual cause remains **unidentified**. Neither a model/runtime difference nor an unarchived original input/prediction-stage difference can be ruled out.

For the strict >0.999 cut, the published DES selection has 1,063 CIDs; the windowed diagnostic has 1,077. They share 1,049; 28 are diagnostic-only and 14 release-only, so 42 of 1,635 memberships differ (2.57%, Jaccard 96.15%). At most 42 individual selection contributions can differ in a further *count on these same DES CIDs*. This does not bound weighted downstream estimates, fitted parameters, or a simulated sample. The [selection-impact record](../../phase2/classification/reconstruction_20260926/selection_impact.json) preserves these counts. At >0.5, 19 DES memberships differ. The gate remains closed.

The most useful original-run assets to resolve the gap are, in order:

1. `DATADES5YRSMP` Pippin output `DES-SN5YR_DES.SNANA.TEXT` and generated `clump.nml`, with SNANA executable revision and command, to establish the exact window.
2. `SNNTESTV19_z_data` data-stage `job.slurm`/CLI and `processed/database.h5` (or per-CID post-window feature tensors), with the exact HEAD/PHOT files and hashes read by that job, to separate preprocessing from inference.
3. Its unrounded `PRED_*.pickle` and exported `predictions.csv`, including job log, to determine whether published four-decimal values came from this weight and data stage.
4. The classification-run `model.pt` checksum and model directory, Pippin/SuperNNova revisions, Python/Torch/Pandas versions, and aggregation/merge log, to identify remaining runtime or artifact substitutions.

The local [reconstruction wrapper](../../scripts/phase2/classification/reconstruct.py) runs the original protocol with `OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 phase2/classification/env/bin/python scripts/phase2/classification/reconstruct.py des` from repository root and enforces the preregistration hash. Its simulation mode enforces the failed DES gate and refuses to run. The [window diagnostic](../../scripts/phase2/classification/diagnose_window.py), [official validator](../../scripts/phase2/classification/run_official_validate.py), and [comparison script](../../scripts/phase2/classification/compare_official_validate.py) are separate reproducible diagnostics. No released labels were used to tune preprocessing, thresholds, or weights; simulation truth Ia is not treated as a high-SNN-probability selection.

## Post-validation residual membership sensitivity

The [declared protocol](../../runs/research_2026_09_26/classifier_residual_sensitivity/protocol.md) restricts the existing 1,020-object frozen observer-score validation to CIDs passing **both** released `PROB_SNNV19 > 0.999` and source-corrected, windowed `pIa > 0.999`. The official `--validate_rnn` output gives the same threshold membership for these CIDs. It admits no reconstruction-only object and changes no fitted coefficient, epoch mask, score template, or primary result. The [script](../../runs/research_2026_09_26/classifier_residual_sensitivity/sensitivity.py), [complete CID membership ledger](../../runs/research_2026_09_26/classifier_residual_sensitivity/membership_and_scores.csv), [result with ordered retained and removed IDs](../../runs/research_2026_09_26/classifier_residual_sensitivity/result.json), and [input/output hash manifest](../../runs/research_2026_09_26/classifier_residual_sensitivity/manifest.json) preserve this as a separate sensitivity.

| Frozen published-mask scores | CIDs | Epochs | Matched product | Information | Fixed log-score gain |
|---|---:|---:|---:|---:|---:|
| Unchanged primary | 1,020 | 39,606 | 192.401284 | 330.336139 | 27.233214 |
| Both classifiers >.999 | 1,006 | 39,148 | 194.542618 | 328.514809 | 30.285213 |
| Removed release-only members | 14 | 458 | −2.141334 | 1.821331 | −3.051999 |

The 14 removed CIDs, in the frozen validation order, are `1701055, 1309288, 1265165, 1287437, 1278602, 1343533, 1632070, 1277945, 1249371, 1370944, 1300470, 1264963, 1879313, 1333442`. The result file lists all 1,006 retained IDs exactly. All **43** discovery-fit CIDs pass both thresholds; their frozen coefficient fit is untouched.

| Field | Retained/primary | Matched product, retained | Information, retained | Fixed gain, retained | Primary fixed gain |
|---|---:|---:|---:|---:|---:|
| C1 | 108/108 | 20.403020 | 28.201794 | 6.302123 | 6.302123 |
| C2 | 107/110 | 16.929795 | 28.502675 | 2.678458 | 1.083840 |
| C3 | 132/134 | 44.476426 | 70.967916 | 8.992468 | 9.045471 |
| E1 | 90/91 | 10.532460 | 26.464135 | −2.699607 | −2.705941 |
| E2 | 91/94 | 12.154793 | 26.043706 | −0.867060 | −1.519366 |
| S1 | 53/55 | 9.979044 | 15.350032 | 2.304028 | 2.092251 |
| S2 | 95/96 | 23.602692 | 29.518177 | 8.843604 | 8.842363 |
| X1 | 92/92 | 3.276257 | 22.107436 | −7.777461 | −7.777461 |
| X2 | 94/95 | 16.632101 | 24.956708 | 4.153748 | 4.446244 |
| X3 | 144/145 | 36.556028 | 56.402231 | 8.354913 | 7.423690 |

The 14 removed objects contributed a net **negative** fixed gain, so the membership gap does not drive the positive frozen validation gain. This is a post-validation, conditional arithmetic check on an intersection cohort, not a new significance calibration, evidence that either classifier is ground truth, or a repair of the failed reproduction gate. Per-object matched products, information and gains were independently checked against the unchanged frozen coefficient and sufficient arrays before summing.
