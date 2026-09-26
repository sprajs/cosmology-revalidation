# Registered forward-comparison robustness continuation (26 September 2026)

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

## Protocol frozen before new sensitivity outcomes

This continuation uses the existing [forward preregistration](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/forward_discrimination/preregistration.json), [handoff](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/phase2/forward_discrimination/HANDOFF.md), and unmodified comparison implementation. The original primary result and outputs remain in `phase2/forward_discrimination/primary/`; new results go only to `runs/research_2026_09_26/forward_robustness/`. The estimand is a retrospective truth-Ia-simulation conditional check of measured, selected DES fit observables. The classifier reconstruction has **not** produced matched simulated SNNV19 probabilities, so no classifier-matched simulation gate is run or assumed. Historical BBC/systematics membership is also unavailable. No result can establish a physical population law, dust mechanism, or cosmology.

Run all seven previously registered feasible sensitivities with the same nine complete forward simulation arms, 300 bootstrap replicates, seed and generated-attempt Poisson pairing used by the original script. The baseline remains the frozen 817 training and 202 held-out common-support DES objects. Each sensitivity recomputes its declared common support before reading fitted-observable outcomes; its cohort may differ from the baseline. Keep original fold-zero CIDs and compare sensitivity contrasts against the original primary only on the **intersection of identical CIDs**, reporting intersection and lost/gained test counts. The runs are:

| Name | Only changed setting |
| --- | --- |
| `all-des` | Include all original DES Hubble-diagram members irrespective of published classifier probability, then same basic cuts/support. This is a contamination/classifier-conditioning stress test, not a matched selection model. |
| `recovered-mask` | Use independently recovered light-curve fit arm; same cuts/support. |
| `strict-cuts` | Add abs(c)<0.2, x1ERR<0.7, PKMJDERR<1, FITPROB>0.01 to real and simulated fits. |
| `bandwidth075` | Redshift kernel width 0.075, truncated at three widths. |
| `bandwidth150` | Redshift kernel width 0.15, truncated at three widths. |
| `reference-qplus05` | Fixed q=+0.5 real magnitude reference; refit seven magnitude offsets on training objects. |
| `reference-qminus1` | Fixed q=-1 real magnitude reference; refit seven magnitude offsets on training objects. |

Within each run, report the primary colour/stretch energy-score and secondary magnitude/colour/stretch energy-score contrasts versus P21 for BS21 and G10, with the registered 300-replicate intervals that resample real train/test objects and common generated attempts. Lower scores are better. For primary-to-sensitivity comparison, subtract the two arm-vs-P21 contrasts per CID on the exact common test intersection, then give the mean and a 10,000-replicate paired-object bootstrap interval (seed 20260926); this interval isolates object variation and does **not** include simulation-sample uncertainty or refitted nuisance variation across the two runs. Report support attrition and field counts. Do not infer robustness from a sensitivity whose common support is too sparse or whose bootstrap fails; preserve the failure and its cause. Keep the generated-attempt index as the simulated denominator key—CIDs can repeat among early generated failures.

A field-provenance discrepancy was discovered after primary outcomes: real fields use the frozen first-PHOT-epoch label, simulated fields use FITRES. It is disclosed separately and is not retroactively called a registered sensitivity. A separate consistency test may later exclude multi-field real objects or derive first-PHOT-epoch simulation fields, but it is not folded into these seven prespecified runs. The later physical audit found that extrapolated low-`R_V` forward truth can imply negative passive-screen extinction at 8000 Å (4,601 historical-law rows among 71,946 written simulation rows). This is a **model-support failure in simulation truth**, not an observed dust or cosmology correction. No extinction curve or low-`R_V` population is silently changed here.

## Completed runs and support

All seven sensitivity runs completed 300 joint bootstrap replicates with the original comparison code, preserving the original primary outputs. New per-run manifests hash their inputs and outputs; `scripts/research_2026_09_26/forward_robustness.py --compare` reverified those hashes and the original primary output hashes before opening scores. The saved [common-CID comparison](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/forward_robustness/comparison/common_cid_comparisons.csv) has 28 BS21/G10 × primary/secondary score checks; [support.json](https://github.com/sprajs/cosmology-revalidation/blob/17487bf659fcbdeeea072221492bac14b04a0a85/runs/research_2026_09_26/forward_robustness/comparison/support.json) gives exact lost/gained CIDs.

| Sensitivity | Train | Test | Test CIDs shared with primary |
| --- | ---: | ---: | ---: |
| Primary baseline | 817 | 202 | 202 |
| All DES regardless of published SNN probability | 1,245 | 314 | 202 |
| Recovered fit mask | 817 | 202 | 202 |
| Strict fitted cuts | 525 | 147 | 147 |
| Redshift bandwidth 0.075 | 759 | 179 | 179 |
| Redshift bandwidth 0.15 | 848 | 213 | 202 |
| Reference q=+0.5 | 817 | 202 | 202 |
| Reference q=-1 | 817 | 202 | 202 |

The table below gives G10 minus P21 **mean** energy-score changes over each sensitivity's full held-out common-support cohort, with its registered 300-replicate 95% interval. Lower is better. These intervals include finite simulation-attempt resampling and refitting the train-only magnitude offsets; they remain conditional on each run's fixed support and do not adjust for multiple diagnostics. The BS21 intervals include zero in all seven sensitivities for both scores; its complete values are in the comparison CSV.

| Sensitivity | Colour/stretch energy Δ [95%] | Magnitude/colour/stretch energy Δ [95%] |
| --- | ---: | ---: |
| Primary baseline | -0.0253 [-0.0765,+0.0110] | -0.0377 [-0.0978,+0.0004] |
| All DES | -0.0304 [-0.0754,+0.0006] | -0.0293 [-0.0808,+0.0010] |
| Recovered mask | -0.0257 [-0.0771,+0.0102] | -0.0382 [-0.0970,-0.0013] |
| Strict cuts | -0.0049 [-0.0550,+0.0280] | -0.0142 [-0.0685,+0.0232] |
| Bandwidth 0.075 | -0.0275 [-0.0811,+0.0052] | -0.0426 [-0.1070,-0.0044] |
| Bandwidth 0.15 | -0.0180 [-0.0600,+0.0063] | -0.0278 [-0.0765,-0.0023] |
| Reference q=+0.5 | -0.0253 [-0.0765,+0.0110] | -0.0377 [-0.0978,+0.0004] |
| Reference q=-1 | -0.0253 [-0.0765,+0.0110] | -0.0377 [-0.0978,+0.0005] |

On exactly shared test CIDs, strict cuts make the G10-minus-P21 colour/stretch contrast **+0.0148** less favorable to G10 than the baseline on those same 147 CIDs; the paired-object bootstrap interval is [+0.0009,+0.0285]. Narrowing the bandwidth changes that contrast by -0.0059 [-0.0118,-0.0004] on 179 common CIDs, while widening changes it by +0.0083 [+0.0005,+0.0164] on 202 common CIDs. These narrower paired-object intervals omit finite-simulation and nuisance-refit uncertainty across runs. They show sensitivity to matching/cuts, not a stable physical ranking. The all-DES colour/stretch value is exactly unchanged on the original 202 CIDs because the same fitted simulated forecasts and conditioning variables are used there; the full-cohort score changes when the additional 112 test objects enter. Reference-q changes are absorbed by train-fitted magnitude offsets to numerical precision in this design, and cannot favor an expansion law.

Several unadjusted secondary-score intervals now fall narrowly below zero, while the strict-cut and all-DES intervals cross zero. The changed support, historical DES tuning of the simulation laws, field-provenance mismatch, unmatched real/SNN classifier membership, and low-`R_V` physical-support failure remain barriers to claiming a preferred population correction. The generated-attempt denominator remains 26,518 per arm and is indexed by `generated_attempt_index`; reused CIDs are never counted as independent attempts.
