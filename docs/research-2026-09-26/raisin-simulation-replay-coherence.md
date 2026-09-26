# Independent archived RAISIN simulation coherence audit

The archived NIR FITRES and LCPLOT products cannot be joined by their printed CID to replay the same fitted simulations. This is a measurement/provenance gate failure, not a measured timing bias or a failure of the SNooPy physical model. No native fits, alternative permutations, or paired timing science scores were run in this audit.

## Direct evidence

The source FITRES has 30,000 rows; LCPLOT has 500 CIDs (1–500), 1,829 accepted data rows and 69,200 model-curve rows. IFIT is always 1 and no rows carry the rejected-data flag −1.

| Same-CID condition | Compatible objects / 500 |
|---|---:|
| Accepted LCPLOT count = FITRES NDOF + 1 | 223 |
| Plotted model clock = FITRES PKMJD, within 0.0056 day | 28 |
| Plotted observing cadence = FITRES SIM_LIBID | 42 |
| Both clock and cadence | 27 |

Every plotted model curve has exactly constant printed `MJD−Tobs` within its CID. The median absolute difference from the same-CID FITRES peak is **52.094 days**, with a maximum **433.004 days**. This cannot be explained by printed precision. Uniform CID offsets from −3 through +3 produce only 9–100 clock matches; there is no simple one-row shift.

All 500 plotted data sequences uniquely match one of the 23 archived SIMLIB J/H cadences. This result is unchanged for time tolerances 0.003, 0.0056 and 0.01 day. The matching requires distinct, same-band source epochs; it allows an accepted subset and makes no assertion about omitted epochs or flags. For example, plotted CID1 has the H cadence of SIMLIB 7, while FITRES CID1 says SIMLIB 3. FITRES CID2 says 7. This identifies a reused cadence, not a unique simulated event.

A separately frozen metadata rule compared each plotted model clock and uniquely identified cadence with all 30,000 FITRES rows. It used no distance, flux, chi-square, or NDOF field. **Every one of the 500 plots has multiple candidates: 12–1,741, median 557.** There are no uniquely identifiable events and no plots without candidates. The original CID survives in only 27 candidate sets. Thus the available metadata neither repairs event identity nor proves that every plotted event is absent from the FITRES product.

The rule was frozen after the same-CID and cadence diagnostics, and after an exploratory clock/count/cadence candidate diagnostic. A discarded, explicitly descriptive per-row-chi-square compatibility check is preserved in `cadence.py` and its ledger; it was not used for the frozen metadata match, selecting a permutation, or an interpretation. In particular, individual plotted contributions are not assumed to sum a correlated native fit objective. This is an integrity diagnostic, not a preregistered physical hypothesis test.

## Source meaning and precision

In pinned SNANA v11_04k `src/snlc_fit.car`, SNLCPLOT takes the printed CID from `SNLC_CCID` (around lines 18180–18188). Data rows come from the native fitted epoch arrays (18241–18274), with rejection flags carried separately. The model clock uses `PKMJD=LCVAL_STORE(IPAR_PEAKMJD)` (18222) and `VMJD=TOBS+PKMJD+MJDOFF` (18363). `GET_FLUX_FITFUN` uses the stored fit parameter vector (6444–6523). The model curve therefore represents a fitted clock, not a freely inferred truth peak or a rest-frame time coordinate.

The text writer `src/sntools_output_text.c` maps data/model/rejected flags to 1/0/−1 (1438–1463) and prints MJD to three decimals and Tobs to two (1572–1581). The 0.0056-day clock tolerance is deliberately generous for this printed-precision audit, not a physical uncertainty. Model Tobs values are on a regular half-step plotting grid. Data-row clock estimates differ from their model clock by at most 0.007 day; those small differences are not treated as a new inconsistency because data times pass through different floating-point arrays and coarser formatting.

The native degree-of-freedom relation is accepted epochs minus free parameters (`snlc_fit.car`, 12811–12821). The archived NIR configuration fixes shape, AV, RV and peak and fits DLMAG, supporting the one-free-parameter count check. The timing and cadence contradictions independently establish the failed join without relying solely on that count assumption. LCPLOT does not preserve a full per-object fit parameter or truth header: the other same-CID shape/dust/redshift/truth assignments cannot be certified from it once event identity fails.

## Provenance and replay interpretation

Both compressed products reproduce their declared Git blob hashes at author commit `b888214a5cbae38ac0bf886488ce7734f5b77e87` in the same nominal directory:

- [Archived FITRES](https://github.com/djones1040/RAISIN_cosmo/blob/b888214a5cbae38ac0bf886488ce7734f5b77e87/output/fit_nir/DES_RAISIN_NIR_SIM/DES_RAISIN_SIM/FITOPT000.FITRES.gz), blob `7dd16e0978e2be4f09ebdc14698aaeca902df26e`.
- [Archived LCPLOT](https://github.com/djones1040/RAISIN_cosmo/blob/b888214a5cbae38ac0bf886488ce7734f5b77e87/output/fit_nir/DES_RAISIN_NIR_SIM/DES_RAISIN_SIM/FITOPT000.LCPLOT.gz), blob `d81d64d3f1f395ad9dbbab71837d89f4811320af`.

The saved SUBMIT.INFO records submission on 2021-11-11, four version jobs, one split per task, nominal FITOPT000, cleanup, and no linked FITOPT000 jobs. This is not a hash-linked execution record for these two files. The available master simulation input names a repeat/seed scheme (`10 123459`), but its execution and event-ID mapping to this product pair are not established. No per-file historical execution logs, unmerged matched plot/table pair, or matched raw HEAD/PHOT were found in the reviewed handoff. Whether different executions, reused identifiers, a merge/renumbering step, or another process caused the incompatibility remains unresolved; no source bug or particular mechanism is established.

Sol's author-grid replay correctly fails its baseline gate: 7/8 FITRES outputs survive, and the archived fit values are not recovered even though replayed epoch tokens generally close. This audit explains why reproducing the rounded plotted flux/error rows is insufficient: it does not validate the attached same-CID redshift, peak and fit metadata. The failed replay must remain preserved. Changing a model, loosening tolerances, or assigning compatible FITRES rows by their distance/objective would not repair this provenance gate.

## Minimum input that enables the timing experiment

Obtain the original simulation HEAD/PHOT pair (or a lossless equivalent) linked to **one** archived NIR and joint optical+NIR fit execution. Required fields are unique generated-event identity, full-precision epoch flux/error/time/band and flags, heliocentric redshift, peak metadata, truth fields and SIM_LIBID. Preserve the generation version/seed and any CID randomization/merge map, and provide the actual fit input, logs, binary/model/KCOR identities, and corresponding accepted-mask exports or a reproducible mask rule.

With those inputs, first reproduce the nominal fixed-peak NIR result and joint-fit result for the same events. Only then compare NIR fits using the original near-truth peak versus the independently recovered joint-fit peak on unchanged simulated observations. The present files do not support a paired timing-bias estimate. A new explicitly generated simulation could test a declared timing mechanism, but would be a separate experiment rather than replay of these archived events.

## Reproducibility

Owned evidence is `runs/research_2026_09_26/astra_design/raisin_simulation_replay_review/`:

- `audit.py`, `coherence-500.csv`, `coherence-result.json`: independent all-row count/clock audit.
- `cadence.py`, `cadence-500.csv`, `cadence-result.json`: source cadence audit and preserved exploratory diagnostics.
- `metadata-match-protocol.json`, `metadata_match.py`, `metadata-candidates-500.csv`, `metadata-match-result.json`: frozen metadata-only association rule and all candidates.
- `verify.py`, `verification.json`: separate Decimal parsing and exhaustive same-band subset matching; no imports from the main audit or collaborator executors; verifies compressed SHA-256 and Git blob identities.
- `manifest.json`: report, source and evidence hashes.

All original files and failed replay outputs remain unchanged. No model/data corrections or cosmological inference follow from this audit.
