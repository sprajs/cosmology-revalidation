# Are omitted negative epochs represented in RAISIN bias simulations?

The archived DES RAISIN simulation cadence reproduces the released data's
positive-epoch membership. It does **not** include the negative DIFFIMG epochs
omitted from those products. The available simulation/fitting configuration
does not explicitly impose a positive-flux cut on newly simulated measurements.
This is a concrete reproduction question for the bias correction; it is not
yet a quantified failure of the published distance correction.

This follow-up uses the sign audit's frozen original-DIFFIMG ledger and the
author repository at commit `b888214a5cbae38ac0bf886488ce7734f5b77e87`.
All acquired source files match their expected Git blobs. The
[cadence protocol](../../runs/research_2026_09_26/raisin_sign_source/cadence-protocol.json)
was saved before matching the simulation epochs, SHA256
`17373c8cf36a010e7956e80e2796dde1b2e2b920edc08c5a70b4c1100ee57b59`.
The data sign finding was already known. No observed flux residual, population
fit or cosmology outcome was optimized here.

For each original row, the rule matches the same physical object's `FIELD`,
optical band and MJD within .00055 day in `DES_RAISIN.simlib`. Multiple matches
are retained explicitly. This matches the input's timestamp precision; it
does not assume a calendar match establishes flux-calibration identity.

| Original DIFFIMG rows for 17 objects | Original positive | Positive in RAISIN / SIMLIB | Original negative | Negative in RAISIN / SIMLIB |
|---|---:|---:|---:|---:|
| All epochs | 5,614 | 5,108 / 5,108 | 4,229 | 0 / 0 |
| Source optical phase −7 to +45 rest days | 811 | 811 / 811 | 56 | 0 / 0 |
| M20 design phase −10 to +40 rest days | 836 | 836 / 836 | 47 | 0 / 0 |
| More than 180 observer days from peak | 3,951 | 3,445 / 3,445 | 3,760 | 0 / 0 |

SIMLIB presence equals RAISIN retention for every one of the 9,843 original
rows. Two original positive rows have multiple time matches in the simulation
table, including one inside both stated fitting windows; no arbitrary match
was selected. These are row counts, not independent exposure or object counts.
The full original row identity, flags and mapping are preserved in
[cadence-matches.csv](../../runs/research_2026_09_26/raisin_sign_source/cadence-matches.csv).

The archived `genSimlib.py` builds observing libraries from an object's
existing MJD/filter list. It imports sky, PSF and zero-point information from
survey simulation libraries. Thus inherited omission from that list can be
carried into a SIMLIB. The exact private input path used by that script is not
execution-linked to the public product, so the source code alone is not proof
of which command caused the omission.

The available nominal generator input uses photostatistical noise smearing
(`SMEARFLAG_FLUX: 1`), a per-object SIMLIB peak/redshift, and source-population
and event-detection cuts. Its master file adds `OIR.J19` scatter. Those event
cuts are different from conditioning every retained epoch on positive measured
flux. The translated and legacy optical simulation fit files leave
`EPCUT_SNRMIN` empty. The released data fit also leaves that setting empty.
Additionally, their phase ranges differ: simulation optical −15 to +45,
released data optical −7 to +45 rest days. These source settings are recorded
as available configurations, not certified historical execution settings.

A simulation that repeats only the surviving observing times but draws new
Gaussian noise does not reproduce selection on the original noise sign. A
correct comparison must either reproduce the actual sign-dependent reduction
for every simulated realization, or restore a justified set of signed data
and model its measurement likelihood. Using a truncated/folded likelihood
without establishing the data transformation would also be unjustified.

The next gate is therefore a source-matched nominal-versus-signed refit,
followed by a coupled simulation that applies the same epoch retention to
every generated realization and recomputes its correction. Sample-level
selection, optical timing passed into NIR, and the covariance-generation gap
remain separate dependencies before cosmology can be updated.

Inputs and exact acquisition records are under
[raisin_sign_source](../../runs/research_2026_09_26/raisin_sign_source/), including
`simulation-acquisition.json` and `simulation-cadence-acquisition.json`.
[Executable comparison](../../scripts/research_2026_09_26/raisin_simulation_cadence.py)
and [results](../../runs/research_2026_09_26/raisin_sign_source/cadence-result.json)
preserve the immutable ledger snapshot and all output hashes. None of the
released photometry, source code or simulation inputs was altered.

An [independent source reparse](raisin-cadence-error-bridge-review.md) now
verifies all 9,843 rows directly from the original FITS, all 23 SIMLIB blocks
and the released products, with zero membership or row-ledger disagreements.
It also identifies an exact numerical error-conversion bridge on 3,133 unique
DES16 author-to-product positive matches:
`sigma_product = round(sigma_author * 1.086 * 0.4*ln(10), 3)`.
This is a magnitude-error round-trip inference; the executed optical conversion
script is still not identified. The earlier 3,132 count additionally required
a unique match to the later DIFFIMG release. The extra product-only row is
DES16S2afz/z at MJD 58030.383, whose later DIFFIMG counterpart has changed
timing and photometry and is not substituted.

The [signed pre-explosion check](raisin-signed-baseline.md) separately tests
the actual author's baseline measurements. Its pooled within-band scatter is
close to the quoted error scale, while group and lag structure remain. This
does not justify assuming a complete independent-Gaussian near-peak likelihood.
