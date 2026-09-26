# Research paused at the user's request — 2026-09-26

The user asked to finish current thinking and stop gracefully. No new experiment
is authorized by this checkpoint itself. Resume only when the user asks. The
research goal is **paused, not complete**. All working-tree changes, executed
sources, unsuccessful gates and original data remain preserved; nothing was
published or committed by this continuation.

## Scientific conclusion at the pause

The data do not justify either “all biases are absent” or “new cosmology has
been established.” There are reproducible processing defects, omitted signed
measurements, calibration ambiguities and a transferable optical residual
pattern. Several apparently large distance responses are not uniquely
identified once numerical instability and alternative intrinsic SEDs are
allowed. A defensible population correction and corrected cosmology remain
downstream of those unresolved issues.

The [E00–E17 status map](programme-status.md) locates all main-path research
directions; the [research index](README.md) links completed bounded experiments.
The [measurement-to-cosmology ledger](measurement-to-cosmology-ledger.md)
distinguishes already-applied corrections from conditional responses and
unresolved calibration provenance. These supersede neither the dated phase-1
report nor the separate physics audit.

High-information findings from the continuation:

- A fixed observer-band residual pattern transfers from 43 discovery SNe to
  1,020 held-out SNe. A broad intrinsic-SED alternative can mimic that residual
  pattern while implying a very different distance response. Residual fit alone
  cannot determine the cosmological correction; arbitrary gray luminosity
  evolution remains exactly degenerate with distance.
- RAISIN's optical product omits negative signed measurements present in its
  author DIFFIMG ancestors. Its simulated cadence inherits that membership.
  The ten-object signed-refit study exposes weak template branches rather than
  a unique omitted-epoch correction. Historical and modern native fits also
  show starting-point sensitivity.
- A source-linked host-mass systematic fits a threshold at 10.44 dex but applies
  it at 10. Its bounded repair restores a missing uncertainty direction; it
  does not change nominal distances or reproduce the full covariance pipeline.
- A WIRC-J/RC1-versus-RC2 filter assignment mismatch has a stable conditional
  all-42 distance response: mean **+0.000483432 mag**. That small mean does not
  bound physical passband, training or population uncertainty.
- Historical HST count-rate nonlinearity inclusion is still unlinked because
  the upstream aperture-magnitude catalog/producer is absent. Ordinary CALWF3
  nonlinearity and the released wavelength-dependent systematic are different
  operations. No generic 0.04-mag correction was imposed.

## Most recent completed work

### Timing simulation

The eight-draw prospective engineering experiment passes generation, physical
support, genuine displaced starts, nine NIR jobs, exact null identity and
independent objective/covariance reviews. Estimated versus true timing changes
the NIR distance by **+0.000353654 mag on average**, with **MC SE 0.001017481 mag**.
This is an unresolved mean response in a conditional single-cadence engineering
sample, not a survey bias bound. See
[the complete record](raisin-prospective-native-engineering.md).

The separate fixed 64-draw experiment retains 7,488 measurements, including
390 negatives, and passes generation. Its initial joint-fit checker stops on
native final-iteration retries for CIDs 22 and 64. Source diagnosis finds real
MINUIT nonconvergence returns which are later overwritten by the ordinary
EXIT return; exported fit flags do not reveal them.

An additive occurrence-aware checker preserves every retry and the original
scientific thresholds. All three remaining joint diagnostics (9 iterations,
actual ±2-day starts) complete and satisfy the cross-run D/peak checks. Internal
warnings remain: nominal/9/minus/plus counts **20/14/10/21**. Native activity is
**69.294204825 s**. These agreements do not establish peak stationarity; **no
64-draw NIR stage or timing-effect summary has been run**.

Exact restart location: `phase2/pt64/occurrence-diagnostics/` (a short symlink into
the immutable prospective run tree), especially `diagnostic-result.json`,
`stationarity-plan.md` and the preserved `joint-failure.json` one directory up.
An output-only MINUIT status/trial-mean patch is prepared under
`phase2/pt64/minuit-trace-design/`; it has not been built or executed. Its next
gate is absent/disabled/enabled numerical identity before using captured means
to independently check both sides of the fitted peak.

### HST calibration and detector controls

Private, pinned CALWF3 3.7.3 exactly reproduces all SCI/ERR/SAMP/TIME pixels for
two current archive science exposures. Archive DQ adds AstroDrizzle bit 4096;
the whole-product equality gate remains failed, with extra archive distortion
products also recorded. This links audited source equations to the current
images, not to historical RAISIN photometry or validated physical noise. See
[native replay](raisin-calwf3-native-replay.md).

The fixed fourteen-dark RAW experiment is complete: seven disjoint primary
pairs plus one overlapping NORMAL sensitivity pair, all six fixed time-weight
powers, reference-only masks, 256 blocks and 247 eligible apertures. Independent
physical-uint16 re-extraction verifies all matrices and contractions. RAW total
moments include unflagged events and detector-state variation; they cannot
isolate the suspected factor-two read-variance term. Eleven inputs retain
indeterminate timing/telemetry flags. See the
[RAW report](../../runs/research_2026_09_26/raisin_hst_pixel_feasibility/dark_ramp_execution/dark-results.md)
and `dark_ramp_execution/root-raw-review/result.json`.

A new first-eight-read native dark control retains original pixel bytes and
time/quality metadata. Both existing eight-read science controls reproduce
their RAW and native FLT files byte-for-byte. The NORMAL pair `idbx41onq` and
`idbx43p7q` was then processed with its original dark switches and nominal
times, in **2.433 s**. All fixed-mask output pixels are finite with positive
quoted pair variance. No output-flag or residual selection was added.

The descriptive ratio of repeat squared difference to summed quoted variance
is **4.201565** across all fixed-mask pixels. Quadrants B/C/A/D give
**0.669390 / 0.616820 / 14.723955 / 0.799958**. The same **247 predefined
aperture-annulus contrasts give 0.547440**. These are distinct spatial
estimands. Independent influence review identifies one pixel, raw (506,945),
contributing **82.8641%** of the total squared difference. It carries native
DQ=32 in the later exposure and has zero weight in all 247 predefined apertures.
Four distinct pixels carrying DQ=32 in either exposure contribute
**8,305.450 of 9,866.900** total squared difference. All primary results remain
unchanged; these pixels were not dropped to obtain a preferred noise scale.
Independent extraction reproduces all five regions, 256 blocks and 247
apertures; independent circle quadrature agrees within 3.94e-15 DN/s.
This single conditional dark pair does not justify rescaling historical
supernova errors.

Exact run: `runs/research_2026_09_26/raisin_hst_pixel_feasibility/calwf3_dark_prefix/`.
Key files are `score-result.json`, `score-state-ledger.json`, the frozen scorer,
`native-DQ-bit-counts.csv`, and both exact null reviews. Dark outputs are
**DN/nominal-second**, whereas flattened science products are electrons/second.
The native formula uses mean gain and adds a fixed dark term even with DARKCORR
omitted. Influence/category/read-history ledgers are in
`calwf3_dark_prefix/independent-score-review/`. Internal selected powers,
fit intervals and variance-term accounting are still needed; no full
fourteen-prefix extension has been launched.

### Conditional dust prediction

The [new signed optical BayeSN design](bayesn-signed-optical-pilot.md) can proceed
on ten complete DES16 ancestors without pretending to reproduce survey
selection. It fits optical measurements using optical-trigger timing and
compares held-out NIR predictions under LCDM versus proper broad distance
priors. Both arms retain the trained intrinsic-SED assumptions. No observed
posterior or held-out NIR score has been computed.

The two redshift-endpoint objects' adapter is built: 15/14 signed optical rows,
6/5 separately sealed NIR rows, corrected released redshift/MW metadata and
the exact released DES J passband. Optical forward, dynamic timing,
derivative, gray-scaling and toy proper-distance-integral checks pass.
The optical 1200→2400 wavelength refinement is **0.001074 quoted sigma**.
These checks do not yet certify J/H forward predictions or actual BayeSN
posterior marginalization. The expanded synthetic-only runtime benchmark
finished before stopping: 80 warmup plus 20 discarded draws per arm took
11.89 s (LCDM) and 7.69 s (broad), with median 63 and maximum 511 NUTS steps.
Its crude single-chain 1,000+1,000 estimates are 238/154 s, not convergence
forecasts. The actual normalized broad-D logpdf agrees with independent
convolution quadrature within 6.75e-14; the normalization integral is one.
No full recovery chains or observed inference were launched. See `checkpoint.json`.

Run: `runs/research_2026_09_26/bayesn_signed_optical_pilot/`; report:
[forward engineering](bayesn-signed-optical-engineering.md). Next: metadata-only
J/H forward gates and an independent static-phase comparison, then the fixed
synthetic recovery suite. Do not open sealed
observed NIR outcomes before optical chains and diagnostics are frozen.

## Resume order

1. Read agent pause checkpoints and verify no old processes remain.
2. Finish the bounded detector influence/variance accounting and two-dimensional
   timing convergence diagnostics; retain their failed primary gates.
3. Prioritize the BayeSN conditional optical-to-NIR prediction because it can
   distinguish observable consequences of distance-prior-driven dust inference.
   Pass synthetic/numerical and chain diagnostics first.
4. Only propagate physically supported, identified corrections through training,
   selection/BBC and shared covariance before fitting corrected cosmology.

No evidence at this pause supports claiming either a bias-free supernova
population or a discovered failure of standard cosmology.
