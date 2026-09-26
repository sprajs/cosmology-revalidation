# RAISIN timing assets: source-selected historical recovery

The missing raw simulation files were not located in the bounded public search. However, per-file repository history recovers a coherent **November 2021 NIR FITRES + joint FITRES + rounded NIR LCPLOT** set. The correct-source eight-object baseline replay now passes its unchanged engineering gates. This does not establish lossless reproduction or a timing-bias result.

## What changed in the archive

The latest author repository stores nominal NIR FITRES last changed in commit `83927ccbb6634172b1120463a483d911dcce9170` (2022-05-02, “new simulations”) beside LCPLOT last changed in `aaa709ead7a7339d56a2a4604d78991329d9368f` (2021-11-11, “fits”). The latest snapshot's same-CID inconsistencies remain valid, as documented in the [independent coherence audit](raisin-simulation-replay-coherence.md). The file history provides a source-based way to choose a historical snapshot without selecting rows by fitted distance or chi-square.

At the LCPLOT's last-modifying commit, the recovered NIR table has 30,000 rows and the joint optical+NIR table 29,995. All 500 archived plots now match their same-CID NIR accepted-row count, peak clock and source SIM_LIBID cadence. Maximum printed clock discrepancy is 0.0005 day. Across all 29,995 common NIR/joint CIDs, 15 `SIM_*` fields plus PKMJDINI, zHEL, zCMB, zHD and MWEBV match exactly. This is metadata coherence, not a full original epoch/noise identity certificate.

| Acquired asset at `aaa709e` | Compressed bytes | Git blob | SHA-256 |
|---|---:|---|---|
| [Nominal NIR FITRES](https://raw.githubusercontent.com/djones1040/RAISIN_cosmo/aaa709ead7a7339d56a2a4604d78991329d9368f/output/fit_nir/DES_RAISIN_NIR_SIM/DES_RAISIN_SIM/FITOPT000.FITRES.gz) | 2,116,981 | `43fa274f059b05b65b7f34af34c4eae2030cffa8` | `c0e6de446c2d4b2179a6570c0a1f166b5384612101fbe5002938a73b8daa3b44` |
| [Joint FITRES](https://raw.githubusercontent.com/djones1040/RAISIN_cosmo/aaa709ead7a7339d56a2a4604d78991329d9368f/output/fit_all/DES_RAISIN_OPTNIR_SIM/DES_RAISIN_SIM/FITOPT000.FITRES.gz) | 3,541,482 | `343479b52f839ae77dbb3254aa21eeee417bc6ca` | `2a1fb2553ce1d3f34c9617bcb83228044d708ae991964a8e59f5d47eb6fbb4e8` |

The historical LCPLOT blob is the already acquired `d81d64d3f1f395ad9dbbab71837d89f4811320af`; it was not downloaded again. The NIR/joint NML files, B18 grid and info, DES KCOR, generation input and SIMLIB have identical Git blobs in the 2021 and latest snapshots. The local model/KCOR/info/NML bytes were also verified directly against the historical Git blob identities.

Root's independently written `scripts/research_2026_09_26/verify_raisin_historical_pair.py` reproduces the compressed blob checks, all 20 shared metadata fields, and all 500 clock/count/cadence checks. Evidence: `runs/research_2026_09_26/raisin_historical_pair_root_review/result.json`.

## Public search and remaining absence

The live author repository has one public branch, `master`, at the same pinned latest commit, no tags and no GitHub release assets. The release repository has one branch and tags v1.0/v1.1, also unchanged from the archived snapshots. Complete, nontruncated trees were inspected:

| Tree | Entries | File blobs | HEAD/PHOT pair located |
|---|---:|---:|---|
| Author `aaa709e` (2021-11-11) | 29,793 | 29,515 | No |
| Author `b888214` (latest) | 29,884 | 29,596 | No |
| Data release v1.0 `318fbc4` | 513 | 488 | No |
| Data release v1.1 `a383c4b` | 516 | 490 | No |

Simulation-related binary candidates are calibration files, not generated photometry. The numerous `biascor_stretch` files are FITRES/LCLIST/LCPLOT outputs, not a lossless original epoch archive. The only file under `SIMLOGS_DES_SPEC` is a master input; per-job generation logs and generated HEAD/PHOT are absent in the inspected trees. The author `.gitignore` excludes simulation logs, compressed files and split fit job directories; this can explain missing working products but does not prove that any particular file existed or was intentionally withheld. The first 100 general commits and complete per-file histories for the two nominal products were retrieved. This is not an exhaustive proof that no other public historical URL exists.

The paper links its online release to [Zenodo 6349657](https://doi.org/10.5281/zenodo.6349657). Its file metadata identifies the v1.0 observed-data GitHub release ZIP (36,141,358 bytes; MD5 `4610a7cbec8d67139167d3a2b42d60b5`). The latest linked version, [6581744](https://doi.org/10.5281/zenodo.6581744), is v1.1 (36,067,671 bytes; MD5 `6cf26573cdc89dfbcd745c44e701b6c7`). Neither record lists a separate simulation asset. Both archives exceed the 30 MB acquisition cap and duplicate already inspected release trees; metadata only was acquired. The [release README](https://github.com/djones1040/RAISIN_DataRelease) describes observed photometry, distances, covariance and fitting assets. The [paper's simulation appendix](https://arxiv.org/html/2201.07801v2#A3) describes models and observing-condition libraries; their availability is not a promise that the original generated epoch files are archived.

The precise missing inputs remain full-precision simulated HEAD/PHOT (or equivalent) with original per-epoch flags and timing, the original heliocentric truth redshift, full accepted/rejected cadence, and generation/merge/seed records linked to this fit execution. The recovered printed tables cannot restore information discarded by LCPLOT formatting. In particular, the joint FITRES supplies a coherent alternative fitted peak, but its original optical epoch data are not recovered here.

## Exact historical source and private build

The recovered 2021 tables explicitly name **SNANA v11_04d**, whereas the later FITRES names v11_04k. Official tag `v11_04d` resolves through annotated tag object `f8c25e0de3c68f818406974467c9736d44d65784` to commit [`10ec91297e4482d593cb5d3d055b10d4aa915071`](https://github.com/RickKessler/SNANA/tree/10ec91297e4482d593cb5d3d055b10d4aa915071).

The source archive is 4,439,172 bytes, SHA-256 `f7758dd94e228f17378f081bd32639064e49d903b061859749552b6353202a35`. Every one of its 201 source files matches the official recursive Git tree blob. Untouched source and a separate build copy are retained under `astra_design/raisin_timing_assets/snana_v11_04d/`.

The private fitter compiled in 19.47 seconds using the existing local GCC/GFortran16, CFITSIO and project GSL runtime. Its only source-file change is upstream-generated disabling of unavailable ROOT/HBOOK output flags. There are no scientific source changes, shared binary swaps or global package changes. Binary SHA-256 is `d42c96c18984434e4e9ddc8bfd80f7cc3a18c4b64ede0467f98116719dcc5b51`. Dependencies resolve, but the original 2021 compiler/binary environment has not been reproduced.

## Frozen baseline-only amendment

The original eight numeric CIDs are retained. The coherent 2021 metadata and original rounded accepted J/H rows produce 16 prepared input files: baseline and byte-identical copy for each object. All 32 epoch flux/error/time tokens are unchanged. Fixed peak uses the archived PKMJDINI token and the same native float32 conversion; conversion values are recorded separately. The historical MW default remains O'Donnell94, and `./snoopy.B18` explicitly names the verified grid rather than allowing a fallback model.

Original engineering criteria are unchanged: all objects present and finite; fixed-parameter and float32 peak closure; archived NDOF and accepted `(band, printed MJD)` multiset; outgoing flux/error print bounds; byte-identical copy outputs; distance difference ≤0.001 mag; data-chi-square difference ≤max(0.01, 0.001×archived data chi-square). These are engineering caps, not a derived precision theorem. The rounding screen and any timing intervention are separate later gates; passing the baseline alone cannot certify them.

`baseline_2021/protocol.json` SHA-256: `399dc9562d8c19af2e18b153cd887f1a576ec2ace9d04c96aa266c56fcc57185`. The resource-only v2 preserves v1 files and limits each native process to 60 seconds, the pair to 120 seconds, and the original pilot to 600 seconds while carrying prior execution time. The prepared SNANA directory is 119 characters, below its 160-character limit; output prefixes remain short and relative. The released baseline and byte-copy jobs completed in 1.58 seconds combined. All eight pass: maximum absolute distance difference from the coherent archived table is 0.000060 mag and maximum data-chi-square difference 0.003750. All 32 accepted epoch flux/error tokens, mask multiplicities and fixed metadata close; the two runs have identical science FITRES rows and LCPLOT bytes. Root independently reproduced these checks using a Decimal parser without importing the runner or gate: `runs/research_2026_09_26/raisin_historical_pair_root_review/native-baseline/result.json`. No paired timing result follows from this engineering pass.

If the rounded replay still fails, `prospective-fallback-design.json` specifies a separate small source-based simulation: eight engineering and 64 new draws at metadata-selected SIMLIB11, fixed nominal stretch/dust, identical optical+NIR latent/noise draws across branches, and true-peak versus estimated-peak NIR fitting. It isolates a conditional timing mechanism and preserves selection/noise limitations; it is not a replacement for historical reproduction or a population/cosmology correction. No generation is executed by this design.


## Truth-header limitation and output-only state gate

The rounded-input reconstruction did not achieve a complete truth-header round trip. Its `SIM_DLMAG` header uses the FITRES name, whereas this historical text reader expects `SIM_DLMU`; the source also has a `SIM_REDSHIFT_CMD` reader key while its writer emits `SIM_REDSHIFT_CMB`. Native output therefore contains −9 for those truth quantities. These headers were not changed to obtain the baseline pass. Fitting uses the independently checked measured redshift and fixed model coordinates; any future simulated-truth comparison must use the preserved archived truth table and an explicit separate mapping. This is not evidence for a physical data or paper defect.

A separate private source copy ports the already reviewed full-state output patch and diagnostic distance-start hook to v11_04d with zero fuzz. The untouched source and original executable remain unchanged. The instrumented build took 19.64 seconds; its binary hash is `58a3d41917361aa738bb66a71f5028aebb11c34e9bb3af82fb7992c522688c8f`. An absent or explicitly zero hook reproduces baseline science FITRES and LCPLOT bytes exactly. Actual post-initialization starts ±0.2 mag change the final unprinted distance by at most 4.77×10⁻⁹ mag. All 96 iteration states reconstruct their native total objective from data quadratic, exported prior and sigma terms within 1.78×10⁻¹⁵. Maximum final-iteration distance change is 2.07×10⁻⁵ mag; maximum fixed-final-covariance amplitude optimum shift is 6.46×10⁻¹⁰ mag. These checks pass for the fixed-shape/dust first-eight problem; they do not validate a free multi-parameter population likelihood. Native runtime was 3.11 seconds.

Evidence is `instrumentation_2021/{state-protocol.json,state-freeze.json,state-result.json,state-arrays.npz}`. The native objective and covariance belong to their respective iteration states; chi-square values across changing covariance states are not interchangeable likelihoods.

## Source-defined rounding screen, frozen but not yet executed

Historical source stores each full-precision epoch in `R8EP_MJD`, but the LCPLOT path passes through `R4EP_ALL/FIT(JEP_MJD)` before promoting it to REAL*8 and printing three decimals. Compiled `MJDOFF=0.0`: the float32 value is the **full MJD**, not a precision-preserving small offset. This explicitly amends the earlier decimal-half-unit description without overwriting it.

For all 32 epochs, the printed token identifies exactly one float32 state. Its inverse full-precision time cell has width 0.00390625 day. Flux/error and redshift/MW header uncertainty instead use the compatible float32 states themselves, since the likelihood consumes those same stored states. Each of eight fixed peak tokens identifies a single float32 state and thus has no nonzero native peak-rounding perturbation. All 32 source-SIMLIB timestamps are uniquely compatible with the derived cells; this is a labelled candidate-input sensitivity, not proof of lossless epoch or flag recovery.

`rounding_2021/protocol.json` freezes 240 endpoint probes over 120 active input coordinates, eight renamed-CID zero controls, 32 deterministic gradient-directed corner checks and eight source-compatible timestamp controls. Original 0.001-mag and data-chi-square engineering caps remain unchanged. The first-order envelope and checked nonlinear remainders are a sensitivity screen, not an exhaustive or probabilistic rounding bound. Execution requires the separate frozen runner/input manifest and native-worker release; no timing intervention is included.

## Evidence

All new files are under `runs/research_2026_09_26/astra_design/raisin_timing_assets/`: public metadata and history responses, acquisition manifests, exact historical FITRES pair, all-500 metadata ledger, source/build verification, frozen baseline scripts/inputs, and prospective fallback design. Source history retrieval and historical table selection occurred before any comparison of numerical historical distance or chi-square outcomes. The earlier incompatible-pair audit and all failed native attempts remain untouched.
