# LCDM baseline reference audit

This is a bounded numerical experiment using exact DESI DR2 Gaussian BAO
compression and the installed Irreducible native SDK. It tests a chosen,
conditional **massless-radiation** flat matter/Lambda model. It fits no
parameters and does not reproduce the full Planck 2018 six-parameter model or
its posterior. The Prospector handoff is a design/source audit; the runnable
variant and the full-model blockers must remain distinct.

## Model and measurement map

| Quantity | Declared value and role |
| --- | --- |
| H0 | 67.4 km/s/Mpc; chosen parameter control |
| Omega_m | 0.315; all matter pressureless |
| Omega_b | 0.049; subset of matter |
| Omega_r | 0.000092; all radiation massless |
| Omega_gamma | 0.0000545; subset of radiation |
| Omega_Lambda | 1-Omega_m-Omega_r = 0.684908 |
| z_drag | 1059; explicitly supplied synthetic-control origin |
| h | 0.674; derived parameter convention |
| physical omega_b, omega_c | Omega_b*h² = 0.022259524; (Omega_m-Omega_b)*h² = 0.120837416 |

These values are chosen controls near common reference values, not a Planck
chain sample, posterior mean or new measurement. Radiation/photon fractions
are independently supplied; no exact Tcmb/Neff conversion is asserted. All
pressureless nonbaryonic matter is in omega_c here, with no massive-neutrino
component. Full Planck 2018 normally uses physical baryon/cold-matter densities,
0.06 eV massive-neutrino evolution, thermal predicted drag and CMB perturbations.
Those closures are unavailable in this Irreducible model and remain blocked.
Primordial amplitude/tilt, optical depth, theta and massive-neutrino fields are
structurally rejected by this packet's closed typed model, rather than ignored.

The two original mean/covariance files are identified in
[experiment.json](experiment.json) and the historical input manifest. Their
exact release is `CobayaSampler/bao_data` commit
`bb0c1c9009dc76d1391300e169e8df38fd1096db`, `desi_bao_dr2/`.
The 13 mean rows retain their exact redshift/DM/DH/DV tag order; the full 13x13
covariance uses those same axes. Both unique paths and hashes are required.
These are released fitted distance summaries, not raw galaxy observations.
Release selection, reconstruction, fiducial calibration and compression
assumptions remain conditional. Unknown cross-probe overlap is not independence;
this experiment makes no joint SN/BAO claim. No raw assets are redistributed.
A “gold dataset” means exact versioned bytes, axis identities and reference
checks, not cosmological truth or a guaranteed valid compression for every model.

## Calculation and frozen comparisons

Shared equations, quadrature, ratios and the normalized density are exclusively
Irreducible's `early_late.hpp` and `bao_conditional.hpp` implementations.
The library retains one full SPD Gaussian factor and computes

    log p(y|model,C) = -1/2 [r^T C^-1 r + log det C + 13 log(2 pi)].

The density uses the product of dimensionless ratio coordinates. Numerical
error is not added to observational C. The thin C++ consumer provides only
bounded typed transport, source declarations and output/status handling;
`lcdm-baseline.native` is this experiment's controller operation, not an
implemented Irreducible CLI/ABI operation.

Representative z = 0, 0.1, 0.7, 2.33, 5 checks E, D_M and D_L. All 13 BAO ratios
and quadratic/log determinant/normalization/log density are compared. A second
point doubles H0 with fixed fractions and supplied drag: E and ratios/density
are invariant, distances halve. This changes physical densities and any thermal
mapping. It is a numerical scaling/identifiability control only; a fixed physical
omega_b/omega_c/Tcmb path is outside this packet.

The frozen external allocations are 1e-9 Mpc + 2e-11*|reference| for distances,
2e-11 relative for E, 1e-11 + 5e-11*|reference| for ratios and absolute 1e-8 for
each density component. Independent reference refinement may occupy at most 5%
of each allocation. The separately enforced projection estimate is <=1e-8.
No row is dropped, covariance jitter added or comparison budget loosened.

[reference.py](reference.py) is reference-only: 60-digit standard-library Decimal,
Newton-derived roots and weights of Legendre P8, composite GL8 directly in z and
u=sqrt(a), and explicit dense scalar LDLT. P8 recurrence/derivative and its
weight formula are reconstructed in source; four rational root guesses are
refined to a declared residual. Pi uses the convergent Machin arctangent identity.
Inputs represent the native binary64 values. Panel counts 128 and 256 provide
the independent refinement test. Permanent reference-only analytic controls
check all GL8 polynomial moments through degree15, a known correlated2x2
LDLT determinant/quadratic and radiation-only distance/ruler limits. This differs from production's adaptive
log-redshift/scale-factor rule and Cholesky. Defined units/c are shared conventions,
not independent evidence for physics. No external numerical code/assets were
copied; original code is BSD-3-Clause and Python's standard library retains PSF
terms. The reference is never a production physics or inference route.

Native policies retain distance/sound absolute1e-12 Mpc and relative2e-14,
ratio absolute1e-14 and relative2e-13. Requested masks, per-point callback/depth
limits, global callback allowance4,000,000, point/query/string/native-payload
limits, arithmetic and factor/solve sensitivity are explicit in
[consumer.cpp](consumer.cpp). The native result reports its actual compiled
producer tolerances. That policy is distinct from comparison and projection
acceptance.

## Run and identities

Build/install Irreducible separately at clean commit
`db4765838fc404a489a0115ee69713f2ea6cd2f4`. This packet pins that SDK's exact
manifest, archive and CLI bytes in [request.json](request.json), so a different
legitimate build also needs a reviewed identity update. From the Reproducible
repository root:

```sh
uv sync --frozen
uv run python experiments/lcdm-baseline/controller.py \
  --engine-source ../irreducible \
  --sdk ../irreducible/evidence/project-review/science/next-operators-20261001/final-main/sdk \
  --name fresh-reviewed-attempt
```

### A fresh clone and a new build identity

The SDK path above is ignored local evidence from this machine, not a packaged
release or durable public archive. A source revision alone cannot reconstruct
its exact compiler/library/archive bytes. A fresh clone can instead build a
**new explicitly identified variant** of the same frozen scientific experiment.
Use a clean detached engine checkout and a fresh build/install directory:

```sh
git clone https://github.com/sprajs/irreducible.git ../irreducible-lcdm
cd ../irreducible-lcdm
git switch --detach db4765838fc404a489a0115ee69713f2ea6cd2f4
rustup component add rustfmt
python3 -m venv .build-tools
.build-tools/bin/python -m pip install cmake==4.1.3
cargo fetch --locked
python3 tools/build.py --profile release --jobs 1
.build-tools/bin/ctest --test-dir build/native-release --output-on-failure
python3 tools/check_install.py --profile release
.build-tools/bin/cmake --install build/native-release --prefix ../reproducible/.work/lcdm-sdk
mkdir -p ../reproducible/.work/lcdm-sdk/bin
cp target/release/irred ../reproducible/.work/lcdm-sdk/bin/irred
cp build/build-manifest-release.json ../reproducible/.work/lcdm-sdk/build-manifest.json
git status --porcelain
cd ../reproducible
uv sync --frozen
uv run python scripts/restore_inputs.py --group bao
uv run python scripts/restore_inputs.py --group bao --verify-only
sha256sum .work/lcdm-sdk/build-manifest.json .work/lcdm-sdk/lib/libirred_core.a .work/lcdm-sdk/bin/irred
```

Prepare dependencies according to Irreducible's pinned getting-started guide;
a different host/toolchain needs its own numerical qualification. The source
status command must be empty. Do not reuse an existing SDK/build directory.
Inspect the new manifest/compiler/flags/source hashes and build ID and deliberately
review new `request.json` engine identity values: manifest SHA-256, archive SHA-256,
CLI SHA-256 and build ID, retaining the same clean source revision. A new
request hash also needs a newly reviewed Prospector consumer handoff when that
origin is used. Keep the model, exact data/axes, references, external allocations
and projection gates unchanged; document the new build identity in the history
and use a new result name. The controller intentionally rejects unreviewed
replacement bytes. Run it with `--engine-source ../irreducible-lcdm --sdk
.work/lcdm-sdk` only after that explicit review. The original local receipt and
its original pins remain immutable; a new build is not that exact historical
binary reproduction.

The controller admits only this native experiment; the generic single-request
runner directs native packets here. It creates a new ignored attempt directory
before admission, retains failures and never overwrites completed attempts.
It checks bounded typed inputs, exact two input roles/order/hashes and parse the same hash-admitted bytes, clean pinned
engine source and all 194 manifest source hashes, the complete installed public
header inventory, installed archive/CLI/compiler/standard-library hashes and actual
CLI discovery build. It hashes SDK contents and packet/controller/reference,
request/origin and source manifests before and after execution, including failed
attempts. It compiles only this fixed consumer with recorded command/source and
executable hashes. Output files and child processes have byte, address-space and
time limits; reference work also runs as a bounded child. Completed/failed files
are made read-only. Full request/transport, discovery, stdout/stderr, independent
references, comparisons and a terminal record remain under ignored `results/`.

Successful numerical checks leave inference unassessed and interpretation
conditional on this massless variant. Agent orchestration identities belong in
separate launch evidence; a CLI run does not infer a model name. A durable
archive is still needed before publication of full scientific results; local
ignored receipts are not a preservation service.

## History

2026-10-01: original strict producer ratio relative2e-14 (absolute1e-14) refused
all requested ratios with `conditioning_budget_exceeded`; the ruler and provider
batch were admitted. The original source hashes, compiled binary, policies, 22,206
callbacks, partial stdout and stderr remain in `initial-frozen-control`,
`diagnosed-frozen-control` and `detailed-frozen-control` below the local result
store. The scalar arithmetic floors and numerator/ruler interval propagation
exceed that ratio allowance. After independent source review, only the producer
ratio relative allowance changed to2e-13, still far below the unchanged external
ratio allocation. No physics, data axes, reference, comparison or projection gate
changed. Every final run executes the original under-budget case as a negative
control and requires that same conditioning refusal.

2026-10-01: `feasible-ratio-policy-control` passed 64 named comparisons using
real exact 13-row DESI DR2 bytes. Quadratic28.070790034811825,
logdet-37.007781378241134, normalization23.892401863321492 and
logdensity-7.477705259946092. Independent quadratic difference4.20e-14 and
logdensity difference2.15e-14 were within absolute1e-8. Largest external-budget
fraction4.70e-6; reference-refinement fraction3.27e-25. Projection estimate1.51e-10
passed the unchanged1e-8 gate; aggregate successful callbacks58,208. This is a
parameter-point numerical result, not a best fit or posterior reproduction.
Subsequent provenance hardening removes inferred agent names and requires exact
unique input paths and full SDK header inventory; earlier receipts remain intact.
These initial dirty-source attempts recorded source hashes but did not retain
complete original source bytes; they are diagnostic history, not independently
reconstructible committed-source runs. Final runs require committed source and
retain exact source snapshots.
