# Full LCDM external reference

This experiment asks whether one explicit full ΛCDM state can produce lensed
CMB spectra, matter power, distances and a predicted drag scale, then feed a
released primary CMB likelihood. CLASS supplies the physics. Reproducible owns
the parameter mapping, ordered input/output adapters, bounded orchestration and
account. This is a reference target for native development; the existing
[native baseline](../lcdm-baseline/README.md) retains its historical identities
and qualification boundaries.

[reference.json](reference.json) fixes CLASS v3.3.0, its original example point,
explicit physical densities/species/recombination/reionization/primordial state,
and four small prediction cases. The source example is conditioned on a joint
Planck primary-plus-lensing result. It is an illustrative starting point here,
not an asserted primary-only best fit or posterior draw. Helium is fixed to the
source point; no new BBN-law agreement is claimed.

The cases evaluate the anchor, a source-supplied precision variant and
`n_s ± 0.005` with every other coordinate fixed. Retain the original CLASS tables
and ordered multipoles. Its dimensionless `D_l` is converted to `C_l` and µK²
for the released likelihood; the linear and Halofit `P(k,z)` products remain
separate. Background H and distance units, printed table precision, six-decimal
drag output and interpolation limits are recorded explicitly. Precision-case
differences are empirical diagnostics, not a certified error bound.

The dedicated command, after a clean pinned CLASS build and an admitted runtime
configuration, is:

```sh
python -B experiments/lcdm-reference/run.py \
  --config .work/lcdm-reference-20261003/config.json \
  --attempt results/lcdm-reference/<fresh-attempt> \
  --deadline-utc <root-grant-deadline-with-UTC-offset>
```

The runtime configuration binds the actual source manifest, compiled executable
and build receipt. Fresh attempt directories retain configurations, raw output,
failures and derived summaries under ignored local storage. This controller does
not add a route to the single-request Irreducible packet runner. External source,
likelihood data and runtime installations remain local and are not redistributed.

The first primary likelihood target is explicitly **Plik-lite TTTEEE + Commander
low-ℓ TT + SimAll low-ℓ EE**, using official clik and the exact released products.
Plik-lite marginalizes foreground nuisance coordinates and retains calibration;
the runtime names/order and all fixed nuisance values must be explicit. Lensed
theory is required. This target excludes the separate lensing likelihood.
Likelihood values and separately declared prior terms remain distinguishable.
The small spectral-index scan is conditional, not a six-parameter fit or a
normalized posterior. DESI BAO and any source-supported SN scores are displayed
separately until their lineage, covariance and prior contracts permit a joint
target.

2026-10-03: all four prediction cases completed in 26.76 seconds using a clean
CLASS build from `0ceb7a9a4c1e444ef5d5d56a8328a0640be91b18`, executable SHA256
`7361da65a9bb037bdc6346dba6408add16114ac4d25c75d8d4df477dd49077d9`.
Each produced lensed TT/EE/TE at ℓ=2…2508, separate linear/Halofit matter spectra
and the same predicted drag ruler, **147.054261 Mpc** at z=1059.928342. Eleven
focused unit/axis/input checks passed. The first attempt completed the CLASS
anchor but failed to read its files: CLASS appends an underscore to the requested
root. The corrected attempt has a new controller identity; both attempts and all
original outputs remain in ignored local storage.

The precision variant changes TT/EE/TE by at most 0.461901/0.0102618/0.0369078 µK²
in `D_l`; these empirical differences have no certified error bound. The linear
matter comparison reaches about 0.0023133 fractional difference while
interpolating onto a different k grid, with one anchor row outside the overlap;
this combines solver and interpolation effects. Background and printed drag
values agree between these two cases. The ns points change the spectra while
holding all other coordinates fixed. No official likelihood, BAO or SN score,
observational fit or native CMB qualification is earned yet. Full local attempt
record SHA256 is `8eb4f479207f53f7889f29ecb395dbdc4304a19d3404623f030b98f3a3340eee`;
a durable publication archive has not been created.

Later on 2026-10-03, the official ESA PLC 3.01 C/Fortran library was built from
the original release archive. The four retained prediction products were scored
in one process, reusing the three initialized likelihood objects without rerunning
CLASS. All three released test-vector checks passed the previously fixed
`1e-6` criterion in raw log-likelihood units; their printed differences were
−1.07424e−9 (Commander), −4.1778e−8 (SimAll) and −3.86228e−9 (Plik-lite).
These checks establish operational repeatability, not a prediction error bound.

| Case | Commander TT | SimAll EE | Plik-lite TTTEEE | Raw log-likelihood sum | Change from anchor |
| --- | ---: | ---: | ---: | ---: | ---: |
| Anchor | −11.627405 | −198.026730 | −292.335511 | −501.989646 | 0 |
| Changed precision | −11.627793 | −198.026730 | −292.426326 | −502.080849 | −0.091203 |
| `n_s − 0.005` | −12.207825 | −198.079388 | −291.832381 | −502.119593 | −0.129947 |
| `n_s + 0.005` | −11.132910 | −198.021347 | −296.083242 | −505.237499 | −3.247853 |

The precision shift is empirical and uncertified; it is material beside the
smaller conditional scan change. This is three fixed spectral-index points,
not a fitted optimum, posterior or uncertainty interval. Raw clik constants
are retained, so the sum is not labelled a goodness-of-fit chi-square.
The actual getters require TT ℓ=0…29 for Commander, EE ℓ=0…29 for SimAll,
and TT/EE/TE ℓ=0…2508 for Plik-lite, followed by `A_planck` in each vector.
These are buffer contracts; the high-ℓ analysis cuts must not be inferred from
the maximum buffer length. ℓ=0,1 inputs are explicit zero placeholders.

`A_planck=1` is fixed. One separately declared relative Gaussian penalty
`−0.5*((A_planck−1)/0.0025)^2` contributes zero to every case. The reviewed
product metadata and generic initializer contain no top-level prior or default.
Component-internal terms and the mapping from this coordinate to the paper's
`y_P` remain unexamined; the external term is a declared construction, not a
certified complete released prior or a normalized posterior. No calibration
uncertainty was measured. BAO/SN combination, inference, observational
interpretation and native full CMB qualification remain open.

The library SHA256 is
`1988ba6c3b4488bbbb0db2841b93ca0223eb441ef3112a6711613079879dbfc5`;
the original scoring receipt is
`a969b48e37355ee8c98a62897eea610e07472e9018d9f36af0605b5d95dabd98`.
All 58 admitted file identities agreed before/after execution. The retained
failed build attempts record Python/Waf compatibility, a GCC configuration
probe correction and build-option/link corrections; likelihood algorithms were
unchanged. Original third-party archives and products remain local.

The retained-product command is [score.py](score.py). Its closed runtime
configuration binds the actual library, resolved Python executable, controller
sources, completed prediction receipt, four product hashes and an explicit flat
inventory of released data and dependencies. Use the getter-derived input order,
the fixed selfcheck criterion and one explicit relative prior. A receipt hash
supplies lineage, not automatic scientific qualification. Execute with an outer
watchdog and capture its exit status and logs, for example:

```sh
timeout --signal=TERM --kill-after=2s 900s \
  python -B experiments/lcdm-reference/score.py \
  --config .work/<admitted-score-config>.json \
  --attempt results/lcdm-reference/<fresh-score-attempt> \
  --deadline-utc <absolute-deadline-with-UTC-offset>
```

The process enforces a stricter 2 GiB address limit and 180-second native phases;
a terminated phase retains its started record without inventing a returned
score. [official_clik.py](official_clik.py) transports the original C API with
managed native errors. The four scored inputs were checked against the actual
SimAll table support before evaluation. Each subsequent attempt needs its own
admitted input, controller and compiled-runtime identities.

[plot.py](plot.py) renders the four retained prediction products after checking
their identities. It uses the emitted sample grids and keeps signed TE and
linear/Halofit matter power distinct. Generated PNG/SVG files remain in ignored
local storage; the source has been rendered and visually inspected.
