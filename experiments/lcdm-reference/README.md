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
