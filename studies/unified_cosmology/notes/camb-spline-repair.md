# Isolating the CAMB nonlinear-lensing spline repair

The installed CAMB 1.6.6 package matches the public PyPI wheel byte for byte: all 41 wheel payload files, including the compiled library and its bundled runtime libraries, agree. Its `RECORD` also verifies. The public source distribution matches all 65 common CAMB/Fortran source files in the tagged release and 22 bundled files in its pinned `forutils` submodule. These establish package and source provenance; the wheel does not contain a complete reproducible-build attestation.

The released source contains the array-indexing defect repaired by [upstream commit 56a95f78](https://github.com/cmbant/CAMB/commit/56a95f78fdd72709c5af6b668141ec0c617777c3). A time array starts at index 0, but physical knots start at 1. Passing the whole array to an explicit-shape spline associates the first ordinate with element 0 and displaces the remaining time coordinates. The repair passes the explicit `1:n` section. Our separate Fortran test confirms this association with a deliberately poisoned element 0. It does **not** demonstrate that the installed calculation actually encountered a NaN or establish its likelihood error.

Two isolated packages have now been built from the identical 1.6.6 source distribution, with the same GCC 16.2.1 compiler and flags. Their only source difference is that argument slice. Both contain exactly the same 31 Python/data files as the public wheel, including the high-l extrapolation and BBN tables. Their environments read other dependencies from the unchanged original environment. The installed wheel was compiled with GCC 14.2.1, so the numerical comparison must retain three arms: original wheel, unpatched local control, and repaired local source. A difference between the last two isolates the code repair under the common local build; the first comparison measures packaging/compiler effects.

The original make-only build receipt is preserved. An explicit second stage copies the high-l template, as upstream `setup.py` normally does, and verifies the full Python/data package. A first preparation attempt also remains archived: it failed while converting a symlinked interpreter path, before a plan or physical evaluation. The corrected preparation preserves the lexical virtual-environment launch path.

The [prepared audit](../results/inference/camb-spline-build-preparation.json) freshly requalifies the original posterior and thermal-reviewed numerical screen. It fixes original screen indices 0 and 31, accuracies 1 and 2, and the three implementations: **12 proposed native likelihood evaluations, none executed here**. Physical parameters, data, calibration priors, full component likelihoods and numerical settings remain those of the corresponding original calculation. The plan includes separate original-versus-control, control-versus-repair, and accuracy2-versus1 comparisons. Raw spectra, actual finalized settings, derived quantities and signed component shifts must be retained. This is a numerical implementation audit, not a posterior update; no existing convergence or statistical gate is relaxed.

Reproduce the build and preparation in a fresh ignored directory:

```bash
python studies/unified_cosmology/code/external_probes/camb_spline_build.py --work .work/unified-cosmology/camb-spline-audit-replay
python studies/unified_cosmology/code/external_probes/camb_spline_finalize.py --work .work/unified-cosmology/camb-spline-audit-replay
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 CLIPY_NOJAX=1 \
.work/unified-cosmology/external-probes/.modern-venv/bin/python \
studies/unified_cosmology/code/inference/camb_spline_audit.py \
  --screen studies/unified_cosmology/results/inference/native-posterior-precision-modern-lcdm-independence.json \
  --review studies/unified_cosmology/results/inference/native-precision-thermal-review-modern-lcdm-independence.json \
  --build-record .work/unified-cosmology/camb-spline-audit-replay/build-record.json \
  --work .work/unified-cosmology/camb-spline-audit-replay/prepared
```

The default compiler uses the already acquired local toolchain; `--compiler` accepts a separately verified compiler. No system package, active environment, original chain or frozen scientific source is changed. [Broader v2 audit](camb-v2-source-audit.md) explains why simply upgrading CAMB would mix this repair with other numerical and physical-model changes.

Work stopped at the user's request to conclude from existing evidence. The unvalidated execution prototype and its design were archived in ignored local storage; they are not part of the public executable research. The [completed preparation](../results/inference/camb-spline-build-preparation.json) and [incomplete execution status with archive hashes](../results/inference/camb-spline-executor-status.json) remain public. No execution validation, independent executor review or physical comparison was completed. All twelve proposed physical evaluations remain an explicit gap, and no numerical or cosmological effect of the repair is claimed.
