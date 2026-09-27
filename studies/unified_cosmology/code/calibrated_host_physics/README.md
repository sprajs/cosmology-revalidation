# Conditional stellar ages from native host spectra

Read the [scientific results and limits](../../notes/calibrated-host-physics.md). The input is all 55 independently associated DESI host spectra, not a catalogue of age medians. This module produces conditional formed-mass age compatibility sets, not a host-age posterior, supernova correction or cosmological likelihood.

First recover the [calibrated host spectra](../calibrated_hosts/README.md) and the pinned stellar data using the acquisition step in the earlier [physical-age package](../../../host_ages/code/physical_ages/README.md). Its existing low-resolution environment and generated libraries remain untouched. This module checks the complete original stellar-data inventory, then compiles its own high-resolution binding. Merely installing a default `fsps` wheel would select the wrong spectral library.

From the repository root:

```bash
.venv/bin/python studies/unified_cosmology/code/calibrated_host_physics/acquire.py --build-environment
.venv/bin/python studies/unified_cosmology/code/calibrated_host_physics/acquire_resolution.py
export SPS_HOME="$PWD/.work/physical-ages/fsps"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
.work/unified-cosmology/calibrated-host-physics/.venv/bin/python studies/unified_cosmology/code/calibrated_host_physics/build_library.py
.work/unified-cosmology/calibrated-host-physics/.venv/bin/python studies/unified_cosmology/code/calibrated_host_physics/build_library.py --refined
.work/unified-cosmology/calibrated-host-physics/.venv/bin/python studies/unified_cosmology/code/calibrated_host_physics/fit.py
.work/unified-cosmology/calibrated-host-physics/.venv/bin/python studies/unified_cosmology/code/calibrated_host_physics/fit.py --refined
.work/unified-cosmology/calibrated-host-physics/.venv/bin/python studies/unified_cosmology/code/calibrated_host_physics/independent_review.py
.work/unified-cosmology/calibrated-host-physics/.venv/bin/python studies/unified_cosmology/code/calibrated_host_physics/report.py
.work/unified-cosmology/calibrated-host-physics/.venv/bin/python studies/unified_cosmology/code/calibrated_host_physics/validate.py
```

`acquire.py` uses an available `FC`/`gfortran`. On this Arch host, no Fortran frontend was installed; the executed route extracted an official signed, hash-pinned compiler package locally and used the matching existing GCC support files. It did not install system packages. A different host should provide a suitable `FC`. Build flags are `-DC3K_LR=0 -DC3K_HR=1`; the acquisition command verifies that the active library is `c3k_hr`. [Dependencies](requirements-lock.txt), [acquisition evidence](../../results/calibrated_host_physics/acquisition.json).

The main environment is unchanged. Raw spectra, the source reference, stellar libraries, fit fixtures and per-case arrays remain under ignored `.work/`. Committed results contain compact summaries, per-host intervals, hashes and validation records. `--targetids` and `--scenarios` on `fit.py` make an explicitly named subset result; they do not overwrite the full result. Refinement uses the eight preselected IDs in [implementation-design.json](implementation-design.json).

`model.problem(targetid, Library(), scenario)` returns `(A, y, C, ages, mass, parameters, metadata)`. The five observed bands and model columns use the same native camera, pixel and emission masks. `C` includes overlap between bands. Stellar columns have `mass=1`; free emission columns have `mass=0` and zero age numerator. The single shared stellar normalization preserves formed-mass ratios. `solve(A,y,C,ages,mass)` returns the best objective, compatibility status, lower/upper age extrema, primal witnesses and numerical diagnostics. The generated fixtures preserve the full matrices and coefficients for independent equation checks.

These finite families do not supply an empirically calibrated age prior or constrain a grey supernova luminosity shift. Statistical compatibility, astrophysical model adequacy, population selection and supernova standardization remain distinct questions.
