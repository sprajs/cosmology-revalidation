# Physical survey selection and correction response

This experiment regenerates simulated supernova photons through DES cadence, noise, pipeline detection, host-redshift selection, native SALT fitting, reconstructed SNNV19 signal selection, and a correction learned on an independent simulation pool. It asks whether an imposed age signal is necessarily erased by these operations. The W22 ages are **mocked progenitor delays**, not measured ages. The experiment does not identify a real cosmological correction.

The results and limitations are in [the scientific note](../../notes/survey-physics-results.md). Small manifests and response covariances are in `../../results/survey_physics/`. Downloaded inputs, binaries, complete failed attempts, generated light curves, classifier outputs, bootstrap draws and native logs remain under the ignored `.work/survey-physics/` directory.

## Inputs and build

`recover.py --extract` downloads the checksum-verified [November 2025 SNDATA_ROOT release](https://zenodo.org/records/17591282), pinned DES Git-tree metadata, and the [PLAsTiCC late-time extrapolation input](https://zenodo.org/records/6672739). No released models, downloadable simulations or classifier weights are committed here. `acquisition_audit.py` checks the 250 original Dovekie LFS identities against their Git blobs and the real LFS batch endpoint. The alternate bundle restores executable ingredients; it does not restore those original mock byte streams.

The native source is [SNANA commit 886408a4e171896db5eaa97e735a655f50cec2db](https://github.com/RickKessler/SNANA/tree/886408a4e171896db5eaa97e735a655f50cec2db). The executed local build was recovered from the repository's preserved research archive, with only the existing output-precision/Hessian instrumentation patch. The same patch is versioned at `studies/light_curve_fitting/code/official-inputs/snana-output-only.patch`.

For a fresh installation, clone that pinned SNANA revision into `.work/survey-physics/SNANA`, apply that output-only patch, then run its supported `autoreconf -fi`, `./configure`, and `make -C src -j4 snlc_sim snlc_fit SALT2mu` build. It requires a C/Fortran toolchain, CFITSIO, GSL, and the native build dependencies documented by SNANA. This does not install or alter system packages. `SURVEY_PHYSICS_LIBRARY_PATH` can select an existing dependency-library directory; the executed host defaults to the archived sysroot. `SURVEY_PHYSICS_ARCHIVE` overrides the location of that preserved archive. Build/source/output hashes are recorded in the validation manifest.

For the two narrowly scoped native sensitivities, copy the build to `.work/survey-physics/SNANA-legacy-eff`, relocate its generated Makefiles to that copy, and apply `legacy_efficiency.patch` and `bbc_bounds.patch`. Rebuild `snlc_sim` and `SALT2mu` there. The first makes legacy post-interpolation probability bounding explicit; the second moves an erroneous post-write capacity check to the correct pre-write index check. Neither patch is used to generate the primary campaign's photons.

The classifier uses [SuperNNova revision fcf8584b64974ef7a238eac718e01be4ed637a1d](https://github.com/supernnova/SuperNNova/tree/fcf8584b64974ef7a238eac718e01be4ed637a1d), Python 3.10 and PyTorch 1.13.1 CPU, with model and normalization files from the official bundle. On this host `classify_all.py` locates the preserved environment. On another host, provide `--python` and `--pythonpath` pointing to an environment containing that revision and its dependencies. The model consumes measured peak timing and measured redshift; generated true ages, type, peak and redshift are not classifier inputs. Passing batch/single checks does not certify equivalence to original published probabilities.

## Reproduction

From the repository root, with its Python dependencies installed. Set `SURVEY_CLASSIFIER_PYTHON` to the absolute Python 3.10 executable containing PyTorch 1.13.1 CPU, and `SURVEY_CLASSIFIER_SOURCE` to the checked-out SuperNNova revision above. The explicit overrides below prevent a fresh checkout from silently relying on this machine's archived environment:

```bash
.venv/bin/python studies/host_ages/code/survey_physics/recover.py --extract
.venv/bin/python studies/host_ages/code/survey_physics/acquisition_audit.py
.venv/bin/python studies/host_ages/code/survey_physics/prepare.py
.venv/bin/python studies/host_ages/code/survey_physics/efficiency.py
.venv/bin/python studies/host_ages/code/survey_physics/run_native.py --workers 4
.venv/bin/python studies/host_ages/code/survey_physics/classify_all.py --python "${SURVEY_CLASSIFIER_PYTHON:?Set classifier Python path}" --pythonpath "${SURVEY_CLASSIFIER_SOURCE:?Set pinned SuperNNova source path}"
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/survey_physics/analyze.py
.venv/bin/python studies/host_ages/code/survey_physics/prepare_positive_dust.py
.venv/bin/python studies/host_ages/code/survey_physics/run_native.py --campaign positive-dust-campaign.json --workers 4
.venv/bin/python studies/host_ages/code/survey_physics/classify_all.py --campaign positive-dust-campaign.json --python "${SURVEY_CLASSIFIER_PYTHON:?Set classifier Python path}" --pythonpath "${SURVEY_CLASSIFIER_SOURCE:?Set pinned SuperNNova source path}"
OPENBLAS_NUM_THREADS=1 .venv/bin/python studies/host_ages/code/survey_physics/analyze.py --positive-dust
.venv/bin/python studies/host_ages/code/survey_physics/native_checks.py legacy
.venv/bin/python studies/host_ages/code/survey_physics/native_checks.py multistart
.venv/bin/python studies/host_ages/code/survey_physics/native_checks.py bbc
.venv/bin/python studies/host_ages/code/survey_physics/positive_bbc.py
.venv/bin/python studies/host_ages/code/survey_physics/validate.py
```

A completed native run is reused only when its frozen input hash matches. Classifier outputs likewise belong to those particular fit tables. To deliberately change a model or seed, create a new campaign/version name rather than reuse existing output names. `prepare.py` and `prepare_positive_dust.py` overwrite generated manifests; the committed manifests describe the completed campaign, including the documented engineering amendments, and should be preserved when reproducing in a new checkout.

## Interpretation safeguards

The original P23 low-RV tail includes invalid negative optical extinction; the `positive_extinction` analysis is only a post-selection subgroup sensitivity. The separate fixed-RV=3.1 campaign changes the population **before** photon generation and supports the physical conclusion within its declared dust family. Neither population is an empirical determination of the true host dust distribution. The released configuration sets host-mass errors to zero, so the correction sees noiseless simulated host mass; real host SED uncertainty is absent.

The independent 40-neighbor predictor subtracts one fitted standardization and one training residual estimate. It does not apply BBC on top. Its uncertainty resamples cadence LIBID blocks, coupling paired arms and independently resampling training and evaluation pools. It conditions on shared released cadence, host library, surface and classifier; it is not an independent-survey calibration. Native production-dimensional BBC is tested separately and fails the available training support. The successful native 1D BBC benchmark is a different estimator.

Five fixed redshift bins and 40-neighbor radius/occupancy gates prevent unsupported extrapolation. Bootstrap covariance explicitly reports the number of jointly supported draws. Failed support remains null. A nonzero nominal closure age association measures remaining model misspecification; it must be separated from the *change* caused by an injected term. The frozen nominal-correction response and the age-refitted correction response are different estimands.
