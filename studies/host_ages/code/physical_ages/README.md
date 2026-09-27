# Physical age bounds from measured light

Read the [scientific interpretation](../../notes/physical-ages-results.md) before treating these compatibility regions as age posteriors. This is a new calculation, not an exact reproduction of C25/R19 age fitting.

First acquire the [galaxy spectroscopy inputs](../galaxy_validation/README.md) and run the [host transport acquisition](../../notes/host-transport-results.md). The required generated interfaces are `.work/galaxy-validation/host-spectrum-age-ledger.csv`, `.work/host-transport/roman-photometry.csv`, `.work/host-transport/des-deep-photometry.csv`, and `.work/host-transport/filters/des-deep-responses.npz`.

From the repository root:

```bash
uv venv --python 3.12.14 .work/physical-ages/.venv
uv pip install --python .work/physical-ages/.venv/bin/python -r studies/host_ages/code/physical_ages/requirements.txt
export OPENBLAS_NUM_THREADS=1
.work/physical-ages/.venv/bin/python studies/host_ages/code/physical_ages/acquire.py
.work/physical-ages/.venv/bin/python studies/host_ages/code/physical_ages/build_library.py
.work/physical-ages/.venv/bin/python studies/host_ages/code/physical_ages/build_library.py --refined
.work/physical-ages/.venv/bin/python studies/host_ages/code/physical_ages/validate.py
.work/physical-ages/.venv/bin/python studies/host_ages/code/physical_ages/fit.py --sample sdss --workers 4
.work/physical-ages/.venv/bin/python studies/host_ages/code/physical_ages/fit.py --sample roman --workers 4
.work/physical-ages/.venv/bin/python studies/host_ages/code/physical_ages/fit.py --sample des --workers 4
.work/physical-ages/.venv/bin/python studies/host_ages/code/physical_ages/sensitivity.py
.work/physical-ages/.venv/bin/python studies/host_ages/code/physical_ages/injections.py
.work/physical-ages/.venv/bin/python studies/host_ages/code/physical_ages/solver_audit.py --prepare
.venv/bin/python studies/host_ages/code/physical_ages/figures.py
```

The main repository environment is unchanged. The physical-age environment deliberately pins its own dependencies; the FSPS wheel contains compiled code, while `acquire.py` retrieves pinned stellar tables and foreground maps. Raw downloads and generated arrays belong below `.work/physical-ages/`. `--limit` is an engineering pilot and is not the full sample. The figure command uses the main environment's pinned Matplotlib.

The optional CUDA calculation requires a working NVIDIA driver. It uses isolated runtime wheels and does not install or modify a system driver:

```bash
uv venv --python 3.12.14 .work/physical-ages/.gpu-venv
uv pip install --python .work/physical-ages/.gpu-venv/bin/python -r studies/host_ages/code/physical_ages/requirements-gpu.txt
.work/physical-ages/.venv/bin/python studies/host_ages/code/physical_ages/gpu_benchmark.py --prepare
.work/physical-ages/.gpu-venv/bin/python studies/host_ages/code/physical_ages/gpu_benchmark.py
```

To independently reproduce the solver-version comparison, create a separate environment with NumPy 2.4.3 and SciPy 1.18.1 and run `solver_audit.py` against the already prepared fixture. The old solver failure is retained as a failed numerical check, not erased or counted as a successful fit. The primary fit uses bounded-variable least squares, not the failing NNLS implementation.

`design.json` preserves the initial model choices and dated pre-fit amendments. `model.py` contains the equations; `fit.py` preserves signs, units, aperture selection, failed models and solver certificates; `sensitivity.py` and `injections.py` expose finite-grid and physical-model dependence. Per-object CSVs and bulk stellar spectra can be regenerated and are excluded from Git.
