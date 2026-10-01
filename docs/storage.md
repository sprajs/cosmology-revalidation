# Data, results and history

Git stores experiment descriptions, configuration, orchestration/visualization
source and small provenance manifests. `data/`, `downloads/`, `results/`, `runs/`,
`simulation/`, `simulations/`, `notebooks/` and `.work/` are ignored. The repository
checker rejects tracked generated products and caps the public tree at 1 MiB
and each packet at eight source files. Change a budget only with an explicit
review of the storage design. `.gitignore` alone cannot untrack
previously committed files. Never use `git add -f` for those stores.

New input declarations name immutable source versions/URLs, bytes, SHA-256,
scientific role, units/axes, calibration, selection and dependence. Acquiring
bytes is not scientific validation. Third-party terms remain applicable even
when downloads are public; original source data are preserved during cleanup.

## Retained historical inputs

[legacy-inputs.json](../sources/legacy-inputs.json) retains 391 frozen identities
from the old campaign, including 34 inputs restored from a local archive.
Original labels and `source_path` fields describe that historical workspace;
they are not current executable paths or a guarantee that an input is raw data.
The old licensing inventory is retained separately in
[legacy-licenses.json](../sources/legacy-licenses.json).

The approximately 22 MiB frozen input archive has been moved to ignored
`data/legacy/frozen-local-inputs.tar.gz`. The recovered local-code archive, its
index and full wider download manifest are also preserved locally under
`data/legacy/`, with their historical Git copies available at the
[previous commit](https://github.com/sprajs/reproducible/tree/1f11997e8935a515cb5f7bb5312c91a38b6c1558/provenance).
No original archive was deleted during this reset. Local ignored storage still
needs a separate backup if relied upon.

In a fresh clone, recover the frozen input archive from that historical snapshot:

```sh
git fetch origin 1f11997e8935a515cb5f7bb5312c91a38b6c1558
mkdir -p data/legacy
git show 1f11997e8935a515cb5f7bb5312c91a38b6c1558:provenance/frozen-local-inputs.tar.gz > data/legacy/frozen-local-inputs.tar.gz
```

Then inspect or restore only the inputs you need:

```sh
python scripts/restore_inputs.py --list
python scripts/restore_inputs.py --group bao --dry-run
python scripts/restore_inputs.py --group bao
python scripts/restore_inputs.py --group bao --verify-only
```

No selection defaults to all 391 historical inputs. Changed existing files,
changed download bytes and an incorrect archive hash are refused. Archive
members are copied individually after checking their expected type/size/hash;
the archive is never blindly extracted. Some public URLs may be unavailable.
Native builds and derived input regeneration need the old documentation and
reviewed new implementations. Restoration does not reproduce an experiment.

## Reset and clone size

The reset removes old generated result trees, Python numerical code, historical
bulk records and the old environment from the active source layout. A short
[historical packet](../experiments/legacy-supernova/README.md) records the scope,
limits and fixed snapshot. The remote Git history was not rewritten; a normal
full clone can still be large. For the current workspace:

```sh
git clone --depth 1 https://github.com/sprajs/reproducible.git
```

Delete disposable run outputs only after recording useful findings and ensuring
published full receipts have an archive. Preserve input originals and unresolved
failures. A scientific result must remain auditable even if the working run store
is later removed.
