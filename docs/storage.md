# Data, results and history

Git stores experiment descriptions, configuration, orchestration/visualization
source and small provenance manifests. `data/`, `downloads/`, `results/`, `runs/`,
`simulation/`, `simulations/`, `notebooks/` and `.work/` are ignored. The repository
checker rejects tracked generated products and caps the public tree at 2 MiB
and each packet at eight source files. Change a budget only with an explicit
review of the storage design. `.gitignore` alone cannot untrack
previously committed files. Never use `git add -f` for those stores.

New input declarations name immutable source versions/URLs, bytes, SHA-256,
scientific role, units/axes, calibration, selection and dependence. Acquiring
bytes is not scientific validation. Third-party terms remain applicable even
when downloads are public; original source data are preserved during cleanup.

## Retained historical inputs

[legacy-inputs.encoded.json](../sources/legacy-inputs.encoded.json) retains 391 frozen identities
from the old campaign, including 34 inputs restored from a local archive.
Original labels and `source_path` fields describe that historical workspace;
they are not current executable paths or a guarantee that an input is raw data.
The old licensing inventory is retained separately in
[legacy-licenses.json](../sources/legacy-licenses.json).
The transport recovers the exact historical JSON bytes, including all 391
identities, field values and order. Older manifest byte hashes remain bound
to their original Git revisions and run records.

2026-10-03 storage review: compacting only manifest whitespace saved 35,611 bytes,
leaving about 70 KiB under the old 1 MiB cap. The new source-only SDSS reader and
the released active46 box consumer need separate controllers, reference/test
source and immutable design contracts while preserving historical packets.
Their existing transport/contract drafts alone total about 85 KiB. The reviewed
2 MiB source cap accommodates these concrete consumers; the eight-file packet
cap and exclusions for generated products and bulk inputs remain enforced.
This storage allowance does not qualify an experiment or permit redistributing
third-party inputs.
Existing request resource limits and run source pins remain unchanged. A rerun
uses its reviewed source revision; the repository allowance does not enlarge
an older controller's source snapshot or numerical/work allocations.

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

## Exact-byte metadata transport

This representation retains the historical 255,902-byte legacy manifest exactly:
SHA256 `7ed81f37e617ce7bb2638f30596165f46ae7ee81ca31788177e90cf57ad5e56f`,
commit `23ebabd606aaceaca469de59c70ec6d7bed87989`, path
`sources/legacy-inputs.json`. The new encoded file and decoder have separate
source identities. A compressed representation is never labelled with the old
raw-file hash. Historical packet source labels such as `sources/legacy-inputs.json`
refer to those decoded original bytes; historical requests, SDKs, snapshots and
receipts remain unchanged. Future attempts record their actual new committed
source and exact current transport/decoder plus decoded-source identities.

The named reader is [metadata_source.py](../scripts/metadata_source.py):
`python -B scripts/metadata_source.py legacy-inputs` writes the exact original JSON
bytes; append `--identity` for all three identities. The restore tool and LCDM
baseline controller both use that same reader. Its closed pinned envelopes admit
only canonical RFC4648 base64, one complete RFC1950 zlib stream without trailing
or unused/unconsumed bytes, exact decoded byte length/SHA, strict finite integer
JSON, unique keys and bounded typed schema. Refusal is explicit before restoration;
there is no fallback to another source or shape. It has no download or physics
operation. Generated decoded inspection files remain ignored, and outputs do not
become observations through decoding.
Invalid inspection arguments refuse before source reads, with a null document
identity and bounded JSON diagnostic; supplied argument text is never echoed.

The source-only [asset-v2 ledger](../sources/actual-data-asset-status.encoded.json)
and [identity companion](../sources/actual-data-asset-status-input-identities.encoded.json)
use the same exact-byte transport. Inspect with `metadata_source.py asset-status`
or `asset-identities`; the [dated provenance account](../sources/actual-data-asset-status-v2.md)
is retained byte-for-byte. Its local paths, old evidence dates and review-pending
wording are historical source metadata, not assertions that those files exist in
fresh clones. The named companion route supplies the public representation of
`fresh_source_identity_ref` without rewriting that original field. All90 asset
rows,15 retained receipt maps,22 source pins, statuses, failures, rights and unknown
covariances remain exact. No primary assets or large receipts are redistributed.

The checker still counts these real encoded source bytes and all decoder,
caller, documentation and test bytes against2MiB. Every packet remains at most
eight source files. This representation changes no cap or scientific policy.
The LCDM baseline's old1MiB source snapshot bound still prevents an attempt from
this complete current tree; use its original reviewed source for historical reruns.
Storage savings alone do not admit NEXT17, a new Box route or all19 execution.
