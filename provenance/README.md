# Provenance

The repository separates current execution paths from immutable historical evidence. Renaming a file does not justify rewriting the hash or path recorded by a past execution.

The [storage guide](../docs/local-storage.md) describes restoration after the
October 2026 cleanup. `downloads.json` expands the frozen download registry;
`frozen-local-inputs.tar.gz` preserves the small inputs without public routes.

| Record | Meaning |
|---|---|
| [inputs.json](inputs.json) | Current input locations, sizes, SHA-256 values, source identities and retrieval URLs |
| [studies.json](studies.json) | Additional original research sources, experimental specifications and scientific notes, with explicit exclusions |
| [study-inputs.json](study-inputs.json) | Historical acquisition identities and download routes for study dependencies |
| [code-origin.json](code-origin.json) | Retained numerical kernels and their original source hashes/symbols |
| [hst-dark-native.json](hst-dark-native.json) | Native calibration provenance for the two frozen dark images |
| [layout.json](layout.json) | Mapping from the earlier repository to the current structure; documented changes to execution paths |
| [edition.json](edition.json) | Current code, environment and documentation inventory |
| [Validation manifest](../validation/manifest.json) | Verified inputs, run outputs and current independent audit code |
| [history](history/) | Immutable manifests from assembly and the September 26 revalidation |
| [licenses.json](licenses.json) | Third-party licensing and attribution inventory |
| [literature.json](literature.json) | Sources for the prospective experimental plan |
| [age-correction-literature.json](age-correction-literature.json) | Primary-source review behind the age-correction experiments |
| [age-correction-execution.json](age-correction-execution.json) | New observational comparisons, controlled recovery, transport studies and their exact validation boundaries |
| [physical-program-execution.json](physical-program-execution.json) | Recovered galaxy spectra and images, physical age bounds, native survey experiments and explicitly unresolved physical gates |
| [physical-program-extension.json](physical-program-extension.json) | Additional nebular observations, infrared source-separation and counterpart measurements, and enlarged native correction simulations |

`path` fields in the input manifest identify current files. `source_path`, source-record paths and historical run configurations describe the workspace at the time the evidence was produced. They may use older names. Those names are preserved as evidence, not presented as active project organization.

## Original research and earlier execution records

The complete pre-reorganization repository is preserved at commit [`17487bf6`](https://github.com/sprajs/reproducible/tree/17487bf659fcbdeeea072221492bac14b04a0a85). Original acquisition, preparation, analysis and diagnostic code, experimental specifications and scientific notes are available in [studies](../studies/README.md). Vendored builds, handoff documents, duplicate source snapshots and bulk generated result sets remain outside the active source tree. No Git history was rewritten.

For example, an original source can be inspected without restoring the old working tree:

```bash
git show 17487bf6:scripts/cosmology/core.py
```

The original local bulk acquisitions and environments were retained outside this checkout. They are not required by the current workflows once the frozen input bundle is present. A Git-only clone still requires the documented public acquisitions and derived input bundle; see [data requirements](../docs/methods/data.md).

Historical manifests in `history/` apply to the source snapshot and keep their original hashes. They should not be used as inventories of the reorganized working tree. Current runs verify their own input and output records through `validation/record.py`; differences between historical and current orchestration code are explicitly recorded in `layout.json`.

The age-correction execution record is the snapshot at commit `93c2026`; its manuscript hashes are historical once the physical extension is added. The later physical-program record and edition inventory bind the updated documentation. Past execution timestamps and hashes are not rewritten to imply that earlier runs used later code.

The first physical-program execution record is the snapshot at commit `ca8cfaab`. Its scientific source and result identities remain unchanged during the subsequent extension; its shared manuscript hashes describe that earlier commit. The separate extension record and current edition inventory bind the later observations, enlarged simulations and revised documentation.
