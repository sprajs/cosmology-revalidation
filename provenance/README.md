# Provenance

The repository separates current execution paths from immutable historical evidence. Renaming a file does not justify rewriting the hash or path recorded by a past execution.

| Record | Meaning |
|---|---|
| [inputs.json](inputs.json) | Current input locations, sizes, SHA-256 values, source identities and retrieval URLs |
| [code-origin.json](code-origin.json) | Retained numerical kernels and their original source hashes/symbols |
| [hst-dark-native.json](hst-dark-native.json) | Native calibration provenance for the two frozen dark images |
| [layout.json](layout.json) | Mapping from the earlier repository to the current structure; documented changes to execution paths |
| [edition.json](edition.json) | Current code, environment and documentation inventory |
| [Validation manifest](../validation/manifest.json) | Verified inputs, run outputs and current independent audit code |
| [history](history/) | Immutable manifests from assembly and the September 26 revalidation |
| [licenses.json](licenses.json) | Third-party licensing and attribution inventory |
| [literature.json](literature.json) | Sources for the prospective experimental plan |

`path` fields in the input manifest identify current files. `source_path`, source-record paths and historical run configurations describe the workspace at the time the evidence was produced. They may use older names. Those names are preserved as evidence, not presented as active project organization.

## Retrieve earlier work

The complete pre-reorganization repository is preserved at commit [`17487bf6`](https://github.com/sprajs/cosmology-revalidation/tree/17487bf659fcbdeeea072221492bac14b04a0a85). The public working tree no longer includes redundant experiments, vendored builds, handoff documents or superseded result sets. No Git history was rewritten.

For example, an original source can be inspected without restoring the old working tree:

```bash
git show 17487bf6:scripts/cosmology/core.py
```

The original local bulk acquisitions and environments were retained outside this checkout. They are not required by the current workflows once the frozen input bundle is present. A Git-only clone still requires the documented public acquisitions and derived input bundle; see [data requirements](../docs/methods/data.md).

Historical manifests in `history/` apply to the source snapshot and keep their original hashes. They should not be used as inventories of the reorganized working tree. Current runs verify their own input and output records through `validation/record.py`; differences between historical and current orchestration code are explicitly recorded in `layout.json`.
