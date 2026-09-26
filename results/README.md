# Results

These are the validated results used in the [manuscript](../README.md), originally executed on 26 September 2026. Directory names identify the scientific role rather than the development session.

| Directory | Contents |
|---|---|
| [baseline](baseline/) | The 19 default workflows, one result set per workflow |
| [alternatives](alternatives/) | Joint CPL, flexible expansion, R19-priority age catalogue, colour-only DES prediction and the imposed age template |
| [robustness](robustness/) | Additional random seeds for four cosmology fits and DES prediction |

Within each analysis, `summary.json` is the scientific result; CSV files contain the associated tables; `run.json` records settings, input/code/output hashes, dependencies and execution status. Large sampler arrays are local and hashed in the run record.

The files were moved without changing their numerical contents or execution records. Some configurations inside old run records therefore retain their original execution-time paths. Those are historical records, not current command recipes: use [configs](../configs/) and the [workflow guide](../docs/workflows.md) for new runs. The [path map](../provenance/layout.json) links earlier locations to the public structure.

Default workflows were rerun and their summaries matched the previous edition exactly. Independent validation and the subsequent structure verification are described in [validation](../validation/README.md). Successful reproduction does not establish physical identification or a complete historical survey reduction.
