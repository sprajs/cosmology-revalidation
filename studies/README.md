# Research studies

The manuscript presents the measured results. These studies retain the original acquisition, preparation, models, statistical tests and scientific findings behind the broader investigation, including negative and unresolved results.

| Topic | Scientific question |
|---|---|
| [Expansion and geometry](expansion/README.md) | Which expansion histories do the distances constrain? |
| [Host ages and brightness residuals](host_ages/README.md) | How do ages and correction definitions change the residual trend? |
| [Dust and colour](dust/README.md) | Where do extinction assumptions become unphysical or non-identifiable? |
| [Calibrated light curves and native fits](light_curve_fitting/README.md) | Do independent flux fits and native uncertainties agree? |
| [Populations, prediction and selection](population_inference/README.md) | Can prediction and population recovery survive selection? |
| [Calibration and spectral degeneracy](calibration/README.md) | Which shared calibration patterns are actually identifiable? |
| [Infrared distances, timing and signed flux](infrared/README.md) | What is measured by optical/infrared differences and timing? |
| [Passbands and photometric lineage](passbands/README.md) | How do filter definitions and conversion choices affect photometry? |
| [Detector repeats and calibration](detectors/README.md) | Which detector combinations are constrained by the repeats? |
| [Spectral models and distance identification](spectral_models/README.md) | Can spectral models identify distance independently of population assumptions? |
| [Physical identities and approximation limits](physical_models/README.md) | Which relations are identities, approximations or empirical models? |
| [Acquisition and source identity](acquisition/README.md) | Can each input be recovered with the intended identity? |

## Reading and validation status

Use the topic summaries first, then the scientific notes for assumptions and interpretation. Notes preserve the evidence from their recorded investigation; old result links point to that pinned historical evidence. Handoffs, project status logs and superseded planning documents are not the reading path.

The main package (`lib`, `workflows`, `research.py`) contains the calculations tested for the manuscript. The additional sources here extend coverage to acquisition, native instrumentation, classifiers, physical checks and experiments that remain inconclusive. They do **not** all have the same validation status. Native SNANA/CALWF3 builds, missing historical assets and unfinished BayeSN or selection inference remain explicit requirements or limitations.

## Execution

Run the manuscript analyses with the [reproduction guide](../docs/workflows.md). To inspect and prepare the additional studies:

```bash
.venv/bin/python studies/manage.py verify
.venv/bin/python studies/manage.py list --topic detectors
.venv/bin/python studies/manage.py prepare --name my-study
```

Preparation copies authored sources and experimental specifications into `.work/my-study/`, using the original relative layout so sibling imports and generated-file dependencies retain their meaning. The public source stays organized by scientific topic. Embedded absolute workspace roots are relocated only in the execution copy, and both hashes are recorded in `preparation.json`. It does not execute code, download inputs, rewrite historical validation claims, or create a functioning native environment automatically.

Restore a required input by its recorded historical path:

```bash
.venv/bin/python studies/manage.py fetch --name my-study --path phase2/official/portable_pilot/assets/salt3_template_0.dat.gz --dry-run
.venv/bin/python studies/manage.py fetch --name my-study --path phase2/official/portable_pilot/assets/salt3_template_0.dat.gz
.venv/bin/python studies/manage.py fetch --name my-study --path phase2/official/portable_pilot/assets/salt3_template_1.dat.gz
.venv/bin/python .work/my-study/scripts/physics_audit/luminosity_check.py
```

The luminosity example needs only the two named model inputs and the primary Python environment. Larger studies require their own acquisition/preparation steps. `list` supplies the exact historical execution path for each source. Run commands from the prepared workspace when their working directory matters, using an absolute interpreter path. Existing workspace destinations and mismatched input hashes are rejected. Inputs are copied, not symlinked back to the frozen bundle.

The [historical Python lock](environment/uv.lock) preserves the earlier base environment, not every optional native or classifier dependency. The current locked environment covers the manuscript package. Additional requirements are visible in study imports and acquisition/build scripts; do not interpret a syntax check as successful execution of those packages.

## What belongs in Git

Keep original source, instrumentation patches, experiment-defining specifications, scientific explanations, source identities and compact manuscript evidence. Keep downloaded code, data, model weights, environments, native build trees, chains and bulk intermediate outputs in ignored local storage. Identical source snapshots share one retained file; superseded execution snapshots remain identified in Git history.

[Source inventory](../provenance/studies.json) maps original locations to the topic structure and records exclusions. [Input identities](../provenance/study-inputs.json) preserve acquisition URLs and hashes. Some excluded vendor files have an explicitly labelled fallback URL into the pinned historical repository; they are fetched into the workspace, never recommitted as authored code. Original acquisition scripts retain the upstream retrieval procedures. Historical source-hash assertions may require their original unrelocated environment; no such assertion is silently updated to make a study pass.
