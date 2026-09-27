# Survey inputs and likelihood interfaces

The [results and reproduction instructions](../../notes/survey-selection.md) describe the observations, source releases, scientific limits and ordered commands. Run commands from the repository root. Downloaded sources, native binaries, expanded arrays and generated event tables belong under `.work/unified-cosmology/survey-selection`; only authored code and compact checksummed results are tracked.

- `acquire*.py`, `assets.py`, `des3yr.py`: pinned public sources and native simulation assets.
- `dovekie.py`: 1,820-row released distance interface; the supplied arrays store precision, so subsets are extracted from its inverse covariance.
- `des3yr_likelihood.py`: separate 20-bin author likelihood representing 329 historical SNe. It is an alternative to Dovekie, with overlapping physical events.
- `observations.py`, `calibration.py`, `ozdes.py`: exact candidate, calibration and host-spectrum joins.
- `dataprep.py`, `classifier_closure.py`, `classifier_full.py`, `diagnose_classifier.py`: measured timing and observed classifier reproduction, retaining the initial failed header-timing diagnostic.
- `prior_classifier_audit.py`: read-only verification of the earlier physical mock campaign's measured timing inputs.
- `des3yr_checks.py`, `review_geometry.py`, `validate.py`: support, bounded independent geometry and provenance checks.

The ordinary analysis scripts use the repository Python environment (`numpy`, `scipy`, `pandas`, `astropy`, `requests`). Native timing additionally requires a built SNANA tree and its compiler/runtime dependencies. Classifier inference requires the pinned SuperNNova source, Python 3.10 and PyTorch 1.13.1 CPU described in the [native recovery instructions](../../../host_ages/code/survey_physics/README.md). These dependencies are explicitly supplied; the scripts do not install packages systemwide.

The distance NPZ interfaces share keys `CID`, `zHD`, `zHEL`, `MU`, `covariance` and `precision`, with exact row ordering. A free global magnitude offset is required. Never add another per-object diagonal error, contamination probability or published systematic covariance to these total covariances. Do not combine overlapping compilations as independent datasets.

`geometry-review.json` applies to explicitly hashed source snapshots. Regenerating that review evaluates the current inference code and records new hashes; a later inference edit does not become reviewed merely because survey validation still passes. Original simulation/classification campaign inputs are required for the historical audit. The final survey check is not a certificate of a complete, independently retrained physical survey likelihood.
