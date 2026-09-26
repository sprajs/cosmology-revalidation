# Classification reconstruction status

> Supporting research note. Scientific findings and limitations are preserved from the recorded investigation; historical result links identify the evidence available at that time. See the [study guide](../../README.md) for current execution status.

Paused at the user's request before inference. Public nominal SNNV19 weights, normalizations and configuration are available. The DES-recommended SuperNNova commit `fcf8584b64974ef7a238eac718e01be4ed637a1d` has been checked out in an owned project directory, and an isolated CPU Python environment is installed with verified imports. Exact historical-environment equivalence and probability reproduction are **unverified**.

No classifier results, simulated acceptance fractions, or probability-calibration claims have been produced. Detailed source findings, input locations, limitations, coordination points and ordered resume commands are in `phase2/classification/HANDOFF.md`. Input/source/environment hashes are in `phase2/classification/handoff-manifest.json`.

The next required step is to preregister probability-reproduction and selection tests, then execute a bounded official-code inference check on raw calibrated flux. Published raw-data preprocessing settings must be restored from the data job, since saved training-stage cli_args alone omit redshift override and observation rejection flags.
