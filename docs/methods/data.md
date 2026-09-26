# Data, extraction and reproducibility

The bundle is an explicit set of scientific inputs, not a mirror of the original workspace. `provenance/inputs.json` records each file's SHA-256, byte count, input group, origin and available public URLs. Original paths are provenance labels only; the code resolves operational paths relative to the repository root.

## Starting points

| Input | Scientific level | Workflows |
|---|---|---|
| Pantheon+ table and STAT+SYS covariance | Released corrected apparent magnitudes and covariance | `cosmology`, `ages` |
| DESI DR2 BAO means/covariance | Released compressed distance measurements | `cosmology`, `bao-shape` |
| Chung supplementary ZIP | Published host-age table literals | `ages` |
| DES SALT3/calibration assets and 12 fixed objectives | Released model plus locally exported calibrated-flux objectives | `des-flux` |
| DES fitted parameter/covariance arrays | Locally derived pre-BBC fit summaries with fixed cohort/folds | `des-predictors` |
| Projected shared-mode matrices | Locally derived flux-residual sufficient statistics and response vectors | `calibration` |
| RAISIN release tables/photometry | Released fits, corrected distances, covariance exports and calibrated photometry | `raisin` |
| Ten signed author precursor files | Public author photometry before positive-only product conversion | `signed-baseline` |
| Source-matched 2021 simulation FITRES pair | Archived author simulations and fitted outputs | `timing` |
| CSP DR3 archive and passbands | Original table literals, physical-filter throughput and reference spectra | `csp-lineage`, `csp-passbands` |
| Eight HST FLTs, PAM and flat references | Current calibrated detector products; explicit frozen aperture design | `hst-geometry`, `hst-repeat`, `hst-flat` |
| Fourteen dark RAWs, BPIXTAB and CCDTAB | Original unsigned detector reads and reference-only mask/gain inputs | `hst-dark-raw` |
| Two first-eight-read dark FLTs | Frozen local CALWF3 3.7.3 derivatives, retaining original dark switches | `hst-dark-calibrated` |

Locally derived DES inputs are part of this bundle. Their calculations are repeatable from that stated level. Their presence does not establish a fresh raw-data reconstruction of the original native fitter, masks or projected-mode construction. No released distance is presented as a directly measured luminosity.

The [dark workflows](hst-dark.md) rebuild detector masks from the two reference tables. Raw moments are re-extracted from all fourteen original RAWs. The calibrated comparison instead starts at the two explicitly labelled native FLT derivatives; it does not rebuild or run CALWF3. The raw cohort and calibrated sensitivity pair overlap, and their results are not independent replications.

## Verification and acquisition

```bash
uv run --frozen python research.py verify
uv run --frozen python research.py fetch --dry-run
uv run --frozen python research.py fetch --group bao
```

`verify` checks every bundled input. Each workflow also verifies its own input groups before execution. `fetch` restores missing public files only when their byte count and SHA-256 match the frozen record. Existing mismatched files are preserved and reported. A newer release or regenerated HST product is not accepted as an equivalent substitute. A hash mismatch or unavailable URL is an unresolved acquisition, not permission to change the input.

Files without a public URL must be restored from the distributed bundle or regenerated through their documented upstream procedure. This includes local DES objective/matrix exports and the original supplementary archive if its historical publisher link is unavailable. The downloader does not claim that a source-code-only checkout can recreate every bundled intermediate.

The age workflow re-extracts its two table files from the ZIP and verifies exact byte identity. The CSP workflow reads only regular archive members directly; it does not extract arbitrary archive paths onto the filesystem. FITRES and photometry parsers require declared columns, consistent row widths and unambiguous identifiers. Covariances retain the original table ordering.

## Output contracts

A run records numerical configuration, seed, lock hash, import-time source hashes, input hashes and output hashes. A failed run remains marked `failed`; a sampler can also finish while being explicitly invalid for posterior reporting. The DES predictor workflow stops before held-out scoring if its four-chain diagnostics fail. Input integrity and numerical convergence do not validate an empirical population or selection assumption.

Cross-platform floating-point order can change last digits. Compare scientific quantities and numerical gates, not compressed NPZ bytes or execution timestamps. Original file bytes are kept immutable. No source download or analysis command contacts authors, submits results, or publishes the research.
