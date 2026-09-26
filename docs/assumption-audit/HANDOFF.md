# Handoff: cross-model dust, age and uncertainty audit

Saved 21 September 2026, Europe/London. Workspace: `/home/szymon/Documents/ChatGPT/supernova`.

## Start here

Read [the synthesis](README.md), then the relevant detailed report: [host ages/Pantheon](pantheon.md), [DES/Dovekie](des.md), [foreground/covariance](foreground-and-covariance.md), or [BAO/CMB](bao.md). This task audited assumptions outside another chat's SALT/SNM correction work. It produced evidence and reproducible diagnostics, not a replacement cosmological analysis.

The completed record is bound by [final-manifest.json](../../runs/assumption_audit/final-manifest.json) and [final-verification.json](../../runs/assumption_audit/final-verification.json). Check their timestamps and hashes before treating later workspace edits as part of this audit. The integrity check is not an independent proof of the astrophysical claims.

A local checkpoint of the small deliverables and preserved code is saved as `runs/assumption_audit/checkpoints/assumption-audit-2026-09-21.tar.gz`, with a SHA-256 sidecar. It excludes bulk maps, original archives and posterior chains; keep those local files too. These audit additions have not been committed: do not discard untracked files with `git clean`.

## Verified findings to retain

1. **Original R19 host-age conversion is inconsistent with the fitted FSPS star-formation history.** MC-Age integrates `K+s*u`, but its spectral generator uses `K*(1+s*u)`. Exact pinned source, archived columns, independent analytic integration and native-normalization quadrature support this. Original global archive: 103 hosts, 105,059,947 draws. On valid draws, the mean/median change across host medians is +0.88175/+0.41962 Gyr; 29 hosts change by more than 1 Gyr. These deterministic transformations are not validated new stellar ages. Revised C25 executable code and posteriors remain unavailable: do not claim this defect is verified in C25.
2. **Original retained chains also contain invalid posterior support.** 1,288,021 draws (1.226%) violate the stated prior, affecting 86 hosts; 5,186 have negative tau. Excluding invalid draws alone moves medians by at most 0.07256 Gyr. CID20048 has 1,019,947 rows; the other 102 have 1,020,000. This 53-row shortfall is recorded, not repaired. C25 has 102 R19 hosts and omits original CID15459. Do not confuse those two different denominator issues.
3. **The host dust model is conditional and differs from the printed equation.** Both birth-cloud and diffuse attenuation act on young stars in executed FSPS; fixed geometry, informative dust/metallicity priors, cosmological age support and a diagonal five-band likelihood restrict the age uncertainty. Foreground-corrected errors remain unchanged. Original joint draws explicitly show dust/metallicity/age posterior correlations. Host/SN foreground error cross-covariance has not been measured.
4. **Ancillary covariance products require care.** Pantheon+ MWEBV/MWCOLORLAW group files do not share the final STATONLY baseline. DES CAL_SALT3 already includes CALSPEC. DES VPEC's difference matrix is not an independent positive-semidefinite covariance, although the full released covariance is positive definite; an illustrative completion has small cosmological effect. Main local fits use the released totals correctly.
5. **Spatial foreground variation is not a global scale.** Actual CSFD-minus-SFD samples leave 99.1% of their centered variation unexplained by a constant plus one global reddening scale (98.4% on the stricter reliable subset). This is not a missing percentage of SN variance. No conversion to standardized distances using a constant extinction coefficient is justified.
6. **Shared correction uncertainty matters.** In the fixed C14-template test, correct shared-slope propagation gives q0 uncertainty 0.0967, versus 0.0743 for an incorrect diagonal replacement. Existing shared-slope code is correct under its fixed-template/independent-prior assumptions. Importance-weight ESS is not chain ESS; chi-square values across different covariances are not evidence comparisons without normalization.
7. **BAO implementation checks pass; its extrapolation has assumptions.** Observable order, Ly-alpha ordering, within-bin covariance and likelihood agree independently. BAO-only q0 changes from about +0.035 to +0.30 when wa's lower bound moves from −3 to −10. Matched long SN+BAO runs agree (−0.4360 versus −0.4350); published-scale cross-bin systematics shift joint q0 by less than 0.0004. Common dust-field cross-probe errors remain unmeasured. CMB references include foreground nuisances and joint ACT–Planck lensing covariance; no CMB-map refit was done.

The rare original age-quadrature error is separate from the SFH mismatch: the worst valid draw differs by 0.53236 Gyr, and executing the original unsplit quadrature reproduces that archived outlier. Typical old-column reconstruction is much tighter; see the per-host error columns rather than claiming every archived age was reproduced exactly.

Prior-valid does not mean converged or correctly weighted: invalid-draw filtering cannot repair a posterior. The independent final review at `runs/assumption_audit/des/r19-prior-peer-review.json` checked aggregates, exact-source boundaries/examples, the raw short member and independent quadrature; it did not independently reparse all 105 million rows.

## Files and provenance

- This task's reports: `docs/assumption-audit/`.
- This task's scripts: `scripts/assumption_audit/`.
- Numerical outputs/manifests/reviews: `runs/assumption_audit/`.
- Full original global archive: `runs/assumption_audit/pantheon/sources/rose2019-campbellG.tar.gz`, MD5 `3420a1f080231d3ba01b7e126cbf1039`, 1,998,731,699 bytes; original source is Zenodo record 3875482. Preserve the archive and validation-circle archive for reproduction.
- Authoritative all-draw results: `pantheon/r19-global-validated-age-{audit.csv,summary.json,manifest.json}` and `pantheon/r19-global-prior-violations.json`, under the output directory. The first stream lives under `pantheon/preliminary/` and contains diagnostic nonfinite results; it is superseded, not silently sanitized.
- Source and archive acquisition records are under `pantheon/sources/`; exact historical FSPS source copies are there too. MC-Age's inspected source is `sources/repos/benjaminrose__mc-age/calculateAge.py`, revision `92713be96a89da991fe53bffcc596a5c0942fc37`.
- Foreground maps: `data/dust/csfd-v2/`; source paper/download record: `sources/updates/2026-09-21-assumptions/` and the foreground manifest. Downloaded maps, source archives and chains are local/ignored; a Git-only copy is insufficient to reproduce everything.
- The age figure is `runs/assumption_audit/pantheon/r19-age-audit.png`; it explicitly conditions on valid draws and is not the revised C25 catalogue.
- DES and BAO producing-code snapshots and source/input/output digests are retained with their manifests. The final manifest inventories reports, scripts and small result files.

## Concurrent work: do not revert or claim it

Other chats modified top-level README/docs, `scripts/cosmology/core.py`, `scripts/standardization/dust_identifiability.py`, and added `phase2/`, physics-audit and salt-dust-audit work. Those are not this task's changes. No commit, push, publication, or author contact was performed here.

A concurrent core edit changed only qbins branches/helper. The original core SHA-256 is `18cc923c532e96791492ae2f72c4c8a8f8fc9cdfface9dfdff8f29cd4178eecc`; the inspected later version is `83f852a3d51e851ed9030ee8f9e5b6c74fae5893bb5c12891a3fe1a8eaa1bc8c`. Both exact files are preserved under `runs/assumption_audit/bao/provenance/`. `bao_verify.py` verifies whole-module AST equivalence after removing precisely those qbins additions. Two long runs originally hashed core at completion, so their manifest contains the later hash: retain that timing limitation; do not rewrite history. Later sampler code captures producing bytes at startup.

## Recheck the saved record

Run from the workspace using its existing `.venv`. These checks do not repeat the expensive scientific chains:

```bash
OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/assumption_audit/des_validate.py
.venv/bin/python scripts/assumption_audit/bao_verify.py
.venv/bin/python scripts/assumption_audit/freeze_audit.py
```

The last command deliberately regenerates the final inventory, so preserve a previous record first if comparing historical state. The relevant detailed reports name the scientific reproduction commands. In particular `pantheon_age_prior_validation.py` streams all 105 million archived draws and is not a quick startup check. `pantheon_age_summary.py` regenerates the crosswalk/figure from the completed results. Do not rerun long chains merely to recover context.

`map_dust.py` used `astropy-healpix==1.1.2` in `/tmp/supernova-assumption-deps`, with `PYTHONPATH` pointing there; it did not change the project environment or lockfile. If that temporary dependency is gone, reinstall it to an isolated target before rerunning the map diagnostic. Large source files are referenced by hashes, not embedded in the final small-artifact inventory.

## Next work, if requested

1. Recover the exact revised C25 age-conversion implementation and joint SFH posterior arrays. Compare with the demonstrated mismatch before applying any transformation to C25 or regressing revised SN residuals.
2. Refit host SEDs under alternative dust laws/geometry and a shared foreground nuisance; retain joint metallicity, dust and SFH posteriors. Propagate through matched SN corrections and selection before claiming a new age–luminosity relation or cosmological result.
3. Resolve the Pantheon grouping-file statistical baseline and DES VPEC generation/full-precision provenance. Test common PS1/Foundation passband uncertainty in the lower-level correction task.
4. Estimate common foreground responses of SN, host and BAO estimators before introducing a cross-probe covariance. An arbitrary correlation coefficient is not a measured uncertainty.

These are identified future scientific dependencies, not unfinished archive processing. All original-global draws have been audited. No pending long sampling run is needed to use the saved findings.
