# Historical supernova and cosmology work

**Historical context. Computation reset on 1 October 2026.**

The old repository tested released distances, flux fits, host ages, infrared
timing, survey selection and detector uncertainty. Its source, manuscript and
execution records are preserved at
[1f11997e8935a515cb5f7bb5312c91a38b6c1558](https://github.com/sprajs/reproducible/tree/1f11997e8935a515cb5f7bb5312c91a38b6c1558).
The original notes remain evidence about those executions. They do not qualify
the replacement software, and their commands are no longer active entry points.

## What should guide new experiments

| Investigation | Historical finding | Boundary to preserve |
| --- | --- | --- |
| Released-distance expansion | Flat ΛCDM recovered the released Pantheon+ matter-density inference closely; the broader BAO inequality test was inconclusive. | Shared data and supplied covariance; conditional numerical agreement is not independent observation. |
| Shared supernova–BAO–CMB fit | The manuscript reported CPL H₀ = 67.460 ± 0.565 km/s/Mpc and q₀ = −0.3310 ± 0.0614. | Higher-accuracy likelihood checks did not converge to the stated tolerances; those tight errors were provisional. |
| Host ages and corrections | Reversing exported corrections changed residual trends; imposed age templates changed expansion estimates. | No empirically identified extra age/dust correction overturned acceleration. Galaxy ages are not progenitor delays. |
| DES flux and prediction | Mean-flux fits agreed, while one local timing Hessian disagreed with profile-supported uncertainty. | A bounded object-level uncertainty defect, not a general failure of a survey release. |
| Infrared timing and covariance | Fitted times repeated initializers; systematic shift vectors did not fully reconstruct covariance. | Initializer recovery is not timing precision. Missing historical transformations remain a reconstruction gap. |
| Detector variance | One pixel contributed 82.9% of a dark-image squared-difference statistic. | Pixel and aperture weights differ; no detector-wide error multiplier follows. |

These are summaries of the old manuscript, not calculations performed during
the redesign. Exact assumptions and results are in that fixed Git snapshot.
Unsuccessful tests and unresolved identification remain part of the history.

## Inputs and reconstruction

[391 frozen input identities](../../sources/legacy-inputs.json) retain original
paths, SHA-256, sizes, source labels and known download URLs. Some are derived
input summaries. [Third-party notices](../../sources/legacy-licenses.json) retain
the existing unresolved redistribution checks.

The 34 inputs without verified public routes are preserved locally in ignored
`data/legacy/frozen-local-inputs.tar.gz`. A fresh clone must recover that exact
archive from the fixed historical commit or a separately identified archive;
the manifest records its hash. Original recovered code is also retained locally
under ignored `data/legacy/`. See [storage](../../docs/storage.md).

Rebuild one scientific question at a time after Prospector review and the needed
Irreducible capabilities. The new background example supplies no supernova,
host-age, detector, survey or CMB inference. No automatic replay of this historical
campaign exists in the new repository.
