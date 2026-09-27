# Unite: source availability and integration limit

The [versioned Unite paper](https://arxiv.org/abs/2609.05053v2), revised 9 September 2026, describes 2,884 likely supernovae obtained by combining and reanalyzing Pantheon+ and DES-SN5YR, including Dovekie updates. Its data-availability statement says the distances, likelihood, host masses and photometry will become public on acceptance. The revision comments explicitly record removal of inactive data links.

The bounded retrieval on 27 September recovered both versions' HTML and PDF and the v2 source archive. The earlier HTML and commented v2 TeX retain two specific release destinations:

- [RyanCamo/UNITE](https://github.com/RyanCamo/UNITE), also written with different capitalization;
- [jasonlee17/SN-Unite_Host_Galaxy_Masses](https://github.com/jasonlee17/SN-Unite_Host_Galaxy_Masses).

Both exact public GitHub repository API requests returned **404**. This does not distinguish private repositories from absent ones. No author-linked Zenodo deposit was found in the retrieved sources. The v2 source archive has **32 entries** (including one directory), comprising manuscript, bibliography, styles and figures; it does not contain the ordered distance vector or covariance. The v1 source request returned 406, while its HTML and PDF were available. PDF extraction warnings and a wrapped-URL artifact are recorded explicitly.

These observations do not establish that the data are unavailable everywhere. They establish that this source-linked acquisition did not recover a usable joint likelihood. We therefore have not incorporated Unite into the measurement. Its published parameter posteriors and plotted Hubble diagram have not been used as substitute observations.

When the release becomes accessible, an integration needs the ordered identifiers, redshifts and corrected distances; matching statistical and systematic covariance or exact likelihood; normalization and anchoring conventions; and release/selection provenance. Its extensive overlap with Pantheon+ and Dovekie requires a **replacement supernova analysis**, not multiplication as an additional independent sample. No cross-sample independence is inferred from a different catalogue name.

The [machine-readable receipt](../results/external_probes/unite-sources.json) pins each response, downloaded source, extracted PDF text and acquisition code. Its conclusions are scoped to the recorded retrieval. All model-call counts are zero. Downloaded files remain ignored.

The [acquisition script](../code/external_probes/unite_sources.py) requires fresh cache and report paths and refuses to overwrite prior evidence. Reproduce the check from the repository root, choosing a new suffix if the paths already exist:

```bash
.work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/external_probes/unite_sources.py \
  --cache .work/unified-cosmology/unite/recheck-v1 \
  --output .work/unified-cosmology/unite/recheck-v1.json
```

A later successful release retrieval requires an actual data/likelihood audit; this inventory does not automatically certify newly returned assets.
