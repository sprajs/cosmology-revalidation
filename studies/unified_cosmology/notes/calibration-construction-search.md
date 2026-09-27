# The unresolved shared-host covariance construction

The near doubling of four shared-host variances is **absent from the published Cepheid-only distance errors**. Those errors agree, within their printed rounding, with the independently reconstructed SN-free host factor. The exact step that produces the larger terms in the Pantheon+ supernova covariance remains unidentified. This is a difference between released products, not a demonstrated programming error or permission to halve a covariance.

Riess et al. Table 6, column **b**, explicitly excludes every supernova in every host. Its adjacent errors were extracted from, and visually checked against, page 30 of [arXiv:2112.04510v3](https://arxiv.org/pdf/2112.04510v3). The alternative column **a**, which excludes supernovae only in the selected host, is not used here.

| Host | Printed column-b σμ, mag | SN-free matrix σμ, mag | Released STATONLY covariance between distinct sibling SNe, mag² | Covariance / printed σμ², allowing ±0.0005 mag rounding |
|---|---:|---:|---:|---:|
| N1448 | 0.037 | 0.03665535 | 0.00265025 | 1.8846–1.9893 |
| N3147 | 0.165 | 0.16536372 | 0.05463299 | 1.9946–2.0189 |
| N5468 | 0.074 | 0.07414204 | 0.01095859 | 1.9744–2.0285 |
| N5643 | 0.052 | 0.05217550 | 0.00538959 | 1.9554–2.0321 |

All four reconstructed σμ values fall inside the printed rounding intervals. Every cross-measurement cell between distinct sibling events has the listed STATONLY value. The exact event names, original row indices and rounding calculations are in [the table audit](../results/distance_ladder/calibration-table6-review.json). The released distance table contains `CEPH_DIST`, but no `CEPH_DISTERR`; its covariance entries are complete SN-pair covariances, not separately certified Cepheid components. Consequently, agreement with the printed errors does not establish the ancestry of every term in the matrix. The [independent source review](calibration-source-review.md) discusses the single Cepheid covariance in the paper's likelihood equation and the distinction between siblings and duplicate observations of one SN.

## What the bounded search recovered

The filename inventory covers the current checkout, including ignored scientific inputs, and the entire local research archive parent. A separate content search reads only small relevant author code and acquisition/provenance files. It also inspects selected compact tar member lists without extracting them. No exact `v6_9_duplicate_cid.cov` file was found. The [supplemental archive search](calibration-archive-headers.md) closes all five deferred large-archive header searches: 49,051 headers across 6.45 GB yield no missing file or new construction recipe. Nested archives were not recursively opened; the extracted SNDATA_ROOT trees were included in the filename inventory. Complete, nontruncated public trees were inspected for the pinned Pantheon+ release, project website, Pippin and SNANA. The exact-filename GitHub code search returned only the already known `PPLUS.yml` reference. Search scope, counts, input hashes, exclusions and query responses are preserved in [the original search record](../results/distance_ladder/calibration-input-search.json) and its separately bound supplement.

The [released PPLUS.yml](https://github.com/PantheonPlusSH0ES/DataRelease/blob/c447f0fea703fcd0fff57de5000947b5ca81286b/Pantheon%2B_Data/7_PIPPIN/PPLUS.yml#L4126) specifies an external `DUP_SIGINT` covariance at scale 1 and includes it in `NOSYS`. It does not identify that file as a Cepheid covariance. Its calibrator list is an event list, not a matrix of host-distance errors.

The recovered [May 2022 Pippin wrapper](https://github.com/dessn/Pippin/blob/b457d7e870090b926e8488e0a3fa56c7e1a62b3e/pippin/create_cov.py) forwards `EXTRA_COVS` and `CALIBRATORS` to the covariance builder. The [June 2022 SNANA builder](https://github.com/RickKessler/SNANA/blob/cb6011e7d95070bede2ec2355c2a8a216b3517cc/util/create_covariance.py#L630) reads external `MU_COV` values indexed by ordered `CID+IDSURVEY` pairs. An isolated synthetic execution confirms that repeated ordered entries overwrite rather than add; the reader neither symmetrizes those entries nor constructs a host covariance. Its calibrator handling suppresses peculiar-velocity/redshift systematic terms. These generic files contain no Cepheid host-to-lightcurve embedding recipe, and their dates do not certify the production versions used for the release. The current builder has later changes and is not substituted for the historical source.

The search also found an author teaching notebook with a simplified two-rung distance ladder. It explicitly omits SNe and does not supply the missing production construction. No private directories or author communications were accessed.

The remaining inputs are concrete: the exact external duplicate-scatter table, the particular SN-free Cepheid mean/covariance used in the release, and the code that maps host covariance onto the repeated SN measurement rows and adds it to the SN covariance. Without them, the search does not distinguish a different covariance version, an additional shared term, or an implementation problem. No new cosmology was fitted and no released data, covariance or active target was changed.

## Reproduction

First restore the pinned distance-ladder and calibration-interface inputs using their acquisition commands. Then run:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/distance_ladder/calibration_input_search.py --acquire

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/distance_ladder/calibration_table6_review.py
```

The search defaults to `Path.home()/.local/share/cosmology-revalidation/archive`; `--archive` selects another restored archive. Optional `--search-public` refreshes the bounded public code queries using `gh`, retaining changed earlier response snapshots. A later filesystem or search index can change search counts; source downloads themselves are hash-pinned. The table comparison uses fixed released input bytes and does not depend on those changing search counts.
