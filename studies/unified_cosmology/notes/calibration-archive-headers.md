# Completion of the deferred archive search

All five large archives deferred by the [earlier covariance-input search](calibration-construction-search.md) have now been checked. Their **49,051 member headers contain no `v6_9_duplicate_cid.cov`**. No new Cepheid covariance construction recipe was found in the candidate text members. The covariance discrepancy remains unresolved; this search does not justify changing the released matrix.

| Archive | Member headers | Candidate entries |
|---|---:|---:|
| SNDATA_ROOT 2024-07-03 | 14,233 | 0 |
| Pantheon+ DataRelease `c447f0f` | 3,224 | 4 |
| Pantheon+ DataRelease `7fc6805` | 2,799 | 1 |
| SNDATA_ROOT 2025-11-12 | 14,374 | 0 |
| SNDATA_ROOT 2026-04-10 | 14,421 | 0 |

The five candidates are the two known cosmology likelihood readers, `PPLUS.yml`, and the already known host-matrix pathname in both source archives. The three text files were inspected in memory. Their hashes match the previously inspected readers and configuration; `PPLUS.yml` again references the absent external file. The two host-matrix entries have 133-byte archive payloads; no numerical covariance was inferred from these pathname matches. The separately restored and validated host matrix remains the input to the existing distance-ladder analysis.

The pass read and SHA-256 hashed **6,452,240,918 compressed bytes**, sequentially with one worker. All five gzip streams passed their trailer checks, all tar header traversals completed, and no per-archive errors occurred. Input sizes, modification times and inode identities remained unchanged during inspection. No files were extracted, no inputs changed, and no network requests or cosmology calculations were made. A small synthetic control separately verified the streaming hash against direct hashing, exact-name detection, and rejection of a corrupted gzip CRC.

Full archive hashes, member-ledger hashes, candidate names and inspected text excerpts are in [the supplemental report](../results/distance_ladder/calibration-archive-headers.json). The earlier search report is retained unchanged and hash-bound. This closes its five explicit header deferrals. Nested archive payloads were counted but were not recursively decompressed; candidate content inspection was restricted to regular UTF-8 text members no larger than 2 MiB. Private production inputs remain outside the search.

To repeat the same read-only pass on the restored archives:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .work/unified-cosmology/external-probes/.modern-venv/bin/python \
  studies/unified_cosmology/code/distance_ladder/calibration_archive_headers.py
```

`--archive` selects another archive root. Complete member ledgers are written under ignored `.work/unified-cosmology/calibration-archive-headers`; the supplemental script and compact report are retained in the repository.
