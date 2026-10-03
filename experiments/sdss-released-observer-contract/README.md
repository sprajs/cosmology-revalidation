# Released SDSS source lineage

This packet retains the released SDSS row identities, their HEAD/update/override
pointers and one original photometry slice. Its dedicated reader performs
serialization and joins only. The physical experiment stays blocked: event,
frame, reduction, exposure, calibration, selection and dependence are unresolved.
No engine request, physical calculation, fit or inference is executed here.

The immutable [candidate](candidate.json) is Prospector's original
`designs/candidate-sdss-released-observer-contract.json`, revision
`4b9a9f4fc93f98b7d1ea1bffe105338162220897`, SHA256
`954d3e2369f21099dbe7fa298e7ee040bab0c5c90ca101f192b79292ff31f5fe`.
This packet narrows it to source lineage and omits its conditional observer and
passband calculation. That calculation has its own
[historical controller](../sn-observer-passband/README.md).
The later source-gate contract is referenced at Prospector
`9abdf9e6856d78b84cd05d1a4fdf7c4365d2b1d9`,
`register/contracts/standard-model-source-gates-v1.json`, SHA256
`ca0ed970215f151d372da450dec7d1f0f24bd307b09d53e6df98c30a30122597`.
It supplies versioned gates without replacing the original candidate.

## Frozen slice

[input-contract.json](input-contract.json) pins nine original assets, schemas,
byte lengths/hashes, bounds and unresolved semantics.
[experiment.json](experiment.json) records source URLs and physical blockers.
The reader parses the same bounded byte objects whose hashes were admitted.
Original assets and decoded outputs stay local; redistribution permission has
not been verified. The reader does not download inputs or extract archives.

All 321 IDSURVEY1 final-table rows retain their original order, HEAD record,
update record and four override rows, including host/group flags, raw literals,
padding and binary record hashes. HEAD pointers account for 35,181 PHOT points.
The chosen CID 6057 is final row 690, HEAD row 2864 and update row 737 (zero based).
Its one-based PHOT interval 313798–313912 contains 115 records. The HEAD-named
`SDSS_allCandidates+BOSS_PHOT.FITS.gz` and renamed
`JLA2014_SDSS_DS17_PHOT.FITS.gz` retain distinct source identities. Their selected
records agree, including 15 signed negative FLUXCAL values. Neither signed flux
nor a source flag supplies a detector nondetection law.

Ten RA rows fail the unchanged printed half-unit comparison: CIDs 1241, 19953,
19940, 12856, 15461, 19913, 20376, 3452, 2789 and 18617. They remain in the ledger.
Differences use an explicit 128-digit Decimal context with rounding traps;
source binary32/binary64 fields and decimal literals remain unchanged. No frame
correction or tolerance slack is introduced.

FITS TUNIT cards are present but blank; EXPTIME is absent. Physical flux/error
units, clock standard, exposure duration/area, luminosity/template, detector law
and cross-epoch covariance remain null. Earlier HEAD, updated redshift/velocity
records and final fitted summaries have distinct ancestry. Unique CID strings
do not establish physical-event independence. The final covariance row order
is preserved, but this reader does not read or subset a covariance.

## Dedicated invocation

From clean committed source, using a verified cache of the nine exact basenames:

```sh
python experiments/sdss-released-observer-contract/reader.py \
  --input-root /path/to/verified/source-cache --attempt source-lineage-v1
```

The reader writes a fresh bounded component under ignored
`results/sdss-released-observer-contract/`. Existing attempts, symlink parents
and escaping paths refuse. It verifies all five packet files against committed
Git blobs and records source revision, Python executable/version, packet and
nine original-file identities before and after execution. Drift revokes
structural completion while retaining checked lineage. Output files become
read-only; permissions do not replace an immutable archive.

Compressed input admission is 96 MiB, HEAD inflation 4 MiB and each PHOT inflation
128 MiB. Headers, text rows/lines and chosen records are bounded. Gzip streams
are drained in bounded chunks through CRC/trailer/EOF. Output is capped at 4 MiB,
with a 64 KiB terminal receipt. Paths and diagnostics retain full-text hashes
with bounded displays; receipt overflow preserves a bounded failure summary and
checked lineage. No physical gate is promoted by structural completion.

## Findings

2026-10-03: the retained ignored prototype v4 regenerated all 321 joins and 115
records, preserved ten RA failures and 15 signed negative fluxes, and matched
the earlier source ledger. Its 15 structural fixtures passed. The selected
uncompressed slice SHA256 was
`d864e65796609a88a6eb1e363bb6d38efc6f37acee4d0508f0591bf769282f41`.
This prototype has its own source/contract identity; its result does not execute
or qualify the present packet.

The present reader's 26 structural/adversarial fixtures passed, including gzip
integrity, exact schemas, raw string storage, unchanged coordinate bounds,
fresh paths, input/runtime drift and terminal-failure preservation. A fresh
committed source attempt and deliberate final review are required before this
packet's own structural finding is recorded. Physical exposure and observational
qualification require separately reviewed source and dependence contracts.
