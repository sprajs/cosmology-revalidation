The original NIR stage stopped before any NIR native process. Its strict null check correctly detected a real serialization mutation: Astropy rewrote 712 spaces as NUL bytes across SUBSURVEY, SNID, IAUC, SIM_MODEL_NAME and SIM_TYPE_NAME. All 103 non-string columns were unchanged. Scalar character representations concealed the padding difference; the full table/column comparison did not.

The additive adapter copies the source file bytes and patches only PEAKMJD after verifying the actual BINTABLE format, big-endian R4 dtype, row stride, data offset, CID coverage, absence of scaling and checksum cards, and original decoded values. The current table has data offset 31,680, row stride 524 and PEAKMJD offset 316. The allowed target comprises exactly 32 bytes across eight rows. Tests show the null file is byte-for-byte identical to its source; the fitted-peak fixture changes 10 of those 32 bytes and no others. All nonpeak columns also match when reread, including original string padding. Native PHOT remains the same symlinked file.

Both failed derived files remain in place, hashed in preserved-failure-manifest.json. No paths were moved. The new wrapper calls the unchanged frozen NIR routine with only adapt_head replaced and uses new derived-start-only-bytes directories. It does not change any fitting, support, convergence, pure-peak or full-null comparison. It carries 26.785140010 seconds into the same 120-second native cap. Preparation and fixture tests used zero native calls.

Release fields are `protocol_sha256`, `freeze_sha256` and `stages: ["nir"]`. Execute only after root release:

```sh
phase2/env-official/bin/python phase2/pte/restricted-peak-engineering/start-only-hook/head-byte-adapter/resume_nir.py --release /absolute/path/to/root-release.json
```

Protocol SHA: 9a16e0bc9ab437c94a094aefc28d3cda59b2c7eeabb49d986f6e97f06a6ce316. Freeze SHA: edfaae20f174bdb10545684e575824360bcf55f8549bbc5ac8130d7f209125f2. No native execution occurred during this preparation.
