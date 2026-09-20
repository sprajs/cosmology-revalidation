# Supernova acceleration dispute: research collection

Prepared 20 September 2026. This workspace gathers literature, public data, source code, provenance and access gaps for a future replication. **No cosmological fits, age regressions, posterior comparisons or tests of the authors' conclusions have been run.** File checks and catalogue counts are acquisition checks only.

The likely exchange is **Son et al. (2025) → Wiseman et al. (2026, Southampton) → Chung et al. (2026)**. It concerns whether supernova age evolution changes the evidence for **accelerating expansion**, particularly acceleration today. It does not claim that the universe never expanded. DES supernova measurements and DESI baryon acoustic oscillation measurements are separate inputs.

Start with these documents:

| Document | Purpose |
|---|---|
| [Research map](docs/research-map.md) | Papers, chronology, claims, definitions and the actual disagreements |
| [Data inventory](docs/data-inventory.md) | What is local, what each release contains, and which versions belong together |
| [Access gaps](docs/open-gaps.md) | Private inputs, unreleased products, failed retrievals and material not mirrored |
| [Popular coverage](docs/media-context.md) | The two relevant Sabine Hossenfelder videos, institutional accounts and journalism |
| [Experimental preparation](docs/replication-boundary.md) | A staged future protocol, without choosing a side or starting analysis |
| [Literature catalogue](catalog/literature.md) | Links to every indexed paper and its local PDF(s) |

The collection contains **39 PDFs** (including alternate versions and the correction), **14 pinned repository snapshots**, and about **11 GB** of local material including original archives and extracted copies. See [collection totals](catalog/collection-summary.json). The core papers, including the published Southampton response and published Yonsei reply, are in `papers/pdf/`; searchable extractions are in `papers/text/`. Supporting papers are included for data provenance and methods; they have not all received a line-by-line methodological review. The [reading ledger](catalog/reading-ledger.json) distinguishes these levels.

Data reside in `data/` and commit-pinned release snapshots under `sources/repos/`. Original repository archives remain in `sources/archives/`. Author code is preserved as evidence; it has not been executed or installed. Downloaded chains are the authors' pre-existing products, not results generated here.

The [acquisition log](catalog/acquisition.jsonl) records source URLs, UTC retrieval times, sizes, hashes and failed attempts. The [repository inventory](catalog/repository_inventory.json) records exact commits and unresolved Git LFS pointers. The [integrity report](catalog/integrity.json) covers downloads and compressed containers. The [table inventory](catalog/table_inventory.json) records observed schemas and counts. Two truncated downloads were quarantined under `sources/failed/`; usable replacements or alternative formats are separate.

Collection helpers use the Python standard library:

```bash
python3 scripts/collect.py papers 2510.13121 2601.13785 2605.21586
python3 scripts/index_collection.py
python3 scripts/inspect_tables.py
```

Existing downloads are deliberately not refreshed in place. For a later update, create a new dated collection or explicitly version the acquisition paths. Consult the inventory before assuming a whole repository is usable: some public simulation objects could not be retrieved, and exact paper-specific analysis configurations remain missing.

Bulk source material is excluded from Git by default. Keep its authors' licenses and citation requirements; this local research collection does not assign a new license to their work.
