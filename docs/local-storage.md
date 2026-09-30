# Restoring data after local cleanup

The October 1, 2026 cleanup removes bulk downloads, simulation/inference outputs,
native build trees, environments and caches. Scientific notes, current source,
designs, compact published results and provenance remain in Git. Removed analysis
arrays require fresh execution; downloading inputs does not reproduce or
scientifically qualify an analysis. The separate Irreducible project will provide
the future computational software; this repository retains the historical research.

## Restore frozen inputs

Python 3.12 or later is sufficient; the downloader uses only the standard library:

```bash
python tools/restore_data.py --list
python tools/restore_data.py --group core --dry-run
python tools/restore_data.py --group core
python tools/restore_data.py --group core --verify-only
```

`core` restores all 391 inputs used by `research.py`. The expanded
[download manifest](../provenance/downloads.json) records 3,478 inputs with their
actual URLs, byte counts and SHA-256 identities. Use `--path` for individual files,
`--group host-transport` for a collection, or `--destination /tmp/input-check` for
an isolated restoration. Downloads go below ignored `data/` and `.work/`.

Different bytes are refused, changed existing files are preserved, responses are
size-bounded and verified files are published atomically. URLs may disappear or
serve changed releases; those are acquisition failures. Representative URLs were
downloaded again during cleanup, not every URL. Every recorded identity was checked
against its existing local file before removal.

Thirty-four frozen inputs have no verified public URL. These include local DES
objectives, projected calibration statistics, aperture designs and two calibrated
dark images. Their exact bytes are retained in the approximately 22 MiB
[frozen input archive](../provenance/frozen-local-inputs.tar.gz), used automatically
by the downloader. This small reproducibility exception prevents substituting newly
generated intermediates for the historical inputs behind published records.

Restore the scientific environment only when needed:

```bash
uv sync --frozen --python 3.12.14 --extra predictors
.venv/bin/python research.py verify
```

## Additional acquisitions and native software

The expanded manifest restores individual downloads, including archives. Use the
existing study scripts for extraction, pinned software, builds and preparation.
These procedures can be expensive and produce new timestamps or numerical bytes.

| Collection | Acquisition and preparation |
|---|---|
| Original studies | `studies/manage.py prepare --name NAME` and its `fetch` command; [study guide](../studies/README.md) |
| Literature | `--group literature`; regenerate PDF text with `pdftotext -layout` |
| SDSS and MaNGA | [Galaxy validation](../studies/host_ages/code/galaxy_validation/README.md) |
| Roman, DES deep-field and HELP images | `studies/host_ages/code/host_transport/acquire.py`, `deep.py`, `infrared.py`; [host transport](../studies/host_ages/notes/host-transport-results.md) |
| Infrared catalogues and beams | [Catalogue observations](../studies/host_ages/code/infrared_photometry/README.md), [beam measurements](../studies/host_ages/code/infrared_resolution/README.md) |
| Stellar libraries and foreground maps | `studies/host_ages/code/physical_ages/acquire.py`; [physical ages](../studies/host_ages/code/physical_ages/README.md) |
| Native survey assets | `studies/host_ages/code/survey_physics/recover.py --extract`; [survey instructions](../studies/host_ages/notes/survey-physics-results.md), [correction support](../studies/host_ages/notes/survey-bbc-support-results.md) |
| Dovekie, DES3YR, classifier, SNDATA | `studies/unified_cosmology/code/survey_selection/acquire.py --bundle`, `acquire_followup.py`; [instructions](../studies/unified_cosmology/code/survey_selection/README.md) |
| Planck, ACT, SPT and author chains | `studies/unified_cosmology/code/external_probes/acquire.py`, `restore_modern.py`, `author_config.py`, `author_acquire.py`; [pinned environments](../studies/unified_cosmology/code/external_probes/README.md) |
| FrankenBlast, YSE and OzDES | `studies/unified_cosmology/code/host_likelihood/acquire.py --posterior-archive`, `ozdes_acquire.py`; [instructions](../studies/unified_cosmology/code/host_likelihood/README.md) |
| DESI calibrated host spectra | [Calibrated hosts](../studies/unified_cosmology/code/calibrated_hosts/README.md) |
| High-resolution stellar response | [Calibrated host physics](../studies/unified_cosmology/code/calibrated_host_physics/README.md) |
| Historical infrared assets | `studies/infrared/code/restore_bayesn_heldout.py`, acquisition scripts in `studies/infrared/code/experiments/raisin_timing_assets/`; [infrared guide](../studies/infrared/README.md) |
| Astronomy software catalogue | `docs/astronomy-software-survey-2026-09-30/collect_ascl_index.py`, `verify_sources.py`; [report](astronomy-software-survey-2026-09-30/REPORT.md). A live registry refresh differs from the historical snapshot. |

## Previously ignored sources and validation

[Recovered source identities](../provenance/recovered-local-code.json) map 1,011
working-tree source paths to 464 unique members of the approximately 4 MiB
[source snapshot archive](../provenance/recovered-local-code.tar.gz). This preserves
diagnostic scripts, old variants, unexecuted drafts, curation tools and native source
snapshots. Clean upstream Git clones had no local changes; their acquisition scripts
retain release pins. These archived sources are historical records, not newly
validated production entrypoints. Inspect them inside an ignored workspace:

```bash
mkdir -p .work/source-recovery
tar -xzf provenance/recovered-local-code.tar.gz -C .work/source-recovery
python validation/restore_data_checks.py
python studies/manage.py verify
```

The [cleanup record](../provenance/storage-cleanup-2026-10-01.json) records storage
scope and checks. Historical GitHub history is preserved. If using a shallow clone,
`git fetch --unshallow origin` restores the multi-gigabyte history locally.
