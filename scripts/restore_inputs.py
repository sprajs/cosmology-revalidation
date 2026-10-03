#!/usr/bin/env python3
"""Restore frozen inputs without installing the scientific Python environment.

Examples: --list; --dry-run; --group bao; --verify-only.
The default selects all 391 historical inputs. Downloads do not run experiments.
Some frozen inputs are derived summaries, not original observations.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tarfile
import tempfile
import urllib.request

from metadata_source import MetadataSourceError, read_document

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def destination(root, relative):
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError(f"Unsafe relative path: {relative}")
    target = root / path
    if not target.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Destination leaves root: {relative}")
    return target


def copy_bounded(source, target, expected):
    count = 0
    while block := source.read(1024 * 1024):
        count += len(block)
        if count > expected:
            raise ValueError("Response exceeds frozen byte count")
        target.write(block)
    if count != expected:
        raise ValueError(f"Expected {expected} bytes; received {count}")


def restore(row, root, bundle, bundle_record):
    target = destination(root, row["path"])
    if target.exists():
        if not target.is_file() or target.stat().st_size != row["bytes"] or sha(target) != row["sha256"]:
            raise ValueError("Existing file differs; preserved without replacement")
        return "already-verified"
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix=target.name + ".", suffix=".partial", dir=target.parent)
    os.close(fd)
    temporary = Path(name)
    try:
        if row["restore"] == "frozen-bundle":
            if sha(bundle) != bundle_record["sha256"]:
                raise ValueError("Frozen local-input archive identity changed")
            with tarfile.open(bundle, "r:gz") as archive:
                member = archive.getmember(row["path"])
                if not member.isfile() or member.size != row["bytes"]:
                    raise ValueError("Archive member differs from frozen record")
                with archive.extractfile(member) as source, temporary.open("wb") as stream:
                    copy_bounded(source, stream, row["bytes"])
        else:
            errors = []
            for url in row["urls"]:
                try:
                    if not url.startswith("https://"):
                        raise ValueError("Only HTTPS download URLs are supported")
                    request = urllib.request.Request(url, headers={"User-Agent": "cosmology-revalidation/0.1"})
                    with urllib.request.urlopen(request, timeout=90) as source, temporary.open("wb") as stream:
                        copy_bounded(source, stream, row["bytes"])
                    if sha(temporary) != row["sha256"]:
                        raise ValueError("Downloaded hash differs from frozen input")
                    break
                except (OSError, ValueError) as exc:
                    errors.append(f"{url}: {exc}")
            else:
                raise ValueError("; ".join(errors) or "No verified download URL")
        if sha(temporary) != row["sha256"]:
            raise ValueError("Restored hash differs from frozen input")
        # Publish without overwriting a file created concurrently.
        target.hardlink_to(temporary)
        return "restored"
    finally:
        temporary.unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", help="Use --list to see groups; 'all' selects every recorded download")
    parser.add_argument("--path", action="append", help="Select exact manifest paths; repeat to select several")
    parser.add_argument("--destination", type=Path, default=ROOT, help="Restoration root, default repository")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    parser.add_argument("--list", action="store_true")
    args = parser.parse_args()
    try:
        manifest, _, source_identity = read_document('legacy-inputs', ROOT)
    except MetadataSourceError as exc:
        print(json.dumps({'source_admission': 'refused', 'document_id': exc.document,
                          'stage': exc.stage, 'error': str(exc)[:256],
                          'manifest_source': exc.identities}, sort_keys=True))
        raise SystemExit(1)
    if args.list:
        groups = {}
        for row in manifest["files"]:
            info = groups.setdefault(row["group"], {"files": 0, "bytes": 0})
            info["files"] += 1
            info["bytes"] += row["bytes"]
        print(json.dumps(groups, indent=2))
        return
    group = args.group or (None if args.path else "all")
    rows = [r for r in manifest["files"] if (group in (None, "all") or r["group"] == group)
            and (not args.path or r["path"] in args.path)]
    if not rows or (args.path and set(args.path) != {r["path"] for r in rows}):
        parser.error("No input matched, or one of the selected paths/groups is unknown")
    failures = []
    statuses = {}
    for row in rows:
        try:
            target = destination(args.destination, row["path"])
            if args.dry_run:
                print(json.dumps({"path": row["path"], "bytes": row["bytes"], "restore": row["restore"],
                                  "urls": row.get("urls", []), "present": target.is_file()}))
                continue
            if args.verify_only:
                if not target.is_file() or target.stat().st_size != row["bytes"] or sha(target) != row["sha256"]:
                    raise ValueError("Input missing or differs from frozen identity")
                status = "verified"
            else:
                status = restore(row, args.destination, ROOT / manifest["frozen_bundle"]["path"], manifest["frozen_bundle"])
            statuses[status] = statuses.get(status, 0) + 1
            print(f"{status}: {row['path']}", flush=True)
        except (OSError, ValueError, KeyError, tarfile.TarError) as exc:
            failures.append({"path": row["path"], "error": str(exc)})
    print(json.dumps({"selected": len(rows), "statuses": statuses, "failures": failures,
                      "manifest_source": source_identity}, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
