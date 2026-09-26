#!/usr/bin/env python3
"""Inspect research sources or prepare an isolated workspace for original studies.

Preparation never runs experiments. Native builds and missing upstream inputs
remain explicit prerequisites; a source inventory is not a reproduction claim.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
HISTORICAL_ROOT = "/home/szymon/Documents/ChatGPT/supernova"


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def read(name):
    return json.loads((ROOT / "provenance" / name).read_text())


def within(root, relative):
    path = Path(relative)
    if path.is_absolute() or ".." in path.parts or not path.parts:
        raise ValueError(f"Expected a relative path: {relative}")
    destination = root / path
    if not destination.resolve().is_relative_to(root.resolve()):
        raise ValueError(f"Path leaves workspace: {relative}")
    return destination


def workspace(name):
    if Path(name).name != name or name in {".", "..", ""}:
        raise ValueError("Workspace name must be one directory name")
    return within(ROOT, f".work/{name}")


def retained():
    return [r for r in read("studies.json")["files"] if "path" in r]


def verify():
    unique = {}
    for record in retained():
        path = within(ROOT, record["path"])
        if digest(path) != record["sha256"]:
            raise ValueError(f"Source identity changed: {record['path']}")
        unique[record["path"]] = record["sha256"]
    return {"source_files_verified": len(unique), "historical_locations": len(retained())}


def prepare(name):
    verify()
    out = workspace(name)
    out.mkdir(parents=True, exist_ok=False)
    records = []
    for record in retained():
        # Notes are reader-facing scientific accounts. Their old paths and hashes
        # remain historical evidence, not current execution inputs.
        if record["disposition"] == "retained-scientific-note":
            continue
        source = within(ROOT, record["path"])
        target = within(out, record["original_path"])
        target.parent.mkdir(parents=True, exist_ok=True)
        content = source.read_bytes()
        relocated = content.replace(HISTORICAL_ROOT.encode(), str(out).encode())
        with target.open("xb") as stream:
            stream.write(relocated)
        target.chmod(source.stat().st_mode & 0o777)
        records.append({"source": record["path"], "target": record["original_path"],
                        "source_sha256": record["sha256"], "workspace_sha256": digest(target),
                        "absolute_root_relocated": relocated != content})
    result = {"scope": "Source preparation only; no scientific result or native build validated.",
              "source_manifest_sha256": digest(ROOT / "provenance/studies.json"),
              "files": records}
    (out / "preparation.json").write_text(json.dumps(result, indent=2) + "\n")
    return {"workspace": str(out), "source_locations_prepared": len(records),
            "absolute_paths_relocated": sum(r["absolute_root_relocated"] for r in records)}


def inputs():
    result = {r["path"]: r for r in read("study-inputs.json")["files"]}
    for r in read("inputs.json")["files"]:
        # Prefer the verified current input bundle when an original path occurs
        # in both manifests. Never substitute a similarly named release.
        result[r["source_path"]] = dict(r, local_path=r["path"], path=r["source_path"])
    return result


def fetch(name, path, dry_run):
    out = workspace(name)
    if not (out / "preparation.json").is_file():
        raise ValueError("Prepare the workspace first")
    item = inputs().get(path)
    if item is None:
        raise ValueError("No frozen input identity for that path; use its study acquisition procedure")
    target = within(out, path)
    if target.exists():
        if digest(target) != item["sha256"]:
            raise ValueError("Existing input differs from the recorded identity")
        return {"path": path, "status": "already-verified"}
    local = within(ROOT, item["local_path"]) if "local_path" in item else None
    url = item.get("url") or next(iter(item.get("urls", [])), None)
    if not url:
        url = item.get("retrieval", {}).get("url")
    if dry_run:
        return {"path": path, "local_copy_available": bool(local and local.is_file()), "url": url,
                "sha256": item["sha256"]}
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=target.parent, prefix=target.name + ".", suffix=".partial", delete=False) as handle:
        temporary = Path(handle.name)
    try:
        if local and local.is_file():
            # Copy, never symlink: an old program must not modify frozen inputs.
            with local.open("rb") as source, temporary.open("wb") as destination:
                shutil.copyfileobj(source, destination)
        elif url and url.startswith("https://"):
            request = urllib.request.Request(url, headers={"User-Agent": "Cosmology-Revalidation/1.0"})
            with urllib.request.urlopen(request, timeout=60) as source, temporary.open("wb") as destination:
                shutil.copyfileobj(source, destination)
        else:
            raise ValueError("Input requires upstream preparation; no verified local copy or download URL")
        if digest(temporary) != item["sha256"]:
            raise ValueError("Input hash mismatch; original identity is required")
        # Atomically publish verified bytes without replacing another output.
        target.hardlink_to(temporary)
    finally:
        temporary.unlink(missing_ok=True)
    return {"path": path, "status": "restored", "sha256": digest(target)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    listing = sub.add_parser("list", help="List source paths and their historical execution locations")
    listing.add_argument("--topic")
    sub.add_parser("verify", help="Verify all retained source, specification and note identities")
    staging = sub.add_parser("prepare", help="Copy sources into .work/NAME; do not run them")
    staging.add_argument("--name", required=True)
    fetching = sub.add_parser("fetch", help="Restore one hash-verified input into a prepared workspace")
    fetching.add_argument("--name", required=True)
    fetching.add_argument("--path", required=True, help="Historical input path from a study specification")
    fetching.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    try:
        if args.command == "verify":
            answer = verify()
        elif args.command == "prepare":
            answer = prepare(args.name)
        elif args.command == "fetch":
            answer = fetch(args.name, args.path, args.dry_run)
        else:
            answer = [r for r in retained() if not args.topic or r["path"].startswith(f"studies/{args.topic}/")]
        print(json.dumps(answer, indent=2))
    except (OSError, ValueError) as exc:
        parser.exit(1, f"{exc}\n")


if __name__ == "__main__":
    main()
