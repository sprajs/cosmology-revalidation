#!/usr/bin/env python3
"""Restore the pinned method descriptions; numerical inputs use prior acquisition."""
from pathlib import Path
import hashlib,json,urllib.request
ROOT=Path(__file__).resolve().parents[4]
def sha(path):
    with path.open("rb") as stream:return hashlib.file_digest(stream,"sha256").hexdigest()
def main():
    record=json.loads((ROOT/"studies/host_ages/results/nebular_dust/sources.json").read_text())
    for source in record["files"]:
        path=ROOT/source["path"]
        if not path.exists():
            path.parent.mkdir(parents=True,exist_ok=True)
            request=urllib.request.Request(source["url"],headers={"User-Agent":"cosmology-revalidation public measurement-method audit"})
            data=urllib.request.urlopen(request,timeout=60).read()
            if hashlib.sha256(data).hexdigest()!=source["sha256"]:
                raise RuntimeError("Upstream method-page bytes changed; inspect and preserve old record: "+source["url"])
            path.write_bytes(data)
        assert sha(path)==source["sha256"] and path.stat().st_size==source["bytes"]
    design=json.loads(Path(__file__).with_name("design.json").read_text())
    for filename,digest in design["inputs"].items():
        path=ROOT/filename
        if not path.exists():raise RuntimeError("Run the galaxy_validation reproduction first: "+filename)
        assert sha(path)==digest
    print("Pinned input and source-description identities verified.")
if __name__=="__main__":main()
