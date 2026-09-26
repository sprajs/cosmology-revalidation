"""Derive immutable row-level lineage handoff from frozen sign-audit ledgers."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def put(path, rows):
    assert rows
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def main():
    lineage = read(OUT / "author-lineage.csv")
    original = read(OUT / "original-retention.csv")
    by_original = {x["original_fits_global_row"]: x for x in original}
    phot = {x["CID"]: x for x in json.loads((ROOT / "runs/research_2026_09_26/raisin_differential/photometry-inventory.json").read_text())}
    tree = json.loads((ROOT / "runs/research_2026_09_26/raisin_differential/mass-threshold/code-tree.json").read_text())
    blob = {x["path"]: x["sha"] for x in tree["tree"]}
    aliases = json.loads((OUT / "protocol.json").read_text())["des_aliases"]
    products = []
    for cid in sorted(aliases):
        p = ROOT / phot[cid]["path"]
        data = p.read_bytes()
        git_blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        author_path = f"data/Photometry/DES_RAISIN/{cid}.snana.dat"
        products.append(dict(CID=cid, release_path=phot[cid]["path"], release_sha256=digest(p), author_git_path=author_path, author_git_blob_sha1=blob.get(author_path, ""), release_git_blob_sha1=git_blob, byte_identical_to_author=str(git_blob == blob.get(author_path))))
    assert len(products) == 17 and all(x["byte_identical_to_author"] == "True" for x in products)
    put(OUT / "author-product-blob-ledger.csv", products)

    candidates = []
    for x in lineage:
        if x["raw_sign"] != "negative":
            continue
        o = by_original.get(x["diffimg_global_row"])
        candidates.append(dict(CID=x["CID"], original_CID=x["original_CID"], source_path=x["author_raw_path"], source_line=x["author_raw_line"], source_observation_index=x["author_raw_index"], band=x["band"], MJD=x["MJD"], FLUXCAL=x["raw_flux"], FLUXCALERR=x["raw_error"], PHOTFLAG=x["raw_PHOTFLAG"], rest_phase=x["rest_phase"], historical_fit_window=x["fit_window"], M20_window=x["m20_window"], released_product_match_count=x["product_candidate_count"], later_diffimg_match_count=x["diffimg_candidate_count"], later_diffimg_global_fits_row=x["diffimg_global_row"], later_diffimg_flux=x["diffimg_flux"], later_diffimg_error=o["original_error"] if o else "", later_diffimg_PHOTFLAG=o["PHOTFLAG"] if o else "", later_diffimg_IMGNUM=o["IMGNUM"] if o else "", later_diffimg_sign=o["original_sign"] if o else "", source_units="FLUXCAL, SNANA reference magnitude 27.5"))
    put(OUT / "author-negative-candidates.csv", candidates)

    object_gates = []
    groups = defaultdict(list)
    for x in lineage:
        groups[x["CID"]].append(x)
    for cid in sorted(groups):
        rows = groups[cid]
        pos = [x for x in rows if x["raw_sign"] == "positive"]
        neg = [x for x in rows if x["raw_sign"] == "negative"]
        unique = [x for x in pos if x["product_candidate_count"] == "1" and x["diffimg_candidate_count"] == "1"]
        def close(x, key, tol):
            return abs(float(x[key])) <= tol
        object_gates.append(dict(CID=cid, year="2016" if cid.startswith("DES16") else "2015", raw_positive=len(pos), raw_negative=len(neg), positive_product_same_time=sum(x["product_candidate_count"] != "0" for x in pos), positive_unique_both=len(unique), positive_product_flux_within_0p001=sum(close(x,"product_flux_delta",.001) for x in unique), positive_product_error_within_0p001=sum(close(x,"product_error_delta",.001) for x in unique), positive_product_error_within_0p005=sum(close(x,"product_error_delta",.005) for x in unique), positive_diffimg_flux_within_0p001=sum(close(x,"diffimg_flux_delta",.001) for x in unique), negative_product_same_time=sum(x["product_candidate_count"] != "0" for x in neg), negative_later_diffimg_unique=sum(x["diffimg_candidate_count"] == "1" for x in neg), negative_later_diffimg_same_sign=sum(x["diffimg_candidate_count"] == "1" and float(x["diffimg_flux"]) < 0 for x in neg), negative_historical_fit_window=sum(x["fit_window"] == "True" for x in neg), negative_fit_window_product_same_time=sum(x["fit_window"] == "True" and x["product_candidate_count"] != "0" for x in neg)))
    put(OUT / "lineage-object-gates.csv", object_gates)
    result = dict(inputs_sha256={p.name:digest(p) for p in [OUT/"protocol.json",OUT/"phase-window-amendment.json",OUT/"author-acquisition.json",OUT/"author-lineage.csv",OUT/"original-retention.csv",ROOT/"runs/research_2026_09_26/raisin_differential/mass-threshold/code-tree.json"]}, counts=dict(author_products_identical=sum(x["byte_identical_to_author"]=="True" for x in products),negative_candidates=len(candidates),des16_negative_candidates=sum(x["CID"].startswith("DES16") for x in candidates),des16_negative_historical_fit_window=sum(x["CID"].startswith("DES16") and x["historical_fit_window"]=="True" for x in candidates),des15_negative_candidates=sum(x["CID"].startswith("DES15") for x in candidates)), caveat="DES16 positive flux rows agree to product rounding, but their errors do not all agree to 0.001; 2015 precursor positives have changed flux and error. Native refit restoration using author-negative errors therefore requires an explicit conditional bridge, not an exact byte-level restoration.")
    (OUT/"lineage-handoff.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__ == "__main__":
    main()
