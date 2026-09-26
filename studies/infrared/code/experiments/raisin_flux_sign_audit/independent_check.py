"""Independent arithmetic and provenance check of the sign-audit ledgers."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]


def rows(name):
    with (OUT/name).open(newline="") as f:
        return list(csv.DictReader(f))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    inventory=rows("inventory-files.csv")
    assert len(inventory)==117
    assert sum(int(x["rows"]) for x in inventory)==23007
    assert sum(int(x["negative"]) for x in inventory)==0
    assert sum(int(x["zero"]) for x in inventory)==0
    assert sum(int(x["positive"]) for x in inventory)==23007
    assert Counter(x["survey"] for x in inventory)=={"CSP":77,"DES":19,"PS1MD":21}

    old=rows("original-retention.csv")
    assert len(old)==9843
    assert sum(x["original_sign"]=="negative" for x in old)==4229
    assert sum(x["original_sign"]=="negative" and x["retained"]=="True" for x in old)==0
    assert sum(x["original_sign"]=="positive" and x["retained"]=="True" for x in old)==5108
    win=[x for x in old if x["original_optical_fit_window"]=="True"]
    assert len(win)==867
    assert Counter((x["original_sign"],x["retained"]) for x in win)=={("negative","False"):56,("positive","True"):811}

    lineage=rows("author-lineage.csv")
    d16=[x for x in lineage if x["CID"].startswith("DES16")]
    assert len(d16)==5378
    assert Counter(x["raw_sign"] for x in d16)=={"positive":3135,"negative":2243}
    assert all(x["product_candidate_count"]!="0" for x in d16 if x["raw_sign"]=="positive")
    assert all(x["product_candidate_count"]=="0" for x in d16 if x["raw_sign"]=="negative")
    product=rows("author-product-blob-ledger.csv")
    assert len(product)==17 and all(x["byte_identical_to_author"]=="True" for x in product)
    for x in product:
        assert sha(ROOT/x["release_path"])==x["release_sha256"]
    candidates=rows("author-negative-candidates.csv")
    assert len(candidates)==2286
    assert sum(x["CID"].startswith("DES16") and x["historical_fit_window"]=="True" for x in candidates)==33
    result={"status":"PASS","checks":["117-file sign inventory","original DIFFIMG sign/phase retention","DES16 signed author precursor","17 author product Git blob identities","negative candidate ledger"],"source_hashes":{name:sha(OUT/name) for name in ["inventory-files.csv","original-retention.csv","author-lineage.csv","author-product-blob-ledger.csv","author-negative-candidates.csv"]}}
    (OUT/"independent-check.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))

if __name__=="__main__":
    main()
