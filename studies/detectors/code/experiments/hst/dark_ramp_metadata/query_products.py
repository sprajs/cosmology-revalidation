"""MAST product metadata for date-ranked candidate dark observations only."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import time

from query_caom import request

HERE=Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=json.loads((HERE/"protocol.json").read_text())
    gate=json.loads((HERE/"products-protocol.json").read_text())
    assert gate["caom_result_sha256"]==sha(HERE/"caom-results.json")
    caom=json.loads((HERE/"caom-results.json").read_text())
    elapsed=caom["elapsed_seconds"]
    transferred=caom["bytes_transferred"]
    started=time.monotonic()
    out={}
    for visit,tiers in caom["visits"].items():
        mid=p["science_visit_midpoints_mjd"][visit]
        records=[]
        for tier in tiers:
            records.extend(tier["rows"])
        # Multiple tiers, if present, retain one CAOM row per observation.
        unique={r["obs_id"]:r for r in records}
        ranked=sorted((r for r in unique.values() if float(r["t_exptime"])>=300),
                      key=lambda r:(abs(float(r["t_min"])-mid),r["obs_id"]))
        selected=ranked[:gate["max_roots_per_visit"]]
        payload={"service":"Mast.Caom.Products","format":"json",
                 "params":{"obsid":",".join(str(r["obsid"]) for r in selected)}}
        result,rec=request(payload,visit+"_products",p["limits"],started,transferred)
        transferred+=rec["bytes"]
        out[visit]={"eligible_caom_exptime_count":len(ranked),
                    "candidate_roots":selected,"products":result.get("data",[]),
                    "query":rec}
    dest=HERE/"product-results.json"
    dest.write_text(json.dumps({"visits":out,"total_metadata_bytes":transferred,
                                "caom_seconds":elapsed,"product_seconds":time.monotonic()-started},indent=2)+"\n")
    print(json.dumps({"visits":{v:{"ranked":len(z["candidate_roots"]),
                                   "products":len(z["products"])} for v,z in out.items()},
                      "total_metadata_bytes":transferred,"sha256":sha(dest)},indent=2))


if __name__=="__main__":
    main()
