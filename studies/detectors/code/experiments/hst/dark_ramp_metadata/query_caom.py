"""Official MAST CAOM metadata discovery; never downloads FITS products."""
from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
URL = "https://mast.stsci.edu/api/v0/invoke"


def sha(data):
    return hashlib.sha256(data).hexdigest()


def request(payload, name, limit, started, transferred):
    if time.monotonic() - started > limit["wall_seconds"]:
        raise TimeoutError("300s metadata cap")
    encoded = urlencode({"request":json.dumps(payload,separators=(",",":"))}).encode()
    req = Request(URL, data=encoded, headers={"Content-type":"application/x-www-form-urlencoded",
                                               "Accept":"text/plain", "User-agent":"codex-metadata-audit/1"})
    with urlopen(req, timeout=35) as response:
        cap = limit["metadata_total_bytes"] - transferred
        data = response.read(cap+1)
        if len(data) > cap:
            raise ValueError("20MB metadata transfer cap")
        code = response.status
    path = HERE / "queries" / f"{name}.json"
    path.write_bytes(data)
    result = json.loads(data)
    if result.get("status") != "COMPLETE":
        raise ValueError(f"MAST status {result.get('status')}: see {path}")
    return result, {"name":name,"url":URL,"payload":payload,"status_code":code,
                    "response_path":str(path),"response_sha256":sha(data),"bytes":len(data)}


def filters(protocol, mid, window, target):
    f=[{"paramName":k,"values":v} for k,v in protocol["caom_filters"].items()
       if k not in ("target_name_freeText",)]
    f.append({"paramName":"t_min","values":[{"min":mid-window,"max":mid+window}]})
    if target:
        f.append({"paramName":"target_name","values":[],
                  "freeText":protocol["caom_filters"]["target_name_freeText"]})
    return f


def count_from(response):
    data=response.get("data",[])
    if len(data)!=1:
        raise ValueError("unexpected MAST count response schema")
    values=list(data[0].values())
    if len(values)!=1:
        raise ValueError("unexpected MAST count columns")
    return int(values[0])


def main():
    protocol=json.loads((HERE/"protocol.json").read_text())
    HERE.joinpath("queries").mkdir(exist_ok=True)
    started=time.monotonic();transferred=0;log=[];out={}
    for visit,mid in protocol["science_visit_midpoints_mjd"].items():
        tiers=[]
        for window in protocol["tiered_time_windows_days"]:
            for target in (True,False):
                base=filters(protocol,mid,window,target)
                label=f"{visit}_{window}d_{'dark_target' if target else 'fallback'}"
                count_payload={"service":"Mast.Caom.Filtered","format":"json",
                               "params":{"columns":"COUNT_BIG(*)","filters":base,"obstype":"all"}}
                count,rec=request(count_payload,label+"_count",protocol["limits"],started,transferred)
                log.append(rec);transferred+=rec["bytes"]
                n=count_from(count)
                if n == 0 and target:
                    continue
                rows=[];page=1
                while len(rows)<n:
                    payload={"service":"Mast.Caom.Filtered","format":"json",
                             "params":{"columns":"*","filters":base,"obstype":"all"},
                             "pagesize":500,"page":page}
                    answer,rec=request(payload,label+f"_page{page}",protocol["limits"],started,transferred)
                    log.append(rec);transferred+=rec["bytes"]
                    got=answer.get("data",[])
                    if not got:
                        break
                    rows.extend(got)
                    page+=1
                    if len(rows)>n+500:
                        raise ValueError("pagination exceeded count")
                tiers.append({"window_days":window,"target_filter":target,"count":n,
                              "retrieved_rows":len(rows),"rows":rows})
                break
            # Eligibility requires FITS-header inspection later; carry first
            # window's observations, then expand if fewer than four distinct
            # candidate roots even exist at CAOM level.
            if tiers and len({r.get("obs_id") for r in tiers[-1]["rows"]})>=4:
                break
        out[visit]=tiers
    saved=HERE/"caom-results.json"
    saved.write_text(json.dumps({"visits":out,"requests":log,"elapsed_seconds":time.monotonic()-started,
                                 "bytes_transferred":transferred},indent=2)+"\n")
    print(json.dumps({"visits":{k:[{"window_days":t["window_days"],"count":t["count"],
                                   "retrieved_rows":t["retrieved_rows"]} for t in v]
                                 for k,v in out.items()},
                      "bytes_transferred":transferred,"elapsed_seconds":time.monotonic()-started,
                      "result_sha256":sha(saved.read_bytes())},indent=2))


if __name__=="__main__":
    main()
