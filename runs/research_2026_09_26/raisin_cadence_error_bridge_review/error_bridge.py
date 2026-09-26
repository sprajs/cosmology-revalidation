"""Direct author-raw/product numerical error-conversion check, no fitted factor."""
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
SIGN=ROOT/"runs/research_2026_09_26/raisin_flux_sign_audit"
F=1.086*0.4*math.log(10)


def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def read_snana(path):
    fields=None;rows=[]
    for line_no,line in enumerate(path.read_text().splitlines(),1):
        if line.startswith("VARLIST:"):
            fields=line.split()[1:]
        elif line.startswith("OBS:"):
            assert fields is not None
            data=line.split()[1:]
            assert len(data)==len(fields),(path,line_no)
            row=dict(zip(fields,data));row["source_line"]=line_no
            rows.append(row)
    return rows


def q(values,p):
    values=sorted(values)
    return values[round(p*(len(values)-1))]


def write_csv(path,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def summary(rows):
    unique=[r for r in rows if r["product_candidate_count"]==1 and r["raw_sign"]=="positive"]
    if not unique:return None
    diffs=[abs(r["error_delta"]) for r in unique]
    ratios=[r["error_ratio"] for r in unique]
    res=[abs(r["factor_unrounded_residual"]) for r in unique]
    return dict(n=len(unique),plain_round_exact=sum(r["plain_round_exact"] for r in unique),factor_round_exact=sum(r["factor_round_exact"] for r in unique),flux_round_exact=sum(r["flux_round_exact"] for r in unique),error_abs_median=q(diffs,.5),error_abs_p90=q(diffs,.9),error_abs_max=max(diffs),error_ratio_median=q(ratios,.5),error_ratio_p90=q(ratios,.9),factor_unrounded_abs_max=max(res))


def main():
    protocol=json.loads((OUT/"protocol.json").read_text())
    amendment=json.loads((OUT/"error-formula-amendment.json").read_text())
    assert amendment["parent_protocol_sha256"]==sha(OUT/"protocol.json")
    for rel,h in amendment["source_sha256"].items():assert sha(ROOT/rel)==h
    acquired=json.loads((SIGN/"author-acquisition.json").read_text())
    assert acquired["plan_sha256"]==sha(SIGN/"author-acquisition-plan.json")
    for record in acquired["records"]:
        assert sha(SIGN/"author-source"/record["path"])==record["sha256"]
    plan=json.loads((SIGN/"author-acquisition-plan.json").read_text())
    inventory={r["CID"]:r for r in json.loads((ROOT/"runs/research_2026_09_26/raisin_differential/photometry-inventory.json").read_text())}
    rows=[]
    for entry in plan["entries"]:
        if entry["role"]!="raw_DES_ancestor":continue
        cid=entry["raisin_CID"]
        raw=read_snana(SIGN/"author-source"/entry["path"])
        released=read_snana(ROOT/inventory[cid]["path"])
        released=[x for x in released if x["FLT"] in "griz"]
        for i,x in enumerate(raw):
            band=x.get("FLT",x.get("BAND"))
            if band not in "griz":continue
            mjd=float(x["MJD"])
            cand=[y for y in released if y["FLT"]==band and abs(float(y["MJD"])-mjd)<=.00055]
            rawflux=float(x["FLUXCAL"]);rawerr=float(x["FLUXCALERR"])
            assert rawerr>0
            y=cand[0] if len(cand)==1 else None
            err=float(y["FLUXCALERR"]) if y else None
            flux=float(y["FLUXCAL"]) if y else None
            rows.append(dict(CID=cid,source_path=entry["path"],source_line=x["source_line"],source_index=i,band=band,MJD=mjd,
                raw_sign="negative" if rawflux<0 else "positive" if rawflux>0 else "zero",raw_flux=rawflux,raw_error=rawerr,
                PHOTFLAG=int(x.get("PHOTFLAG",-1)),ZPFLUX=x.get("ZPFLUX",""),PSF=x.get("PSF",""),SKYSIG=x.get("SKYSIG",""),GAIN=x.get("GAIN",""),
                product_candidate_count=len(cand),product_line=y["source_line"] if y else "",product_flux=flux if y else "",product_error=err if y else "",
                error_delta=(err-rawerr) if y else "",error_ratio=(err/rawerr) if y else "",
                predicted_error=round(rawerr*F,3),predicted_plain_round=round(rawerr,3),
                factor_unrounded_residual=(err-rawerr*F) if y else "",
                factor_round_exact=(err==round(rawerr*F,3)) if y else False,
                plain_round_exact=(err==round(rawerr,3)) if y else False,
                flux_round_exact=(flux==round(rawflux,3)) if y else False))
    write_csv(OUT/"raw-product-error-rows.csv",rows)
    groups={}
    for key,make in (("object",lambda r:r["CID"]),("band",lambda r:r["band"]),("PHOTFLAG",lambda r:str(r["PHOTFLAG"]))):
        grouped=defaultdict(list)
        for r in rows:
            if r["CID"].startswith("DES16"):grouped[make(r)].append(r)
        groups[key]={label:summary(value) for label,value in sorted(grouped.items())}
    d16=[r for r in rows if r["CID"].startswith("DES16")]
    d15=[r for r in rows if r["CID"].startswith("DES15")]
    result=dict(status="PASS",factor=F,formula="round(raw_error * 1.086 * 0.4 * ln(10),3)",source_mechanism="Numerically exact; author optical execution source not located.",
        des16_rows=len(d16),des16_positive=sum(r["raw_sign"]=="positive" for r in d16),des16_negative=sum(r["raw_sign"]=="negative" for r in d16),
        des16_unique_positive=summary(d16),des16_positive_ambiguous=sum(r["raw_sign"]=="positive" and r["product_candidate_count"]>1 for r in d16),
        des16_negative_in_product=sum(r["raw_sign"]=="negative" and r["product_candidate_count"]>0 for r in d16),des15_unique_positive=summary(d15),
        grouped=groups,inputs_sha256={"protocol.json":sha(OUT/"protocol.json"),"error-formula-amendment.json":sha(OUT/"error-formula-amendment.json"),"author-acquisition-plan.json":sha(SIGN/"author-acquisition-plan.json")},output_sha256=sha(OUT/"raw-product-error-rows.csv"))
    (OUT/"error-bridge-result.json").write_text(json.dumps(result,indent=2)+"\n")
    assert result["des16_unique_positive"]["factor_round_exact"]==3133
    print(json.dumps({k:v for k,v in result.items() if k!="grouped"},indent=2))

if __name__=="__main__":main()
