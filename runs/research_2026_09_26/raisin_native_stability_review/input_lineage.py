"""Independent SNANA input-row provenance and multiplicity check."""
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
NATIVE=ROOT/"runs/research_2026_09_26/astra_design/raisin_signed_refit"


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def obs(path):
    lines=path.read_text().splitlines();cols=None;found=[];nobs=None
    for line in lines:
        p=line.split()
        if not p:continue
        if p[0]=="VARLIST:":cols=p[1:]
        elif p[0]=="NOBS:":nobs=int(p[1])
        elif p[0]=="OBS:":
            assert cols and len(cols)==len(p)-1
            d=dict(zip(cols,p[1:]));d["original_line"]=line
            found.append(d)
    assert nobs==len(found)
    return found


def numeric_key(r):
    return (float(r["MJD"]),r.get("FLT",r.get("BAND")),float(r["FLUXCAL"]),float(r["FLUXCALERR"]))


def main():
    cohort=list(csv.DictReader((NATIVE/"cohort.csv").open()))
    ledger=list(csv.DictReader((NATIVE/"source-row-ledger.csv").open()))
    result=[]
    for c in cohort:
        cid=c["CID"];raw=obs(ROOT/c["raw_path"]);release=obs(ROOT/c["released_path"])
        arm={a:obs(NATIVE/"data"/f"RSR_{a}"/f"{cid}.snana.dat") for a in "RABH"}
        rr=NATIVE/"data/RSR_R"/f"{cid}.snana.dat"
        rc=NATIVE/"data/RSR_Rcopy"/f"{cid}.snana.dat"
        assert sha(rr)==sha(ROOT/c["released_path"])==sha(rc)
        rawopt=[x for x in raw if x.get("FLT",x.get("BAND")) in "griz"]
        pos=[x for x in rawopt if float(x["FLUXCAL"])>0];neg=[x for x in rawopt if float(x["FLUXCAL"])<=0]
        assert len(pos)==int(c["author_positive"]) and len(neg)==int(c["author_nonpositive"])
        ra=[x for x in arm["A"] if x["FLT"] in "griz"]
        rb=[x for x in arm["B"] if x["FLT"] in "griz"]
        assert Counter(map(numeric_key,ra))==Counter(map(numeric_key,pos))
        assert Counter(map(numeric_key,rb))==Counter(map(numeric_key,rawopt))
        assert arm["H"][:len(release)]==release
        assert Counter(map(numeric_key,arm["H"][len(release):]))==Counter(map(numeric_key,neg))
        assert [x["original_line"] for x in arm["R"]]==[x["original_line"] for x in release]
        assert Counter(map(numeric_key,[x for x in arm["A"] if x["FLT"] not in "griz"]))==Counter(map(numeric_key,[x for x in release if x["FLT"] not in "griz"]))
        source=[x for x in ledger if x["CID"]==cid]
        assert len(source)==len(rawopt)
        assert Counter((float(x["MJD"]),x["band"],float(x["FLUXCAL"]),float(x["FLUXCALERR"])) for x in source)==Counter(map(numeric_key,rawopt))
        result.append(dict(CID=cid,release_optical=sum(x["FLT"] in "griz" for x in release),author_positive=len(pos),author_nonpositive=len(neg),R_optical=sum(x["FLT"] in "griz" for x in arm["R"]),A_optical=len(ra),B_optical=len(rb),H_optical=sum(x["FLT"] in "griz" for x in arm["H"]),R_byte_identical=True,Rcopy_byte_identical=True,A_positive_multiplicity_exact=True,B_signed_multiplicity_exact=True,H_release_prefix_exact=True,H_negative_multiplicity_exact=True))
    with (OUT/"input-lineage.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
    summary=dict(status="PASS",objects=len(result),release_optical=sum(r["release_optical"] for r in result),author_positive=sum(r["author_positive"] for r in result),author_nonpositive=sum(r["author_nonpositive"] for r in result),output_sha256=sha(OUT/"input-lineage.csv"))
    (OUT/"input-lineage-result.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))

if __name__=="__main__":main()
