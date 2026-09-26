"""Exact Git-blob author DES raw-to-product sign lineage, fixed aliases."""
import csv
import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("audit",OUT/"executed_source.py")
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    plan=json.loads((OUT/"author-acquisition-plan.json").read_text())
    acquired=json.loads((OUT/"author-acquisition.json").read_text())
    assert acquired["plan_sha256"]==sha(OUT/"author-acquisition-plan.json")
    for record in acquired["records"]:
        assert sha(OUT/"author-source"/record["path"])==record["sha256"]
    aliases=json.loads((OUT/"protocol.json").read_text())["des_aliases"]
    phot={x["CID"]:x for x in json.loads((ROOT/"runs/research_2026_09_26/raisin_differential/photometry-inventory.json").read_text())}
    original=list(csv.DictReader((OUT/"original-retention.csv").open()))
    allrows=[]; object_summary=[]
    for entry in plan["entries"]:
        if entry["role"]!="raw_DES_ancestor":continue
        cid=entry["raisin_CID"]
        assert aliases[cid]==entry["original_CID"]
        _,raw=audit.read_raisin(OUT/"author-source"/entry["path"])
        product_head,product=audit.read_raisin(ROOT/phot[cid]["path"])
        product=[r for r in product if r["FLT"] in "griz"]
        peak=float(product_head["PEAKMJD"].split()[0]);redshift=float(product_head["REDSHIFT_HELIO"].split()[0])
        prior=[r for r in original if r["raisin_CID"]==cid]
        local=[]
        for row in raw:
            band=row.get("FLT",row.get("BAND"))
            if band not in "griz":continue
            mjd=float(row["MJD"]);flux=float(row["FLUXCAL"]);error=float(row["FLUXCALERR"])
            p=[r for r in product if r["FLT"]==band and abs(float(r["MJD"])-mjd)<=.00055]
            d=[r for r in prior if r["band"]==band and abs(float(r["MJD"])-mjd)<=.00055]
            phase=(mjd-peak)/(1+redshift)
            local.append(dict(CID=cid,original_CID=entry["original_CID"],author_raw_path=entry["path"],author_raw_index=row["observation_index"],author_raw_line=row["line_number"],band=band,MJD=mjd,raw_flux=flux,raw_error=error,raw_sign=audit.sign(flux),raw_PHOTFLAG=row.get("PHOTFLAG",""),
                rest_phase=phase,fit_window=-7<=phase<=45,m20_window=-10<=phase<=40,product_candidate_count=len(p),diffimg_candidate_count=len(d),
                product_index=p[0]["observation_index"] if len(p)==1 else "",product_flux=float(p[0]["FLUXCAL"]) if len(p)==1 else "",product_error=float(p[0]["FLUXCALERR"]) if len(p)==1 else "",
                product_flux_delta=float(p[0]["FLUXCAL"])-flux if len(p)==1 else "",product_error_delta=float(p[0]["FLUXCALERR"])-error if len(p)==1 else "",
                diffimg_global_row=d[0]["original_fits_global_row"] if len(d)==1 else "",diffimg_flux=float(d[0]["original_flux"]) if len(d)==1 else "",diffimg_flux_delta=float(d[0]["original_flux"])-flux if len(d)==1 else ""))
        allrows.extend(local)
        object_summary.append(dict(CID=cid,raw_rows=len(local),raw_negative=sum(x["raw_sign"]=="negative" for x in local),raw_positive=sum(x["raw_sign"]=="positive" for x in local),negative_in_product=sum(x["raw_sign"]=="negative" and x["product_candidate_count"]>0 for x in local),positive_in_product=sum(x["raw_sign"]=="positive" and x["product_candidate_count"]>0 for x in local),raw_to_diffimg_rows=sum(x["diffimg_candidate_count"]>0 for x in local),fitwindow_negative=sum(x["fit_window"] and x["raw_sign"]=="negative" for x in local),fitwindow_negative_in_product=sum(x["fit_window"] and x["raw_sign"]=="negative" and x["product_candidate_count"]>0 for x in local),product_optical_rows=len(product)))
    audit.put_csv(OUT/"author-lineage.csv",allrows,list(allrows[0]))
    audit.put_csv(OUT/"author-object-summary.csv",object_summary,list(object_summary[0]))
    summary={
        "raw_rows":len(allrows),"raw_negative":sum(r["raw_sign"]=="negative" for r in allrows),"raw_positive":sum(r["raw_sign"]=="positive" for r in allrows),
        "negative_in_product":sum(r["raw_sign"]=="negative" and r["product_candidate_count"]>0 for r in allrows),
        "positive_in_product":sum(r["raw_sign"]=="positive" and r["product_candidate_count"]>0 for r in allrows),
        "raw_to_diffimg_rows":sum(r["diffimg_candidate_count"]>0 for r in allrows),
        "fitwindow_negative":sum(r["fit_window"] and r["raw_sign"]=="negative" for r in allrows),
        "fitwindow_negative_in_product":sum(r["fit_window"] and r["raw_sign"]=="negative" and r["product_candidate_count"]>0 for r in allrows),
        "ambiguous_product_rows":sum(r["product_candidate_count"]>1 for r in allrows),
        "raw_positive_matched_product_rounding_count":sum(r["raw_sign"]=="positive" and r["product_candidate_count"]==1 and abs(r["product_flux_delta"])<=.001 for r in allrows),
    }
    product_same={}
    for cid in ["DES15E2mhy","DES16S1agd"]:
        src=OUT/"author-source"/f"data/Photometry/DES_RAISIN/{cid}.snana.dat"
        rel=ROOT/phot[cid]["path"]
        product_same[cid]=sha(src)==sha(rel)
    result={"author_acquisition_sha256":sha(OUT/"author-acquisition.json"),"summary":summary,"author_product_byte_identical_to_release":product_same,"note":"Exact Git blob raw precursor matches by band/MJD, while 2015 raw files have limited in-season rows. This does not identify the transformation script or imply same DIFFIMG revision for every unmatched row."}
    (OUT/"author-lineage-result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":main()
