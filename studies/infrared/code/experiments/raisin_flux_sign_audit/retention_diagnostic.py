"""Source-row DIFFIMG sign/flag/phase retention ledger; no noise fit."""
import csv
import hashlib
import importlib.util
import json
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("audit", OUT / "executed_source.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    amendment = json.loads((OUT / "phase-window-amendment.json").read_text())
    assert amendment["original_protocol_sha256"] == sha(OUT / "protocol.json")
    for rel, digest in amendment["nml_sha256"].items():
        assert sha(ROOT / rel) == digest
    protocol = json.loads((OUT / "protocol.json").read_text())
    prior = json.loads((ROOT / "runs/research_2026_09_26/raisin_differential/photometry-inventory.json").read_text())
    pmap = {row["CID"]: row for row in prior}
    aliases = protocol["des_aliases"]
    original, heads = audit.load_original("DIFFIMG", ROOT / "data/des-diffimg/DES-SN5YR_DIFFIMG_HEAD.FITS.gz", ROOT / "data/des-diffimg/DES-SN5YR_DIFFIMG_PHOT.FITS.gz", set(aliases.values()))
    all_rows = []
    for rcid, ocid in sorted(aliases.items()):
        header, raisin = audit.read_raisin(ROOT / pmap[rcid]["path"])
        peak = float(header["PEAKMJD"].split()[0])
        redshift = float(header["REDSHIFT_HELIO"].split()[0])
        ropt = [r for r in raisin if r["FLT"] in "griz"]
        for old in original[ocid]:
            if old["band"] not in "griz":
                continue
            phase = (old["MJD"] - peak) / (1 + redshift)
            matches = [r for r in ropt if r["FLT"] == old["band"] and abs(float(r["MJD"]) - old["MJD"]) <= .00055]
            nearest = min(matches, key=lambda r: abs(float(r["MJD"]) - old["MJD"])) if matches else None
            all_rows.append(dict(raisin_CID=rcid,original_CID=ocid,original_head_index=old["head_index"],original_fits_global_row=old["fits_global_row"],original_local_row=old["local_row"],
                band=old["band"],MJD=old["MJD"],original_flux=old["flux"],original_error=old["error"],original_sign=audit.sign(old["flux"]),PHOTFLAG=old["PHOTFLAG"],IMGNUM=old["IMGNUM"],
                rest_phase=phase,original_optical_fit_window=(-7<=phase<=45),m20_window=(-10<=phase<=40),offseason=abs(old["MJD"]-peak)>180,
                raisin_match_count=len(matches),retained=bool(matches),raisin_index=nearest["observation_index"] if nearest else "",raisin_line=nearest["line_number"] if nearest else "",
                raisin_flux=float(nearest["FLUXCAL"]) if nearest else "",raisin_error=float(nearest["FLUXCALERR"]) if nearest else "",
                flux_delta=float(nearest["FLUXCAL"])-old["flux"] if nearest else "",error_delta=float(nearest["FLUXCALERR"])-old["error"] if nearest else ""))
    audit.put_csv(OUT / "original-retention.csv", all_rows, list(all_rows[0]))
    summary = {}
    for window in ["all", "original_optical_fit_window", "m20_window", "offseason"]:
        subset = all_rows if window == "all" else [r for r in all_rows if r[window]]
        summary[window] = {
            "rows": len(subset),
            "by_sign": {s: {"total": sum(r["original_sign"]==s for r in subset), "retained": sum(r["original_sign"]==s and r["retained"] for r in subset)} for s in ["negative", "zero", "positive", "nonfinite"]},
            "retained": sum(r["retained"] for r in subset),
            "ambiguous": sum(r["raisin_match_count"]>1 for r in subset),
        }
    by_flag = Counter((r["PHOTFLAG"],r["original_sign"],r["retained"]) for r in all_rows)
    flags = [dict(PHOTFLAG=k[0],sign=k[1],retained=k[2],rows=v) for k,v in sorted(by_flag.items())]
    audit.put_csv(OUT / "flag-retention.csv", flags, ["PHOTFLAG","sign","retained","rows"])
    per_object=[]
    for rcid in sorted(aliases):
        rows=[r for r in all_rows if r["raisin_CID"]==rcid]
        per_object.append(dict(CID=rcid,original_ID=aliases[rcid],original_optical=len(rows),original_negative=sum(r["original_sign"]=="negative" for r in rows),original_positive=sum(r["original_sign"]=="positive" for r in rows),negative_retained=sum(r["original_sign"]=="negative" and r["retained"] for r in rows),positive_retained=sum(r["original_sign"]=="positive" and r["retained"] for r in rows),
            fitwindow_negative=sum(r["original_optical_fit_window"] and r["original_sign"]=="negative" for r in rows),fitwindow_negative_retained=sum(r["original_optical_fit_window"] and r["original_sign"]=="negative" and r["retained"] for r in rows)))
    audit.put_csv(OUT / "per-object-retention.csv",per_object,list(per_object[0]))
    matched=[r for r in all_rows if r["retained"]]
    result=dict(protocol_sha256=sha(OUT/"protocol.json"),phase_window_amendment_sha256=sha(OUT/"phase-window-amendment.json"),head_hits=heads,summary=summary,
        positive_retained_flux_abs_delta_quantiles={str(q):float(sorted(abs(r["flux_delta"]) for r in matched)[round(q*(len(matched)-1))]) for q in [0,.5,.9,.99,1]},
        positive_retained_error_abs_delta_quantiles={str(q):float(sorted(abs(r["error_delta"]) for r in matched)[round(q*(len(matched)-1))]) for q in [0,.5,.9,.99,1]},
        matched_flux_decimal_rounding_count=sum(abs(r["flux_delta"])<=.001 for r in matched),matched_error_decimal_rounding_count=sum(abs(r["error_delta"])<=.001 for r in matched),
        note="Exact same-band/time candidate retention in the released DIFFIMG revision; no historic author-pipeline cause inferred from this alone.")
    (OUT/"retention-result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result["summary"],indent=2))


if __name__ == "__main__":
    main()
