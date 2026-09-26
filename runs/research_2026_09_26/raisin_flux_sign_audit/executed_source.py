"""Frozen RAISIN flux-sign inventory and source-row DES epoch matching."""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/research_2026_09_26/raisin_flux_sign_audit"
PROTOCOL = OUT / "protocol.json"
PRIOR = ROOT / "runs/research_2026_09_26/raisin_differential"
SMP = ROOT / "sources/repos/des-science__DES-SN5YR@1.3/0_DATA/DES-SN5YR_DES"
DIFF = ROOT / "data/des-diffimg"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def put_csv(path, rows, columns):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def read_raisin(path):
    header = {}
    columns = None
    rows = []
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        if not line or line.startswith("#") or ":" not in line:
            continue
        key, value = line.split(":", 1)
        value = value.strip()
        if key == "VARLIST":
            columns = value.split()
        elif key == "OBS":
            assert columns is not None
            values = value.split()
            assert len(values) == len(columns)
            rows.append(dict(zip(columns, values), line_number=lineno, observation_index=len(rows)))
        elif key not in header:
            header[key] = value
    assert len(rows) == int(header["NOBS"])
    return header, rows


def phase_bin(phase):
    if phase < -60:
        return "lt_minus60"
    if phase < -15:
        return "minus60_to_minus15"
    if phase <= 45:
        return "minus15_to_45"
    if phase <= 100:
        return "45_to_100"
    return "gt_100"


def sign(flux):
    return "nonfinite" if not math.isfinite(flux) else ("negative" if flux < 0 else "zero" if flux == 0 else "positive")


def inventory(phot):
    files = []
    aggregates = defaultdict(Counter)
    negative_rows = []
    for obj in phot:
        path = ROOT / obj["path"]
        head, rows = read_raisin(path)
        assert head["SNID"].split()[0] == obj["CID"]
        assert head["SURVEY"].split()[0] == obj["survey"]
        peak = float(head["PEAKMJD"].split()[0])
        redshift = float(head["REDSHIFT_HELIO"].split()[0])
        assert redshift > -1
        counts = Counter()
        for row in rows:
            flux = float(row["FLUXCAL"])
            error = float(row["FLUXCALERR"])
            band = row["FLT"]
            mjd = float(row["MJD"])
            phase = (mjd - peak) / (1 + redshift)
            bucket = phase_bin(phase)
            label = sign(flux)
            offseason = abs(mjd - peak) > 180
            counts[label] += 1
            counts["offseason"] += int(offseason)
            counts["nonpositive_error"] += int(not math.isfinite(error) or error <= 0)
            for level, value in [("survey", obj["survey"]), ("band", f'{obj["survey"]}/{band}'), ("phase", f'{obj["survey"]}/{bucket}'), ("survey_band_phase", f'{obj["survey"]}/{band}/{bucket}')]:
                aggregates[(level, value)][label] += 1
                aggregates[(level, value)]["offseason"] += int(offseason)
                aggregates[(level, value)]["nonpositive_error"] += int(not math.isfinite(error) or error <= 0)
            if label == "negative":
                negative_rows.append(dict(CID=obj["CID"],survey=obj["survey"],band=band,MJD=mjd,flux=flux,error=error,phase=phase,line_number=row["line_number"]))
        assert counts["negative"] == obj["negative_flux_count"]
        files.append(dict(CID=obj["CID"],survey=obj["survey"],path=obj["path"],rows=len(rows),**{k:counts[k] for k in ["negative","zero","positive","nonfinite","offseason","nonpositive_error"]}))
    put_csv(OUT / "inventory-files.csv", files, list(files[0]))
    summary = [dict(level=level,value=value,**{k:c[k] for k in ["negative","zero","positive","nonfinite","offseason","nonpositive_error"]}) for (level,value),c in sorted(aggregates.items())]
    put_csv(OUT / "inventory-strata.csv", summary, list(summary[0]))
    put_csv(OUT / "negative-rows.csv", negative_rows, ["CID","survey","band","MJD","flux","error","phase","line_number"])
    return files, summary


def as_text(value):
    return value.decode().strip() if isinstance(value, bytes) else str(value).strip()


def load_original(label, headpath, photpath, aliases):
    by_cid = defaultdict(list)
    with fits.open(headpath, memmap=False) as hd:
        head = hd[1].data
        for i in range(len(head)):
            cid = as_text(head["SNID"][i])
            if cid in aliases:
                by_cid[cid].append((i,int(head["PTROBS_MIN"][i]),int(head["PTROBS_MAX"][i]),as_text(head["IAUC"][i])))
    selected = {}
    with fits.open(photpath, memmap=False) as hd:
        phot = hd[1].data
        for cid in sorted(aliases):
            hits = by_cid[cid]
            selected[cid] = []
            for head_index,start,end,iauc in hits:
                assert start >= 1 and end >= start
                for global_index in range(start-1,end):
                    mjd = float(phot["MJD"][global_index])
                    if not math.isfinite(mjd) or mjd < 0:
                        continue
                    selected[cid].append(dict(source=label,CID=cid,head_index=head_index,head_IAUC=iauc,
                        fits_global_row=global_index,local_row=global_index-(start-1),
                        MJD=mjd,band=as_text(phot["BAND"][global_index]),
                        flux=float(phot["FLUXCAL"][global_index]),error=float(phot["FLUXCALERR"][global_index]),
                        PHOTFLAG=int(phot["PHOTFLAG"][global_index]),IMGNUM=int(phot["IMGNUM"][global_index])))
    return selected, {cid:len(by_cid[cid]) for cid in sorted(aliases)}


def match(phot,aliases):
    pmap = {x["CID"]:x for x in phot}
    original = {}
    headcounts = {}
    for source, hp, pp in [
        ("SMP",SMP/"DES-SN5YR_DES_HEAD.FITS.gz",SMP/"DES-SN5YR_DES_PHOT.FITS.gz"),
        ("DIFFIMG",DIFF/"DES-SN5YR_DIFFIMG_HEAD.FITS.gz",DIFF/"DES-SN5YR_DIFFIMG_PHOT.FITS.gz")]:
        original[source],headcounts[source] = load_original(source,hp,pp,set(aliases.values()))
    candidates = []
    per_row = []
    for raisin_cid,original_cid in sorted(aliases.items()):
        _, rows = read_raisin(ROOT / pmap[raisin_cid]["path"])
        for raisin in rows:
            band = raisin["FLT"]
            if band not in "griz":
                continue
            time = float(raisin["MJD"])
            rflux = float(raisin["FLUXCAL"])
            rerr = float(raisin["FLUXCALERR"])
            for source in ["SMP","DIFFIMG"]:
                nearby = [r for r in original[source][original_cid] if r["band"]==band and abs(r["MJD"]-time)<=.01]
                primary = [r for r in nearby if abs(r["MJD"]-time)<=.00055]
                per_row.append(dict(raisin_CID=raisin_cid,original_CID=original_cid,raisin_observation_index=raisin["observation_index"],raisin_line=raisin["line_number"],source=source,band=band,raisin_MJD=time,raisin_flux=rflux,raisin_error=rerr,primary_candidate_count=len(primary),secondary_candidate_count=len(nearby),primary_original_negative_count=sum(r["flux"]<0 for r in primary),secondary_original_negative_count=sum(r["flux"]<0 for r in nearby)))
                for old in nearby:
                    dt=abs(old["MJD"]-time)
                    ferr=abs(old["error"]-rerr)
                    signed_error=abs(old["flux"]-rflux)
                    abs_error=abs(abs(old["flux"])-rflux)
                    candidates.append(dict(raisin_CID=raisin_cid,original_CID=original_cid,raisin_observation_index=raisin["observation_index"],raisin_line=raisin["line_number"],source=source,band=band,raisin_MJD=time,raisin_flux=rflux,raisin_error=rerr,
                        original_head_index=old["head_index"],original_head_IAUC=old["head_IAUC"],original_fits_global_row=old["fits_global_row"],original_local_row=old["local_row"],original_MJD=old["MJD"],original_flux=old["flux"],original_error=old["error"],original_PHOTFLAG=old["PHOTFLAG"],original_IMGNUM=old["IMGNUM"],
                        dt_days=dt,primary_time=dt<=.00055,flux_signed_delta=signed_error,flux_absolute_delta=abs_error,error_delta=ferr,
                        signed_decimal_match=dt<=.00055 and signed_error<=.001 and ferr<=.001,
                        absolute_decimal_match=dt<=.00055 and abs_error<=.001 and ferr<=.001,
                        positive_only_decimal_match=dt<=.00055 and old["flux"]>0 and signed_error<=.001 and ferr<=.001))
    put_csv(OUT/"match-candidates.csv",candidates,list(candidates[0]))
    put_csv(OUT/"match-rows.csv",per_row,list(per_row[0]))
    return candidates,per_row,headcounts,original


def main():
    protocol = json.loads(PROTOCOL.read_text())
    assert protocol["status"] == "before_current_sign_inventory_and_epoch_matching"
    for rel,digest in protocol["inputs_sha256"].items():
        assert sha(ROOT/rel)==digest,rel
    assert not (OUT/"result.json").exists()
    phot=json.loads((PRIOR/"photometry-inventory.json").read_text())
    aliases=protocol["des_aliases"]
    files,strata=inventory(phot)
    candidates,rows,headcounts,original=match(phot,aliases)
    totals={source:dict(raisin_optical_rows=sum(r["source"]==source for r in rows),
        primary_any=sum(r["source"]==source and r["primary_candidate_count"]>0 for r in rows),
        primary_unique=sum(r["source"]==source and r["primary_candidate_count"]==1 for r in rows),
        primary_ambiguous=sum(r["source"]==source and r["primary_candidate_count"]>1 for r in rows),
        secondary_any=sum(r["source"]==source and r["secondary_candidate_count"]>0 for r in rows),
        primary_negative_candidate_rows=sum(r["source"]==source and r["primary_original_negative_count"]>0 for r in rows),
        signed_decimal_match_rows=len({(c["raisin_CID"],c["raisin_observation_index"]) for c in candidates if c["source"]==source and c["signed_decimal_match"]}),
        absolute_decimal_match_rows=len({(c["raisin_CID"],c["raisin_observation_index"]) for c in candidates if c["source"]==source and c["absolute_decimal_match"]}),
        positive_only_decimal_match_rows=len({(c["raisin_CID"],c["raisin_observation_index"]) for c in candidates if c["source"]==source and c["positive_only_decimal_match"]}),
        original_head_hits=headcounts[source],original_rows=sum(len(x) for x in original[source].values())) for source in ["SMP","DIFFIMG"]}
    inventory_totals={k:sum(x[k] for x in files) for k in ["rows","negative","zero","positive","nonfinite","offseason","nonpositive_error"]}
    result=dict(protocol_sha256=sha(PROTOCOL),physical_files=len(files),inventory_totals=inventory_totals,
        survey_inventory={x["value"]:{k:x[k] for k in ["negative","zero","positive","nonfinite","offseason","nonpositive_error"]} for x in strata if x["level"]=="survey"},
        des_aliases=len(aliases),match_by_source=totals,
        interpretation="Descriptive integrity and time/identity candidates only; no arbitrary scale fit, noise calibration or cosmology inference.")
    (OUT/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    (OUT/"executed_source.py").write_bytes(Path(__file__).read_bytes())
    manifest={p.name:sha(p) for p in [PROTOCOL,OUT/"executed_source.py",OUT/"inventory-files.csv",OUT/"inventory-strata.csv",OUT/"negative-rows.csv",OUT/"match-candidates.csv",OUT/"match-rows.csv",OUT/"result.json"]}
    (OUT/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="match_by_source"},indent=2))
    for source in totals:
        print(source,{k:v for k,v in totals[source].items() if k!="original_head_hits"})


if __name__ == "__main__":
    main()
