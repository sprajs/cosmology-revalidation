"""Independent FITS/SIMLIB/RAISIN reparse for the fixed 17 DES aliases."""
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

from astropy.io import fits

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
SIGN = ROOT / "runs/research_2026_09_26/raisin_flux_sign_audit"
SOURCE = ROOT / "runs/research_2026_09_26/raisin_sign_source"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def decode(x):
    return x.decode().strip() if isinstance(x, bytes) else str(x).strip()


def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_product(path):
    header = {}
    rows = []
    for line_no, line in enumerate(path.read_text().splitlines(), 1):
        if line.startswith("OBS:"):
            parts = line.split()
            rows.append((parts[2], float(parts[1]), line_no))
        elif ":" in line and not line.startswith("#"):
            k, v = line.split(":", 1)
            header[k.strip()] = v.strip()
    return header, rows


def main():
    protocol = json.loads((OUT / "protocol.json").read_text())
    assert all(sha(ROOT / p) == h for p, h in protocol["inputs_sha256"].items())
    aliases = json.loads((SIGN / "protocol.json").read_text())["des_aliases"]
    inventory = {x["CID"]: x for x in json.loads((ROOT / "runs/research_2026_09_26/raisin_differential/photometry-inventory.json").read_text())}
    product = {}
    for cid in aliases:
        product[cid] = read_product(ROOT / inventory[cid]["path"])

    sim = defaultdict(list)
    lib = field = None
    total_sim = 0
    libs = set()
    for line_no, line in enumerate((SOURCE / "sim/simlibs/DES_RAISIN.simlib").read_text().splitlines(), 1):
        z = line.split()
        if not z:
            continue
        if z[0] == "LIBID:":
            lib = int(z[1]); field = None; libs.add(lib)
        elif z[0] == "FIELD:":
            field = z[1]
        elif z[0] == "S:":
            assert field and lib is not None
            sim[(field, z[3])].append((float(z[1]), lib, line_no))
            total_sim += 1
    assert total_sim == 6877 and len(libs) == 23

    # HEAD is indexed by original CID, without consulting the prior row ledger.
    heads = defaultdict(list)
    with fits.open(ROOT / "data/des-diffimg/DES-SN5YR_DIFFIMG_HEAD.FITS.gz", memmap=False) as hdu:
        data = hdu[1].data
        for i in range(len(data)):
            cid = decode(data["SNID"][i])
            if cid in aliases.values():
                heads[cid].append((i, int(data["PTROBS_MIN"][i]), int(data["PTROBS_MAX"][i])))
    assert len(heads) == 17 and all(len(v) == 1 for v in heads.values())

    rows = []
    with fits.open(ROOT / "data/des-diffimg/DES-SN5YR_DIFFIMG_PHOT.FITS.gz", memmap=False) as hdu:
        data = hdu[1].data
        for cid in sorted(aliases):
            original_id = aliases[cid]
            head_index, start, end = heads[original_id][0]
            header, product_rows = product[cid]
            peak = float(header["PEAKMJD"].split()[0])
            zhel = float(header["REDSHIFT_HELIO"].split()[0])
            for global_row in range(start - 1, end):
                mjd = float(data["MJD"][global_row])
                if mjd < 0:
                    continue
                band = decode(data["BAND"][global_row])
                if band not in "griz":
                    continue
                flux = float(data["FLUXCAL"][global_row])
                error = float(data["FLUXCALERR"][global_row])
                phase = (mjd - peak) / (1 + zhel)
                pmatches = [(t,ln) for b,t,ln in product_rows if b == band and abs(t - mjd) <= .00055]
                smatches = [(t,li,ln) for t,li,ln in sim[(cid,band)] if abs(t-mjd) <= .00055]
                rows.append(dict(CID=cid, original_ID=original_id, HEAD_index=head_index, FITS_global_row=global_row,
                    MJD=mjd, band=band, flux=flux, error=error, PHOTFLAG=int(data["PHOTFLAG"][global_row]), IMGNUM=int(data["IMGNUM"][global_row]),
                    sign="negative" if flux<0 else "positive" if flux>0 else "zero", rest_phase=phase,
                    fit_window=-7<=phase<=45, m20_window=-10<=phase<=40, offseason=abs(mjd-peak)>180,
                    product_matches=len(pmatches), product_lines=";".join(str(x[1]) for x in pmatches),
                    simlib_matches=len(smatches), simlib_lines=";".join(str(x[2]) for x in smatches),
                    simlib_LIBIDs=";".join(str(x[1]) for x in smatches)))
    assert len(rows) == 9843
    write_csv(OUT / "independent-cadence-rows.csv", rows)

    root_rows = read_csv(SOURCE / "cadence-matches.csv")
    by_key = {(x["raisin_CID"], int(x["original_fits_global_row"])): x for x in root_rows}
    assert len(by_key) == len(root_rows) == len(rows)
    mismatches = []
    for row in rows:
        key = row["CID"], row["FITS_global_row"]
        expected = by_key[key]
        equal = (row["band"] == expected["band"] and abs(row["MJD"]-float(expected["MJD"]))<1e-7
            and abs(row["flux"]-float(expected["original_flux"]))<1e-5 and row["PHOTFLAG"]==int(expected["PHOTFLAG"])
            and row["IMGNUM"]==int(expected["IMGNUM"])
            and row["product_matches"] == int(expected["raisin_match_count"])
            and row["simlib_matches"] == int(expected["simlib_matches"])
            and row["simlib_lines"] == expected["simlib_line"]
            and row["simlib_LIBIDs"] == expected["simlib_LIBID"]
            and row["fit_window"] == (expected["original_optical_fit_window"]=="True")
            and row["m20_window"] == (expected["m20_window"]=="True")
            and row["offseason"] == (expected["offseason"]=="True"))
        if not equal:
            mismatches.append(key)
    strata = {}
    for name, fn in (("all",lambda x:True),("fit_window",lambda x:x["fit_window"]),("m20_window",lambda x:x["m20_window"]),("offseason",lambda x:x["offseason"])):
        sub=[x for x in rows if fn(x)]
        strata[name]={str(k):v for k,v in Counter((x["sign"],x["product_matches"]>0,x["simlib_matches"]>0) for x in sub).items()}
    result = dict(status="PASS" if not mismatches else "FAIL", rows=len(rows), simlib_rows=total_sim, simlib_LIBIDs=len(libs), aliases=len(aliases), mismatches=mismatches,
        membership_discordance=sum((x["product_matches"]>0)!=(x["simlib_matches"]>0) for x in rows),
        all_negative=sum(x["sign"]=="negative" for x in rows), all_negative_simlib=sum(x["sign"]=="negative" and x["simlib_matches"]>0 for x in rows),
        fit_negative=sum(x["fit_window"] and x["sign"]=="negative" for x in rows), fit_negative_simlib=sum(x["fit_window"] and x["sign"]=="negative" and x["simlib_matches"]>0 for x in rows),
        fit_positive=sum(x["fit_window"] and x["sign"]=="positive" for x in rows), fit_positive_simlib=sum(x["fit_window"] and x["sign"]=="positive" and x["simlib_matches"]>0 for x in rows),
        ambiguous=[dict(CID=x["CID"], FITS_global_row=x["FITS_global_row"], product_matches=x["product_matches"], simlib_matches=x["simlib_matches"], simlib_lines=x["simlib_lines"]) for x in rows if x["simlib_matches"]>1],
        strata=strata, output_sha256=sha(OUT/"independent-cadence-rows.csv"))
    (OUT/"independent-cadence-result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k != "strata"},indent=2))
    assert not mismatches and result["membership_discordance"]==0

if __name__=="__main__":
    main()
