"""Independent read-only parser/check for the historical RAISIN fit batch."""
from collections import Counter, defaultdict
from decimal import Decimal
from pathlib import Path
import csv
import gzip
import hashlib
import json
import math

import numpy as np

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "runs/research_2026_09_26/raisin_historical_native"
OUT = BASE / "independent-review"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def records(path, marker):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as f:
        names = None
        result = []
        for line in f:
            if line.startswith("VARNAMES:"):
                names = line.split()[1:]
            elif line.startswith(marker + ":"):
                words = line.split()[1:]
                assert names is not None and len(words) == len(names), (path, line)
                result.append(dict(zip(names, words)))
    return result


def half_print_unit(s):
    d = Decimal(s)
    return float(Decimal(5).scaleb(d.as_tuple().exponent - 1))


def covariance(row, rounded_upper=False):
    errkeys = ("STRETCHERR", "AVERR", "DLMAGERR")
    covkeys = ((0, 1, "COV_STRETCH_AV"), (0, 2, "COV_STRETCH_DLMAG"),
               (1, 2, "COV_AV_DLMAG"))
    c = np.zeros((3, 3))
    for i, key in enumerate(errkeys):
        x = float(row[key])
        if rounded_upper:
            x += half_print_unit(row[key])
        c[i, i] = x*x
    for i, j, key in covkeys:
        x = float(row[key])
        if rounded_upper:
            x = math.copysign(max(0., abs(x)-half_print_unit(row[key])), x)
        c[i, j] = c[j, i] = x
    return c


def matrix_status(row):
    c = covariance(row)
    u = covariance(row, True)
    eig = np.linalg.eigvalsh(c)
    ueig = np.linalg.eigvalsh(u)
    corr = [abs(c[i,j])/math.sqrt(c[i,i]*c[j,j]) for i,j in ((0,1),(0,2),(1,2))]
    bound = np.zeros((3,3))
    for i,key in enumerate(("STRETCHERR","AVERR","DLMAGERR")):
        x=float(row[key]);h=half_print_unit(row[key]);bound[i,i]=2*x*h+h*h
    for i,j,key in ((0,1,"COV_STRETCH_AV"),(0,2,"COV_STRETCH_DLMAG"),(1,2,"COV_AV_DLMAG")):
        bound[i,j]=bound[j,i]=half_print_unit(row[key])
    # Weyl bound: every rounded matrix changes lambda_min by <= operator norm <= Frobenius norm.
    return float(eig[0]), float(ueig[0]), max(corr), float(np.linalg.norm(bound,ord="fro"))


def main():
    protocol_path = BASE / "fit-protocol-shortprefix.json"
    protocol = json.loads(protocol_path.read_text())
    summary = json.loads((BASE / "summary.json").read_text())
    hash_failures = []
    for listing in (protocol["inputs_sha256"], protocol["parent_inputs_sha256"],
                    summary["verified_output_hashes"]):
        for name, h in listing.items():
            p = ROOT / name
            if not p.exists() or digest(p) != h:
                hash_failures.append(name)
    author_path = ROOT / "runs/research_2026_09_26/astra_design/raisin_signed_refit/author-optical-FITOPT000.FITRES.gz"
    author = {r["CID"]: r for r in records(author_path, "SN")}
    entries, masks, detail = {}, {}, []
    statuses = Counter()
    nplot = 0
    for job in protocol["jobs"]:
        folder = (ROOT / job["nml"]).parent
        x = json.loads((folder / "execution.json").read_text())
        assert x["returncode"] == 0 and x["protocol_sha256"] == digest(protocol_path)
        assert x["log_sha256"] == digest(folder / "fit.log")
        rows = records(folder / "fit.FITRES.TEXT", "SN")
        plots = records(folder / "fit.LCPLOT.TEXT", "OBS")
        if job["stage"] == "full":
            assert sorted(row["CID"] for row in rows) == sorted(protocol["cohort"])
        else:
            assert [row["CID"] for row in rows] == [protocol["cohort"][0]]
        assert all(row["CID"] in protocol["cohort"] for row in plots)
        statuses[(job["stage"], str(job["law"]), "jobs")] += 1
        nplot += len(plots)
        if job["stage"] != "full":
            continue
        for row in rows:
            cid = row["CID"]
            accepted = [r for r in plots if r["CID"] == cid and int(r["DATAFLAG"]) == 1]
            free = 4 if job["timing"] == "free" else 3
            assert len(accepted) == int(row["NDOF"]) + free
            key = (int(job["law"]), job["timing"], job["arm"], float(job["start"]), cid)
            assert key not in entries
            entries[key] = row
            masks[key] = Counter((r["BAND"],r["MJD"],r["FLUXCAL"]) for r in accepted)
            mineig, upper_eig, corr, rounding_bound = matrix_status(row)
            detail.append({"CID":cid,"law":int(job["law"]),"timing":job["timing"],
                           "arm":job["arm"],"start":float(job["start"]),
                           "min_eigenvalue":mineig,"rounding_favorable_min_eigenvalue":upper_eig,
                           "rounding_perturbation_bound":rounding_bound,
                           "max_abs_pair_correlation":corr,"ERRFLAG_FIT":int(row["ERRFLAG_FIT"]),
                           "accepted":len(accepted),
                           "accepted_negative":sum(float(r["FLUXCAL"])<0 for r in accepted)})
    assert len(entries) == 270
    cases = []
    for law in (94,99):
        for timing in ("free","fixed"):
            for arm in ("R","A","B","H"):
                for cid in protocol["cohort"]:
                    keys = [(law,timing,arm,start,cid) for start in (1.,.85,1.15)]
                    if not all(key in entries for key in keys):
                        continue
                    vals = [float(entries[key]["DLMAG"]) for key in keys]
                    cases.append({"CID":cid,"law":law,"timing":timing,"arm":arm,
                                  "spread":max(vals)-min(vals),
                                  "same_mask":masks[keys[0]]==masks[keys[1]]==masks[keys[2]]})
    primary = [c for c in cases if c["law"]==94]
    parent_cases=list(csv.DictReader((BASE/"start-spreads.csv").open()))
    parent_case_index={(int(x["law"]),x["timing"],x["arm"],x["CID"]):x for x in parent_cases}
    assert len(parent_case_index)==len(cases)
    case_differences=[]
    for c in cases:
        pc=parent_case_index[(c["law"],c["timing"],c["arm"],c["CID"])]
        if abs(c["spread"]-float(pc["DLMAG_range"]))>1e-12 or c["same_mask"]!=(pc["same_accepted_mask"]=="True"):
            case_differences.append((c["law"],c["timing"],c["arm"],c["CID"]))
    parent_ledger = list(csv.DictReader((BASE/"fit-ledger.csv").open()))
    parent_index = {(int(r["law"]),r["timing"],r["arm"],float(r["start"]),r["CID"]):r for r in parent_ledger}
    assert len(parent_index)==len(entries)
    ledger_differences=[]
    for key,row in entries.items():
        parent=parent_index[key]
        for name in ("DLMAG","AV","STRETCH","PKMJD","FITCHI2","NDOF","RV","RVERR","PKMJDINI","PKMJDERR"):
            if float(row[name])!=float(parent[name]):
                ledger_differences.append((key,name))
    nominal = []
    for cid in protocol["cohort"]:
        old=entries[(94,"free","R",1.,cid)]
        observed = float(old["DLMAG"])
        nominal.append({"CID":cid,"historical":observed,"author":float(author[cid]["DLMAG"]),
                        "delta":observed-float(author[cid]["DLMAG"]),
                        "delta_av":float(old["AV"])-float(author[cid]["AV"]),
                        "delta_stretch":float(old["STRETCH"])-float(author[cid]["STRETCH"]),
                        "delta_peak":float(old["PKMJD"])-float(author[cid]["PKMJD"]),
                        "delta_data_chi2":float(old["FITCHI2"])-float(author[cid]["FITCHI2"]),
                        "delta_ndof":int(old["NDOF"])-int(author[cid]["NDOF"])})
    parent_baseline={r["CID"]:r for r in csv.DictReader((BASE/"author-and-version-baseline.csv").open())}
    baseline_differences=[]
    for item in nominal:
        cid=item["CID"]
        for law in (94,99):
            r=entries[(law,"free","R",1.,cid)]
            for field in ("DLMAG","AV","STRETCH","PKMJD","FITCHI2","NDOF"):
                computed=float(r[field])-float(author[cid][field])
                parent=float(parent_baseline[cid][f"old{law}_minus_author_{field}"])
                if abs(computed-parent)>1e-12:
                    baseline_differences.append((cid,law,field))
    pair_differences=[]
    for p in csv.DictReader((BASE/"paired-start-responses.csv").open()):
        left,right=p["contrast"].split("-")
        for field in ("DLMAG","AV","STRETCH","PKMJD"):
            l=float(entries[(94,p["timing"],left,float(p["start"]),p["CID"])][field])
            r=float(entries[(94,p["timing"],right,float(p["start"]),p["CID"])][field])
            if abs((l-r)-float(p[field]))>1e-12:
                pair_differences.append((p["CID"],p["timing"],p["start"],p["contrast"],field))
    failures = [x for x in detail if x["min_eigenvalue"]<=0]
    robust = [x for x in failures if x["rounding_favorable_min_eigenvalue"]<=0]
    rigorous_rounding_failures=[x for x in failures if -x["min_eigenvalue"]>x["rounding_perturbation_bound"]]
    output = {"protocol_sha256":digest(protocol_path),
              "summary_sha256":digest(BASE/"summary.json"),
              "source_sha256":digest(Path(__file__)),
              "hash_files_checked":sum(len(x) for x in (protocol["inputs_sha256"],protocol["parent_inputs_sha256"],summary["verified_output_hashes"])),
              "hash_failures":hash_failures,"job_counts":{str(k):v for k,v in statuses.items()},
              "parent_ledger_numeric_differences":ledger_differences,
              "parent_case_differences":case_differences,
              "parent_baseline_differences":baseline_differences,
              "parent_pair_differences":pair_differences,
              "full_records":len(entries),"plot_rows_all_jobs":nplot,
              "errflag_zero":sum(x["ERRFLAG_FIT"]==0 for x in detail),
              "rv_value_1p518_count":sum(float(r["RV"])==1.518 for r in entries.values()),
              "rv_error_zero_count":sum(float(r["RVERR"])==0 for r in entries.values()),
              "fixed_peak_equal_seed_count":sum(float(r["PKMJD"])==float(r["PKMJDINI"]) for k,r in entries.items() if k[1]=="fixed"),
              "cases":len(primary),"same_masks":sum(x["same_mask"] for x in primary),
              "spreads_gt_0p001":sum(x["spread"]>.001 for x in primary),
              "spreads_gt_0p01":sum(x["spread"]>.01 for x in primary),
              "max_spread":max(x["spread"] for x in primary),
              "max_spread_case":max(primary,key=lambda x:x["spread"]),
              "covariance_nonpositive":len(failures),"covariance_nonpositive_beyond_print_rounding":len(robust),
              "covariance_nonpositive_proven_by_weyl_rounding_bound":len(rigorous_rounding_failures),
              "covariance_min_eigenvalue_min":min(x["min_eigenvalue"] for x in failures),
              "max_abs_pair_correlation":max(x["max_abs_pair_correlation"] for x in detail),
              "author_abs_max":max(abs(x["delta"]) for x in nominal),
              "author_abs_median":float(np.median([abs(x["delta"]) for x in nominal])),
              "nominal_deltas":nominal,
              "engineering_fitres_byte_equal":(BASE/"fits-shortprefix/engineering/94/free/R/1.0/fit.FITRES.TEXT").read_bytes()==(BASE/"fits-shortprefix/engineering_copy/94/free/Rcopy/1.0/fit.FITRES.TEXT").read_bytes(),
              "engineering_fitres_rows_equal":records(BASE/"fits-shortprefix/engineering/94/free/R/1.0/fit.FITRES.TEXT","SN")==records(BASE/"fits-shortprefix/engineering_copy/94/free/Rcopy/1.0/fit.FITRES.TEXT","SN"),
              "engineering_lcplot_byte_equal":(BASE/"fits-shortprefix/engineering/94/free/R/1.0/fit.LCPLOT.TEXT").read_bytes()==(BASE/"fits-shortprefix/engineering_copy/94/free/Rcopy/1.0/fit.LCPLOT.TEXT").read_bytes()}
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"result.json").write_text(json.dumps(output,indent=2)+"\n")
    with (OUT/"covariance-diagnostics.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=detail[0].keys());w.writeheader();w.writerows(detail)
    with (OUT/"independent-case-spreads.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cases[0].keys());w.writeheader();w.writerows(cases)
    print(json.dumps({k:v for k,v in output.items() if k!="nominal_deltas"},indent=2))


if __name__=="__main__":
    main()
