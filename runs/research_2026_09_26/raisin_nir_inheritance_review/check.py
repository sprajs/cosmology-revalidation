"""Read-only exact-CID provenance join for ten released DES NIR fits."""
from pathlib import Path
import csv
import hashlib
import json
import re
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
REL = ROOT / "sources/repos/djones1040__RAISIN_DataRelease@a383c4b"
SRC = ROOT / "sources/repos/RickKessler__SNANA@v11_04k/src"
SIM = ROOT / "runs/research_2026_09_26/raisin_sign_source"
COHORT = ROOT / "runs/research_2026_09_26/raisin_historical_native/fit-protocol-shortprefix.json"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def fitres(path):
    names = None
    out = {}
    for line in path.read_text().splitlines():
        if line.startswith("VARNAMES:"):
            names = line.split()[1:]
        elif line.startswith("SN:"):
            words = line.split()[1:]
            assert names and len(words)==len(names), (path,line)
            row = dict(zip(names,words))
            assert row["CID"] not in out
            out[row["CID"]] = row
    return out


def simlib(path):
    result = {}
    libid = None
    peak = None
    for no,line in enumerate(path.read_text().splitlines(),1):
        if line.startswith("LIBID:"):
            libid=line.split()[1];peak=None
        if line.startswith("REDSHIFT:") and "PEAKMJD:" in line:
            peak=line.split("PEAKMJD:",1)[1].split()[0]
        if line.startswith("FIELD:"):
            name=line.split()[1]
            result.setdefault(name,[]).append({"libid":libid,"peak":peak,"line":no})
    return result


def main():
    cohort=json.loads(COHORT.read_text())["cohort"]
    assert len(cohort)==len(set(cohort))==10
    dirs={b:REL/f"distances/w/{b}_dist/RAISIN_combined_FITOPT000.FITRES"
          for b in ("nir","optical","opticalnir")}
    tables={b:fitres(p) for b,p in dirs.items()}
    optical_shift_paths={n:REL/f"distances/w/nir_dist/RAISIN_combined_FITOPT{n:03d}.FITRES" for n in range(19,23)}
    optical_shift_tables={n:fitres(p) for n,p in optical_shift_paths.items()}
    simpath=SIM/"sim/simlibs/DES_RAISIN.simlib"
    sims=simlib(simpath)
    author_nir_path=OUT/"author/output/fit_nir/DES_RAISIN.FITRES.TEXT"
    author_combined_path=OUT/"author/DES_RAISIN_optnir.FITRES.TEXT"
    author_nir=fitres(author_nir_path)
    author_combined=fitres(author_combined_path)
    t0path=OUT/"author/data/raisin_t0.txt"
    t0={}
    for line in t0path.read_text().splitlines():
        words=line.split()
        if words and words[0].startswith("DES16"):
            assert words[0] not in t0
            t0[words[0]]=words[1]
    inputs=[COHORT,simpath,SRC/"snana.car",SRC/"snlc_fit.car",REL/"README.md",
            SIM/"raisin_cosmo/genSimlib.py",SIM/"sim/inputs/DES/sim_DES_SNOOPY.input",
            SIM/"fit/sim/DES_RAISIN_optnir.nml",REL/"lcfitting/REFAC_DES_RAISIN_nir_sys.nml",
            REL/"lcfitting/REFAC_DES_RAISIN_optical_sys.nml",
            REL/"lcfitting/REFAC_DES_RAISIN_opticalnir_sys.nml",*dirs.values(),
            *optical_shift_paths.values(),author_nir_path,author_combined_path,t0path,OUT/"author-acquisition.json",
            OUT/"source-v11_03c/acquisition.json",OUT/"source-v11_03c/snana.car",OUT/"source-v11_03c/snlc_fit.car"]
    output=[];problems=[]
    for cid in cohort:
        phot=REL/f"photometry/RAISIN/DES_RAISIN/{cid}.snana.dat"
        inputs.append(phot)
        peaks=[line.split(":",1)[1].split()[0] for line in phot.read_text().splitlines() if line.startswith("PEAKMJD:")]
        row={"CID":cid,"phot_peak_entries":"|".join(peaks),"phot_peak_entry_count":len(peaks),
             "simlib_entries":json.dumps(sims.get(cid,[]),sort_keys=True),
             "author_t0":t0.get(cid,"")}
        row["phot_peak_expected_float32_print"]=(f"{np.float32(float(peaks[0])):.4f}" if peaks else "")
        if len(set(peaks))!=1 or not peaks:problems.append({"cid":cid,"problem":"missing/conflicting photometry peak entries","values":peaks})
        if len(sims.get(cid,[]))!=1:problems.append({"cid":cid,"problem":"SIMLIB field join not unique","matches":sims.get(cid,[])})
        if cid not in t0:problems.append({"cid":cid,"problem":"missing author t0"})
        elif len(sims.get(cid,[]))==1 and sims[cid][0]["peak"]!=f"{float(t0[cid]):.1f}":
            problems.append({"cid":cid,"problem":"SIMLIB peak differs from author t0 formatted to one decimal"})
        for name,table in (("author_nir",author_nir),("author_combined",author_combined)):
            if cid not in table:
                problems.append({"cid":cid,"problem":f"missing {name} raw FITRES"});continue
            for key in ("PKMJDINI","PKMJD","PKMJDERR","STRETCH","STRETCHERR","AV","AVERR","NDOF","ERRFLAG_FIT"):
                row[f"{name}_{key}"]=table[cid][key]
        if cid in author_nir and peaks:
            if row["author_nir_PKMJDINI"]!=row["phot_peak_expected_float32_print"] or row["author_nir_PKMJD"]!=row["phot_peak_expected_float32_print"]:
                problems.append({"cid":cid,"problem":"author raw NIR fixed peak differs from float32 header print"})
        for branch in dirs:
            if cid not in tables[branch]:
                problems.append({"cid":cid,"problem":f"missing {branch} FITRES"});continue
            fit=tables[branch][cid]
            for key in ("PKMJD","PKMJDERR","STRETCH","STRETCHERR","AV","AVERR","RV","RVERR","NDOF","CUTFLAG_SNANA"):
                row[f"{branch}_{key}"]=fit[key]
        if cid in tables["nir"] and peaks:
            row["nir_minus_header_day"]=float(row["nir_PKMJD"])-float(peaks[0])
            if any(float(row[f"nir_{key}"])!=want for key,want in (("PKMJDERR",0.),("STRETCH",1.),("STRETCHERR",0.),("AV",0.),("AVERR",0.))):
                problems.append({"cid":cid,"problem":"NIR fixed-parameter output mismatch"})
        output.append(row)
    OUT.mkdir(exist_ok=True)
    with (OUT/"des16-peak-inheritance.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=output[0].keys());writer.writeheader();writer.writerows(output)
    result={"plan_sha256":sha(OUT/"review-plan.json"),"script_sha256":sha(Path(__file__)),
            "source_hashes":{str(p.relative_to(ROOT)):sha(p) for p in inputs},
            "n_cid":len(output),"n_problems":len(problems),"problems":problems,
            "max_abs_nir_minus_header_day":max(abs(r["nir_minus_header_day"]) for r in output),
            "nir_fixed_parameter_rows":sum(float(r["nir_PKMJDERR"])==0 and float(r["nir_STRETCH"])==1 and float(r["nir_STRETCHERR"])==0 and float(r["nir_AV"])==0 and float(r["nir_AVERR"])==0 for r in output),
            "simlib_peak_equal_header_count":sum(len(sims.get(r["CID"],[]))==1 and float(sims[r["CID"]][0]["peak"])==float(r["phot_peak_entries"].split("|")[0]) for r in output),
            "author_t0_simlib_closure_count":sum(len(sims.get(r["CID"],[]))==1 and sims[r["CID"]][0]["peak"]==f"{float(r['author_t0']):.1f}" for r in output),
            "author_raw_nir_header_float32_closure_count":sum(r["author_nir_PKMJDINI"]==r["phot_peak_expected_float32_print"] and r["author_nir_PKMJD"]==r["phot_peak_expected_float32_print"] for r in output),
            "author_raw_nir_to_release_peak_equal_count":sum(float(r["author_nir_PKMJD"])==float(r["nir_PKMJD"]) for r in output),
            "released_nir_des_optical_zp_fitopts_equal_nominal_count":{
                str(n):sum(cid in optical_shift_tables[n] and all(optical_shift_tables[n][cid][k]==tables["nir"][cid][k] for k in ("PKMJD","STRETCH","AV","DLMAG","DLMAG_biascor","mures","MASS_CORR")) for cid in cohort)
                for n in optical_shift_paths},
            "ledger_sha256":sha(OUT/"des16-peak-inheritance.csv")}
    (OUT/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k!="source_hashes"},indent=2))


if __name__=="__main__":main()
