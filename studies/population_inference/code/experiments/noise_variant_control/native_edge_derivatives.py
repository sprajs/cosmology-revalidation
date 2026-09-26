"""Output-only fixed-coordinate native derivative probes on seven selected edges."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import subprocess
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
EXE=ROOT/"phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe"
PY=ROOT/"phase2/env-official/bin/python"
EXPORT=ROOT/"scripts/research_2026_09_26/export_flux_objectives.py"
SELECT={"P21_rho000":(1294,7349,20101),
        "P21_rho090":(1294,7349),
        "P21_noisetrue120":(7349,20101)}
PLAN_SHA="583c2b783b66a4fff53b216670d3828f6b33e6875af5ab255b5d5eda57148dab"


def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for part in iter(lambda:f.read(1048576),b""):
            h.update(part)
    return h.hexdigest()


def plan_gate():
    assert sha(HERE/"native-edge-plan.md")==PLAN_SHA
    assert sha(EXE)==json.loads((HERE/"P21_rho000/prepared.json").read_text())["binary_sha256"]
    import pandas as pd
    q=pd.read_csv(HERE/"object_scores/selection-ledger.csv").set_index(["arm","CID"])
    for arm,ids in SELECT.items():
        flag=pd.read_csv(HERE/arm/"mean-edge-ledger.csv")
        for cid in ids:
            assert cid in set(flag.CID) and bool(q.loc[(arm,cid),"selected_archived_quality_classifier"])


def key(t,b,f,e):
    return (float(t),str(b).strip(),float(f),float(e))


def make_masked_input(out):
    from collections import Counter
    ledger=[]
    for arm,ids in SELECT.items():
        version=f"NE_{arm}"
        dest=out/"inputs"/version
        dest.mkdir(parents=True)
        stem=f"PH2_pilot02_{arm}"
        source=ROOT/"phase2/literature/simulations/outputs"/stem
        with fits.open(source/f"{stem}_HEAD.FITS") as hh, fits.open(source/f"{stem}_PHOT.FITS") as pp:
            photo=pp[1].data
            flags=photo["PHOTFLAG"].copy()
            photo["PHOTFLAG"][:]=np.bitwise_or(flags,1<<29)
            index={int(row["SNID"]):row for row in hh[1].data}
            for cid in ids:
                h=index[cid]
                lo=int(h["PTROBS_MIN"])-1;hi=int(h["PTROBS_MAX"])
                pools=defaultdict(deque)
                for i in range(lo,hi):
                    raw=photo[i]
                    pools[key(raw["MJD"],raw["BAND"],raw["FLUXCAL"],raw["FLUXCALERR"])].append(i)
                with np.load(HERE/arm/"objectives"/f"objective_{cid}.npz",allow_pickle=False) as obj:
                    req=Counter(key(*row) for row in zip(obj["MJD"],obj["band"],obj["data_flux"],obj["data_fluxerr"]))
                    assert not any(len(pools[k])>1 and n<len(pools[k]) for k,n in req.items())
                    selected=[]
                    for row in zip(obj["MJD"],obj["band"],obj["data_flux"],obj["data_fluxerr"]):
                        k=key(*row);assert pools[k]
                        rawidx=pools[k].popleft()
                        photo["PHOTFLAG"][rawidx]=int(flags[rawidx]) & ~(1<<29)
                        selected.append(rawidx)
                assert len(selected)==len(set(selected))
                ledger.append(dict(arm=arm,CID=cid,accepted_raw_zero_based=selected,
                                   original_accepted_epochs=len(selected)))
            hh.writeto(dest/f"{version}_HEAD.FITS")
            pp.writeto(dest/f"{version}_PHOT.FITS")
        # The saved HEAD PHOTFILE keyword still names the original basename.
        os.link(dest/f"{version}_PHOT.FITS",dest/f"{stem}_PHOT.FITS")
        (dest/f"{version}.LIST").write_text(f"{version}_HEAD.FITS\n")
        (dest/f"{version}.README").write_text((source/f"{stem}.README").read_text())
    (out/"mask-ledger.json").write_text(json.dumps(ledger,indent=2)+"\n")
    return ledger


def fixed_nml(text,out,arm,seed,prefix):
    text=re.sub(r"(?m)^\s*PRIVATE_DATA_PATH\s*=.*$",f"    PRIVATE_DATA_PATH = '{HERE/'native_edge_full'/'inputs'}'",text)
    text=re.sub(r"(?m)^\s*VERSION_PHOTOMETRY\s*=.*$",f"    VERSION_PHOTOMETRY = 'NE_{arm}'",text)
    text=re.sub(r"(?m)^\s*OPT_SNCID_LIST\s*=.*$","    OPT_SNCID_LIST = 3",text)
    for name,value in (("CUTWIN_TREST","-999.,999."),("FITWIN_TREST","-999.,999."),("DELCHI2_REJECT","1.0E9")):
        text=re.sub(r"(?m)^\s*"+name+r"\s*=.*$",f"    {name} = {value}",text)
    text=text.replace("&SNLCINP","&SNLCINP\n    PHOTFLAG_MSKREJ = 536870912",1)
    text=re.sub(r"(?m)^\s*SNCID_LIST_FILE\s*=.*$",f"    SNCID_LIST_FILE = '{seed}'",text)
    text=re.sub(r"(?m)^\s*TEXTFILE_PREFIX\s*=.*$",f"    TEXTFILE_PREFIX = '{prefix}'",text)
    text=text.replace("&FITINP","&FITINP\n    LFIXPAR_ALL = T",1)
    assert "LFIXPAR_ALL = T" in text and "OPT_MWCOLORLAW = -99" in text
    return text


def job(out,arm,label,coordinate,ids):
    path=out/arm/label
    path.mkdir(parents=True)
    seed=path/"seed.FITRES"
    lines=["VARNAMES: CID PKMJD x0 x1 c"]
    for cid in ids:
        with np.load(HERE/arm/"objectives"/f"objective_{cid}.npz",allow_pickle=False) as f:
            pars=f["parameters_x0_x1_c_t0"].copy()
        if coordinate:
            j,sign,delta=coordinate
            pars[j]+=sign*delta
        x0,x1,c,t0=pars
        assert -5<x1<5
        lines.append(f"SN: {cid} {t0:.17g} {x0:.17g} {x1:.17g} {c:.17g}")
    seed.write_text("\n".join(lines)+"\n")
    baseline=(HERE/arm/"fit.nml").read_text()
    nml=fixed_nml(baseline,out,arm,seed,path/"fit")
    (path/"fit.nml").write_text(nml)
    (path/"cids.txt").write_text("\n".join(map(str,ids))+"\n")
    return path


def run_one(path):
    env=os.environ.copy()
    env.update(SNANA_DIR=str(ROOT/"phase2/official/build/SNANA-audit-v3"),
               SNDATA_ROOT=str(ROOT/"phase2/official/inputs/SNDATA_ROOT"),
               LD_LIBRARY_PATH=str(ROOT/"phase2/official/build/sysroot/usr/lib"),
               OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1",MKL_NUM_THREADS="1")
    with (path/"fit.log").open("x") as log:
        fit=subprocess.run([str(EXE),str(path/"fit.nml")],cwd=path,env=env,stdout=log,stderr=subprocess.STDOUT)
    assert fit.returncode==0,path
    exp=subprocess.run([str(PY),str(EXPORT),"--log",str(path/"fit.log"),"--output",str(path/"objectives"),
                        "--expected-cids",str(path/"cids.txt")],cwd=ROOT,env=env,capture_output=True,text=True)
    (path/"export.log").write_text(exp.stdout+exp.stderr)
    assert exp.returncode==0,(path,exp.stderr)
    return str(path.relative_to(HERE))


def main(half):
    plan_gate()
    out=HERE/("native_edge_half" if half else "native_edge_full")
    out.mkdir(exist_ok=False)
    if half:
        assert (HERE/"native_edge_full/mask-ledger.json").exists()
        (out/"mask-ledger.json").write_bytes((HERE/"native_edge_full/mask-ledger.json").read_bytes())
    else:
        make_masked_input(out)
    factor=.5 if half else 1.
    steps={1:.001*factor,2:.0001*factor,3:.01*factor}
    protocol={"selection":SELECT,"steps_x1_c_t0":steps,"original_free_fit_C_retained_for_projection":True,
              "derivative_only_mask": "PHOTFLAG bit29 freezes exact original accepted raw-PHOT rows; widened phase and disabled clipping",
              "plan_sha256":PLAN_SHA,"source_sha256":sha(Path(__file__)),"binary_sha256":sha(EXE),
              "original_score_manifest_sha256":sha(HERE/"object_scores/manifest.json"),
              "fixed_input_mask_sha256":sha(out/"mask-ledger.json")}
    (out/"protocol.json").write_text(json.dumps(protocol,indent=2)+"\n")
    jobs=[]
    for arm,ids in SELECT.items():
        if not half:
            jobs.append(job(out,arm,"baseline",None,ids))
        for j,h in steps.items():
            for sign in (-1,1):
                jobs.append(job(out,arm,f"p{j}_{'plus' if sign>0 else 'minus'}",(j,sign,h),ids))
    # Each job writes only its private path; archived input and modified masked
    # PHOT files are read-only after preparation. No shared writable scratch.
    with ThreadPoolExecutor(max_workers=2) as pool:
        for completed in pool.map(run_one,jobs):
            print("finished",completed,flush=True)
    result={"jobs":len(jobs),"object_exports":sum(len(SELECT[arm]) for arm in SELECT)*(7 if not half else 6),
            "all_exported":True,"source_sha256":sha(Path(__file__)),"plan_sha256":PLAN_SHA}
    (out/"run-gate.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--half",action="store_true")
    main(ap.parse_args().half)
