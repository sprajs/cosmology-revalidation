"""Independent native FITRES/LCPLOT/input/multistart consistency review."""
import bisect
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parents[3]
OUT=Path(__file__).resolve().parent
NATIVE=ROOT/"runs/research_2026_09_26/astra_design/raisin_signed_refit"


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def csvrows(path):
    with path.open(newline="") as f:return list(csv.DictReader(f))


def table(path,prefix):
    columns=None;rows=[]
    for line in path.read_text().splitlines():
        x=line.split()
        if not x:continue
        if x[0]=="VARNAMES:":columns=x[1:]
        elif x[0]==prefix:
            assert columns and len(x)-1==len(columns),(path,len(x),len(columns))
            rows.append(dict(zip(columns,x[1:])))
    return rows


def input_observations(path):
    fields=None;out=[]
    for ln,line in enumerate(path.read_text().splitlines(),1):
        x=line.split()
        if not x:continue
        if x[0]=="VARLIST:":fields=x[1:]
        elif x[0]=="OBS:":
            assert len(x)-1==len(fields)
            d=dict(zip(fields,x[1:]))
            if d["FLT"] in "griz":out.append(dict(line=ln,band=d["FLT"],mjd=float(d["MJD"]),flux=float(d["FLUXCAL"]),error=float(d["FLUXCALERR"])))
    return out


def close(a,b,issues,label,rtol=1e-10,atol=1e-10):
    if not math.isclose(float(a),float(b),rel_tol=rtol,abs_tol=atol):issues.append((label,float(a),float(b)))


def match_native(rows,source,issues,tag):
    # Match all observed LCPLOT rows, accepted or rejected, by band/date and
    # printed flux/error. Exact count and one-to-one matching are required.
    if len(rows)!=len(source):
        issues.append((tag,"data row count",len(rows),len(source)))
        return dict(matched=0,negative=0,maxdt=None)
    groups=defaultdict(list)
    for i,r in enumerate(source):groups[r["band"]].append((r["mjd"],i))
    for v in groups.values():v.sort()
    used=set();maxdt=0.;ambig=0
    for y in rows:
        band=y["BAND"];mjd=float(y["MJD"]);flux=float(y["FLUXCAL"]);error=float(y["FLUXCAL_ERR"])
        candidates=groups[band];dates=[z[0] for z in candidates]
        lo=bisect.bisect_left(dates,mjd-.0021);hi=bisect.bisect_right(dates,mjd+.0021)
        possible=[]
        for t,i in candidates[lo:hi]:
            r=source[i]
            ftol=max(1e-8,abs(r["flux"])*5.1e-5)
            etol=max(1e-8,abs(r["error"])*5.1e-5)
            if abs(r["flux"]-flux)<ftol and abs(r["error"]-error)<etol and i not in used:
                possible.append((abs(t-mjd),i))
        if len(possible)>1:ambig+=1
        if not possible:
            issues.append((tag,"unmatched LCPLOT",y["CID"],band,mjd,flux,error))
            continue
        dt,i=min(possible);used.add(i);maxdt=max(maxdt,dt)
        if (flux<0)!=(source[i]["flux"]<0):issues.append((tag,"sign mismatch",i))
    if len(used)!=len(source):issues.append((tag,"source rows unmatched",len(source)-len(used)))
    return dict(matched=len(used),ambiguous_candidates=ambig,maxdt=maxdt)


def main():
    plan=json.loads((OUT/"review-plan.json").read_text())
    assert all(sha(ROOT/k)==h for k,h in plan["inputs_sha256"].items())
    manifest=json.loads((NATIVE/"input-manifest.json").read_text())
    manifest_hash=sha(NATIVE/"input-manifest.json")
    assert manifest_hash=="35a9c9e375898cf644a2d354d5c106a3abdfd1ec759a7ce7d0a23557d38ec80a"
    assert all(sha(ROOT/k)==h for k,h in manifest["inputs_sha256"].items())
    stability=json.loads((NATIVE/"stability-protocol.json").read_text())
    assert stability["primary_manifest_sha256"]==manifest_hash
    assert sha(NATIVE/"stability_controls.py")==stability["source_sha256"]
    assert all(sha(ROOT/j["nml"])==j["sha256"] for j in stability["jobs"])
    cohort=csvrows(NATIVE/"cohort.csv");cids=[x["CID"] for x in cohort]
    assert cids==manifest["cohort"]
    stored={ (r["CID"],r["timing"],r["arm"]):r for r in csvrows(NATIVE/"multistart-comparison.csv")}
    stored_pairs={ (r["CID"],r["timing"],r["comparison"]):r for r in csvrows(NATIVE/"paired-start-sensitivity.csv")}
    assert len(stored)==80 and len(stored_pairs)==40
    inputs={(arm,cid):input_observations(NATIVE/"data"/f"RSR_{arm}"/f"{cid}.snana.dat") for arm in "RABH" for cid in cids}
    counts={arm:sum(len(inputs[arm,cid]) for cid in cids) for arm in "RABH"}
    assert counts=={"R":3135,"A":3135,"B":5378,"H":5378},counts
    # Preserve exact release rows in R, and A/B/H multiplicities from generated text.
    assert all(len(inputs["B",cid])==len(inputs["A",cid])+int(next(x for x in cohort if x["CID"]==cid)["author_nonpositive"]) for cid in cids)
    allfit={};allmask={};checks=[];issues=[]
    for start in (1.0,.85,1.15):
        for timing in ("free","fixed"):
            for arm in "RABH":
                base=NATIVE/"fits"/("full" if start==1.0 else "multistart")
                path=base/timing/arm if start==1.0 else base/str(start)/timing/arm
                execution=json.loads((path/"execution.json").read_text())
                if execution["returncode"]!=0:issues.append(("nonzero native returncode",str(path)))
                if sha(path/"fit.log")!=execution["log_sha256"]:issues.append(("native log hash",str(path)))
                reference=(manifest_hash if start==1.0 else sha(NATIVE/"stability-protocol.json"))
                refkey="manifest_sha256" if start==1.0 else "protocol_sha256"
                if execution[refkey]!=reference:issues.append(("execution protocol link",str(path)))
                fitpath=path/"fit.FITRES.TEXT"
                if not fitpath.exists():fitpath=path/"fit"
                fits=table(fitpath,"SN:");plots=table(path/"fit.LCPLOT.TEXT","OBS:")
                if len(fits)!=10:issues.append(("FITRES row count",str(path),len(fits)))
                by_plot=defaultdict(list)
                for x in plots:
                    if x["DATAFLAG"]!="0":by_plot[x["CID"]].append(x)
                cov={}
                for line in (path/"fit.log").read_text().splitlines():
                    z=line.split()
                    if z and z[0]=="PHASE2_HESSIAN:":
                        values=np.array([float(v) for v in z[2:]])
                        assert len(values)==20
                        cov[z[1]]=values[-16:].reshape(4,4)
                for f in fits:
                    cid=f["CID"];key=start,timing,arm,cid
                    if cid not in cids:issues.append(("unexpected CID",key))
                    if int(f["ERRFLAG_FIT"])!=0:issues.append(("ERRFLAG_FIT",key,f["ERRFLAG_FIT"]))
                    close(f["RV"],1.518,issues,(key,"RV"),atol=1e-5)
                    close(f["RVERR"],0,issues,(key,"RVERR"))
                    peak=float(f["PKMJD"]);ini=float(f["PKMJDINI"])
                    if timing=="fixed":
                        close(peak,ini,issues,(key,"fixed peak"),atol=1e-5)
                        close(f["PKMJDERR"],0,issues,(key,"fixed peak error"))
                    elif float(f["PKMJDERR"])<=0:issues.append(("free peak error",key))
                    if not (.7<=float(f["STRETCH"])<=1.3):issues.append(("stretch grid",key))
                    matrix=cov.get(cid)
                    if matrix is None:issues.append(("missing covariance",key))
                    else:
                        M=(matrix+matrix.T)/2
                        minimum=np.linalg.eigvalsh(M[:3,:3] if timing=="fixed" else M).min()
                        if minimum<=0:issues.append(("nonpositive native covariance",key,float(minimum)))
                    data=by_plot[cid];src=inputs[arm,cid]
                    closure=match_native(data,src,issues,key)
                    accepted=[y for y in data if y["DATAFLAG"]=="1"]
                    nfree=4 if timing=="free" else 3
                    close(len(accepted),float(f["NDOF"])+nfree,issues,(key,"accepted NDOF"))
                    mask=sorted((y["BAND"],y["MJD"],y["FLUXCAL"],y["FLUXCAL_ERR"]) for y in accepted)
                    allmask[key]=mask;allfit[key]=f
                    checks.append(dict(start=start,timing=timing,arm=arm,CID=cid,ERRFLAG_FIT=int(f["ERRFLAG_FIT"]),N_input_optical=len(src),N_data_export=len(data),N_accepted=len(accepted),N_accepted_negative=sum(float(y["FLUXCAL"])<0 for y in accepted),N_source_negative=sum(x["flux"]<0 for x in src),match_count=closure["matched"],max_date_error_day=closure["maxdt"],covariance_min_eigenvalue=float(minimum) if matrix is not None else "",DLMAG=float(f["DLMAG"]),AV=float(f["AV"]),STRETCH=float(f["STRETCH"]),PKMJD=peak,PKMJDINI=ini,RV=float(f["RV"]),FITCHI2=float(f["FITCHI2"])))
    if len(allfit)!=240:issues.append(("fit total",len(allfit)))
    summary=[];different=[]
    for cid in cids:
        for timing in ("free","fixed"):
            for arm in "RABH":
                fs=[allfit[s,timing,arm,cid] for s in (1.0,.85,1.15)]
                masks=[allmask[s,timing,arm,cid] for s in (1.0,.85,1.15)]
                same=masks[0]==masks[1]==masks[2]
                if not same:
                    different.append(dict(CID=cid,timing=timing,arm=arm,accepted_counts=[len(m) for m in masks],start1_vs_085=list((Counter(masks[0])-Counter(masks[1])).elements()),start085_vs_1=list((Counter(masks[1])-Counter(masks[0])).elements())))
                key=cid,timing,arm;expected=stored[key]
                if same!=(expected["same_accepted_mask"]=="True"):issues.append(("mask summary",key))
                changes={name:np.ptp([float(f[name]) for f in fs]) for name in ("DLMAG","AV","STRETCH","PKMJD")}
                chi=[float(f["FITCHI2"]) for f in fs]
                if not same and expected["comparable_objective"]!="False":issues.append(("different mask objective labelled comparable",key))
                for name,val in [(n+"_range",changes[n]) for n in changes]+[("objective_range",max(chi)-min(chi)),("objective_improvement",chi[0]-min(chi))]:close(val,expected[name],issues,(key,name))
                best=(1.0,.85,1.15)[int(np.argmin(chi))]
                close(best,expected["best_start_by_objective"],issues,(key,"best_by_data_FITCHI2"))
                summary.append(dict(CID=cid,timing=timing,arm=arm,same_mask=same,DLMAG_range=float(changes["DLMAG"]),AV_range=float(changes["AV"]),STRETCH_range=float(changes["STRETCH"]),PKMJD_range=float(changes["PKMJD"]),FITCHI2_range=max(chi)-min(chi),FITCHI2_improvement_initial=chi[0]-min(chi)))
    if len(summary)!=80:issues.append(("case count",len(summary)))
    for cid in cids:
        for timing in ("free","fixed"):
            for arm0,arm1,label in (("A","B","B_minus_A"),("R","H","H_minus_R")):
                vals=[float(allfit[s,timing,arm1,cid]["DLMAG"])-float(allfit[s,timing,arm0,cid]["DLMAG"]) for s in (1.0,.85,1.15)]
                expected=stored_pairs[cid,timing,label]
                for field,val in [("initial_delta",vals[0]),("start085_delta",vals[1]),("start115_delta",vals[2]),("minimum_delta",min(vals)),("maximum_delta",max(vals)),("range_delta",max(vals)-min(vals))]:close(val,expected[field],issues,(cid,timing,label,field))
    maxrow=max(summary,key=lambda r:r["DLMAG_range"])
    assert (maxrow["CID"],maxrow["timing"],maxrow["arm"])==("DES16C1cim","free","H")
    result=dict(status="PASS" if not issues else "FAIL",fits=len(allfit),comparisons=len(summary),same_masks=sum(r["same_mask"] for r in summary),different_masks=different,
        DLMAG_range_over_0p001=sum(r["DLMAG_range"]>.001 for r in summary),DLMAG_range_over_0p01=sum(r["DLMAG_range"]>.01 for r in summary),max_DLMAG_range=maxrow,
        max_same_mask_data_FITCHI2_improvement=max(r["FITCHI2_improvement_initial"] for r in summary if r["same_mask"]),
        input_optical_counts=counts,negative_exported_all_arms=sum(r["N_source_negative"] for r in checks),accepted_negative_all_fits=sum(r["N_accepted_negative"] for r in checks),
        source_hashes_closed=True,issues=issues[:100],issue_count=len(issues))
    (OUT/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    for name,rows in (("fit-checks.csv",checks),("comparison-checks.csv",summary)):
        with (OUT/name).open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    print(json.dumps({k:v for k,v in result.items() if k not in ("different_masks","issues")},indent=2))
    assert not issues

if __name__=="__main__":main()
