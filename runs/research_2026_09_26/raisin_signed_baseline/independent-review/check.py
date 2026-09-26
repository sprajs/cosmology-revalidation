"""Independent scalar-sum and unordered-pair audit of frozen signed baseline."""
import csv
import hashlib
import json
import math
import statistics
from collections import Counter, defaultdict
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
OUT=Path(__file__).resolve().parent
BASE=OUT.parent


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(newline="") as f:return list(csv.DictReader(f))


def raw_rows(path):
    columns=None;result=[]
    for lineno,line in enumerate(path.read_text().splitlines(),1):
        parts=line.split()
        if not parts:continue
        if parts[0]=="VARLIST:":columns=parts[1:]
        if parts[0]!="OBS:":continue
        assert columns and len(parts)-1==len(columns)
        d=dict(zip(columns,parts[1:]))
        b=d.get("FLT",d.get("BAND"))
        if b not in "griz":continue
        f=float(d["FLUXCAL"]);s=float(d["FLUXCALERR"]);t=float(d["MJD"])
        assert all(map(math.isfinite,(f,s,t))) and s>0
        result.append(dict(line=lineno,MJD=t,band=b,flux=f,error=s,flag=d.get("PHOTFLAG","")))
    return result


def close(a,b,label,issues,rtol=1e-11,atol=1e-11):
    if not math.isclose(float(a),float(b),rel_tol=rtol,abs_tol=atol):issues.append((label,a,b))


def pair_bin(lag):
    for lo,hi in [(0,.5),(.5,7),(7,30),(30,180),(180,math.inf)]:
        if lo<=lag<hi:return lo,None if math.isinf(hi) else hi
    raise AssertionError(lag)


def main():
    plan=json.loads((OUT/"review-plan.json").read_text())
    assert all(sha(ROOT/k)==v for k,v in plan["inputs_sha256"].items())
    protocol=json.loads((BASE/"protocol.json").read_text())
    assert sha(BASE/"protocol.json")=="fbaac7dca06161c53b6c180af4281ddbe1bcf98e91b81efdb73a9d1e92085d3d"
    assert all(sha(ROOT/k)==v for k,v in protocol["source_hashes"].items())
    root_result=json.loads((BASE/"result.json").read_text())
    assert root_result["protocol_sha256"]==sha(BASE/"protocol.json")
    assert sha(ROOT/"scripts/research_2026_09_26/raisin_signed_baseline.py")==sha(BASE/"executed-source.py")
    assert all(sha(BASE/k)==v for k,v in root_result["output_sha256"].items())
    cases=read_csv(ROOT/"runs/research_2026_09_26/astra_design/raisin_signed_refit/cohort.csv")
    assert [c["CID"] for c in cases]==protocol["cohort"]
    root_groups={(int(r["cut_days"]),r["CID"],r["band"]):r for r in read_csv(BASE/"baseline-groups.csv")}
    root_pairs={(int(r["cut_days"]),r["CID"],r["band"],float(r["lag_lower"]),r["lag_upper"],r["same_integer_MJD"]):r for r in read_csv(BASE/"pair-summaries.csv")}
    root_flags={(int(r["cut_days"]),r["CID"],r["band"],r["PHOTFLAG"]):r for r in read_csv(BASE/"flag-strata.csv")}
    root_ledger={(int(r["cut_days"]),r["CID"],int(r["source_line"])):r for r in read_csv(BASE/"baseline-source-rows.csv")}
    issues=[];count={};seen_groups=set();seen_pairs=set();seen_flags=set();seen_ledger=set();groups=[]
    for cut in (180,365):
        n_epoch=n_neg=0
        for case in cases:
            cid=case["CID"];peak=float(case["peak_header"])
            selected=[x for x in raw_rows(ROOT/case["raw_path"]) if x["MJD"]<peak-cut]
            n_epoch+=len(selected);n_neg+=sum(x["flux"]<0 for x in selected)
            for row in selected:
                key=(cut,cid,row["line"]);stored=root_ledger.get(key)
                if stored is None:issues.append(("missing source row",key));continue
                seen_ledger.add(key)
                for name,stored_name in [("MJD","MJD"),("flux","flux"),("error","error")]:close(row[name],stored[stored_name],(key,name),issues)
                if row["band"]!=stored["band"] or row["flag"]!=stored["PHOTFLAG"]:issues.append(("source metadata",key))
            for band in "griz":
                rr=[x for x in selected if x["band"]==band];n=len(rr);assert n>2
                w=[1/(x["error"]*x["error"]) for x in rr]
                W=sum(w);mean=sum(a["flux"]*wi for a,wi in zip(rr,w))/W;se=math.sqrt(1/W)
                resid=[(x["flux"]-mean)/x["error"] for x in rr]
                Q=sum(v*v for v in resid);dof=n-1
                key=cut,cid,band;stored=root_groups.get(key)
                if stored is None:issues.append(("missing group",key));continue
                seen_groups.add(key)
                for name,val in [("N",n),("negative",sum(x["flux"]<0 for x in rr)),("baseline_FLUXCAL",mean),("baseline_conditional_diag_SE",se),("baseline_conditional_diag_z",mean/se),("Q_after_intercept",Q),("dof",dof),("Q_per_dof",Q/dof)]:close(val,stored[name],(key,name),issues)
                groups.append(dict(cut=cut,CID=cid,band=band,N=n,negative=sum(x["flux"]<0 for x in rr),baseline=mean,SE=se,Q=Q,dof=dof,Q_per_dof=Q/dof))
                flag_groups=defaultdict(list)
                for i,x in enumerate(rr):flag_groups[x["flag"]].append(i)
                for flag,indices in flag_groups.items():
                    fk=cut,cid,band,flag;stored_flag=root_flags.get(fk)
                    if stored_flag is None:issues.append(("missing flag stratum",fk));continue
                    seen_flags.add(fk)
                    close(len(indices),stored_flag["N"],(fk,"N"),issues)
                    close(sum(rr[i]["flux"]<0 for i in indices),stored_flag["negative"],(fk,"negative"),issues)
                    fover=[rr[i]["flux"]/rr[i]["error"] for i in indices]
                    close(statistics.median(fover),stored_flag["median_flux_over_error"],(fk,"median_flux_over_error"),issues)
                    close(sum(resid[i]**2 for i in indices)/len(indices),stored_flag["mean_projected_residual_squared"],(fk,"mean_projected_residual_squared"),issues)
                pb=defaultdict(list)
                for i in range(n):
                    a=rr[i]
                    for j in range(i+1,n):
                        b=rr[j];lag=abs(a["MJD"]-b["MJD"]);lo,hi=pair_bin(lag)
                        night=math.floor(a["MJD"])==math.floor(b["MJD"])
                        normdiff2=(a["flux"]-b["flux"])**2/(a["error"]**2+b["error"]**2)
                        # For distinct indices, Cov[r_i,r_j] = -1/(sigma_i sigma_j W).
                        excess=resid[i]*resid[j]+1/(a["error"]*b["error"]*W)
                        pb[(lo,"" if hi is None else str(hi),str(night))].append((normdiff2,excess))
                for (lo,hi,night),vals in pb.items():
                    pk=cut,cid,band,lo,hi,night;stored_pair=root_pairs.get(pk)
                    if stored_pair is None:issues.append(("missing pair bin",pk));continue
                    seen_pairs.add(pk)
                    close(len(vals),stored_pair["N_pairs"],(pk,"N_pairs"),issues)
                    close(sum(v[0] for v in vals)/len(vals),stored_pair["mean_normalized_difference_squared"],(pk,"mean_normalized_difference_squared"),issues)
                    close(sum(v[1] for v in vals)/len(vals),stored_pair["mean_projected_product_excess"],(pk,"mean_projected_product_excess"),issues)
        count[cut]=dict(epochs=n_epoch,negative=n_neg,groups=sum(g["cut"]==cut for g in groups))
        original=next(x for x in root_result["summary"] if x["cut_days"]==cut)
        for name,val in [("N_epochs",n_epoch),("N_negative",n_neg),("N_groups",count[cut]["groups"]),
                         ("pooled_Q_per_dof",sum(g["Q"] for g in groups if g["cut"]==cut)/sum(g["dof"] for g in groups if g["cut"]==cut))]:close(val,original[name],(cut,name),issues)
    for label,seen,stored in [("groups",seen_groups,root_groups),("pairs",seen_pairs,root_pairs),("flags",seen_flags,root_flags),("source",seen_ledger,root_ledger)]:
        if seen!=set(stored):issues.append(("key closure",label,len(seen),len(stored),len(set(stored)-seen)))
    assert count=={180:dict(epochs=3476,negative=1736,groups=40),365:dict(epochs=2856,negative=1428,groups=40)},count
    result=dict(status="PASS" if not issues else "FAIL",protocol_sha256=sha(BASE/"protocol.json"),root_source_sha256=sha(BASE/"executed-source.py"),root_output_hashes_pass=True,
        counts=count,keys=dict(groups=len(seen_groups),pairs=len(seen_pairs),flag_strata=len(seen_flags),source_rows=len(seen_ledger)),issues=issues[:100],issue_count=len(issues))
    (OUT/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    with (OUT/"groups.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(groups[0]));w.writeheader();w.writerows(groups)
    print(json.dumps(result,indent=2))
    assert not issues

if __name__=="__main__":main()
