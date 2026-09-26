#!/usr/bin/env python3
"""Audit every original R19 global draw against the stated age-model prior.

An explicit follow-up to the preserved first pass. Do not silently drop invalid
archive draws: retain full counts, first examples, and original full-vs-valid
quantiles. FSPS-consistent ages are reported conditional on valid draws.
"""
from pathlib import Path
import hashlib
import io
import json
import re
import tarfile
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from pantheon_age_semantics import ROOT, OUT, SRC, moments


def sha(p):
    with p.open("rb") as f:return hashlib.file_digest(f,"sha256").hexdigest()


def violations(d,ca):
    z,dust,tau,start,turn,phi,c=d[:,:7].T
    return dict(nonfinite_parameter=~np.isfinite(d[:,:7]).all(axis=1),
        logZ_bound=~((z>-2.5)&(z<.5)),dust_bound=~((dust>=0)&(dust<=.9)),
        tau_bound=~((tau>.1)&(tau<10)),
        start_bound=~((start>.5)&(start<turn-2)),
        transition_bound=~((turn>2.5)&(turn<=ca)),
        phi_bound=~((phi>-1.520838)&(phi<1.520838)),c_bound=~((c>-45)&(c<-5)))


def main():
    original_script=OUT/"provenance/pantheon_age_archive.first-pass.py"
    dependencies=[Path(__file__),Path(__file__).with_name("pantheon_age_semantics.py")]
    execution_hashes={str(p.relative_to(ROOT)):sha(p) for p in dependencies}
    archive=SRC/"rose2019-campbellG.tar.gz"
    photo=ROOT/"sources/repos/benjaminrose__mc-age/data/campbell_global.tsv"
    ph=pd.read_csv(photo,sep="\t").set_index("SNID")
    metadata=json.loads((SRC/"zenodo-3875482.json").read_text())
    entry=next(x for x in metadata["files"] if x["key"]=="campbellG.tar.gz")
    with archive.open("rb") as f:md5=hashlib.file_digest(f,"md5").hexdigest()
    assert archive.stat().st_size==entry["size"] and entry["checksum"]=="md5:"+md5
    cosmo=[FlatLambdaCDM(H0=70,Om0=.27,Tcmb0=tc) for tc in [0,2.725]]
    rows=[];details=[]
    with tarfile.open(archive,mode="r|gz") as bundle:
        for member in bundle:
            if not member.isfile() or not member.name.endswith("_chain.tsv"):continue
            cid=int(re.search(r"SN(\d+)_",member.name).group(1));redshift=float(ph.loc[cid,"redshift"])
            raw=bundle.extractfile(member).read()
            member_sha=hashlib.sha256(raw).hexdigest()
            d=pd.read_csv(io.BytesIO(raw),sep="\t",comment="#",header=None).to_numpy();del raw
            trials=[]
            for co in cosmo:
                ca=float(co.age(redshift).value)
                bad=violations(d,ca);valid=~np.logical_or.reduce(list(bad.values()))
                paired=valid&np.isfinite(d[:,7])
                tau,start,turn,phi=d[paired,2:6].T
                old,_=moments(tau,start,turn,np.tan(phi),ca,False)
                assert len(old)>0 and np.isfinite(old).all()
                trials.append((float(np.median(abs(old-d[paired,7]))),ca,bad,paired,old))
            chosen=int(np.argmin([x[0] for x in trials]));mad,ca,bad,valid,old=trials[chosen]
            tau,start,turn,phi=d[valid,2:6].T
            new,_=moments(tau,start,turn,np.tan(phi),ca,True)
            assert np.isfinite(new).all()
            original=d[valid,7];delta=new-original;error=abs(old-original)
            qold=np.quantile(original,[.16,.5,.84]);qnew=np.quantile(new,[.16,.5,.84])
            finite_archive=np.isfinite(d[:,7]);qall=np.quantile(d[finite_archive,7],[.16,.5,.84])
            cc=np.corrcoef(np.c_[new,original,d[valid,0],d[valid,1]].T)
            invalid=~valid;ids=np.flatnonzero(invalid)
            detail=dict(CID=cid,member=member.name,member_sha256=member_sha,
                draws_total=len(d),draws_used=int(valid.sum()),draws_excluded=int(invalid.sum()),
                archive_nonfinite_age=int((~finite_archive).sum()),
                violations={k:int(v.sum()) for k,v in bad.items()},
                rows_with_negative_tau=int((d[:,2]<=0).sum()),
                example_invalid_rows=[dict(row_zero_based=int(j),parameters=d[j,:7].tolist(),archived_age=float(d[j,7])) for j in ids[:5]])
            details.append(detail)
            rows.append(dict(CID=cid,redshift=redshift,draws_total=len(d),draws_valid=int(valid.sum()),
                draws_excluded=int(invalid.sum()),invalid_fraction=float(invalid.mean()),
                Tcmb0_K=[0.,2.725][chosen],original_full_age_median=qall[1],
                original_valid_age_q16=qold[0],original_valid_age_median=qold[1],original_valid_age_q84=qold[2],
                fsps_valid_age_q16=qnew[0],fsps_valid_age_median=qnew[1],fsps_valid_age_q84=qnew[2],
                median_age_change_valid=qnew[1]-qold[1],median_effect_of_removing_invalid=qold[1]-qall[1],
                mean_age_change_valid=float(delta.mean()),
                fraction_valid_draws_abs_change_gt_1Gyr=float(np.mean(abs(delta)>1)),
                code_archive_abs_error_median=mad,code_archive_abs_error_p99=float(np.quantile(error,.99)),
                code_archive_abs_error_max=float(error.max()),
                fraction_valid_code_error_gt_p01Gyr=float(np.mean(error>.01)),
                original_age_dust_corr=cc[1,3],fsps_age_dust_corr=cc[0,3],
                original_age_logZ_corr=cc[1,2],fsps_age_logZ_corr=cc[0,2],dust_logZ_corr=cc[2,3]))
            pd.DataFrame(rows).to_csv(OUT/"r19-global-validated-age-audit.csv",index=False)
            (OUT/"r19-global-prior-violations.json").write_text(json.dumps(details,indent=2,allow_nan=False)+"\n")
            print(f"{len(rows):3d} CID {cid:5d}: valid {valid.sum()}/{len(d)}; median {qold[1]:.5f} -> {qnew[1]:.5f}; source error {mad:.3g}",flush=True)
    frame=pd.DataFrame(rows);change=frame.median_age_change_valid
    output=dict(hosts=len(frame),draws_total=int(frame.draws_total.sum()),draws_valid=int(frame.draws_valid.sum()),
        draws_excluded=int(frame.draws_excluded.sum()),hosts_with_excluded_draws=int((frame.draws_excluded>0).sum()),
        max_host_invalid_fraction=float(frame.invalid_fraction.max()),
        all_violations={key:sum(x["violations"][key] for x in details) for key in details[0]["violations"]},
        negative_tau_draws=sum(x["rows_with_negative_tau"] for x in details),
        age_change_conditional_on_valid_draws=dict(mean=float(change.mean()),median=float(change.median()),
           min=float(change.min()),max=float(change.max()),n_abs_gt_p1=int((abs(change)>.1).sum()),
           n_abs_gt_p5=int((abs(change)>.5).sum()),n_abs_gt_1=int((abs(change)>1).sum()),
           original_median_below4=int((frame.original_valid_age_median<4).sum()),
           fsps_median_below4=int((frame.fsps_valid_age_median<4).sum()),
           n_crossing_4Gyr=int(((frame.original_valid_age_median<4)!=(frame.fsps_valid_age_median<4)).sum())),
        invalid_draw_exclusion_median_impact_max_Gyr=float(abs(frame.median_effect_of_removing_invalid).max()),
        source_validation=dict(max_host_median_abs_error=float(frame.code_archive_abs_error_median.max()),
           max_host_p99_abs_error=float(frame.code_archive_abs_error_p99.max()),
           max_host_max_abs_error=float(frame.code_archive_abs_error_max.max()),
           Tcmb0_counts=frame.Tcmb0_K.value_counts().to_dict()),
        posterior_correlations={col:dict(min=float(frame[col].min()),median=float(frame[col].median()),
            max=float(frame[col].max()),n_abs_gt_p5=int((abs(frame[col])>.5).sum()))
            for col in ["original_age_dust_corr","fsps_age_dust_corr","original_age_logZ_corr","fsps_age_logZ_corr","dust_logZ_corr"]},
        limitation="All archived draws were audited. Age shifts condition on the explicitly counted subset accepted by the stated MC-Age prior and finite original age. Invalid draws cannot be repaired as posterior samples by merely recalculating ages. No photometry was refit and no claim is made about private C25 posterior values or cosmological impact.")
    (OUT/"r19-global-validated-age-summary.json").write_text(json.dumps(output,indent=2,allow_nan=False)+"\n")
    inputs=[archive,photo,original_script]
    provenance=dict(executed_script_hashes_captured_at_start=execution_hashes,
        input_sha256={str(p.relative_to(ROOT)):sha(p) for p in inputs},archive_verified_md5=md5,
        prerequisite_first_pass="r19-global-posterior-age-audit.csv is diagnostic only, superseded by explicit prior validation",
        versions=dict(numpy=np.__version__,pandas=pd.__version__))
    assert all(sha(ROOT/p)==value for p,value in execution_hashes.items())
    (OUT/"r19-global-validated-age-manifest.json").write_text(json.dumps(provenance,indent=2,allow_nan=False)+"\n")
    print(json.dumps(output,indent=2,allow_nan=False),flush=True)


if __name__=="__main__":main()
