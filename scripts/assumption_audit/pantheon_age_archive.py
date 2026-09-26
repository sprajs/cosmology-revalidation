#!/usr/bin/env python3
"""Recompute all original R19 global posterior ages with the fitted FSPS SFH.

This changes no photometry or parameter posterior, estimates no age-HR slope,
and does not substitute these ages for the private revised C25 posteriors.
"""
from pathlib import Path
import hashlib
import json
import tarfile
import re
import io
import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
from pantheon_age_semantics import ROOT, OUT, SRC, moments


def digest(path, algorithm="sha256"):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def main():
    archive=SRC/"rose2019-campbellG.tar.gz"
    metadata=json.loads((SRC/"zenodo-3875482.json").read_text())
    entry=next(x for x in metadata["files"] if x["key"]=="campbellG.tar.gz")
    assert archive.stat().st_size==entry["size"]
    md5=digest(archive,"md5")
    assert "md5:"+md5==entry["checksum"]
    photometry=ROOT/"sources/repos/benjaminrose__mc-age/data/campbell_global.tsv"
    ph=pd.read_csv(photometry,sep="\t").set_index("SNID")
    rows=[]
    cosmo=[FlatLambdaCDM(H0=70,Om0=.27,Tcmb0=tc) for tc in [0.,2.725]]
    with tarfile.open(archive,mode="r|gz") as bundle:
        for member in bundle:
            if not member.isfile() or not member.name.endswith("_chain.tsv"):
                continue
            cid=int(re.search(r"SN(\d+)_",member.name).group(1))
            z=float(ph.loc[cid,"redshift"])
            # pandas probes seekability; tarfile's streaming reader lacks that
            # method, so buffer this one member only (never the whole archive).
            payload=bundle.extractfile(member).read()
            d=pd.read_csv(io.BytesIO(payload),sep="\t",comment="#",header=None).to_numpy()
            del payload
            tau,start,transition,phi=d[:,2:6].T
            slope=np.tan(phi)
            trials=[]
            for cosmology in cosmo:
                ca=float(cosmology.age(z).value)
                old,_=moments(tau,start,transition,slope,ca,False)
                trials.append((float(np.median(abs(old-d[:,7]))),ca,old))
            choice=int(np.argmin([x[0] for x in trials]))
            mad,ca,old=trials[choice]
            new,_=moments(tau,start,transition,slope,ca,True)
            change=new-d[:,7]
            error=abs(old-d[:,7])
            corr=np.corrcoef(np.c_[new,d[:,7],d[:,0],d[:,1]].T)
            qold=np.quantile(d[:,7],[.16,.5,.84])
            qnew=np.quantile(new,[.16,.5,.84])
            rows.append(dict(CID=cid,member=member.name,draws=len(d),redshift=z,
                Tcmb0_K=[0.,2.725][choice],
                original_age_q16=qold[0],original_age_median=qold[1],original_age_q84=qold[2],
                fsps_age_q16=qnew[0],fsps_age_median=qnew[1],fsps_age_q84=qnew[2],
                median_age_change=qnew[1]-qold[1],mean_age_change=float(change.mean()),
                fraction_draws_abs_change_gt_1Gyr=float(np.mean(abs(change)>1)),
                code_archive_abs_error_median=mad,code_archive_abs_error_p99=float(np.quantile(error,.99)),
                code_archive_abs_error_max=float(error.max()),fraction_code_archive_error_gt_p01Gyr=float(np.mean(error>.01)),
                original_age_dust_corr=corr[1,3],fsps_age_dust_corr=corr[0,3],
                original_age_logZ_corr=corr[1,2],fsps_age_logZ_corr=corr[0,2],
                dust_logZ_corr=corr[2,3]))
            pd.DataFrame(rows).to_csv(OUT/"r19-global-posterior-age-audit.csv",index=False)
            print(f"{len(rows):3d} CID {cid:5d}: {qold[1]:.4f} -> {qnew[1]:.4f} Gyr; source median error {mad:.3g}",flush=True)
    frame=pd.DataFrame(rows)
    change=frame.median_age_change
    result=dict(hosts=len(rows),total_draws=int(frame.draws.sum()),
        retained_chain_policy="All archived rows; no undocumented burn-in selection or fresh posterior fit",
        age_estimand="formed-mass-weighted mean age of FSPS sfh=5, integrated for each archived parameter draw",
        change_in_host_posterior_medians=dict(mean=float(change.mean()),median=float(change.median()),
           min=float(change.min()),max=float(change.max()),n_abs_gt_p1=int((abs(change)>.1).sum()),
           n_abs_gt_p5=int((abs(change)>.5).sum()),n_abs_gt_1=int((abs(change)>1).sum()),
           quantiles=change.quantile([.025,.16,.5,.84,.975]).to_dict()),
        source_validation=dict(max_host_median_abs_error=float(frame.code_archive_abs_error_median.max()),
           max_host_p99_abs_error=float(frame.code_archive_abs_error_p99.max()),
           Tcmb0_counts=frame.Tcmb0_K.value_counts().to_dict(),
           max_fraction_draws_error_gt_p01=float(frame.fraction_code_archive_error_gt_p01Gyr.max())),
        posterior_correlations={col:dict(min=float(frame[col].min()),median=float(frame[col].median()),
                                        max=float(frame[col].max()),n_abs_gt_p5=int((abs(frame[col])>.5).sum()))
                                for col in ["original_age_dust_corr","fsps_age_dust_corr","original_age_logZ_corr","fsps_age_logZ_corr","dust_logZ_corr"]},
        limitation="Reprocesses original R19 posterior draws under historical FSPS semantics. This is not a new photometric fit, an official corrected catalogue, a revised C25 posterior, or a cosmological effect estimate.")
    (OUT/"r19-global-posterior-age-summary.json").write_text(json.dumps(result,indent=2)+"\n")
    inputs=[archive,photometry,Path(__file__),Path(__file__).with_name("pantheon_age_semantics.py")]
    manifest=dict(inputs={str(p.relative_to(ROOT)):digest(p) for p in inputs},
        zenodo_record="https://zenodo.org/records/3875482",archive_md5_verified=md5,
        fsps_revision="ae31b2f63d865354ce944e5c22eba6e93e01e67d")
    (OUT/"r19-global-posterior-age-manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(json.dumps(result,indent=2),flush=True)


if __name__=="__main__":main()
