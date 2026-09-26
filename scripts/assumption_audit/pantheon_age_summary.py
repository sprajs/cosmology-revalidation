#!/usr/bin/env python3
"""Summarize the completed, finite original-R19 age audit only."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/"runs/assumption_audit/pantheon"


def main():
    table=OUT/"r19-global-validated-age-audit.csv"
    summary=OUT/"r19-global-validated-age-summary.json"
    details=OUT/"r19-global-prior-violations.json"
    d=pd.read_csv(table);s=json.loads(summary.read_text());v=json.loads(details.read_text())
    assert len(d)==s["hosts"]==len(v)==103 and d.CID.nunique()==len(d)
    assert np.isfinite(d.select_dtypes(include="number").to_numpy()).all()
    photo=ROOT/"sources/repos/benjaminrose__mc-age/data/campbell_global.tsv"
    original=OUT/"sources/rose2019-ages-campbellG.tsv"
    revised=ROOT/"data/derived/age_signal/all_age_rows.csv"
    ph=pd.read_csv(photo,sep="\t");old=pd.read_csv(original,sep="\t",comment="#",header=None)
    rev=pd.read_csv(revised);rev=rev[rev.age_source=="R19"]
    ids=set(d.CID);rids=set(rev.CID.astype(int))
    cross=dict(original_archive_hosts=len(d),original_archive_unique_CIDs=d.CID.nunique(),
        original_photometry_rows=len(ph),original_summary_rows=len(old),revised_C25_R19_rows=len(rev),
        original_sources_same_CID_set=ids==set(ph.SNID)==set(old[0].astype(int)),
        original_absent_from_C25=sorted(ids-rids),C25_absent_from_original=sorted(rids-ids),
        note="All original global hosts retained. CID15459 absence in C25 is not explained by these public tables.")
    (OUT/"r19-sample-crosswalk.json").write_text(json.dumps(cross,indent=2,allow_nan=False)+"\n")
    with plt.rc_context({"font.size":10,"axes.spines.top":False,"axes.spines.right":False}):
        fig,axes=plt.subplots(1,2,figsize=(11,4.8),layout="constrained")
        ax=axes[0]
        ax.errorbar(d.original_valid_age_median,d.fsps_valid_age_median,
            xerr=np.vstack([d.original_valid_age_median-d.original_valid_age_q16,d.original_valid_age_q84-d.original_valid_age_median]),
            yerr=np.vstack([d.fsps_valid_age_median-d.fsps_valid_age_q16,d.fsps_valid_age_q84-d.fsps_valid_age_median]),
            fmt="none",ecolor="#78909c",alpha=.2,linewidth=.65,zorder=0)
        scatter=ax.scatter(d.original_valid_age_median,d.fsps_valid_age_median,c=100*d.invalid_fraction,
            cmap="viridis",s=28,edgecolor="white",linewidth=.4)
        ax.plot([0,13],[0,13],color="#455a64",ls="--",lw=1)
        ax.set(xlim=(0,13),ylim=(0,13),xlabel="Original age posterior median (Gyr)",
            ylabel="Median using the fitted FSPS history (Gyr)",title="Same parameter draws; different age integral")
        fig.colorbar(scatter,ax=ax,label="Archived draws outside stated prior (%)",shrink=.85)
        ax=axes[1]
        changes=d.median_age_change_valid.to_numpy()
        bins=np.arange(np.floor(changes.min()*2)/2, np.ceil(changes.max()*2)/2+.51,.5)
        ax.hist(changes,bins=bins,color="#1565a0",edgecolor="white")
        ax.axvline(0,color="#455a64",lw=1)
        ax.set(xlabel="Change in host age posterior median (Gyr)",ylabel="Number of hosts",
            title=f"Original R19 global sample: {len(d)} hosts")
        fig.suptitle("Audit of archived host ages",fontsize=15)
        fig.supxlabel("Ages shown conditional on stated-prior-valid draws. No photometry refit; not the revised C25 catalogue.",fontsize=9)
        fig.savefig(OUT/"r19-age-audit.png",dpi=180)
        plt.close(fig)
    inputs=[Path(__file__),table,summary,details,photo,original,revised]
    manifest={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}
    (OUT/"r19-age-summary-manifest.json").write_text(json.dumps(manifest,indent=2,allow_nan=False)+"\n")
    print(json.dumps(cross,indent=2))


if __name__=="__main__":main()
