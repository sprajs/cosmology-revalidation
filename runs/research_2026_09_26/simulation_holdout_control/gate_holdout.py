"""Native objective/raw-PHOT closure with explicit 1000-object LCPLOT cap."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DESIGN = ROOT / "runs/research_2026_09_26/astra_design/simulation_holdout_design"
BASE = ROOT / "scripts/research_2026_09_26/simulation_residual_control.py"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main(arm):
    assert arm in ("P21","G10")
    path = HERE/arm
    assert not (path/"provenance-ledger.csv").exists()
    spec=importlib.util.spec_from_file_location("pilot_match_rows",BASE)
    native=importlib.util.module_from_spec(spec);assert spec.loader;spec.loader.exec_module(native)
    frozen=json.loads((DESIGN/"manifest.json").read_text())
    for rel,wanted in frozen["sha256"].items():
        assert digest(ROOT/rel)==wanted,rel
    prep=json.loads((path/"prepared.json").read_text())
    fit=json.loads((path/"fit-gate.json").read_text())
    exp=json.loads((path/"export-gate.json").read_text())
    assert prep["selected_count"] == frozen["counts"][arm]["holdout"]
    assert fit["returncode"] == exp["returncode"] == 0 and fit["graceful"]
    assert exp["max_objective_error"] < 1e-7
    assert exp["objects"] == len(fit["objective_cids"])
    assert len(fit["fitres_cids"])/prep["selected_count"] >= .95
    sys.path.insert(0,str(ROOT/"scripts/phase2/official"))
    from audit_fits import read_fit
    newfit=read_fit(path/"fit.FITRES.TEXT")
    newfit.CID=newfit.CID.astype(int)
    assert newfit.CID.is_unique
    newfit=newfit.set_index("CID")
    assert set(newfit.index)==set(fit["fitres_cids"])
    cohort=pd.read_csv(DESIGN/f"{arm}-cohort.csv").set_index("CID")
    version=f"PH2_pilot02_{arm}"
    base=ROOT/f"phase2/literature/simulations/outputs/{version}"
    with fits.open(base/f"{version}_HEAD.FITS",memmap=True) as hh, fits.open(base/f"{version}_PHOT.FITS",memmap=True) as pp:
        heads=hh[1].data;phot=pp[1].data
        index={int(row["SNID"]):row for row in heads}
        assert len(index)==len(heads)
        raw={}
        for cid in fit["fitres_cids"]:
            h=index[cid]
            lo,hi=int(h["PTROBS_MIN"]),int(h["PTROBS_MAX"])
            q=phot[lo-1:hi].copy()
            assert len(q)==int(h["NOBS"])
            raw[cid]=q
    lc=pd.read_csv(path/"fit.LCPLOT.TEXT",sep=r"\s+",comment="#",header=None)
    assert lc.shape[1]==10
    nml=(path/"fit.nml").read_text()
    log=(path/"fit.log").read_text()
    assert "MXLC_PLOT = 1000" in nml
    assert "SNLCPAK-table: MXLC_PLOT=    1000" in log
    source=(ROOT/"phase2/official/build/SNANA-audit-v3/src/snana.F90").read_text()
    assert "if ( N_SNLC_PLOT < MXLC_PLOT ) THEN" in source
    assert "IF ( MADE_LCPLOT )  N_SNLC_PLOT = N_SNLC_PLOT + 1" in source
    fit_source_path=ROOT/"phase2/official/build/SNANA-audit-v3/src/snlc_fit.F90"
    fit_source=fit_source_path.read_text()
    assert "epoch           = EPLIST_FIT(ifitdata)" in fit_source
    assert "LFITDATA = ifitdata .LE. NFITDATA" in fit_source
    assert "'PHASE2_FLUX:'" in fit_source
    objective_order=[];seen=set()
    for line in log.splitlines():
        if line.startswith("PHASE2_OBJECTIVE:"):
            cid=int(line.split()[1])
            if cid not in seen:
                seen.add(cid);objective_order.append(cid)
    assert set(objective_order)==set(fit["objective_cids"])
    successful_order=[cid for cid in objective_order if cid in set(fit["fitres_cids"])]
    plotted=set(lc.loc[lc[6]==1,0].astype(int))
    assert len(plotted)==min(1000,len(successful_order))
    assert plotted==set(successful_order[:1000]), "sporadic pre-cap or post-cap LCPLOT membership"
    ids=sorted(set(fit["fitres_cids"]) & set(fit["objective_cids"]))
    rows=[]
    for cid in ids:
        with np.load(path/"objectives"/f"objective_{cid}.npz",allow_pickle=False) as obj:
            mjd,band,y,err=(obj[k] for k in ("MJD","band","data_flux","data_fluxerr"))
            assert len(mjd)>=5
            out=np.column_stack([mjd,band,y,err])
            q=raw[cid]
            phot_rows=np.column_stack([q["MJD"].astype(float),np.char.strip(q["BAND"].astype(str)),
                                       q["FLUXCAL"].astype(float),q["FLUXCALERR"].astype(float)])
            dt,df,de=native.match_rows(out,phot_rows,True)
            assert dt<.005 and df<.001 and de<.001,(cid,dt,df,de)
            arch=cohort.loc[cid]
            fitrow=newfit.loc[cid]
            assert str(arch.field)==str(q["FIELD"][0]).strip()
            assert abs(float(fitrow.PKMJDINI)-float(arch.PKMJDINI))<.005
            assert abs(float(fitrow.zHEL)-float(arch.zHEL))<1e-6
            pars=obj["parameters_x0_x1_c_t0"]
            dp=np.array([float(fitrow[k]) for k in ("x0","x1","c","PKMJD")])-pars
            assert abs(dp[0])<2e-8 and abs(dp[1])<1e-5 and abs(dp[2])<1e-5 and abs(dp[3])<.005,(cid,dp)
            C=obj["frozen_flux_covariance"]
            assert C.shape==(len(out),len(out))
            mineig=float(np.linalg.eigvalsh(C).min())
            assert mineig>0
            plot=lc.loc[(lc[0]==cid)&(lc[6]==1)]
            has_plot=len(plot)>0
            if has_plot:
                plot_rows=np.column_stack([plot[1].to_numpy(float),plot[7].to_numpy(str),
                                           plot[4].to_numpy(float),plot[5].to_numpy(float)])
                assert len(plot_rows)==len(out),(cid,len(plot_rows),len(out))
                dt_plot,_,_=native.match_rows(out,plot_rows,False)
                assert dt_plot<.005,(cid,dt_plot)
            else:
                dt_plot=np.nan
            rows.append(dict(CID=cid,arm=arm,law="approx_minus99",epochs=len(out),
                             raw_mjd_max_abs=dt,raw_flux_max_abs=df,raw_error_max_abs=de,
                             lcplot_available=has_plot,lcplot_mjd_max_abs=dt_plot,
                             first_phot_field=str(arch.field),
                             archived_basic_quality=bool(arch.basic_quality_pass),
                             original_pkmjdini=float(arch.PKMJDINI),new_pkmjdini=float(fitrow.PKMJDINI),
                             archived_pkmjd=float(arch.PKMJD),new_pkmjd=float(fitrow.PKMJD),
                             min_epoch_cov_eigenvalue=mineig))
    frame=pd.DataFrame(rows)
    frame.to_csv(path/"provenance-ledger.csv",index=False,float_format="%.17g")
    result={"cohort":"holdout","arm":arm,"law":"approx_minus99",
            "selected":prep["selected_count"],"exported":len(frame),
            "success_fraction":len(frame)/prep["selected_count"],
            "full_pilot_adequate":len(frame)/prep["selected_count"]>=.95,
            "epochs":int(frame.epochs.sum()),
            "max_raw_mjd_difference":float(frame.raw_mjd_max_abs.max()),
            "max_raw_flux_difference":float(frame.raw_flux_max_abs.max()),
            "max_raw_error_difference":float(frame.raw_error_max_abs.max()),
            "lcplot_verified_objects":int(frame.lcplot_available.sum()),
            "lcplot_unavailable_objects":int((~frame.lcplot_available).sum()),
            "max_lcplot_mjd_difference_when_available":float(frame.lcplot_mjd_max_abs.max()),
            "max_pkmjdini_difference":float(np.max(abs(frame.new_pkmjdini-frame.original_pkmjdini))),
            "max_pkmjd_difference":float(np.max(abs(frame.new_pkmjd-frame.archived_pkmjd))),
            "min_epoch_cov_eigenvalue":float(frame.min_epoch_cov_eigenvalue.min()),
            "native_objective_mask_and_raw_phot_closure_all":True,
            "lcplot_cap":1000,
            "lcplot_first_successful_objective_order_only":True,
            "lcplot_source_sha256":digest(ROOT/"phase2/official/build/SNANA-audit-v3/src/snana.F90"),
            "native_accepted_epoch_source_sha256":digest(fit_source_path),
            "gate_limit":"LCPLOT generated only up to MXLC_PLOT=1000; beyond cap native objective records the accepted rows; every row matched exact raw HEAD/PHOT.",
            "source_sha256":digest(Path(__file__)),
            "base_match_rows_sha256":digest(BASE),
            "provenance_ledger_sha256":digest(path/"provenance-ledger.csv")}
    (path/"provenance-gate.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("gate_limit",)},indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("arm",choices=("P21","G10"));args=parser.parse_args();main(args.arm)
