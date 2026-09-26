"""Frozen P21/G10 disjoint holdout and combined fixed-direction score."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DESIGN = ROOT / "runs/research_2026_09_26/astra_design/simulation_holdout_design"
PILOT = ROOT / "runs/research_2026_09_26/simulation_residual_control"
COMMON = ROOT / "runs/research_2026_09_26/common_classifier_residual"
VAL = ROOT / "runs/research_2026_09_26/astra_design/validation1020"
PLAN = HERE / "score-plan.md"
EXPECTED_PLAN = "ee108850f9d64a3689f7f9c50e8667ee6d89cecb56547b1542b1893ee104288f"


def digest(p):
    h = hashlib.sha256()
    with Path(p).open("rb") as stream:
        for part in iter(lambda: stream.read(1048576), b""):
            h.update(part)
    return h.hexdigest()


def frozen_gates():
    assert digest(PLAN) == EXPECTED_PLAN
    design = json.loads((DESIGN / "manifest.json").read_text())
    for rel, expected in design["sha256"].items():
        assert digest(ROOT / rel) == expected, rel
    gate = json.loads((HERE / "inference/inference-gate.json").read_text())
    assert gate["count"] == 3477 and gate["max_single_batch_abs"] <= 1e-6
    assert json.loads((COMMON / "score-independent-verification.json").read_text())["result"] == "passed"
    paths = {}
    for arm in ("P21", "G10"):
        path = HERE / arm
        prep = json.loads((path / "prepared.json").read_text())
        fit = json.loads((path / "fit-gate.json").read_text())
        exp = json.loads((path / "export-gate.json").read_text())
        prov = json.loads((path / "provenance-gate.json").read_text())
        tangent = json.loads((path / "tangent-gate.json").read_text())
        edge = json.loads((path / "mean-edge-gate.json").read_text())
        assert prep["law_option"] == -99 and prep["selected_count"] == design["counts"][arm]["holdout"]
        assert fit["returncode"] == exp["returncode"] == 0 and fit["graceful"]
        assert prov["full_pilot_adequate"] and prov["success_fraction"] >= .95
        assert tangent["all_rank4"] and edge["eligible_objects"] == prov["exported"]
        assert digest(path / "provenance-ledger.csv") == prov["provenance_ledger_sha256"]
        assert digest(path / "mean-edge-ledger.csv") == edge["ledger_sha256"]
        paths[arm] = path
    return paths


def make_object_scores(paths):
    dest = HERE / "object_scores"
    dest.mkdir(exist_ok=False)
    source = ROOT / "scripts/salt_dust_audit/flux_response.py"
    spec = importlib.util.spec_from_file_location("flux_response", source)
    fr = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(fr)
    model, bands, _, zp = fr.build_model()
    magoff = {str(row["Filter Name"])[-1]: float(row["Primary Mag"]) for row in zp}
    with np.load(VAL / "frozen-discovery-coefficients.npz", allow_pickle=False) as frozen:
        coeff = frozen["basis_mean"].astype(float)
    assert coeff.shape == (3,)
    probs = pd.read_csv(HERE / "inference/holdout-probabilities.csv")
    support = pd.read_csv(HERE / "truth-support-ledger.csv")
    records = []
    arrays = []
    sys.path.insert(0, str(ROOT / "scripts/phase2/official"))
    from audit_fits import read_fit
    for arm, path in paths.items():
        cohort = pd.read_csv(DESIGN / f"{arm}-cohort.csv").set_index("CID")
        native = read_fit(path / "fit.FITRES.TEXT")
        native.CID = native.CID.astype(int)
        native = native.set_index("CID")
        ids = pd.read_csv(path / "provenance-ledger.csv").CID.astype(int).tolist()
        assert set(ids) == set(native.index)
        pred = probs.loc[probs.arm == arm].set_index("CID")
        tags = support.loc[support.arm == arm].set_index("CID")
        for j,cid in enumerate(ids):
            with np.load(path / "objectives" / f"objective_{cid}.npz", allow_pickle=False) as obj:
                pars = obj["parameters_x0_x1_c_t0"]
                row = SimpleNamespace(CID=str(cid), x0=float(pars[0]), x1=float(pars[1]),
                                      c=float(pars[2]), PKMJD=float(pars[3]),
                                      zHEL=float(obj["zHEL"][0]))
                points = pd.DataFrame({"MJD":obj["MJD"],"BAND":obj["band"],
                                       "FLUXCAL":obj["data_flux"],"FLUXCALERR":obj["data_fluxerr"]})
                order = points.sort_values(["MJD","BAND"],kind="stable").index.to_numpy()
                _, audit, _, diagnostic = fr.analyze(row,points,{"MWEBV":obj["MWEBV"][0]},
                                                       model,bands,magoff,False)
                f = obj["model_flux"][order]
                y = obj["data_flux"][order]
                b = obj["band"][order]
                C = obj["frozen_flux_covariance"][np.ix_(order,order)]
                J = audit["jacobian_flux"].copy()
                J[:,0] = -fr.K*f
                L = np.linalg.cholesky(C)
                wj = solve_triangular(L,J,lower=True)
                U,s,_ = np.linalg.svd(wj,full_matrices=True)
                assert int(np.sum(s > s[0]*1e-10)) == 4
                Q = U[:,4:]
                residual = Q.T @ solve_triangular(L,y-f,lower=True)
                griz = np.column_stack([-fr.K*f*(b==band) for band in "griz"])
                obs = np.column_stack([griz[:,0]-griz[:,1],griz[:,2]-griz[:,1],griz[:,3]-griz[:,1]])
                T = Q.T @ solve_triangular(L,obs,lower=True)
                u = T.T @ residual
                F = T.T @ T
                a = float(coeff @ u)
                information = float(coeff @ F @ coeff)
                assert information > 0 and np.isfinite(a)
                archived = cohort.loc[cid]
                fitrow = native.loc[cid]
                assert str(archived.field) in ("C1","C2","C3","E1","E2","S1","S2","X1","X2","X3")
                quality = (abs(float(fitrow.x1)) < 3 and abs(float(fitrow.c)) < .3
                           and float(fitrow.x1ERR) < 1 and float(fitrow.PKMJDERR) < 2
                           and float(fitrow.cERR) < 1.5 and float(fitrow.FITPROB) > .001
                           and .025 < float(fitrow.zHD) < 1.2)
                records.append(dict(arm=arm,law="approx_minus99",CID=cid,field=str(archived.field),
                                    zHEL=float(archived.zHEL),zbin=int(archived.zbin),
                                    cell=f"{archived.field}_{int(archived.zbin)}",LIBID=int(archived.SIM_LIBID),
                                    SNRMAX1_archived=float(archived.SNRMAX1),
                                    support_class=str(tags.loc[cid,"support_class"]),
                                    archived_basic_quality=bool(archived.basic_quality_pass),
                                    new_basic_quality=bool(quality),pIa=float(pred.loc[cid,"pIa"]),
                                    gt999=bool(pred.loc[cid,"gt999"]),epochs=len(y),
                                    projected_dimension=len(residual),
                                    projected_chi2=float(residual @ residual),
                                    matched_filter=a,information=information,fixed_gain=a-information/2,
                                    nuisance_condition=float(s[0]/s[-1]),
                                    max_native_independent_mean_sigma=float(np.max(np.abs((audit["flux_model"]-f)/obj["data_fluxerr"][order]))),
                                    independent_derivative_halfstep_error=float(diagnostic["derivative_relative_error"])))
                arrays.append((u,F))
            if (j+1)%250 == 0:
                print(f"scored {arm} {j+1}/{len(ids)}",flush=True)
    frame = pd.DataFrame(records).sort_values(["arm","CID"])
    lookup = {(r["arm"],r["CID"]):x for r,x in zip(records,arrays)}
    frame.to_csv(dest/"holdout-object-scores.csv",index=False,float_format="%.17g")
    np.savez_compressed(dest/"holdout-sufficient-arrays.npz",
                        arm=frame.arm.to_numpy(dtype="U3"),CID=frame.CID.to_numpy(int),
                        u=np.array([lookup[(r.arm,r.CID)][0] for r in frame.itertuples()]),
                        F=np.array([lookup[(r.arm,r.CID)][1] for r in frame.itertuples()]))
    sources = [Path(__file__),PLAN,source,VAL/"frozen-discovery-coefficients.npz",
               HERE/"generator-noise-provenance.json",HERE/"native_design/manifest.json",
               HERE/"inference/holdout-probabilities.csv",HERE/"truth-support-ledger.csv"]
    sources += [DESIGN/"manifest.json",*(DESIGN/f"{arm}-cohort.csv" for arm in ("P21","G10")),
                HERE/"inference/manifest.json",HERE/"inference-source-closure.json",
                HERE/"truth-support-manifest.json",HERE/"baseline-gate-failure.json",
                HERE/"gate_holdout.py",HERE/"source_snapshots/native_prepare_fit.py",
                HERE/"source_snapshots/holdout_inference.py",
                HERE/"source_snapshots/baseline_gate_attempt.py"]
    sources += [p/f for p in paths.values() for f in ("prepared.json","fit-gate.json","export-gate.json",
                    "provenance-gate.json","tangent-gate.json","mean-edge-gate.json",
                    "provenance-ledger.csv","mean-edge-ledger.csv","objectives/manifest.json","fit.FITRES.TEXT")]
    (dest/"manifest.json").write_text(json.dumps({
        "inputs_sha256":{str(p.relative_to(ROOT)):digest(p) for p in sources},
        "outputs_sha256":{str(p.relative_to(ROOT)):digest(p) for p in dest.iterdir() if p.name!="manifest.json"},
        "frozen_coefficients":coeff.tolist(),"scope":"Full native holdout scores; no transport or outcome-selected cut"},
        indent=2)+"\n")
    return {"objects":len(frame),"by_arm":frame.groupby("arm").size().to_dict(),
            "mean_edge_flagged":{arm:json.loads((path/"mean-edge-gate.json").read_text())["flagged_epochs"]
                                 for arm,path in paths.items()}}


def main():
    paths = frozen_gates()
    assert not (HERE/"object_scores").exists()
    print(json.dumps(make_object_scores(paths),indent=2))


if __name__ == "__main__":
    main()
