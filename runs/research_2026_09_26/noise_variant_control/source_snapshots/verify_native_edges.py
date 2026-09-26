"""Independent finite-difference, projection and fixed-weight replacement audit."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_triangular

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
SELECT={"P21_rho000":(1294,7349,20101),"P21_rho090":(1294,7349),
        "P21_noisetrue120":(7349,20101)}


def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for part in iter(lambda:f.read(1048576),b""):
            h.update(part)
    return h.hexdigest()


def main():
    m=json.loads((HERE/"native_edge_sensitivity/manifest.json").read_text())
    for rel,digest in {**m["inputs_sha256"],**m["outputs_sha256"]}.items():
        assert sha(ROOT/rel)==digest,rel
    assert sha(HERE/"native_edge_derivatives.py")==sha(HERE/"source_snapshots/native_edge_derivatives.py")
    assert sha(HERE/"compare_native_edges.py")==sha(HERE/"source_snapshots/compare_native_edges.py")
    for scope,count in (("full",7),("half",6)):
        base=HERE/f"native_edge_{scope}"
        assert json.loads((base/"run-gate.json").read_text())["object_exports"]==count*7
        for arm,ids in SELECT.items():
            labels=[f"p{j}_{sign}" for j in (1,2,3) for sign in ("minus","plus")]
            if scope=="full": labels=["baseline"]+labels
            for label in labels:
                manifest=json.loads((base/arm/label/"objectives/manifest.json").read_text())
                assert manifest["count"]==len(ids)
                assert set(manifest["checks"])==set(map(str,ids))
                assert max(abs(v["objective_residual"]) for v in manifest["checks"].values())<1e-7
    obj=pd.read_csv(HERE/"native_edge_sensitivity/native-edge-object-scores.csv")
    assert len(obj)==7 and not obj.duplicated(["arm","CID"]).any()
    with np.load(HERE/"native_edge_sensitivity/native-jacobians.npz",allow_pickle=False) as js:
        coeff=np.load(ROOT/"runs/research_2026_09_26/astra_design/validation1020/frozen-discovery-coefficients.npz",allow_pickle=False)["basis_mean"]
        # Explicit scalar observer magnitude direction, independent of score helper.
        gauge=np.array([coeff[0],-coeff.sum(),coeff[1],coeff[2]])
        max_M=max_I=max_G=max_Q=0.
        for r in obj.itertuples():
            with np.load(HERE/r.arm/"objectives"/f"objective_{r.CID}.npz",allow_pickle=False) as a:
                y=a["data_flux"];f=a["model_flux"];C=a["frozen_flux_covariance"];b=a["band"]
            L=np.linalg.cholesky(C)
            ry=solve_triangular(L,y-f,lower=True)
            observer=-.4*np.log(10)*f*np.array([gauge["griz".index(str(x))] for x in b])
            ro=solve_triangular(L,observer,lower=True)
            for scale in ("full","half"):
                J=js[f"{r.arm}_{r.CID}_{scale}"]
                U,s,_=np.linalg.svd(solve_triangular(L,J,lower=True),full_matrices=False)
                assert int(np.sum(s>s[0]*1e-10))==4
                projected_y=ry-U@(U.T@ry)
                projected_o=ro-U@(U.T@ro)
                M=float(projected_o@projected_y);I=float(projected_o@projected_o)
                G=M-I/2;Q=float(projected_y@projected_y)
                max_M=max(max_M,abs(M-getattr(r,f"native_{scale}_M")))
                max_I=max(max_I,abs(I-getattr(r,f"native_{scale}_I")))
                max_G=max(max_G,abs(G-getattr(r,f"native_{scale}_G")))
                max_Q=max(max_Q,abs(Q-getattr(r,f"native_{scale}_Q")))
                assert abs(M-getattr(r,f"native_{scale}_M"))<1e-9
                assert abs(I-getattr(r,f"native_{scale}_I"))<1e-9
                assert abs(G-getattr(r,f"native_{scale}_G"))<1e-9
                assert abs(Q-getattr(r,f"native_{scale}_Q"))<1e-9
    support=pd.read_csv(HERE/"transport/supported-object-ledger.csv")
    trans=pd.read_csv(HERE/"native_edge_sensitivity/fixed-transport-sensitivity.csv")
    max_fixed=0.
    for row in trans.itertuples():
        q=support.loc[(support.stage==row.stage)&(support.transport==row.transport)&(support.arm==row.arm)]
        pairs=obj.loc[obj.arm==row.arm].set_index("CID")
        matched=q.loc[q.CID.isin(pairs.index)]
        w=matched.base_weight.to_numpy(float)
        dM=float(w@pairs.loc[matched.CID,"delta_M"].to_numpy(float)) if len(matched) else 0.
        dI=float(w@pairs.loc[matched.CID,"delta_I"].to_numpy(float)) if len(matched) else 0.
        old=pd.read_csv(HERE/"transport/transport-summary.csv")
        ref=old.loc[(old.stage==row.stage)&(old.transport==row.transport)&(old.arm==row.arm)].iloc[0]
        amp=(ref.weighted_M+dM)/(ref.weighted_I+dI)
        max_fixed=max(max_fixed,abs(amp-row.native_J_amplitude))
        assert abs(amp-row.native_J_amplitude)<1e-10
    result={"status":"passed","manifest_hashes":len(m["inputs_sha256"])+len(m["outputs_sha256"]),
            "full_object_exports":49,"half_object_exports":42,"selected_pairs":7,
            "max_independent_native_M_error":max_M,"max_independent_native_I_error":max_I,
            "max_independent_native_G_error":max_G,"max_independent_native_Q_error":max_Q,
            "max_fixed_transport_amplitude_error":max_fixed}
    (HERE/"native_edge_sensitivity/independent-verification.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__=="__main__":
    main()
