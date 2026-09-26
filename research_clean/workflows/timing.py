"""Coherent 2021 simulation-table joins, cuts and timing residuals."""

import numpy as np
from lib.paths import DATA
from lib.records import fitres, write_rows

DEFAULTS = {}


def stats(x):
    x = np.asarray(x, float)
    if not len(x) or not np.isfinite(x).all():
        raise ValueError("Empty/nonfinite statistic input")
    return {
        "n": len(x),
        "mean": float(x.mean()),
        "sd": float(x.std()),
        "q025_median_q975": np.quantile(x, [0.025, 0.5, 0.975]).tolist(),
    }


def run(out, cfg):
    nir = fitres(DATA / "timing/nir.FITRES.gz")
    joint = fitres(DATA / "timing/optnir.FITRES.gz")
    common = nir.index[nir.index.isin(joint.index)]
    if len(nir) != 30000 or len(joint) != 29995 or len(common) != 29995:
        raise ValueError("Source-matched 2021 counts changed")
    n, o = nir.loc[common], joint.loc[common]
    for name in ["SIM_PKMJD", "SIM_DLMAG", "SIM_AV", "SIM_STRETCH", "SIM_RV"]:
        if not np.array_equal(n[name].to_numpy(), o[name].to_numpy()):
            raise ValueError(
                f"Generating truth differs between joined branches: {name}"
            )
    av = o.AV.to_numpy() < 0.3 * o.RV.to_numpy()
    stretch = (o.STRETCH.to_numpy() > 0.75) & (o.STRETCH.to_numpy() < 1.185)
    selected = av & stretch
    selected_ids = set(common[selected])
    bins = []
    for lo, hi in zip([0, 0.2, 0.3, 0.4, 0.5, 0.6], [0.2, 0.3, 0.4, 0.5, 0.6, 1.0]):
        mask = (n.zHD.to_numpy() >= lo) & (n.zHD.to_numpy() < hi)
        good = mask & selected
        bins.append(
            {
                "zlo": lo,
                "zhi": hi,
                "common_rows": int(mask.sum()),
                "selected": int(good.sum()),
                "mean_sim_AV_before": float(n.SIM_AV.to_numpy()[mask].mean())
                if mask.any()
                else None,
                "mean_sim_AV_after": float(n.SIM_AV.to_numpy()[good].mean())
                if good.any()
                else None,
            }
        )
    write_rows(out / "redshift_selection.csv", bins)
    write_rows(
        out / "membership.csv",
        [
            {
                "CID": cid,
                "joint_present": cid in joint.index,
                "passes_AV_stretch": cid in selected_ids,
            }
            for cid in nir.index
        ],
    )
    return {
        "counts": {
            "nir": len(nir),
            "joint": len(joint),
            "missing_joint": len(nir) - len(common),
            "selected": int(selected.sum()),
            "fail_AV_only": int((~av & stretch).sum()),
            "fail_stretch_only": int((av & ~stretch).sum()),
            "fail_both": int((~av & ~stretch).sum()),
        },
        "missing_joint_CIDs": nir.index[~nir.index.isin(joint.index)].tolist(),
        "nir_peak_minus_truth": stats(n.PKMJD - n.SIM_PKMJD),
        "joint_peak_minus_truth": stats(o.PKMJD - o.SIM_PKMJD),
        "nir_fit_minus_initializer": stats(n.PKMJD - n.PKMJDINI),
        "nir_distance_minus_truth": stats(n.DLMAG - n.SIM_DLMAG),
        "selected_nir_peak_minus_truth": stats(
            (n.PKMJD - n.SIM_PKMJD).to_numpy()[selected]
        ),
        "scope": "Descriptive accounting for source-matched archived simulation tables. Near-truth timing is a simulation procedure; it is not a measured real-data timing bias. Selection is explicit and missing rows remain in the ledger.",
    }
