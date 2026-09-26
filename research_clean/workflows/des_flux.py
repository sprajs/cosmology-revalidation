"""Independent fits to fixed exported DES flux objectives."""

import numpy as np
from lib.flux_engine import Engine, read_inputs
from lib.flux_fit import Case
from lib.records import write_json

DEFAULTS = {
    "objects": "all",
    "interpolations": ["sncosmo", "scipy"],
    "two_starts": True,
    "heldout_epochs": True,
}


def run(out, cfg):
    import hashlib

    engine = Engine()
    _, _, native = read_inputs()
    cids = sorted(native) if cfg["objects"] == "all" else list(map(str, cfg["objects"]))
    if not cids or len(set(cids)) != len(cids) or any(c not in native for c in cids):
        raise ValueError("Objects must be unique IDs from the frozen 12-object cohort")
    if not cfg["interpolations"] or not set(cfg["interpolations"]) <= {
        "sncosmo",
        "scipy",
    }:
        raise ValueError("Unknown interpolation")
    full, heldout = [], []
    for cid in cids:
        case = Case(engine, cid)
        for interpolation in cfg["interpolations"]:
            fit, model, jacobian, covariance = case.fit(
                interpolation=interpolation, two_starts=cfg["two_starts"]
            )
            fit["interpolation"] = interpolation
            fit["delta_mB_x1_c_t0"] = (
                (np.array(fit["x"]) - case.x0)
                * np.array([-2.5 / np.log(10), 1.0, 1.0, 1.0])
            ).tolist()
            full.append(fit)
        if cfg["heldout_epochs"]:
            testmask = np.array(
                [
                    int(
                        hashlib.sha256(
                            f"independent-flux-v1|{cid}|{int(np.floor(t + 0.5))}".encode()
                        ).hexdigest(),
                        16,
                    )
                    % 5
                    == 0
                    for t in case.t
                ]
            )
            train, test = np.flatnonzero(~testmask), np.flatnonzero(testmask)
            if not len(train) or not len(test):
                raise ValueError("Frozen split has no training or held-out epochs")
            for family in ["salt", "f99", "phase_colour"]:
                fit, prediction, jacobian, covariance = case.fit(
                    family, train, two_starts=cfg["two_starts"]
                )
                fit["prediction"] = case.predictive(
                    fit, prediction, jacobian, covariance, train, test
                )
                fit["train_rows"] = train.tolist()
                fit["test_rows"] = test.tolist()
                heldout.append(fit)
        write_json(out / "fits.json", full)
        write_json(out / "heldout.json", heldout)
    failures = [
        x
        for x in full
        if not x["success"] or x["boundary"] or x["hessian_covariance"] is None
    ]
    return {
        "objects": cids,
        "fits": len(full),
        "invalid_full_fits": len(failures),
        "heldout_fits": len(heldout),
        "valid_heldout_fits": sum(x["prediction"]["valid"] for x in heldout),
        "max_abs_delta_mB_x1_c_t0": np.max(
            abs(np.array([x["delta_mB_x1_c_t0"] for x in full])), axis=0
        ).tolist(),
        "scope": "Frozen calibrated observations, mask, redshift, calibration and covariance; an independent fitter and integrator. This does not reconstruct detector pixels, certify native error construction, or identify a dust law from a winning predictive score.",
    }
