#!/usr/bin/env python3
"""Parametric calibration of the TITAN/ZTF *working statistical model*.

The observed host summaries/design and selected sample are fixed. This checks
optimization, nominal coverage and conditional power, not physical age recovery,
selection, joint SED inference or astrophysical identifiability.
"""
import argparse, datetime, json
from pathlib import Path
import numpy as np
import pandas as pd
import titan, ztf
from analyse import ROOT, HERE, sha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=300)
    ap.add_argument("--work", type=Path, default=ROOT / ".work/environment-validation")
    args = ap.parse_args()
    source = args.work / "titan-ztf-join.csv"
    q = pd.read_csv(source)
    q = q[q.selected_titan].reset_index(drop=True)
    fold = q.fold.copy()
    hg = q.host_group.copy()
    q = ztf.prepare(q)
    q["fold"] = fold
    q["host_group"] = hg
    controls, bc, _, sig = ztf.fit(
        q, titan.CONTROL, titan.matrix, lambda d, b, n: titan.variance(d, b, n, False)
    )
    names = titan.CONTROL + ["host_age"]
    x = titan.matrix(q, names)
    rng = np.random.default_rng(20260928)
    base = np.r_[bc, 0.0]
    v = titan.variance(q, base, names, False) + sig**2
    # The same noise draws serve the paired zero and nonzero conditional truths.
    noises = rng.normal(size=(args.draws, len(q))) * np.sqrt(v)
    out = {}
    rows = []
    for truth in [0.0, -0.03]:
        true = base.copy()
        true[-1] = truth
        bs = []
        ses = []
        failed = []
        for rep, noise in enumerate(noises):
            mock = q.copy()
            mock["y"] = x @ true + noise
            try:
                fit, _, _, _ = ztf.fit(
                    mock,
                    names,
                    titan.matrix,
                    lambda d, b, n: titan.variance(d, b, n, False),
                )
                b = fit["parameters"]["host_age"]
                se = fit["standard_errors"]["host_age"]
                bs.append(b)
                ses.append(se)
                rows.append({"truth": truth, "replicate": rep, "estimate": b, "se": se})
            except RuntimeError as e:
                failed.append({"replicate": rep, "error": str(e)})
        bs = np.array(bs)
        ses = np.array(ses)
        coverage = float(np.mean(abs(bs - truth) <= 1.96 * ses))
        detect = float(np.mean(abs(bs) > 1.96 * ses))
        out[str(truth)] = {
            "replicates_requested": args.draws,
            "converged": len(bs),
            "failures": failed,
            "mean_estimate": float(bs.mean()),
            "empirical_sd": float(bs.std(ddof=1)),
            "mean_nominal_se": float(ses.mean()),
            "coverage95": coverage,
            "coverage_binomial_mc_se": float(
                np.sqrt(coverage * (1 - coverage) / len(bs))
            ),
            "two_sided_detection_fraction": detect,
            "detection_binomial_mc_se": float(np.sqrt(detect * (1 - detect) / len(bs))),
        }
        print(truth, json.dumps(out[str(truth)]), flush=True)
    pd.DataFrame(rows).to_csv(
        args.work / "titan-working-model-injections.csv", index=False
    )
    result = {
        "created_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(__file__),
        "dependencies_sha256": {
            p.name: sha(p)
            for p in [HERE / "titan.py", HERE / "ztf.py", HERE / "analyse.py"]
        },
        "input_join_sha256": sha(source),
        "seed": 20260928,
        "sample": len(q),
        "injections": out,
        "scope": "Parametric statistical self-calibration only: fixed observed posterior-median predictors, independent Gaussian noise from fitted working likelihood; no hidden age/dust, selection or cross-object calibration perturbation.",
    }
    (
        ROOT / "studies/host_ages/results/environment_validation/titan-calibration.json"
    ).write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
