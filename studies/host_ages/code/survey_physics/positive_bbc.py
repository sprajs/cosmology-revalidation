#!/usr/bin/env python3
"""Native1D BBC on the positive-dust mocks; distinct from the independent predictor."""
import json, re
import numpy as np
from common import *
from native_checks import filter_fitres, command


def main():
    source = WORK / "bbc/redshift_1D/bbc.input"
    pairs = [
        ("nominal", "SPP_TRAIN_NOMINAL", "SPP_EVAL_NOMINAL"),
        ("age_nominal_training", "SPP_TRAIN_NOMINAL", "SPP_EVAL_AGE"),
        ("age_matched_training", "SPP_TRAIN_AGE", "SPP_EVAL_AGE"),
    ]
    result = {
        "scope": "Native opt_biascor=2 redshift-only benchmark with alpha/beta/mass-step fitted on each evaluation arm. Pure Ia, no BEAMS contamination. Native COV is conditional statistical M0DIF covariance, not correction-training Monte Carlo or survey systematics. Separate from the 40-neighbor predictor and not combined with it.",
        "code_sha256": sha(__file__),
        "source_input_sha256": sha(source),
        "arms": {},
    }
    for label, t, e in pairs:
        train = filter_fitres(t)
        ev = filter_fitres(e)
        text = source.read_text()
        for key, value in [("datafile", ev), ("simfile_biascor", train)]:
            text = re.sub(r"(?m)^" + key + r"=.*$", key + "=" + str(value), text)
        folder = WORK / "bbc" / ("positive_" + label)
        folder.mkdir(exist_ok=True)
        path = folder / "bbc.input"
        path.write_text(text)
        r = command(WORK / "SNANA/bin/SALT2mu.exe", [path], folder, "bbc")
        log = (folder / "bbc.log").read_text()
        r["graceful"] = r["returncode"] == 0 and "Done." in log[-1000:]
        r["input_sha256"] = sha(path)
        if r["graceful"]:
            vals = [
                line
                for line in (folder / "bbc.COV").read_text().splitlines()
                if line and not line.startswith("#")
            ]
            n = int(vals[0])
            c = np.array([float(v) for v in vals[1:]]).reshape(n, n)
            assert np.isfinite(c).all() and np.max(abs(c - c.T)) < 1e-8
            r["covariance_dimension"] = n
            r["covariance_min_eigenvalue"] = float(np.linalg.eigvalsh(c).min())
            r["fit_parameters"] = {
                m.group(1): {
                    "value": float(m.group(2)),
                    "conditional_error": float(m.group(3)),
                }
                for m in re.finditer(
                    r"par\s+\d+\s+\(\s*\d+\)\s+(alpha0|beta0|gamma0)\s+([\deE+.-]+)\s+\+/-\s+([\deE+.-]+)",
                    log,
                )
            }
            r["retained_count"] = int(re.search(r"Wrote (\d+) SN", log).group(1))
            r["outputs"] = {
                n: sha(folder / n) for n in ["bbc.FITRES", "bbc.M0DIF", "bbc.COV"]
            }
        else:
            r["failure_tail"] = log[-3000:]
        result["arms"][label] = r
    (RESULTS / "positive-dust-native-bbc.json").write_text(
        json.dumps(result, indent=2) + "\n"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
