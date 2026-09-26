"""Independent direct-matrix validation of the three analytic amplitude arms."""
from pathlib import Path
import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    rows = [json.loads(line) for line in (OUT / "draw-results.jsonl").read_text().splitlines()]
    with np.load(OUT / "draws.npz", allow_pickle=False) as z:
        h = z["mean"]; d = z["diagonal_variance"]; v = z["loading"]
        y = z["flux"]; positive = z["positive"]; archive = z["archived_positive"]
    C = np.diag(d)+np.outer(v,v)
    maximum = {x: 0.0 for x in ("constrained", "unconstrained", "score")}
    worst = {}
    for i, row in enumerate(rows):
        assert row["draw"] == i
        masks = {"signed_gaussian": np.ones(74,bool),
                 "naive_positive_gaussian": positive[i],
                 "archive_mask_gaussian": archive}
        for name, mask in masks.items():
            S = C[np.ix_(mask,mask)]
            x = np.linalg.solve(S, h[mask])
            num = float(x @ y[i,mask])
            den = float(x @ h[mask])
            a = num/den
            estimated = row["estimators"][name]
            gaps = {"constrained": abs(float(np.clip(a,0,4))-estimated["amplitude"]),
                    "unconstrained": abs(a-estimated["unconstrained_amplitude"]),
                    "score": abs(num-den-estimated["score_at_truth"])}
            for key, gap in gaps.items():
                if gap > maximum[key]:
                    maximum[key] = gap; worst[key] = {"draw": i, "estimator": name}
    result = {"pass": all(x < 1e-12 for x in maximum.values()),
              "comparisons": len(rows)*3, "max_abs_gaps": maximum,
              "worst_locations": worst,
              "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in
                                [OUT / "draws.npz", OUT / "draw-results.jsonl", Path(__file__)]}}
    (OUT / "analytic-check.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
