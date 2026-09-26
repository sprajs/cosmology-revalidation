"""Read-only paired diagnostics for the immutable recovery-result stream."""
from pathlib import Path
import hashlib
import json
import math
import statistics

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
STREAM = OUT / "draw-results.jsonl"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mcse(values):
    return statistics.stdev(values)/math.sqrt(len(values)) if len(values) > 1 else None


def main():
    rows = [json.loads(line) for line in STREAM.read_text().splitlines()]
    assert len(rows) == 128 and [r["draw"] for r in rows] == list(range(128))
    names = ("naive_positive_gaussian", "archive_mask_gaussian",
             "joint_censored", "full_pattern_conditioned")
    paired_log = {}
    for name in names:
        valid = [(r["draw"], r["estimators"][name]["delta_log_distance"],
                  r["estimators"]["signed_gaussian"]["delta_log_distance"])
                 for r in rows if r["estimators"][name]["delta_log_distance"] is not None
                 and r["estimators"]["signed_gaussian"]["delta_log_distance"] is not None]
        diffs = [a-b for _,a,b in valid]
        paired_log[name + "_minus_signed"] = {"finite_n": len(valid),
                                               "excluded_draw_ids": sorted(set(range(128))-{i for i,_,_ in valid}),
                                               "mean": statistics.mean(diffs) if diffs else None,
                                               "mcse": mcse(diffs)}
    numerical = {}
    for name in ("joint_censored", "full_pattern_conditioned"):
        fits = [r["estimators"][name] for r in rows]
        numerical[name] = {
            "independent_adaptive_candidate_checks": sum(len(x["candidate_checks"]) for x in fits),
            "max_coarse_fine_loglikelihood_gap": max(x["coarse_fine_loglikelihood_gap"] for x in fits),
            "max_coarse_fine_amplitude_gap": max(x["coarse_fine_amplitude_gap"] for x in fits),
            "max_hermite_adaptive_gap": max(x["max_hermite_vs_adaptive_gap"] for x in fits),
            "max_adaptive_relative_error_estimate": max(x["max_adaptive_relative_error_estimate"] for x in fits),
            "max_score_halfstep_gap": max(x["score_step_gap"] for x in fits),
            "zero_draw_ids": [r["draw"] for r in rows if r["estimators"][name]["boundary"] == "zero"],
            "upper_draw_ids": [r["draw"] for r in rows if r["estimators"][name]["boundary"] == "upper"],
            "coarse_or_fine_tie_draw_ids": [r["draw"] for r in rows if len(r["estimators"][name]["coarse_ties"]) > 1
                                            or len(r["estimators"][name]["fine_ties"]) > 1],
        }
    result = {"scope": "Read-only postprocessing of frozen synthetic recovery fits",
              "paired_log_distance_finite_only": paired_log,
              "numerical": numerical,
              "input_sha256": {str(p.relative_to(ROOT)): sha(p) for p in
                               [STREAM, OUT / "summary.json", OUT / "full-run-state.json",
                                OUT / "execution-protocol.json", Path(__file__)]}}
    (OUT / "postprocess.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"paired_log_distance": paired_log, "numerical": numerical}, indent=2))


if __name__ == "__main__":
    main()
