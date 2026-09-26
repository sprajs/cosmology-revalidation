"""Independent hash and ledger arithmetic checks for common-classifier score."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / "score"


def digest(p):
    h = hashlib.sha256()
    with p.open("rb") as stream:
        for part in iter(lambda: stream.read(1048576), b""):
            h.update(part)
    return h.hexdigest()


def main():
    sources = json.loads((OUT / "manifest.json").read_text())
    for rel, expected in {**sources["inputs_sha256"], **sources["outputs_sha256"]}.items():
        assert digest(ROOT / rel) == expected, rel
    pred = pd.read_csv(HERE / "inference/simulation-probabilities.csv")
    realpred = pd.read_csv(HERE / "inference/real-probabilities.csv")
    old = pd.read_csv(ROOT / "runs/research_2026_09_26/classifier_residual_sensitivity/membership_and_scores.csv")
    ids = set(realpred.loc[realpred.gt999, "CID"])
    assert len(ids) == 1006
    assert ids == set(old.loc[old.retained_both_gt_0p999, "CID"])
    assert pred.groupby("arm").size().to_dict() == {"G10": 256, "P21": 256}
    assert pred.groupby("arm").gt999.sum().to_dict() == {"G10": 151, "P21": 152}
    src = pd.read_csv(ROOT / "runs/research_2026_09_26/simulation_residual_control/score/object-scores.csv")
    src = src.merge(pred[["arm", "CID", "gt999"]], on=["arm", "CID"], validate="many_to_one")
    result = pd.read_csv(OUT / "transport-summary.csv")
    ledger = pd.read_csv(OUT / "supported-object-ledger.csv")
    field = pd.read_csv(OUT / "field-contributions.csv")
    boot = pd.read_csv(OUT / "bootstrap.csv.gz")
    assert len(result) == 40 and len(boot) == 40000
    groupcols = ("arm", "law", "transport", "stage", "support_source")
    for row in result.itertuples():
        mask = np.logical_and.reduce([ledger[c] == getattr(row, c) for c in groupcols])
        weights = ledger.loc[mask]
        assert len(weights) == row.supported_sim and weights.CID.is_unique
        q = src.loc[(src.arm == row.arm) & (src.law == row.law)].set_index("CID")
        q = q.loc[weights.CID]
        w = weights.base_weight.to_numpy(float)
        assert abs(w.sum()-row.supported_real) < 1e-8
        assert abs(w @ q.matched_filter.to_numpy(float)-row.weighted_a) < 1e-10
        assert abs(w @ q.information.to_numpy(float)-row.weighted_I) < 1e-10
        assert abs(w @ q.fixed_gain.to_numpy(float)-row.weighted_G) < 1e-10
        assert abs(row.weighted_a/row.weighted_I-row.transported_amplitude) < 1e-12
        assert abs(row.weighted_G/row.supported_real-row.transported_mean_gain) < 1e-12
        fmask = np.logical_and.reduce([field[c] == getattr(row,c) for c in groupcols])
        f = field.loc[fmask]
        assert f.n.sum() == row.supported_sim
        for name in ("a", "I", "G"):
            assert abs(f["weighted_" + name].sum()-getattr(row,"weighted_" + name)) < 1e-10
        bmask = np.logical_and.reduce([boot[c] == getattr(row,c) for c in groupcols])
        b = boot.loc[bmask]
        assert len(b) == 1000 and (~b.valid).sum() == row.bootstrap_failures
        valid = b.loc[b.valid]
        assert abs(valid.amplitude.quantile(.025)-row.bootstrap_amplitude_ci_low) < 1e-10
        assert abs(valid.amplitude.quantile(.975)-row.bootstrap_amplitude_ci_high) < 1e-10
    report = {"verified_manifest_files": len(sources["inputs_sha256"])+len(sources["outputs_sha256"]),
              "real_intersection": len(ids), "sim_classifier_counts": {"P21": 152, "G10": 151},
              "transport_rows": len(result), "bootstrap_rows": len(boot),
              "result": "passed"}
    (HERE / "score-independent-verification.json").write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
