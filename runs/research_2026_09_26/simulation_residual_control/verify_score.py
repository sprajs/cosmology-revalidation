"""Independent arithmetic and file-integrity checks for the frozen score output."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCORE = HERE / "score"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> None:
    manifest = json.loads((SCORE / "manifest.json").read_text())
    for rel, wanted in {**manifest["inputs_sha256"], **manifest["outputs_sha256"]}.items():
        assert digest(ROOT / rel) == wanted, rel
    scores = pd.read_csv(SCORE / "object-scores.csv")
    summary = pd.read_csv(SCORE / "transport-summary.csv")
    fields = pd.read_csv(SCORE / "field-contributions.csv")
    boot = pd.read_csv(SCORE / "bootstrap.csv.gz")
    with np.load(SCORE / "sufficient-arrays.npz", allow_pickle=False) as arrays:
        assert np.array_equal(arrays["CID"], scores.CID.to_numpy())
        # The archived arm/law arrays have object dtype; never unpickle them.
        # The numeric CID order and manifest-hashed CSV establish row identity.
        c = np.asarray(manifest["frozen_coefficients"], float)
        aa = arrays["u"] @ c
        ii = np.einsum("i,nij,j->n", c, arrays["F"], c)
    assert np.max(abs(aa - scores.matched_filter)) < 1e-12
    assert np.max(abs(ii - scores.information)) < 1e-12
    assert np.max(abs(aa - ii / 2 - scores.fixed_gain)) < 1e-12
    assert scores.groupby(["arm", "law"]).size().eq(255).all()
    assert scores.groupby(["arm", "CID"]).law.nunique().eq(2).all()
    assert len(boot) == 1000 * len(summary)
    checks = []
    for row in summary.itertuples():
        q = scores.loc[(scores.arm == row.arm) & (scores.law == row.law)]
        if row.stage == "archived_quality":
            q = q.loc[q.archived_basic_quality]
        elif row.stage == "new_quality":
            q = q.loc[q.new_basic_quality]
        elif row.stage == "p21_grid_nonnegative" and row.arm == "P21":
            q = q.loc[q.support_class == "nonnegative_declared_grid"]
        assert len(q) == row.selected
        assert abs(q.matched_filter.sum() - row.raw_a) < 1e-10
        assert abs(q.information.sum() - row.raw_I) < 1e-10
        assert abs(q.fixed_gain.sum() - row.raw_G) < 1e-10
        cells = set(row.supported_cell_ids.split("|"))
        qcells = q.cell if row.transport == "field_z" else q.cell + "_" + np.where(q.SNRMAX1_archived < 15, "lt15", "ge15")
        q = q.loc[qcells.isin(cells)].copy()
        q["transport_cell"] = qcells.loc[q.index]
        assert len(q) == row.supported_objects
        assert abs(q.matched_filter.sum() - row.supported_unweighted_a) < 1e-10
        assert abs(q.information.sum() - row.supported_unweighted_I) < 1e-10
        assert abs(q.fixed_gain.sum() - row.supported_unweighted_G) < 1e-10
        f = fields.loc[(fields.arm == row.arm) & (fields.law == row.law)
                       & (fields.stage == row.stage) & (fields.transport == row.transport)]
        assert f.n.sum() == row.supported_objects
        for key in ("raw_a", "raw_I", "raw_G"):
            expected = getattr(row, "supported_unweighted_" + key[-1].upper() if key != "raw_a" else "supported_unweighted_a")
            assert abs(f[key].sum() - expected) < 1e-10
        assert abs(f.weighted_a.sum() / f.weighted_I.sum() - row.transported_amplitude) < 1e-10
        assert abs(f.weighted_G.sum() / row.supported_real_objects - row.transported_mean_gain) < 1e-10
        b = boot.loc[(boot.arm == row.arm) & (boot.law == row.law)
                     & (boot.stage == row.stage) & (boot.transport == row.transport)]
        assert len(b) == 1000
        assert int((~b.valid).sum()) == row.bootstrap_failures
        good = b.loc[b.valid]
        assert abs(good.amplitude.quantile(.025) - row.bootstrap_amplitude_ci_low) < 1e-10
        assert abs(good.amplitude.quantile(.975) - row.bootstrap_amplitude_ci_high) < 1e-10
        checks.append(dict(arm=row.arm, law=row.law, stage=row.stage,
                           transport=row.transport, selected=row.selected,
                           supported=row.supported_objects, real=row.supported_real_objects))
    report = {"input_and_output_hashes_verified": len(manifest["inputs_sha256"]) + len(manifest["outputs_sha256"]),
              "objects": len(scores), "transport_rows": len(summary),
              "bootstrap_rows": len(boot), "max_score_arithmetic_error": float(max(
                  np.max(abs(aa - scores.matched_filter)),
                  np.max(abs(ii - scores.information)),
                  np.max(abs(aa - ii / 2 - scores.fixed_gain)))),
              "rows": checks}
    (HERE / "score-independent-verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "rows"}, indent=2))


if __name__ == "__main__":
    main()
