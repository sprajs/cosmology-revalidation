"""Tag frozen holdout truth support on the predeclared optical grid."""
from pathlib import Path
import hashlib
import json
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DESIGN = ROOT / "runs/research_2026_09_26/astra_design/simulation_holdout_design"
sys.path.insert(0, str(ROOT / "scripts/salt_dust_audit"))
from snana_extinction import SnanaExtinction


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    lib = ROOT / "runs/salt_dust_audit/snana_extinction/libsnana_extinction.so"
    ext = SnanaExtinction(lib)
    waves = np.arange(2800., 8001., 10.)
    for rv in (.7, 1., 3.1, 6.):
        assert np.allclose(ext(waves, rv, option=-99),
                           ext(waves, rv, option=99, historical=True), rtol=0, atol=1e-10)
    rows = []
    for arm in ("P21", "G10"):
        cohort = pd.read_csv(DESIGN / f"{arm}-cohort.csv")
        for q in cohort.itertuples():
            av, rv = float(q.AV), float(q.RV)
            item = dict(arm=arm, CID=int(q.CID), generated_attempt_index=int(q.generated_attempt_index),
                        AV=av, RV=rv)
            if arm == "G10" and av == rv == -9:
                item.update(support_class="no_explicit_host_screen_sentinel",
                            A8000_historical=np.nan, min_grid_A_historical=np.nan)
            elif not np.isfinite(av) or av < 0:
                item.update(support_class="invalid_AV", A8000_historical=np.nan, min_grid_A_historical=np.nan)
            elif av == 0:
                item.update(support_class="zero_AV", A8000_historical=0., min_grid_A_historical=0.)
            elif not np.isfinite(rv) or not .01 <= rv <= 8:
                item.update(support_class="RV_outside_audit_range", A8000_historical=np.nan,
                            min_grid_A_historical=np.nan)
            else:
                vals = ext(waves, rv, av/rv, option=99, historical=True)
                item.update(support_class=("negative_passive_attenuation" if vals.min() < -1e-12
                                           else "nonnegative_declared_grid"),
                            A8000_historical=float(vals[-1]), min_grid_A_historical=float(vals.min()))
            rows.append(item)
    frame = pd.DataFrame(rows)
    assert len(frame) == 3477
    out = HERE / "truth-support-ledger.csv"
    assert not out.exists()
    frame.to_csv(out, index=False, float_format="%.17g")
    report = {"counts": {f"{a}:{s}":int(n) for (a,s),n in frame.groupby(["arm","support_class"]).size().items()},
              "grid_angstrom": [2800,8000,10], "cut": -1e-12,
              "inputs_sha256": {str(p.relative_to(ROOT)):digest(p) for p in
                  [Path(__file__),lib,ROOT/"scripts/salt_dust_audit/snana_extinction.py",
                   DESIGN/"manifest.json",*(DESIGN/f"{arm}-cohort.csv" for arm in ("P21","G10"))]},
              "ledger_sha256":digest(out)}
    (HERE/"truth-support-manifest.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report["counts"],indent=2))


if __name__ == "__main__":
    main()
