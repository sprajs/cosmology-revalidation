"""Execute native-fit gates for the frozen, coupled P21 noise variants.

This module performs no simulation or residual selection.  The frozen design
and the inherited native objective/flux gates define the experiment.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DESIGN = ROOT / "runs/research_2026_09_26/astra_design/noise_variant_design"
PROTO = ROOT / "runs/research_2026_09_26/astra_design/noise-variant-protocol.md"
HANDOFF = DESIGN / "handoff-manifest.json"
OUT = ROOT / "runs/research_2026_09_26/noise_variant_control"
NATIVE_DESIGN = OUT / "native_design"
BASE = ROOT / "scripts/research_2026_09_26/simulation_residual_control.py"
ARMS = ("P21_rho000", "P21_rho090", "P21_noisetrue120")
EXPECTED_PROTO = "ca48b11823ac3d74f32b707421e95d698718c03111d8ffc8258a605ce8c91016"
EXPECTED_HANDOFF = "6670646abca020b65814c0d57c29e56c0226d1783a3a39ad3a7e454c1e6942d2"


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1048576), b""):
            h.update(block)
    return h.hexdigest()


def preflight():
    assert digest(PROTO) == EXPECTED_PROTO
    assert digest(HANDOFF) == EXPECTED_HANDOFF
    frozen = json.loads(HANDOFF.read_text())["sha256"]
    # The design, source audit, and selected cohort are immutable at execution.
    for rel, expected in frozen.items():
        assert digest(ROOT / rel) == expected, rel
    return frozen


def native_module():
    spec = importlib.util.spec_from_file_location("pilot_native", BASE)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader
    spec.loader.exec_module(module)
    module.DESIGN = NATIVE_DESIGN
    module.OUT = OUT
    module.preflight = preflight
    return module


def prepare(arm):
    assert arm in ARMS
    preflight()
    path = OUT / arm
    path.mkdir(parents=True, exist_ok=False)
    NATIVE_DESIGN.mkdir(parents=True, exist_ok=True)
    ids = [int(x) for x in (DESIGN / f"{arm}-cids.txt").read_text().split()]
    assert len(ids) == len(set(ids)) == 256
    q = pd.read_csv(DESIGN / f"{arm}-cohort.csv")
    assert len(q) == 256 and set(q.CID) == set(ids)
    q = q.set_index("CID").loc[ids]
    assert np.isfinite(q[["t0_double", "mB_double", "x1_double", "c_double"]].to_numpy(float)).all()
    assert np.isfinite(q.x0.to_numpy(float)).all() and (q.x0 > 0).all()
    assert (q.generated_attempt_index.to_numpy() == pd.read_csv(DESIGN / "P21-cohort.csv").set_index("CID").loc[ids].generated_attempt_index.to_numpy()).all()
    native_cohort = q.reset_index().copy()
    native_cohort["field_first_phot"] = native_cohort.field
    native_cohort.to_csv(NATIVE_DESIGN / f"{arm}-cohort.csv", index=False, float_format="%.17g")
    seed = path / "measured_initial_values.FITRES"
    x0 = 10 ** ((10.635 - q.mB_double.to_numpy(float)) / 2.5)
    seed.write_text("VARNAMES: CID PKMJD x0 x1 c\n" + "".join(
        f"SN: {cid} {t0:.17g} {amp:.17g} {x1:.17g} {color:.17g}\n"
        for cid, t0, amp, x1, color in zip(ids, q.t0_double, x0, q.x1_double, q.c_double)))
    base = ROOT / "phase2/checkpoint/official-inputs/snana_forward_p21.nml"
    actual = ROOT / "phase2/official/inputs/snana_forward_p21.nml"
    assert digest(base) == digest(actual)
    original = base.read_text()
    assert "OPT_MWCOLORLAW = 99" in original and "LFIXPAR_ALL" not in original
    assert "OPT_SNCID_LIST" not in original and "SNCID_LIST_FILE" not in original
    assert "SNTABLE_LIST      = 'FITRES(text:host)'" in original
    old_version = "VERSION_PHOTOMETRY = 'PH2_pilot02_P21'"
    assert original.count(old_version) == 1
    version = f"PH2_pilot02_{arm}"
    text = original.replace(old_version, f"VERSION_PHOTOMETRY = '{version}'")
    text = text.replace("&SNLCINP", "&SNLCINP\n"
                        "    MXLC_PLOT = 1000\n"
                        "    OPT_SNCID_LIST = 2\n"
                        f"    SNCID_LIST_FILE = '{seed.resolve()}'", 1)
    text = text.replace("SNTABLE_LIST      = 'FITRES(text:host)'",
                        "SNTABLE_LIST      = 'FITRES(text:host) LCPLOT(text:col)'", 1)
    text = re.sub(r"(?m)^    TEXTFILE_PREFIX\s*=.*$",
                  f"    TEXTFILE_PREFIX = '{(path / 'fit').resolve()}'", text, count=1)
    text = text.replace("OPT_MWCOLORLAW = 99", "OPT_MWCOLORLAW = -99", 1)
    assert "OPT_MWCOLORLAW = -99" in text and "OPT_SETPKMJD   = 20" in text
    nml = path / "fit.nml"
    nml.write_text(text)
    info = {"arm": arm, "law": "approx_minus99", "law_option": -99,
            "cohort": "noise_variant_common256", "selected_ids": ids,
            "selected_count": 256, "opt_sncid_list": 2,
            "fitter_parameters_free": True, "truth_coordinates_used": False,
            "generated_attempt_indices": q.generated_attempt_index.astype(int).tolist(),
            "nml_sha256": digest(nml), "seed_sha256": digest(seed),
            "source_sha256": digest(Path(__file__)), "base_runner_sha256": digest(BASE),
            "base_nml_sha256": digest(base), "binary_sha256": digest(ROOT / "phase2/official/build/SNANA-audit-v3/bin/snlc_fit.exe"),
            "native_cohort_sha256": digest(NATIVE_DESIGN / f"{arm}-cohort.csv"),
            "frozen_cohort_sha256": digest(DESIGN / f"{arm}-cohort.csv"),
            "frozen_ids_sha256": digest(DESIGN / f"{arm}-cids.txt")}
    (path / "prepared.json").write_text(json.dumps(info, indent=2) + "\n")
    return {"arm": arm, "ids": 256, "version": version, "nml_sha256": digest(nml)}


def action(name, arm):
    assert arm in ARMS
    preflight()
    path = OUT / arm
    native = native_module()
    return {"fit": native.run_fit, "export": native.export,
            "gate": native.gate, "tangent": native.tangent_gate,
            "mean_edges": native.mean_edges}[name](path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=("prepare", "fit", "export", "gate", "tangent", "mean_edges"))
    ap.add_argument("arm", choices=ARMS)
    args = ap.parse_args()
    print(json.dumps(prepare(args.arm) if args.action == "prepare" else action(args.action, args.arm), indent=2))


if __name__ == "__main__":
    main()
