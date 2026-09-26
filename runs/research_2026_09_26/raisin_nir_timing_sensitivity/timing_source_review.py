"""Read-only provenance check for the available author simulation timing chain."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
TIMING = ROOT / "runs/research_2026_09_26/raisin_simulation_assets/timing"
TREE = ROOT / "runs/research_2026_09_26/raisin_differential/mass-threshold/code-tree.json"
SIM = ROOT / "runs/research_2026_09_26/raisin_sign_source/sim/inputs/DES/sim_DES_SNOOPY.input"
SOURCE = OUT / "SNANA-v11_03c-source/src"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    acq = json.loads((TIMING / "acquisition.json").read_text())
    tree = json.loads(TREE.read_text())
    blobs = {r["path"]: r for r in tree["tree"] if r["type"] == "blob"}
    checks = []
    for row in acq:
        path = TIMING / Path(row["source_path"]).name
        data = path.read_bytes()
        blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
        checks.append({"source_path": row["source_path"], "sha256": sha(path), "git_blob": blob,
                       "closed": sha(path) == row["sha256"] and blob == row["git_blob"] == blobs[row["source_path"]]["sha"]})
    sim_text = SIM.read_text()
    snlc_sim = (SOURCE / "snlc_sim.c").read_text()
    textio = (SOURCE / "sntools_dataformat_text.c").read_text()
    nml_checks = {r["source_path"]: all(k in (TIMING / Path(r["source_path"]).name).read_text()
                                        for k in ("INISTP_SHAPE = 0.0", "INISTP_AV = 0.0",
                                                  "INISTP_PEAKMJD = 0.0", "FILTLIST_FIT = 'JH'"))
                  for r in acq}
    chain = {
        "nominal_sim_uses_SIMLIB_peak": "USE_SIMLIB_PEAKMJD: 1" in sim_text,
        "nominal_sim_sigma_search_peak_0p01": "GENSIGMA_SEARCH_PEAKMJD:  0.01" in sim_text,
        "sim_parser_aliases_search_sigma": '"GENSIGMA_PEAKMJD GENSIGMA_SEARCH_PEAKMJD"' in snlc_sim,
        "sim_gaussian_peak_rule": "PEAKMJD_SMEAR = GENLC.PEAKMJD + smear" in snlc_sim,
        "sim_stores_SEARCH_PEAKMJD": "SNDATA.SEARCH_PEAKMJD = GENLC.PEAKMJD_SMEAR" in snlc_sim,
        "text_writer_uses_SEARCH_PEAKMJD": 'SNDATA.SEARCH_PEAKMJD );' in textio and '"PEAKMJD:' in textio,
        "text_reader_maps_PEAKMJD": 'strcmp(word0,"PEAKMJD:")' in textio and '&SNDATA.SEARCH_PEAKMJD' in textio,
    }
    result = {"pass": all(x["closed"] for x in checks) and all(nml_checks.values()) and all(chain.values()),
              "author_fit_nml_blobs": checks, "author_nml_fixed_JH": nml_checks,
              "available_simulation_chain": chain,
              "qualification": "This available simulated-header path is not shown to have generated the released observed DES16 PEAKMJD headers or the effective published bias correction.",
              "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in
                                [TIMING / "acquisition.json", SIM, SOURCE / "snlc_sim.c",
                                 SOURCE / "sntools_dataformat_text.c", TREE, Path(__file__)]}}
    (OUT / "timing-source-review.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"pass": result["pass"], "chain": chain}, indent=2))


if __name__ == "__main__":
    main()
