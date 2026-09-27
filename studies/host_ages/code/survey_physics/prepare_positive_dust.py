#!/usr/bin/env python3
"""Separate physical-support campaign: intervene in dust before selection, not after."""
import datetime, gzip, json
import numpy as np
from extinction import fitzpatrick99
from common import *
from prepare import replace_key

src = DATA / "models/population_pdf/DES-SN5YR_P23/DES-SN5YR_DES_S3_W22.DAT.gz"
text = gzip.open(src, "rt").read()
start = text.index("VARNAMES: RV ")
end = text.index("VARNAMES: EBV ", start)
pdf = WORK / "inputs/W22-fixed-RV31.DAT"
pdf.write_text(
    "# Derived hypothesis: remove W22 mass-conditioned RV PDF; fixed RV=3.1 is defined in simulation input. Other released PDF blocks retained.\n"
    + text[:start]
    + text[end:]
)
# Optical domain is prescribed, and a wider UV/NIR positivity check is also recorded.
wavelength = np.linspace(1000, 30000, 2901).astype(float)
a = fitzpatrick99(wavelength, 3.1, 3.1)
assert a.min() > 0
jobs = []
for pool, arm, n, seed in [
    ("train", "nominal", 6000, 270921),
    ("eval", "nominal", 3000, 270922),
    ("eval", "age", 3000, 270922),
    ("train", "age", 6000, 270921),
]:
    name = f"SPP_{pool.upper()}_{arm.upper()}"
    text = (WORK / "inputs" / f"SPH_{pool.upper()}_{arm.upper()}.input").read_text()
    for key, val in {
        "GENVERSION": name,
        "GENPREFIX": name,
        "NGENTOT_LC": n,
        "RANSEED": seed,
        "GENPDF_FILE": pdf,
        "GENPEAK_RV": 3.1,
        "GENSIGMA_RV": "0 0",
        "GENRANGE_RV": "3.1 3.1",
    }.items():
        text = replace_key(text, key, val)
    path = WORK / "inputs" / (name + ".input")
    path.write_text(text)
    jobs.append(
        {
            "name": name,
            "pool": pool,
            "arm": arm,
            "attempts": n,
            "seed": seed,
            "input": str(path),
            "input_sha256": sha(path),
        }
    )
r = {
    "frozen_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "reason": "Released low-RV dust has negative optical extinction. This independent physical-support population fixes RV=3.1 before generating all photons, selection and fitting, and retains W22 mocked ages/width plus the released EBV and intrinsic-colour distributions. This is a changed hypothesis, not a repair of empirical RV inference.",
    "code_sha256": sha(__file__),
    "pdf_source_sha256": sha(src),
    "derived_pdf_sha256": sha(pdf),
    "RV": 3.1,
    "positive_extinction_check_wavelength_A": [1000, 30000],
    "min_A_per_EBV1": float(a.min()),
    "jobs": jobs,
    "analysis": "Apply frozen 40-neighbor training-only correction with independent seeds; no post-selection extinction screen. Same quality/classifier/bins/bootstrap as main benchmark; Independent age-injected training refits standardization and correction; paired training seeds match the nominal population.",
}
(RESULTS / "positive-dust-campaign.json").write_text(json.dumps(r, indent=2) + "\n")
