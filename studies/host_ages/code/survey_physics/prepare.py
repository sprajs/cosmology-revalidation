#!/usr/bin/env python3
"""Prepare explicit, source-traceable DES/W22 interventions before selection."""
import datetime
import gzip
import json
import re
from common import ROOT, HERE, WORK, RESULTS, DATA, ARCHIVE, sha


def replace_key(text, key, value):
    pattern = r"(?m)^" + re.escape(key) + r":[^\n]*"
    assert len(re.findall(pattern, text)) <= 1, key
    return (
        re.sub(pattern, key + ": " + str(value), text)
        if re.search(pattern, text)
        else text + "\n" + key + ": " + str(value) + "\n"
    )


def main():
    inp = WORK / "inputs"
    inp.mkdir(exist_ok=True)
    (WORK / "simulations").mkdir(exist_ok=True)
    (WORK / "fits").mkdir(exist_ok=True)
    # Valid grid-node probabilities, documented separately from legacy saturation.
    effsrc = DATA / "models/searcheff/SEARCHEFF_zHOST_DES-SN5YR_OBS.DAT.gz"
    lines = gzip.open(effsrc, "rt").read().splitlines()
    bad = []
    for i, line in enumerate(lines):
        if line.startswith("HOSTEFF:"):
            cells = line.split()
            value = float(cells[3])
            if not 0 <= value <= 1:
                bad.append({"line": i + 1, "original_probability": value})
            cells[3] = f"{min(1, max(0, value)):.6f}"
            lines[i] = " ".join(cells)
    eff = inp / "host-efficiency-W22-bounded.DAT"
    eff.write_text("\n".join(lines) + "\n")
    # W22 SN_age is a synthetic progenitor delay in Gyr, not mean_age (Myr).
    hostsrc = DATA / "simlib/DES/DES-SN5YR_W22.HOSTLIB.gz"
    hostage = inp / "W22-age-brightness.HOSTLIB.gz"
    total = 0
    amin, amax = float("inf"), -float("inf")
    with (
        gzip.open(hostsrc, "rt") as src,
        gzip.open(hostage, "wt", compresslevel=1) as out,
    ):
        cols = None
        for line in src:
            if line.startswith("VARNAMES:"):
                cols = line.split()[1:]
                assert "SNMAGSHIFT" not in cols
                agecol = cols.index("SN_age") + 1
                line = line.rstrip() + " SNMAGSHIFT\n"
            elif line.startswith("GAL:"):
                a = float(line.split()[agecol])
                assert 0 <= a <= 14, a
                amin, amax = min(amin, a), max(amax, a)
                line = line.rstrip() + f" {-0.03 * (a - 3):.8f}\n"
                total += 1
            out.write(line)
    base = DATA / "sample_input_files/DES-SN5YR/base_files/sim"
    model = (base / "sim_ia_salt_des5yr.input").read_text()
    survey = (base / "sim_des5yr_survey.input").read_text().split("SIMGEN_DUMP:")[0]
    model = re.sub(r"(?m)^INPUT_FILE_INCLUDE:.*$", survey, model)
    pdfbase = DATA / "models/population_pdf/DES-SN5YR_P23"
    # Keep exactly W22's width-age map while substituting the released sys1
    # intrinsic-colour/dust/alpha/beta/MAG_OFFSET block as one named alternative.
    w22 = gzip.open(pdfbase / "DES-SN5YR_DES_S3_W22.DAT.gz", "rt").read()
    dust = gzip.open(pdfbase / "DES-SN5YR_DES_S3_P23sys1.DAT.gz", "rt").read()
    w22_width = w22[w22.index("VARNAMES: SALT2x1 ") : w22.index("VARNAMES: SALT2c ")]
    dust_tail = dust[dust.index("VARNAMES: SALT2c ") :].split("VARNAMES: SALT2x1 ")[0]
    dust_pdf = inp / "W22-width-P23sys1-dust.DAT"
    dust_pdf.write_text(
        "# Explicit alternative: W22 age-width map + released P23sys1 colour/dust block.\n"
        + w22_width
        + dust_tail
    )
    columns = "CID GENZ GALZTRUE LIBID RA DEC MWEBV MU PEAKMJD SNRMAX SNRMAX2 SNRMAX3 NOBS TRESTMIN TRESTMAX CUTMASK SIM_EFFMASK GALID LOGMASS_TRUE LOGSFR SN_age mean_age SALT2mB SALT2x1 SALT2c SALT2alpha SALT2beta AV RV MAGSMEAR_COH".split()
    settings = {
        "KCOR_FILE": "$SNDATA_ROOT/kcor/Dovekie/calib_DES-SN5YR_DES.fits.gz",
        "GENMODEL": "$SNDATA_ROOT/models/SALT3/SALT3.DOVEKIE",
        "GENMODEL_EXTRAP_LATETIME": str(
            WORK / "inputs/SNIa_Extrap_LateTime_2expon.TEXT"
        ),
        "GENPDF_FILE": str(pdfbase / "DES-SN5YR_DES_S3_W22.DAT.gz"),
        "HOSTLIB_FILE": str(hostsrc),
        "HOSTLIB_DZTOL": "0.015 0 0",
        "GENRANGE_REDSHIFT": "0.05 1.2",
        "HOSTLIB_STOREVAR": "LOGMASS_TRUE,LOGSFR,SN_age,mean_age",
        "SEARCHEFF_zHOST_FILE": str(eff),
        "OMEGA_MATTER": ".315",
        "OMEGA_LAMBDA": ".685",
        "FORMAT_MASK": "32",
        "PATH_SNDATA_SIM": str(WORK / "simulations"),
    }
    for key, value in settings.items():
        model = replace_key(model, key, value)
    model += "\nSIMGEN_DUMPALL: " + str(len(columns)) + "\n " + " ".join(columns) + "\n"
    jobs = []
    for pool, arm, n, seed in [
        ("train", "nominal", 12000, 270911),
        ("eval", "nominal", 6000, 270912),
        ("eval", "age", 6000, 270912),
        ("eval", "dust", 6000, 270912),
        ("train", "age", 12000, 270911),
    ]:
        name = f"SPH_{pool.upper()}_{arm.upper()}"
        text = model
        for key, value in {
            "GENVERSION": name,
            "GENPREFIX": name,
            "NGENTOT_LC": n,
            "RANSEED": seed,
        }.items():
            text = replace_key(text, key, value)
        if arm == "age":
            text = replace_key(text, "HOSTLIB_FILE", hostage)
            text = replace_key(text, "HOSTLIB_MSKOPT", 262)
            # Bit4 enables the shift but does not automatically load this column.
            # Do not enable bit64: that also loads W22's unused mock Av column.
            text = replace_key(
                text,
                "HOSTLIB_STOREVAR",
                "LOGMASS_TRUE,LOGSFR,SN_age,mean_age,SNMAGSHIFT",
            )
        if arm == "dust":
            text = replace_key(text, "GENPDF_FILE", dust_pdf)
        path = inp / (name + ".input")
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
    sources = [
        hostsrc,
        effsrc,
        pdfbase / "DES-SN5YR_DES_S3_W22.DAT.gz",
        pdfbase / "DES-SN5YR_DES_S3_P23sys1.DAT.gz",
        base / "sim_ia_salt_des5yr.input",
        base / "sim_des5yr_survey.input",
    ]
    result = {
        "prepared_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "code_sha256": sha(__file__),
        "protocol_sha256": sha(HERE / "protocol.json"),
        "campaign_selected_after_engineering_benchmark_before_scientific_contrasts": True,
        "jobs": jobs,
        "mock_host_rows": total,
        "mock_SN_age_Gyr_range": [amin, amax],
        "out_of_probability_range_grid_nodes": bad,
        "inputs": [{"path": str(p), "sha256": sha(p)} for p in sources],
        "derived": [
            {"path": str(p), "sha256": sha(p)} for p in [hostage, eff, dust_pdf]
        ],
        "generated_attempts_total": sum(j["attempts"] for j in jobs),
        "scientific_limits": [
            "W22 host properties and ages are simulated assumptions.",
            "Dust alternative changes the complete released P23sys1 block, not only RV.",
            "No exact-original classifier equivalence or SALT retraining.",
            "Bounded probability interpolation is explicit and requires sensitivity to legacy semantics.",
        ],
    }
    (RESULTS / "campaign.json").write_text(json.dumps(result, indent=2) + "\n")
    print(
        json.dumps(
            {k: v for k, v in result.items() if k not in ["inputs", "derived", "jobs"]},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
