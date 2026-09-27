#!/usr/bin/env python3
"""Recorded native BBC support, bounded-efficiency and optimizer sensitivity checks."""
import argparse
import json
import re
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor
from common import *
from analyze import load, selected, fitres


def command(exe, args, folder, stem):
    folder.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    with (folder / (stem + ".log")).open("w") as out:
        run = subprocess.run(
            [str(exe), *map(str, args)],
            cwd=folder,
            env=native_env(),
            stdout=out,
            stderr=subprocess.STDOUT,
        )
    return {
        "returncode": run.returncode,
        "seconds": time.monotonic() - started,
        "log_sha256": sha(folder / (stem + ".log")),
    }


def legacy():
    folder = WORK / "legacy-efficiency"
    folder.mkdir(exist_ok=True)
    text = (
        (WORK / "inputs/SPH_EVAL_NOMINAL.input")
        .read_text()
        .replace("SPH_EVAL_NOMINAL", "SPH_EVAL_LEGACY")
    )
    text = re.sub(
        r"(?m)^SEARCHEFF_zHOST_FILE:.*$",
        "SEARCHEFF_zHOST_FILE: "
        + str(DATA / "models/searcheff/SEARCHEFF_zHOST_DES-SN5YR_OBS.DAT.gz"),
        text,
    )
    path = folder / "legacy.input"
    path.write_text(text)
    exe = WORK / "SNANA-legacy-eff/bin/snlc_sim.exe"
    r = command(exe, [path], folder, "generation")
    from analyze import dump

    n = dump("SPH_EVAL_NOMINAL")
    l = dump("SPH_EVAL_LEGACY")
    cols = [
        "CID",
        "LIBID",
        "GENZ",
        "PEAKMJD",
        "GALID",
        "SN_age",
        "AV",
        "RV",
        "SALT2x1",
        "SALT2c",
    ]
    common = n.merge(
        l, on=["CID", "LIBID"], suffixes=("_nominal", "_legacy"), validate="one_to_one"
    )
    r.update(
        input_sha256=sha(path),
        executable_sha256=sha(exe),
        compared_occurrences=len(common),
        nominal_attempts=len(n),
        legacy_attempts=len(l),
    )
    r["common_latents_equal"] = {
        c: bool((common[c + "_nominal"] == common[c + "_legacy"]).all())
        for c in cols[2:]
    }
    r["accepted_nominal"] = int(n.FLAG_ACCEPT.sum())
    r["accepted_legacy"] = int(l.FLAG_ACCEPT.sum())
    r["acceptance_disagreements"] = int(
        (common.FLAG_ACCEPT_nominal != common.FLAG_ACCEPT_legacy).sum()
    )
    r["clipped_runtime_count"] = (
        (folder / "generation.log").read_text().count("SURVEY_PHYSICS_EFF_BOUND")
    )
    (RESULTS / "legacy-efficiency.json").write_text(json.dumps(r, indent=2) + "\n")
    print(r)


def filter_fitres(name):
    f, _ = load(name)
    ids = set(selected(f, "primary").index)
    source = WORK / "fits" / name / "fit.FITRES.TEXT"
    dest = WORK / "bbc" / (name + ".FITRES")
    dest.parent.mkdir(exist_ok=True)
    # Preserve released/native column semantics; only declared quality+pIa cohort filtering.
    dest.write_text(
        "\n".join(
            line
            for line in source.read_text().splitlines()
            if not line.startswith("SN:") or line.split()[1] in ids
        )
        + "\n"
    )
    return dest


def bbc():
    train = filter_fitres("SPH_TRAIN_NOMINAL")
    ev = filter_fitres("SPH_EVAL_NOMINAL")
    source = DATA / "sample_input_files/DES-SN5YR/base_files/bbc/BBC_des5yr.input"
    base = source.read_text().split("#END_YAML", 1)[1]
    # Pure Ia benchmark: omit contaminated BEAMS prior and classifier-likelihood factor.
    omitted = [
        "cid_reject_file",
        "varname_pIa",
        "simfile_ccprior",
        "idsurvey_list_probcc0",
    ]
    for key in omitted:
        base = re.sub(r"(?m)^" + key + r"=.*\n", "", base)
    settings = {
        "datafile": str(ev),
        "simfile_biascor": str(train),
        "surveygroup_biascor": "'DES(zbin=0.075)'",
        "u13": "0",
        "p13": "0",
        "zmin": ".05",
        "u2": "1",
    }
    for key, value in settings.items():
        base = re.sub(r"(?m)^" + key + r"=.*$", key + "=" + value, base)
    outcomes = {}
    for label, opt in [
        ("released_4D", 4336),
        ("redshift_1D", 2),
        ("bounds_fixed_4D", 4336),
        ("bounds_fixed_snr50_4D", 4336),
    ]:
        folder = WORK / "bbc" / label
        folder.mkdir(exist_ok=True)
        text = re.sub(r"(?m)^opt_biascor=.*$", "opt_biascor=" + str(opt), base)
        text = re.sub(r"(?m)^prefix=.*$", "prefix=bbc", text)
        if opt == 2:
            # 1D cohort-average benchmark cannot support separate sparse field cells.
            text = re.sub(r"(?m)^fieldGroup_biascor=.*\n", "", text)
        if "snr50" in label:
            text += "\nsnrmin_sigint_biascor=50\n"
        path = folder / "bbc.input"
        path.write_text(text)
        exe = WORK / (
            "SNANA-legacy-eff/bin/SALT2mu.exe"
            if label.startswith("bounds_fixed")
            else "SNANA/bin/SALT2mu.exe"
        )
        r = command(exe, [path], folder, "bbc")
        r["executable_sha256"] = sha(exe)
        if label.startswith("bounds_fixed"):
            r["patch_sha256"] = sha(HERE / "bbc_bounds.patch")
        if "snr50" in label:
            r["changed_gate"] = (
                "Separate sparse-support diagnostic: intrinsic-scatter calibration SNR>50 instead of released60."
            )
        log = (folder / "bbc.log").read_text()
        r.update(
            input_sha256=sha(path),
            reference_source_sha256=sha(source),
            opt_biascor=opt,
            graceful=r["returncode"] == 0
            and "Done." in log[-1000:]
            and "FATAL ERROR ABORT" not in log,
            tail=log[-3500:],
        )
        outcomes[label] = r
    (RESULTS / "native-bbc.json").write_text(json.dumps(outcomes, indent=2) + "\n")
    print(json.dumps(outcomes, indent=2))


def multistart():
    f, _ = load("SPH_EVAL_NOMINAL")
    f = selected(f, "primary").sort_values("zHD")
    subset = f.iloc[__import__("numpy").linspace(0, len(f) - 1, 24).astype(int)]
    base = (WORK / "fits/SPH_EVAL_NOMINAL/fit.nml").read_text()

    def one(sign):
        folder = WORK / "multistart" / ("plus" if sign > 0 else "minus")
        folder.mkdir(parents=True, exist_ok=True)
        initial = folder / "initial.FITRES"
        initial.write_text(
            "VARNAMES: CID PKMJD x0 x1 c\n"
            + "".join(
                f"SN: {i} {row.PKMJD+sign*4:.17g} {row.x0:.17g} {sign*2} {row.c:.17g}\n"
                for i, row in subset.iterrows()
            )
        )
        text = base.replace(
            "&SNLCINP",
            "&SNLCINP\n OPT_SNCID_LIST = 3\n SNCID_LIST_FILE = '" + str(initial) + "'",
        )
        path = folder / "fit.nml"
        path.write_text(text)
        r = command(WORK / "SNANA/bin/snlc_fit.exe", [path], folder, "fit")
        g = fitres(folder / "fit.FITRES.TEXT")
        c = subset.join(g, rsuffix="_alternative", validate="one_to_one")
        r["returned"] = len(g)
        r["requested"] = len(subset)
        r["max_absolute_differences"] = {
            x: float((c[x] - c[x + "_alternative"]).abs().max())
            for x in ["mB", "x1", "c", "PKMJD", "FITCHI2"]
        }
        r["chi2_worse_by_gt01"] = int((c.FITCHI2_alternative - c.FITCHI2 > 0.1).sum())
        r["chi2_better_by_gt01"] = int((c.FITCHI2 - c.FITCHI2_alternative > 0.1).sum())
        c.to_csv(folder / "comparison.csv")
        r["comparison_sha256"] = sha(folder / "comparison.csv")
        r["input_sha256"] = sha(path)
        r["initialization_sha256"] = sha(initial)
        return sign, r

    with ThreadPoolExecutor(max_workers=2) as pool:
        r = dict(pool.map(one, [-1, 1]))
    (RESULTS / "native-multistart.json").write_text(json.dumps(r, indent=2) + "\n")
    print(r)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("check", choices=["bbc", "legacy", "multistart"])
    a = p.parse_args()
    globals()[a.check]()
