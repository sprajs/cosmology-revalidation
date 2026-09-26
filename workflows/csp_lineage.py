"""Rebuild CSP physical-filter lineage from original archive table literals."""

import tarfile
from decimal import Decimal
from collections import Counter, defaultdict
from lib.paths import DATA
from lib.photometry import read_raisin
from lib.records import write_rows

DEFAULTS = {}


def decimal(value):
    return str(Decimal(value).normalize())


def run(out, cfg):
    raw = []
    snpy = []
    nir = {"Y", "Ydw", "J", "Jrc2", "Jdw", "H", "Hdw"}
    merge = {"Jdw": "J", "Hdw": "H"}
    with tarfile.open(DATA / "csp/photometry.tgz") as archive:
        handle = archive.extractfile("DR3/SN_photo.dat")
        if handle is None:
            raise ValueError("Missing raw CSP table")
        for line in handle.read().decode().splitlines():
            parts = line.split()
            if len(parts) != 5:
                raise ValueError("Invalid raw photometry row")
            if parts[1] in nir:
                raw.append((parts[0].lower(), parts[1], *map(decimal, parts[2:])))
        for member in archive:
            if not member.isfile() or not member.name.endswith("_snpy.txt"):
                continue
            lines = archive.extractfile(member).read().decode().splitlines()
            name = lines[0].split()[0].lower()
            band = None
            for line in lines[1:]:
                a = line.split()
                if not a or a[0].startswith("#"):
                    continue
                if a[0] == "filter":
                    band = a[1]
                    continue
                if len(a) != 3:
                    raise ValueError("Invalid SNooPy row")
                if band in nir:
                    snpy.append((name, band, *map(decimal, a)))
    expected = Counter((name, merge.get(b, b), t, m, e) for name, b, t, m, e in raw)
    actual = Counter(snpy)
    physical = defaultdict(set)
    by_object = defaultdict(list)
    for name, b, t, m, e in raw:
        physical[name, merge.get(b, b), t].add(b)
    for row in snpy:
        by_object[row[0]].append(row)
    mapping = {"J": "J", "Jrc2": "j", "Y": "Y", "Ydw": "y", "H": "H"}
    inverse = {v: k for k, v in mapping.items()}
    folder = DATA / "raisin/photometry/RAISIN/CSPDR3_RAISIN"
    objects = []
    ledger = []
    mismatches = []
    for filename in (folder / "CSPDR3_RAISIN.LIST").read_text().splitlines():
        if not filename.strip() or filename.startswith("#"):
            continue
        header, observations = read_raisin(folder / filename)
        cid = header["SNID"].split()[0]
        key = "sn" + cid.lower()
        released = []
        for r in observations:
            if r["FLT"] not in mapping.values():
                continue
            released.append(
                (
                    r["FLT"],
                    decimal(Decimal(r["MJD"]) - 53000),
                    decimal(r["MAG"]),
                    decimal(r["MAGERR"]),
                    decimal(r["FLUXCAL"]),
                    decimal(r["FLUXCALERR"]),
                )
            )
        predicted = []
        for _, b, t, m, e in by_object[key]:
            f = 10 ** (-0.4 * (float(m) - 27.5))
            err = f * (10 ** (0.4 * float(e)) - 1)
            predicted.append(
                (mapping[b], t, m, e, decimal(f"{f:.5e}"), decimal(f"{err:.5e}"))
            )
        pc, rc = Counter(predicted), Counter(released)
        objects.append(
            {
                "CID": cid,
                "predicted_rows": len(predicted),
                "released_rows": len(released),
                "exact_converter_identity": pc == rc,
            }
        )
        for r, count in rc.items():
            bands = physical.get((key, inverse[r[0]], r[1]), set())
            ledger.append(
                {
                    "CID": cid,
                    "released_band": r[0],
                    "time": r[1],
                    "multiplicity": count,
                    "physical_filter": "|".join(sorted(bands)),
                    "unambiguous": len(bands) == 1,
                    "converter_identity": pc[r] >= count,
                }
            )
        for kind, counts in [("predicted_only", pc - rc), ("released_only", rc - pc)]:
            for row, count in counts.items():
                mismatches.append(
                    {
                        "CID": cid,
                        "kind": kind,
                        "row": " ".join(row),
                        "multiplicity": count,
                    }
                )
    write_rows(out / "objects.csv", objects)
    write_rows(out / "physical_filter_ledger.csv", ledger)
    write_rows(
        out / "mismatches.csv", mismatches, ["CID", "kind", "row", "multiplicity"]
    )
    return {
        "raw_NIR_rows": len(raw),
        "snpy_NIR_rows": len(snpy),
        "raw_to_snpy_exact_after_label_merge": expected == actual,
        "raw_unmatched": sum((expected - actual).values()),
        "snpy_unmatched": sum((actual - expected).values()),
        "listed_objects": len(objects),
        "converter_exact_objects": sum(x["exact_converter_identity"] for x in objects),
        "raw_Jdw_rows": sum(r[1] == "Jdw" for r in raw),
        "released_unique_Jdw_rows": sum(
            r["multiplicity"] for r in ledger if r["physical_filter"] == "Jdw"
        ),
        "scope": "Exact printed-value and metadata lineage. Merged labels do not establish an empirical calibration correction; ambiguous, missing and duplicate rows are retained. The physical-filter labels are assigned by identity/time, never by matching brightness.",
    }
