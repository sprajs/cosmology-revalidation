"""Quoted flat-reference variance on the unchanged HST aperture support."""

from pathlib import Path
from collections import Counter
from contextlib import ExitStack
import json
import numpy as np
from astropy.io import fits
from lib.paths import DATA
from lib.records import sha256 as sha

P = DATA / "hst"
DEFAULTS = {}


def run(out, config):
    O = out
    state = json.loads((P / "stage_a.json").read_text())
    assert sha(P / "operators.npz") == state["operators_sha256"]
    with np.load(P / "operators.npz") as z:
        ops = {k: z[k] for k in ("flat", "aperture", "annulus", "pam")}
    offsets = {
        (r["candidate"], r["exposure"]): (r["start"], r["end"])
        for r in state["operator_rows"]
    }
    candidates = [r for r in state["rows"] if r.get("secondary")]
    assert len(candidates) == 170
    results = []
    contrasts = {}
    inputs = [Path(__file__), P / "stage_a.json", P / "operators.npz"]
    for i, e in enumerate(state["exposures"]):
        path = P / "exposures" / e["filename"]
        assert sha(path) == e["sha256"]
        inputs.append(path)
        samp = Counter()
        time_values = []
        with ExitStack() as stack:
            f = stack.enter_context(
                fits.open(path, memmap=True, do_not_scale_image_data=True)
            )
            assert f[0].header["LFLTFILE"] == "N/A"
            refs = []
            for key in ("PFLTFILE", "DFLTFILE"):
                rp = P / "flat_refs" / f[0].header[key].split("$")[-1]
                rf = stack.enter_context(fits.open(rp, memmap=True))
                inputs.append(rp)
                assert rf["SCI"].shape == (1024, 1024)
                assert all(
                    rf["SCI"].header[k] == 1 and f["SCI"].header[k] == 1
                    for k in ("LTM1_1", "LTM2_2")
                )
                # Logical coordinates = LTM*physical + LTV. Both LTM=1;
                # reference pixel coordinate is science coordinate +5.
                dx = rf["SCI"].header["LTV1"] - f["SCI"].header["LTV1"]
                dy = rf["SCI"].header["LTV2"] - f["SCI"].header["LTV2"]
                assert dx == dy == 5
                refs.append(rf)
            for c in candidates:
                lo, hi = offsets[c["id"], i]
                idx = ops["flat"][lo:hi]
                y, err = [f[k].data.ravel()[idx].astype(float) for k in ("SCI", "ERR")]
                dq = f["DQ"].data.ravel()[idx]
                valid = np.isfinite(y) & np.isfinite(err) & (err > 0) & (dq == 0)
                a, b, p = [
                    ops[k][lo:hi].astype(float) for k in ("aperture", "annulus", "pam")
                ]
                k = e["photflam"] / state["exposures"][0]["photflam"]
                w = k * (
                    a * p * valid
                    - np.sum(a * p * valid) / np.sum(b * valid) * b * valid
                )
                iy, ix = idx // 1014 + 5, idx % 1014 + 5
                relative = []
                for rf in refs:
                    q = rf["SCI"].data[iy, ix].astype(float)
                    qe = rf["ERR"].data[iy, ix].astype(float)
                    assert np.isfinite(q[valid]).all() and (q[valid] > 0).all()
                    assert np.isfinite(qe[valid]).all() and (qe[valid] >= 0).all()
                    rr = np.zeros(len(idx))
                    rr[valid] = (qe[valid] / q[valid]) ** 2
                    relative.append(rr)
                term = w[valid] * y[valid]
                rv = (relative[0] + relative[1])[valid]
                flatvar = float(np.sum(term * term * rv))
                v = float(np.sum((w[valid] * err[valid]) ** 2))
                assert flatvar <= v * (1 + 1e-5)
                one = dict(
                    candidate=c["id"],
                    visit=e["visit"],
                    ordinal=e["ordinal"],
                    filename=e["filename"],
                    quoted_V=v,
                    flat_V=flatvar,
                    pflat_V=float(np.sum(term * term * relative[0][valid])),
                    dflat_V=float(np.sum(term * term * relative[1][valid])),
                    flat_fraction=flatvar / v,
                )
                results.append(one)
                contrasts[c["id"], e["visit"], e["ordinal"]] = {
                    **one,
                    "idx": idx[valid],
                    "coefficient": term,
                    "relative_variance": rv,
                }
                sv = f["SAMP"].data.ravel()[idx][valid]
                tv = f["TIME"].data.ravel()[idx][valid]
                samp.update(map(int, sv))
                time_values.extend(map(float, tv))
            (O / (e["filename"] + "-sample-support.json")).write_text(
                json.dumps(
                    {
                        "exposure": e["filename"],
                        "sampled_operator_occurrences": sum(samp.values()),
                        "sample_counts": dict(samp),
                        "time_min_max": [min(time_values), max(time_values)],
                        "header": {
                            key: f[0].header.get(key)
                            for key in (
                                "EXPTIME",
                                "NSAMP",
                                "SAMP_SEQ",
                                "CAL_VER",
                                "DARKFILE",
                                "PFLTFILE",
                                "DFLTFILE",
                                "READNSEA",
                                "READNSEB",
                                "READNSEC",
                                "READNSED",
                            )
                        },
                    },
                    indent=2,
                )
                + "\n"
            )
    pairs = []
    for c in candidates:
        for visit in ("search", "template"):
            a, b = [contrasts[c["id"], visit, j] for j in (2, 4)]
            common, ia, ib = np.intersect1d(a["idx"], b["idx"], return_indices=True)
            assert np.array_equal(
                a["relative_variance"][ia], b["relative_variance"][ib]
            )
            covariance = float(
                np.sum(
                    a["coefficient"][ia]
                    * b["coefficient"][ib]
                    * a["relative_variance"][ia]
                )
            )
            v = a["quoted_V"] + b["quoted_V"]
            fv = a["flat_V"] + b["flat_V"]
            assert abs(covariance) <= np.sqrt(a["flat_V"] * b["flat_V"]) + 1e-12
            pairs.append(
                dict(
                    candidate=c["id"],
                    visit=visit,
                    quoted_V=v,
                    flat_marginal_V=fv,
                    same_detector_pixel_flat_covariance=covariance,
                    V_minus_twice_flat_covariance=v - 2 * covariance,
                    flat_fraction=fv / v,
                    common_valid_detector_pixels=len(common),
                )
            )
    summary = {}
    for visit in ("search", "template"):
        rr = [r for r in pairs if r["visit"] == visit]
        v = sum(r["quoted_V"] for r in rr)
        summary[visit] = {
            "pairs": len(rr),
            "sum_quoted_V": v,
            "sum_flat_marginal_V": sum(r["flat_marginal_V"] for r in rr),
            "pooled_flat_fraction": sum(r["flat_marginal_V"] for r in rr) / v,
            "sum_twice_flat_covariance": 2
            * sum(r["same_detector_pixel_flat_covariance"] for r in rr),
            "pooled_variance_fraction_removed_by_shared_flat": 2
            * sum(r["same_detector_pixel_flat_covariance"] for r in rr)
            / v,
        }
    out = {
        "summary": summary,
        "exposure_rows": results,
        "pair_rows": pairs,
        "interpretation": "Pipeline diagonal flat-error component and same-reference/same-detector-pixel covariance under independent reference pixels. Actual cross-pixel, cross-reference and ramp covariance are not determined; this is not an empirical noise correction.",
    }
    (O / "result.json").write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    (O / "manifest.json").write_text(
        json.dumps({str(q): sha(q) for q in sorted(set(inputs))}, indent=2) + "\n"
    )
    return {"visits": summary, "scope": out["interpretation"]}
