"""Independent read-only audit of fixed CALWF3 flat-reference components."""
from __future__ import annotations

from collections import Counter
from contextlib import ExitStack
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
FLAT = HERE.parent
BASE = FLAT.parents[1]


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def check(a, b, label):
    if abs(a - b) > 1e-10 + 1e-10 * abs(b):
        raise AssertionError(f"{label}: {a} versus {b}")
    return abs(a - b)


def main():
    protocol = json.loads((HERE / "protocol.json").read_text())
    for r in protocol["inputs"]:
        p = Path(r["path"])
        assert sha(p) == r["sha256"] and p.stat().st_size == r["bytes"], p
    state = json.loads((BASE / "pixel_execution/stage_a.json").read_text())
    root = json.loads((FLAT / "result.json").read_text())
    with np.load(BASE / "pixel_execution/operators.npz", allow_pickle=False) as f:
        fields = {k: f[k] for k in ("flat", "aperture", "annulus", "pam")}
    offsets = {(r["candidate"], r["exposure"]): (r["start"],r["end"])
               for r in state["operator_rows"]}
    objects = [r for r in state["rows"] if r.get("secondary")]
    assert len(objects) == 170
    root_exp = {(r["candidate"],r["filename"]):r for r in root["exposure_rows"]}
    root_pairs = {(r["candidate"],r["visit"]):r for r in root["pair_rows"]}
    assert len(root_exp) == 1360 and len(root_pairs) == 340
    saved = {}
    maximum = Counter()
    support = []
    for exposure_index, e in enumerate(state["exposures"]):
        image = Path(e["path"])
        with ExitStack() as stack:
            science = stack.enter_context(fits.open(image, memmap=True, do_not_scale_image_data=True))
            assert science["SCI"].shape == (1014,1014)
            header = science[0].header
            assert header["LFLTFILE"] == "N/A"
            refs = {}
            for kind in ("PFLTFILE", "DFLTFILE"):
                p = FLAT.parent / header[kind].split("$")[-1]
                ref = stack.enter_context(fits.open(p, memmap=True))
                assert ref["SCI"].shape == (1024,1024)
                sh,rh = science["SCI"].header,ref["SCI"].header
                assert (sh["LTM1_1"],sh["LTM2_2"],rh["LTM1_1"],rh["LTM2_2"]) == (1,1,1,1)
                assert (rh["LTV1"]-sh["LTV1"],rh["LTV2"]-sh["LTV2"]) == (5,5)
                refs[kind] = ref
            occurrences = Counter()
            unique_indices = set()
            time_min,time_max = float("inf"),float("-inf")
            for obj in objects:
                cid = obj["id"]
                lo,hi = offsets[(cid,exposure_index)]
                pix = fields["flat"][lo:hi]
                assert len(np.unique(pix)) == len(pix)
                sci = science["SCI"].data.ravel()[pix].astype(np.float64)
                err = science["ERR"].data.ravel()[pix].astype(np.float64)
                dq = science["DQ"].data.ravel()[pix]
                valid = (dq == 0) & np.isfinite(sci) & np.isfinite(err) & (err > 0)
                ap = fields["aperture"][lo:hi].astype(np.float64)
                ann = fields["annulus"][lo:hi].astype(np.float64)
                pam = fields["pam"][lo:hi].astype(np.float64)
                ka = e["photflam"] / state["exposures"][0]["photflam"]
                scale = np.dot(valid, ap*pam) / np.dot(valid, ann)
                weight = ka * valid * (ap*pam - scale*ann)
                quoted = float(np.dot(weight*err, weight*err))
                yref,xref = pix // 1014 + 5, pix % 1014 + 5
                assert xref.min() >= 5 and yref.min() >= 5
                assert xref.max() <= 1018 and yref.max() <= 1018
                tau = {}
                component = {}
                for kind,label in (("PFLTFILE","pflat_V"),("DFLTFILE","dflat_V")):
                    rr = refs[kind]
                    ref_sci = rr["SCI"].data[yref,xref].astype(np.float64)
                    ref_err = rr["ERR"].data[yref,xref].astype(np.float64)
                    assert np.all(np.isfinite(ref_sci[valid])) and np.all(ref_sci[valid] > 0)
                    assert np.all(np.isfinite(ref_err[valid])) and np.all(ref_err[valid] >= 0)
                    tau[kind] = (ref_err[valid]/ref_sci[valid])**2
                    component[label] = float(np.dot((weight[valid]*sci[valid])**2,
                                                    tau[kind]))
                flatvar = component["pflat_V"] + component["dflat_V"]
                expected = root_exp[(cid,e["filename"])]
                for name,value in (("quoted_V",quoted),("flat_V",flatvar),
                                   ("pflat_V",component["pflat_V"]),
                                   ("dflat_V",component["dflat_V"]),
                                   ("flat_fraction",flatvar/quoted)):
                    maximum[name] = max(maximum[name], check(value,expected[name],
                                                               f"{e['filename']} CID{cid} {name}"))
                saved[(cid,e["visit"],e["ordinal"])] = {
                    "pixels":pix[valid].astype(np.int32),
                    "coefficient":(weight[valid]*sci[valid]),
                    "tau2":tau["PFLTFILE"]+tau["DFLTFILE"],
                    "quoted":quoted,"flat":flatvar}
                samp=science["SAMP"].data.ravel()[pix][valid]
                sec=science["TIME"].data.ravel()[pix][valid]
                occurrences.update(int(x) for x in samp)
                unique_indices.update(int(x) for x in pix[valid])
                time_min=min(time_min,float(np.min(sec)))
                time_max=max(time_max,float(np.max(sec)))
            old_support = json.loads((FLAT / (e["filename"]+"-sample-support.json")).read_text())
            assert dict(occurrences) == {int(k):v for k,v in old_support["sample_counts"].items()}
            assert (time_min,time_max) == tuple(old_support["time_min_max"])
            unique=np.fromiter(sorted(unique_indices),dtype=np.int64)
            unique_samp=Counter(int(x) for x in science["SAMP"].data.ravel()[unique])
            unique_time=science["TIME"].data.ravel()[unique]
            support.append({"filename":e["filename"],"visit":e["visit"],
                            "operator_occurrences":sum(occurrences.values()),
                            "occurrence_SAMP":dict(occurrences),
                            "unique_valid_detector_pixels":len(unique),
                            "unique_SAMP":dict(unique_samp),
                            "unique_TIME_min_max":[float(unique_time.min()),float(unique_time.max())],
                            "NSAMP_header":header["NSAMP"],"EXPTIME_header":header["EXPTIME"]})
    pair_sums = {v:Counter() for v in ("search","template")}
    for obj in objects:
        for visit in ("search","template"):
            first,second=[saved[(obj["id"],visit,j)] for j in (2,4)]
            # Independent map-join instead of the original intersect1d route.
            index = {int(p):j for j,p in enumerate(second["pixels"])}
            cross=0.0
            for j,p in enumerate(first["pixels"]):
                k=index.get(int(p))
                if k is not None:
                    check(first["tau2"][j],second["tau2"][k],"shared reference same pixel")
                    cross += first["coefficient"][j]*second["coefficient"][k]*first["tau2"][j]
            v=first["quoted"]+second["quoted"]
            fv=first["flat"]+second["flat"]
            expected=root_pairs[(obj["id"],visit)]
            for name,value in (("quoted_V",v),("flat_marginal_V",fv),
                               ("same_detector_pixel_flat_covariance",cross),
                               ("V_minus_twice_flat_covariance",v-2*cross),
                               ("flat_fraction",fv/v)):
                maximum[name]=max(maximum[name],check(value,expected[name],
                                                        f"{visit} CID{obj['id']} {name}"))
            pair_sums[visit]["quoted"]+=v
            pair_sums[visit]["flat"]+=fv
            pair_sums[visit]["twice_cov"]+=2*cross
    for visit in ("search","template"):
        rs=root["summary"][visit]
        for name,value in (("sum_quoted_V",pair_sums[visit]["quoted"]),
                           ("sum_flat_marginal_V",pair_sums[visit]["flat"]),
                           ("sum_twice_flat_covariance",pair_sums[visit]["twice_cov"]),
                           ("pooled_flat_fraction",pair_sums[visit]["flat"]/pair_sums[visit]["quoted"]),
                           ("pooled_variance_fraction_removed_by_shared_flat",
                            pair_sums[visit]["twice_cov"]/pair_sums[visit]["quoted"])):
            maximum["summary_"+name]=max(maximum["summary_"+name],
                                          check(value,rs[name],visit+" summary "+name))
    result={"status":"PASS independent detector/reference mapping and 1360 exposure/340 pair flat arithmetics",
            "no_science_flux_or_ERR_modified":True,
            "protocol_sha256":sha(HERE/"protocol.json"),
            "root_result_sha256":sha(FLAT/"result.json"),
            "maximum_abs_difference":dict(maximum),
            "support":support,
            "interpretation":"Only same-reference, same-detector-pixel covariance; crosspixel/crossreference correlations and true ramp noise unresolved"}
    (HERE/"result.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],
                      "max_gap":max(maximum.values()),
                      "total_unique_valid_by_exposure":[r["unique_valid_detector_pixels"] for r in support]},indent=2))


if __name__ == "__main__":
    main()
