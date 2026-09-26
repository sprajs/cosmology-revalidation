#!/usr/bin/env python3
"""Independent FITS-level HST checks; writes only validation/reports/hst.json.

Run from any directory with .venv/bin/python. No workflow or lib
module is imported. Frozen result tables are comparators, never input pixels.
The added deletion diagnostics describe influence; they do not change a sample
or provide a detector-noise correction or independent-pixel significance.
"""
from __future__ import annotations
import csv
import hashlib
import argparse
import json
import math
from pathlib import Path
import time
import numpy as np
from astropy.io import fits
from scipy.integrate import quad

BASE = Path(__file__).resolve().parents[1]
DATA, RESULTS = BASE / "data", BASE / "results/baseline"
CHECKS, INPUTS = [], {}


def track(path):
    path = Path(path)
    if str(path.relative_to(BASE)) not in INPUTS:
        INPUTS[str(path.relative_to(BASE))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return path


def load(path):
    return json.loads(track(path).read_text())


def check(name, actual, expected, rtol=1e-10, atol=1e-10):
    a, b = np.asarray(actual), np.asarray(expected)
    passed = a.shape == b.shape and bool(np.allclose(a, b, rtol=rtol, atol=atol))
    CHECKS.append({"name": name, "passed": passed, "rtol": rtol, "atol": atol,
                   "max_absolute_gap": float(np.max(np.abs(a-b))) if a.shape == b.shape else None})
    if not passed:
        raise AssertionError(name)


def csvrows(path):
    with track(path).open() as handle:
        return list(csv.DictReader(handle))


def influence(numerator, denominator, groups):
    """Fixed group deletions and concentration, not a confidence interval."""
    n, d, g = np.asarray(numerator), np.asarray(denominator), np.asarray(groups)
    rows = [{"group": int(k), "n": int(np.sum(g == k)),
             "ratio_with_group_deleted": float((n.sum()-n[g == k].sum())/(d.sum()-d[g == k].sum())),
             "fraction_squared_signal": float(n[g == k].sum()/n.sum())} for k in np.unique(g)]
    return {"full_ratio": float(n.sum()/d.sum()), "deletion_ratio_range":
            [min(r["ratio_with_group_deleted"] for r in rows), max(r["ratio_with_group_deleted"] for r in rows)],
            "largest_group_signal_fraction": max(r["fraction_squared_signal"] for r in rows), "groups": rows}


def science():
    p = DATA / "hst"
    state = load(p / "stage_a.json")
    with np.load(track(p / "operators.npz")) as z:
        ops = {k: z[k] for k in z.files}
    candidates = [r for r in state["rows"] if r.get("secondary")]
    check("geometry frozen support", [len(state["rows"]), len(candidates), sum(bool(r.get("strict")) for r in state["rows"])], [385, 170, 0], 0, 0)
    offsets = {(r["candidate"], r["exposure"]): slice(r["start"], r["end"]) for r in state["operator_rows"]}
    measurement, max_cancel, min_a, min_b = {}, 0., 1., 1.
    units, references = {}, {}
    min_distance = min(math.hypot(a["east"]-b["east"],a["north"]-b["north"]) for i,a in enumerate(candidates) for b in candidates[i+1:])
    for i,e in enumerate(state["exposures"]):
        path = track(p / "exposures" / e["filename"])
        with fits.open(path, memmap=False) as f:
            units[e["filename"]] = {key: f["SCI"].header.get(key) for key in ("BUNIT", "BSCALE", "BZERO")}
            references[e["filename"]] = {key: f[0].header.get(key) for key in ("CAL_VER", "NLINFILE", "PFLTFILE", "DFLTFILE")}
            refs = []
            for key in ("PFLTFILE", "DFLTFILE"):
                name = f[0].header[key].split("$")[-1]
                with fits.open(track(p / "flat_refs" / name), memmap=False) as ff:
                    refs.append((name, ff["SCI"].data.astype(float), ff["ERR"].data.astype(float)))
                    check("flat offset " + e["filename"] + key, [ff["SCI"].header[k]-f["SCI"].header[k] for k in ("LTV1","LTV2")], [5,5], 0,0)
            sci, err, dq = [f[k].data.ravel() for k in ("SCI","ERR","DQ")]
            for c in candidates:
                sl = offsets[c["id"],i]
                indices = ops["flat"][sl]
                a,b,pam = [ops[k][sl].astype(float) for k in ("aperture","annulus","pam")]
                y,s,quality = sci[indices].astype(float), err[indices].astype(float), dq[indices]
                good = np.isfinite(y) & np.isfinite(s) & (s>0) & (quality==0)
                min_a = min(min_a, float((a*pam)[good].sum()/(a*pam).sum()))
                min_b = min(min_b, float(b[good].sum()/b.sum()))
                aa,bb = (a*pam)[good], b[good]
                weights = aa - bb*(aa.sum()/bb.sum())
                weights *= e["photflam"]/state["exposures"][0]["photflam"]
                max_cancel = max(max_cancel, abs(float(weights.sum())))
                # A detector-index keyed coefficient map independently implements
                # shared-reference covariance across dithered exposures.
                ix = indices[good]
                ry,rx = ix//1014+5, ix%1014+5
                rv = sum((rs[ry,rx]/rf[ry,rx])**2 for _,rf,rs in refs)
                coeff = weights*y[good]
                measurement[c["id"],e["visit"],e["ordinal"]] = {
                    "flux": float(np.dot(weights,y[good])), "variance": float(np.dot(weights**2,s[good]**2)),
                    "flat_variance": float(np.dot(coeff**2,rv)),
                    "flat_map": {int(j): (float(co),float(variance)) for j,co,variance in zip(ix,coeff,rv)},
                    "refs": [r[0] for r in refs]}
    old = load(RESULTS/"hst-repeat/summary.json")
    oldflat = load(RESULTS/"hst-flat/summary.json")
    report = {}
    for visit in ("search","template"):
        d,v,fv,cv = [],[],[],[]
        for c in candidates:
            a,b = [measurement[c["id"],visit,k] for k in (2,4)]
            assert a["refs"] == b["refs"]
            d.append(a["flux"]-b["flux"]);v.append(a["variance"]+b["variance"])
            fv.append(a["flat_variance"]+b["flat_variance"])
            common = a["flat_map"].keys() & b["flat_map"].keys()
            cv.append(sum(a["flat_map"][j][0]*b["flat_map"][j][0]*a["flat_map"][j][1] for j in common))
        d,v = np.array(d),np.array(v)
        intercept = float(np.average(d,weights=1/v))
        ratio, centered = float(np.dot(d,d)/v.sum()), float(np.sum((d-intercept)**2/v)/(len(d)-1))
        check(visit+" repeat ratios", [ratio,centered], [old["visits"][visit]["sum_d2_over_sum_V"],old["visits"][visit]["centered_variance_ratio"]])
        flatfrac,shared = float(sum(fv)/v.sum()),float(2*sum(cv)/v.sum())
        check(visit+" flat variance", [flatfrac,shared], [oldflat["visits"][visit]["pooled_flat_fraction"],oldflat["visits"][visit]["pooled_variance_fraction_removed_by_shared_flat"]])
        report[visit] = {"n":len(d),"ratio":ratio,"centered_ratio":centered,"weighted_intercept":intercept,
            "flat_fraction":flatfrac,"shared_flat_fraction":shared,
            "tile_influence":influence(d*d,v,[c["tile"] for c in candidates]),
            "naive_flat_adjusted_ratio":float(np.dot(d,d)/(v.sum()-2*sum(cv))),
            "interpretation":"One repeat pair per visit, with many spatial locations. Tile deletion is descriptive, not repeated-exposure replication."}
    return {"visits":report,"units":units,"references":references,
            "maximum_constant_image_residual":max_cancel,"minimum_aperture_coverage":min_a,
            "minimum_annulus_coverage":min_b,"minimum_site_spacing_arcsec":min_distance,
            "geometry_scope":"Frozen sky operators independently applied; full sky-WCS rebuild is the separate hst-geometry replay. Support uses held-out finite/error/DQ quality, but not held-out brightness amplitudes or sign."}


def dark_support():
    p=DATA/"hst-dark"
    bad=np.zeros((1024,1024),bool)
    for name in ("3562029fi_bpx.fits","3562028ni_bpx.fits"):
        with fits.open(track(p/"references"/name)) as f:
            for r in f["BPIX"].data:
                if str(r["CCDAMP"]).strip() not in ("N/A","ABCD") or int(r["CCDCHIP"]) not in (-999,1) or float(r["CCDGAIN"]) not in (-999.,2.5) or not int(r["VALUE"]):continue
                x,y,n=int(r["PIX1"])+4,int(r["PIX2"])+4,int(r["LENGTH"])
                coords=np.arange(n)
                xx,yy=(x+coords,np.full(n,y)) if int(r["AXIS"])==1 else (np.full(n,x),y+coords)
                good=(xx>=0)&(xx<1024)&(yy>=0)&(yy<1024)
                bad[yy[good],xx[good]]=True
    active=~bad;active[:5]=False;active[1019:]=False;active[:,:5]=False;active[:,1019:]=False
    check("reference-only dark pixels",active.sum(),993750,0,0)
    plan=load(p/"design.json")
    radii=plan["spatial_operator"]["radii_pixels"]
    extent=math.ceil(max(radii)+.5)
    coords=range(-extent,extent+1)
    # Numerical integration uses vertical chords, independent of analytic
    # quadrant primitives in lib.dark_geometry.
    def area(x,y,r):
        lo,hi=max(x-.5,-r),min(x+.5,r)
        if hi<=lo:return 0.
        def length(t):
            chord=math.sqrt(max(0.,r*r-t*t))
            return max(0.,min(y+.5,chord)-max(y-.5,-chord))
        breaks=[v for yy in (y-.5,y+.5) if abs(yy)<r for v in (-math.sqrt(r*r-yy*yy),math.sqrt(r*r-yy*yy)) if lo<v<hi]
        return quad(length,lo,hi,points=breaks,epsabs=1e-11,epsrel=1e-11)[0]
    arrays=[np.array([[area(x,y,r) for x in coords] for y in coords]) for r in radii]
    check("independent circle quadrature areas",[a.sum() for a in arrays],np.pi*np.array(radii)**2,1e-11,1e-9)
    a=arrays[0];b=arrays[2]-arrays[1];sites=[]
    for y in range(32,1024,64):
        for x in range(32,1024,64):
            yy,xx=slice(y-extent,y+extent+1),slice(x-extent,x+extent+1)
            aa,bb=a*active[yy,xx],b*active[yy,xx]
            if aa.sum()/a.sum()>=.9 and bb.sum()/b.sum()>=.75:
                w=aa-bb*aa.sum()/bb.sum();check(f"dark weight cancellation {x},{y}",w.sum(),0,0,1e-10)
                sites.append((x,y,yy,xx,w))
    check("independent dark aperture support",len(sites),247,0,0)
    return active,sites,plan


def calibrated(active,sites):
    p=DATA/"hst-dark/calibrated";arr=[]
    for root in ("idbx41onq","idbx43p7q"):
        with fits.open(track(p/(root+"_flt.fits")),memmap=False) as f:
            assert f["SCI"].header["BUNIT"]=="COUNTS/S"
            arr.append({k:np.pad(f[k].data.astype(float),5) for k in ("SCI","ERR","DQ")})
    d=arr[1]["SCI"]-arr[0]["SCI"];v=arr[1]["ERR"]**2+arr[0]["ERR"]**2
    ratio=float(np.sum(d[active]**2)/np.sum(v[active]))
    ad,av,ag=[],[],[]
    for x,y,yy,xx,w in sites:
        ad.append(float(np.sum(w*d[yy,xx])));av.append(float(np.sum(w*w*v[yy,xx])));ag.append((y//256)*4+x//256)
    ad,av=np.array(ad),np.array(av)
    old=load(RESULTS/"hst-dark-calibrated/summary.json")
    check("calibrated pixel and independent-quadrature aperture ratios",[ratio,np.dot(ad,ad)/av.sum()],[old["regions"][0]["repeat_to_quoted_ratio"],old["apertures"]["ratio"]],1e-9,1e-10)
    flat=np.flatnonzero(active);rank=flat[np.argsort(d.ravel()[flat]**2)[::-1]];total=np.sum(d[active]**2)
    yn,xn=np.unravel_index(rank[:10],active.shape)
    check("largest pixel identity",[xn[0],yn[0]],[506,945],0,0)
    yy,xx=np.indices(active.shape);groups=((yy//64)*16+xx//64)[active]
    pixelinf=influence(d[active]**2,v[active],groups)
    bit32=active & (((arr[0]["DQ"].astype(np.uint16)|arr[1]["DQ"].astype(np.uint16))&32)!=0)
    # Counterfactual deletion recorded explicitly, never used as a replacement.
    peak=rank[0];flagless=active&~bit32
    return {"pixel_ratio":ratio,"aperture_ratio":float(np.dot(ad,ad)/av.sum()),"aperture_n":len(sites),
        "largest_pixel_fraction":float(d.ravel()[peak]**2/total),"top10_pixel_fraction":float(np.sum(d.ravel()[rank[:10]]**2)/total),
        "dq32_pixels":int(bit32.sum()),"dq32_squared_signal_fraction":float(np.sum(d[bit32]**2)/total),
        "counterfactual_without_dq32_ratio":float(np.sum(d[flagless]**2)/np.sum(v[flagless])),
        "pixel_block_influence":pixelinf,"aperture_16_region_influence":influence(ad*ad,av,ag),
        "interpretation":"Counterfactual deletion diagnoses concentration only. DQ32 and squared residuals are outcomes, so neither is used to redefine the frozen mask. Pixel and aperture results are different estimands."}


def raw(active,sites,plan):
    t=np.array(plan["time_operator"]["t_seconds_nominal"]);powers=plan["time_operator"]["powers"]
    X=np.column_stack((np.ones(7),t));hs=[]
    for power in powers:
        w=np.abs(np.linspace(-1,1,7))**power
        hs.append((np.linalg.pinv(X*np.sqrt(w[:,None]))*np.sqrt(w)[None,:])[1])
    h=np.array(hs);check("independent WLS slope constant annihilation",h.sum(axis=1),np.zeros(len(h)),0,1e-15)
    check("independent WLS slope normalization",h@t,np.ones(len(h)),0,1e-12)
    qdef={"B":(slice(5,512),slice(5,512)),"C":(slice(5,512),slice(512,1019)),"A":(slice(512,1019),slice(5,512)),"D":(slice(512,1019),slice(512,1019))}
    with np.load(track(RESULTS/"hst-dark-raw/read-matrix.npz")) as z:frozen={k:z[k] for k in z.files}
    old_ap={(r["pair"],int(r["x"]),int(r["y"]),float(r["power"])):float(r["signed_dn_per_nominal_s"]) for r in csvrows(RESULTS/"hst-dark-raw/apertures.csv")}
    rows=[];maximum_apgap=0.;maximum_qgap=0.;headers={}
    for ip,pair in enumerate(plan["pairs"]):
        readings=[]
        for root in (pair["earlier"],pair["later"]):
            with fits.open(track(DATA/"hst-dark/raw"/(root+"_raw.fits")),memmap=False,uint=True) as f:
                check(root+" read times",[f["SCI",k].header["SAMPTIME"] for k in range(15,8,-1)],t,0,1e-6)
                headers[root]=f[0].header["EXPFLAG"]
                readings.append(np.stack([f["SCI",k].data.astype(np.float64) for k in range(15,8,-1)]))
        d=readings[1]-readings[0];del readings
        for iq,(q,(yy,xx)) in enumerate(qdef.items()):
            v=d[:,yy,xx][:,active[yy,xx]];mean=v.mean(axis=1);centered=v-mean[:,None]
            # Independent centered-covariance decomposition, rather than the
            # original uncentered Gram construction.
            covariance=np.cov(v,bias=True)/2
            gamma=covariance+np.outer(mean,mean)/2
            check(f"raw matrix {pair['pair']} {q}",gamma,frozen["quadrant_gamma"][ip,0,iq],1e-10,1e-8)
            check(f"raw means {pair['pair']} {q}",mean,frozen["quadrant_mean"][ip,0,iq],1e-10,1e-10)
            actual=np.mean((h@v)**2,axis=1)/2
            frommatrix=np.einsum("ij,jk,ik->i",h,gamma,h)
            maximum_qgap=max(maximum_qgap,float(np.max(np.abs(actual-frommatrix))))
            check(f"raw direct matrix identity {pair['pair']} {q}",actual,frommatrix,1e-10,1e-10)
            direct_center=float(h[0]@covariance@h[0]);diag_center=float(np.dot(h[0]**2,np.diag(covariance)))
            rows.append({"pair":pair["pair"],"kind":pair["kind"],"quadrant":q,"primary_half_squared_slope":float(actual[0]),
                "coherent_fraction":float((h[0]@mean)**2/(2*actual[0])),"centered_full_temporal_over_diagonal_moment":direct_center/diag_center})
        for x,y,yy,xx,w in sites:
            slopes=h@np.array([np.sum(one[yy,xx]*w) for one in d])
            expected=[old_ap[pair["pair"],x,y,float(power)] for power in powers]
            maximum_apgap=max(maximum_apgap,float(np.max(np.abs(slopes-expected))))
            check(f"raw independent aperture {pair['pair']} {x},{y}",slopes,expected,1e-8,1e-8)
    return {"pair_quadrants":rows,"max_direct_matrix_gap":maximum_qgap,"max_independent_aperture_gap":maximum_apgap,
        "primary_disjoint_pairs":7,"overlapping_sensitivity_pairs":1,"indeterminate_exposures":sum(v=="INDETERMINATE" for v in headers.values()),
        "temporal_comparison_scope":"Centered spatial moments of paired raw-read differences contain detector events and changes, not isolated electronic read noise. Full/diagonal comparison diagnoses mathematical importance of temporal cross terms, not an ERR correction."}


def main():
    global RESULTS
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", default="baseline")
    args=parser.parse_args()
    RESULTS=BASE/"results"/args.results
    start=time.monotonic()
    report={"schema":1,"audit":"Independent HST FITS and quadrature checks","science":science()}
    active,sites,plan=dark_support()
    report["dark_calibrated"]=calibrated(active,sites)
    report["dark_raw"]=raw(active,sites,plan)
    report["checks"]={"count":len(CHECKS),"passed":sum(c["passed"] for c in CHECKS),"details":CHECKS}
    report["input_sha256"]=INPUTS
    report["script_sha256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    report["elapsed_seconds"]=time.monotonic()-start
    out=BASE/"validation/reports/hst.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"report":str(out),"checks":report["checks"]["passed"],"seconds":report["elapsed_seconds"]}))

if __name__=="__main__":main()
