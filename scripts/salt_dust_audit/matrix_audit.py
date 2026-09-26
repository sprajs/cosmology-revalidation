#!/usr/bin/env python3
"""Audit released covariance geometry without inventing signed dust corrections.

Reproduce: OPENBLAS_NUM_THREADS=2 .venv/bin/python scripts/salt_dust_audit/matrix_audit.py
Inputs are read-only; all outputs live in runs/salt_dust_audit/matrix_audit.
"""
from pathlib import Path
import hashlib
import json
import os
import platform
import subprocess
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import scipy
from scipy.linalg import cho_factor, cho_solve, eigh
from scipy.special import expit

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "runs/salt_dust_audit/matrix_audit"
RELEASES = {
    "original": {"root": "sources/repos/des-science__DES-SN5YR@1.3/4_DISTANCES_COVMAT",
                 "hd": "DES-SN5YR_HD.csv", "meta": "DES-SN5YR_HD+MetaData.csv",
                 "alpha": .16087, "beta": 3.11780, "gamma": .03754},
    "Dovekie": {"root": "sources/repos/des-science__DES-SN5YR/4_DISTANCES_COVMAT",
                "hd": "DES-Dovekie_HD.csv", "meta": "DES-Dovekie_Metadata.csv",
                "alpha": .169, "beta": 3.14, "gamma": .033},
}
INPUTS = set()


def read_table(path):
    INPUTS.add(path)
    if path.read_text().startswith("CID,"):
        d = pd.read_csv(path, dtype={"CID": str})
    else:
        d = pd.read_csv(path, sep=r"\s+", comment="#", dtype={"CID": str})
        d = d.drop(columns="VARNAMES:")
    assert not d[["CID", "IDSURVEY"]].duplicated().any()
    return d


def inv_spd(c):
    ans = cho_solve(cho_factor(c, lower=True), np.eye(len(c)))
    return (ans + ans.T) / 2


def load_matrix(path, n):
    INPUTS.add(path)
    if path.suffix == ".npz":
        with np.load(path) as d:
            assert int(d["nsn"][0]) == n
            p = np.zeros((n, n), dtype=np.float64)
            assert len(d["cov"]) == n * (n + 1) // 2
            p[np.triu_indices(n)] = d["cov"]
        p += np.triu(p, 1).T
        c = inv_spd(p)
        assert np.max(np.abs(p @ c - np.eye(n))) < 1e-8
        return c
    a = np.loadtxt(path)
    assert int(a[0]) == n and a.size == 1 + n*n
    c = a[1:].reshape(n, n)
    assert np.max(np.abs(c-c.T)) < 1e-12
    return c


def mu(z, zh, om=.3, w0=-1., wa=0., nodes=96):
    g, w = np.polynomial.legendre.leggauss(nodes)
    x = z[:, None] * (g+1)/2
    de = (1+x)**(3*(1+w0+wa)) * np.exp(-3*wa*x/(1+x))
    integ = z/2 * np.sum(w / np.sqrt(om*(1+x)**3 + (1-om)*de), axis=1)
    return 5*np.log10(299792.458/70 * (1+zh)*integ) + 25


def design(hd, meta):
    z, zh = hd.zHD.to_numpy(), hd.zHEL.to_numpy()
    columns = {"intercept": np.ones(len(hd))}
    checks = {}
    for name, value in [("om", .3), ("w0", -1.), ("wa", 0.)]:
        def deriv(eps):
            return (mu(z, zh, **{name: value+eps}) - mu(z, zh, **{name: value-eps}))/(2*eps)
        columns["dmu_d"+name] = deriv(1e-4)
        checks[name+"_derivative_step_halving_max"] = float(np.max(abs(deriv(1e-4)-deriv(5e-5))))
    checks["quadrature_doubling_max_mag"] = float(np.max(abs(mu(z, zh)-mu(z, zh, nodes=192))))
    columns.update({"dmu_dalpha": meta.x1.to_numpy(), "dmu_dbeta": -meta.c.to_numpy(),
                    "dmu_dgamma": expit((meta.HOST_LOGMASS.to_numpy()-10)/.001)-.5,
                    "dmu_dbias_scale": -meta.biasCor_mu.to_numpy(),
                    "z_linear": z, "MWEBV_proxy": meta.MWEBV.to_numpy(),
                    "DES_indicator": (hd.IDSURVEY.to_numpy()==10).astype(float)})
    return columns, checks


def correlation_and_center(x, p):
    one = np.ones(len(x))
    x = x - np.outer(one, (one @ p @ x)/(one @ p @ one))
    gram = x.T @ p @ x
    norms = np.sqrt(np.diag(gram))
    return gram / np.outer(norms, norms), x, norms


def matrix_health(c):
    ev = eigh(c, eigvals_only=True, check_finite=False, driver="evr")
    return {"min_eigenvalue_mag2": float(ev[0]), "max_eigenvalue_mag2": float(ev[-1]),
            "condition_number": float(ev[-1]/ev[0]),
            "median_diagonal_sigma_mag": float(np.median(np.sqrt(np.diag(c))))}


def audit_release(label, cfg):
    print("Starting", label, flush=True)
    out = OUT / label
    out.mkdir(parents=True, exist_ok=True)
    base = ROOT / cfg["root"]
    hd = read_table(base/cfg["hd"])
    raw = read_table(base/cfg["meta"])
    keys = pd.MultiIndex.from_frame(hd[["CID", "IDSURVEY"]])
    same_order = keys.equals(pd.MultiIndex.from_frame(raw[["CID", "IDSURVEY"]]))
    meta = raw.set_index(["CID", "IDSURVEY"]).loc[keys].reset_index()
    n = len(hd)
    assert np.max(abs(hd.MU-meta.MU)) < 1e-8
    suffix = ".txt.gz" if label == "original" else ".npz"
    cstat = load_matrix(base/("STATONLY"+suffix), n)
    ctotal = load_matrix(base/("STAT+SYS"+suffix), n)
    if label == "original":
        assert np.max(abs(cstat)) == 0
        cstat += np.diag(hd.MUERR_FINAL.to_numpy()**2)
        ctotal += cstat
    # Statistical matrices are diagonal in both releases. Assert rather than assume.
    assert np.max(abs(cstat-np.diag(np.diag(cstat)))) < 1e-12
    ps = np.diag(1/np.diag(cstat))
    pt = inv_spd(ctotal)
    sig = np.sqrt(np.diag(cstat))
    columns, numeric = design(hd, meta)
    names = list(columns)
    x = np.column_stack(list(columns.values()))
    rows = hd[["CID", "IDSURVEY", "zHD", "zHEL", "MU"]].copy()
    rows.insert(0, "row_index_0based", np.arange(n))
    for col in ["c", "x1", "HOST_LOGMASS", "MWEBV", "biasCor_mu"]:
        rows[col] = meta[col].to_numpy()
    for name in names:
        rows[name] = columns[name]
    rows["stat_sigma_mag"] = sig
    rows["published_color_term_mag"] = -cfg["beta"]*meta.c.to_numpy()
    rows["published_stretch_term_mag"] = cfg["alpha"]*meta.x1.to_numpy()
    rows["published_host_term_from_rounded_mass_mag"] = cfg["gamma"]*columns["dmu_dgamma"]
    rows["published_bias_correction_term_mag"] = -meta.biasCor_mu.to_numpy()
    rows.to_csv(out/"rows_and_correction_jacobian.csv", index=False)

    metrics = {}
    for metric, p in [("stat", ps), ("total", pt)]:
        gram, xc, norms = correlation_and_center(x[:,1:], p)
        pd.DataFrame(gram, index=names[1:], columns=names[1:]).to_csv(out/f"profiled_design_cosines_{metric}.csv")
        ev = np.linalg.eigvalsh(gram)
        metrics[metric] = {"profiled_normalized_design_eigenvalues": ev.tolist(),
                           "condition_number": float(ev[-1]/ev[0]),
                           "variance_inflation_factors": dict(zip(names[1:], np.diag(np.linalg.inv(gram)).tolist()))}

    # Local GLS operators map a positive additive data distance shift into fitted
    # parameters. No external BAO/CMB prior and no dust-prior interpretation.
    operators, operator_info = {}, {}
    families = {"flat_LCDM": [0,1], "flat_wCDM": [0,1,2], "flat_w0waCDM": [0,1,2,3]}
    parameter_names = ["intercept_mag", "Omega_m", "w0", "wa"]
    op_rows = rows[["row_index_0based", "CID", "IDSURVEY"]].copy()
    projection_rows = []
    for metric, p in [("stat", ps), ("total", pt)]:
        for family, idx in families.items():
            key = metric+"_"+family
            h = x[:,idx]
            f = h.T @ p @ h
            a = np.linalg.solve(f, h.T @ p)
            assert np.max(abs(a @ h - np.eye(len(idx)))) < 1e-8
            operators[key] = a
            operator_info[key] = {"parameters": [parameter_names[i] for i in idx], "fisher_condition": float(np.linalg.cond(f)),
                                  "conditional_parameter_covariance": np.linalg.inv(f).tolist(),
                                  "operator_times_design_identity_max_error": float(np.max(abs(a@h-np.eye(len(idx)))))}
            for j,i in enumerate(idx):
                op_rows[key+"__"+parameter_names[i]] = a[j]
                for k,name in enumerate(names[4:], start=4):
                    projection_rows.append({"model_and_metric": key, "parameter": parameter_names[i], "perturbation": name,
                                            "response_per_unit_perturbation": float(a[j] @ x[:,k])})
    op_rows.to_csv(out/"cosmology_response_operator.csv", index=False)
    pd.DataFrame(projection_rows).to_csv(out/"cosmology_response_to_correction_columns.csv", index=False)

    summary = {"n": n, "metadata_original_order_matches_HD": same_order,
               "stat_matrix_vs_published_sigma_max_abs_mag": float(np.max(abs(sig-hd["MUERR_FINAL" if label=="original" else "MUERR"].to_numpy()))),
               "stat": matrix_health(cstat), "total": matrix_health(ctotal), "numeric_checks": numeric,
               "design_collinearity": metrics, "cosmology_response": operator_info, "components": []}
    eigen_rows, mode_rows, response_cov, geometry_rows = [], [], [], []
    diag_rows = rows[["row_index_0based", "CID", "IDSURVEY"]].copy()
    factors = {}
    grouped_calibration = None
    standalone_calspec = None
    cs_sum = np.zeros((n,n))
    whitened_one = 1/sig
    q = whitened_one/np.linalg.norm(whitened_one)
    # Pstat whitened and common-intercept-profiled comparison designs.
    xp = x[:,1:]/sig[:,None]
    xp -= np.outer(q, q@xp)
    xp /= np.linalg.norm(xp, axis=0)

    files = sorted((base/"SingleSYS_CovMatrix").glob("*"+suffix))
    for path in files:
        name = path.name.removesuffix(suffix)
        print(label, name, flush=True)
        cs = load_matrix(path, n)
        if label != "original":
            cs -= cstat
        cs = (cs+cs.T)/2
        cs_sum += cs
        if name in {"CALIBplusSALT3", "CAL_SALT3"}:
            grouped_calibration = cs.copy()
        if name == "CALSPEC":
            standalone_calspec = cs.copy()
        # Raw eigenvalues can be distorted by float32 inversion roundoff for the
        # BEAMS-suppressed objects whose effective errors reach hundreds of mag.
        # Retain the raw spectrum, but compress in the statistical metric.
        ev_raw = eigh(cs, eigvals_only=True, check_finite=False, driver="evr")
        ks = cs/sig[:,None]/sig[None,:]
        ev, vectors = eigh(ks, check_finite=False, driver="evr")
        # This is an explicit numerical compression threshold, not a significance cut.
        tol = max(1e-12, ev[-1]*1e-5)
        keep = np.flatnonzero(abs(ev) > tol)[np.argsort(abs(ev[abs(ev)>tol]))[::-1]]
        signs = np.sign(ev[keep])
        v = vectors[:,keep]
        # Fixed arbitrary orientation, stable against eigensolver sign choices.
        v *= np.where(v[np.argmax(abs(v), axis=0), np.arange(len(keep))] >= 0, 1., -1.)
        bw = v*np.sqrt(abs(ev[keep]))[None,:]
        b = sig[:,None]*bw
        residual = np.linalg.norm(ks-(bw*signs)@bw.T)/np.linalg.norm(ks)
        bp = bw.copy()
        bp -= np.outer(q,q@bp)
        factors[name] = (bp, signs)
        component = {"name": name, "min_eigenvalue_mag2": float(ev_raw[0]),
                     "max_eigenvalue_mag2": float(ev_raw[-1]), "trace_mag2": float(np.trace(cs)),
                     "median_diagonal_sigma_mag": float(np.median(np.sqrt(np.maximum(np.diag(cs),0)))),
                     "max_diagonal_sigma_mag": float(np.sqrt(max(np.diag(cs)))),
                     "compression_threshold_stat_whitened": float(tol), "retained_rank": int(len(keep)),
                     "retained_positive_rank":int(np.sum(signs>0)), "retained_negative_rank":int(np.sum(signs<0)),
                     "compression_relative_stat_whitened_frobenius_residual": float(residual),
                     "negative_eigenvalue_absolute_sum_mag2": float(-ev_raw[ev_raw<0].sum()),
                     "stat_whitened_min_eigenvalue":float(ev[0]),"stat_whitened_max_eigenvalue":float(ev[-1]),
                     "stat_whitened_negative_frobenius_fraction":float(np.linalg.norm(ev[ev<0])/np.linalg.norm(ev)),
                     "stat_whitened_leading_positive_trace_fraction": float(ev[-1]/ev[ev>0].sum()),
                     "standalone_PSD_at_compression_tolerance": bool(np.all(signs>0))}
        summary["components"].append(component)
        diag_rows[name+"__variance_mag2"] = np.diag(cs)
        eigen_rows.extend({"component": name,"index_ascending_0based": i,"raw_eigenvalue_mag2": float(e),
                           "stat_whitened_eigenvalue_dimensionless":float(ev[i])} for i,e in enumerate(ev_raw))
        for k,eidx in enumerate(keep):
            df = rows[["row_index_0based", "CID", "IDSURVEY"]].copy()
            df["component"] = name; df["mode_0based"] = k
            df["eigenvalue_sign"] = int(signs[k])
            df["loading_mag_arbitrary_sign"] = b[:,k]
            mode_rows.append(df)
            bn = bp[:,k]/np.linalg.norm(bp[:,k])
            for j,coln in enumerate(names[1:]):
                geometry_rows.append({"component": name,"mode_0based":k,"eigenvalue_sign":int(signs[k]),"design_column":coln,
                                      "absolute_stat_whitened_profiled_cosine":float(abs(bn@xp[:,j]))})
        for op_name, a in operators.items():
            cv = a@cs@a.T
            # This carries released covariance through a fixed local response,
            # not the refitted cosmological likelihood of the alternative pipeline.
            response_cov.append({"component":name,"model_and_metric":op_name,
                                 "parameters":operator_info[op_name]["parameters"],
                                 "propagated_parameter_covariance":cv.tolist(),
                                 "propagated_parameter_sigma":np.sqrt(np.maximum(np.diag(cv),0)).tolist()})

    csy = ctotal-cstat
    diff = cs_sum-csy
    summary["single_component_additivity"] = {
        "interpretation":"Naive sum includes overlapping calibration group and separate CALSPEC; use grouping_adjusted_additivity below.",
        "sum_single_minus_total_sys_relative_frobenius": float(np.linalg.norm(diff)/np.linalg.norm(csy)),
        "sum_single_minus_total_sys_relative_stat_whitened_frobenius":float(np.linalg.norm(diff/sig[:,None]/sig[None,:])/np.linalg.norm(csy/sig[:,None]/sig[None,:])),
        "sum_single_minus_total_sys_max_abs_mag2":float(np.max(abs(diff))),
        "sum_single_trace_mag2":float(np.trace(cs_sum)),"released_total_sys_trace_mag2":float(np.trace(csy))}
    # The +cal COVOPT substring includes CALSPEC. Confirm its rank contribution
    # and compare disjoint groups, rather than treating all filenames as disjoint.
    remainder = csy-(cs_sum-standalone_calspec)
    rw = remainder/sig[:,None]/sig[None,:]
    ew,vw = eigh(rw,check_finite=False)
    cal_ev = eigh((grouped_calibration-standalone_calspec)/sig[:,None]/sig[None,:],eigvals_only=True)
    total_ev = eigh(csy/sig[:,None]/sig[None,:],eigvals_only=True)
    remain_tol = max(1e-12,total_ev[-1]*1e-5)
    remain_keep = np.flatnonzero(abs(ew)>remain_tol)
    summary["grouping_adjusted_additivity"] = {
        "operation":"Csys_total minus (sum of all individual files minus separately listed CALSPEC)",
        "relative_stat_whitened_frobenius":float(np.linalg.norm(rw)/np.linalg.norm(csy/sig[:,None]/sig[None,:])),
        "calibration_minus_CALSPEC_positive_rank":int(np.sum(cal_ev>cal_ev[-1]*1e-5)),
        "calibration_minus_CALSPEC_negative_rank":int(np.sum(cal_ev< -cal_ev[-1]*1e-5)),
        "remainder_min_stat_whitened_eigenvalue":float(ew[0]),
        "remainder_max_stat_whitened_eigenvalue":float(ew[-1]),
        "remainder_resolved_positive_rank_at_1e_minus5_total_max":int(np.sum(ew>remain_tol)),
        "remainder_resolved_negative_rank_at_1e_minus5_total_max":int(np.sum(ew< -remain_tol)),
        "remainder_trace_mag2":float(np.trace(remainder)),
        "remainder_identity":"Unassigned. Original README lists SIGINT_MODEL but file absent; numerical resemblance is not proof of identity." if label=="original" else "Consistent with finite-precision residue after removing CALSPEC overlap."}
    if len(remain_keep):
        rem_rows=[]
        for k,i in enumerate(remain_keep):
            vec=vw[:,i]*np.sqrt(abs(ew[i]))*sig
            vec*=1 if vec[np.argmax(abs(vec))]>=0 else -1
            rd=rows[["row_index_0based","CID","IDSURVEY"]].copy()
            rd["mode_0based"]=k;rd["eigenvalue_sign"]=int(np.sign(ew[i]))
            rd["loading_mag_arbitrary_sign"]=vec
            rem_rows.append(rd)
        pd.concat(rem_rows,ignore_index=True).to_csv(out/"unassigned_grouping_remainder_modes.csv",index=False)
    gram_rows = []
    for a, (ba,sa) in factors.items():
        for b, (bb,sb) in factors.items():
            ab = float(np.sum((ba.T@bb)**2*np.outer(sa,sb)))
            aa = float(np.sum((ba.T@ba)**2*np.outer(sa,sa)))
            bbn = float(np.sum((bb.T@bb)**2*np.outer(sb,sb)))
            gram_rows.append({"component_a":a,"component_b":b,
                              "stat_whitened_profiled_covariance_frobenius_cosine":ab/np.sqrt(aa*bbn)})
    pd.DataFrame(eigen_rows).to_csv(out/"systematic_eigenspectra.csv.gz",index=False)
    pd.concat(mode_rows,ignore_index=True).to_csv(out/"systematic_spectral_modes.csv.gz",index=False)
    pd.DataFrame(geometry_rows).to_csv(out/"mode_design_collinearity.csv",index=False)
    pd.DataFrame(gram_rows).to_csv(out/"component_covariance_overlap.csv",index=False)
    diag_rows.to_csv(out/"systematic_diagonal_variances.csv",index=False)
    (out/"propagated_parameter_covariances.json").write_text(json.dumps(response_cov,indent=2)+"\n")
    (out/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    return summary


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    summary = {label:audit_release(label,cfg) for label,cfg in RELEASES.items()}
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    primary_source_paths = [
        "sources/repos/RickKessler__SNANA/util/create_covariance.py",
        "papers/text/2511.07517v3.txt",
        "papers/text/2401.02945v2.txt",
    ]
    for cfg in RELEASES.values():
        primary_source_paths += [cfg["root"]+"/README.md",cfg["root"]+"/SingleSYS_CovMatrix/README.md"]
        primary_source_paths.append(str(Path(cfg["root"]).parent/"7_PIPPIN_FILES/D5yr_analysis.yml"))
    INPUTS.update(ROOT/p for p in primary_source_paths)
    manifest = {"created_utc":datetime.now(timezone.utc).isoformat(),
                "git_head":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
                "script_sha256":sha(Path(__file__)),
                "inputs":{str(p.relative_to(ROOT)):sha(p) for p in sorted(INPUTS)},
                "outputs":{str(p.relative_to(ROOT)):sha(p) for p in sorted(OUT.rglob("*")) if p.is_file() and p.name not in {"manifest.json", "verification.json"}},
                "environment":{"python":platform.python_version(),"numpy":np.__version__,"pandas":pd.__version__,
                               "scipy":scipy.__version__,"OPENBLAS_NUM_THREADS":os.environ.get("OPENBLAS_NUM_THREADS")},
                "configuration":{"cosmology_linearization":{"Omega_m":.3,"w0":-1.,"wa":0.,"flat":True},
                                 "common_intercept":"profiled in collinearity; fitted in response operators",
                                 "eigenmode_compression":"statistically whiten; retain absolute eigenvalues > max(1e-12,1e-5 * largest positive eigenvalue), including negative eigenvalue sign",
                                 "mode_sign":"largest absolute loading positive; not physical perturbation direction",
                                 "statistical_metric":"released STATONLY covariance; avoid inclusion of systematic under test in geometry"}}
    (OUT/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print("Complete:",OUT,flush=True)


if __name__ == "__main__":
    main()
