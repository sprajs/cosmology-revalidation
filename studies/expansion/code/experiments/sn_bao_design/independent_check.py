"""Outcome-free numerical review of the frozen SN/BAO feasibility design.

Only CID, zHD and zHEL are requested from the SN table. No magnitude is loaded
or scored. GLS coefficients use whitened QR rather than the primary normal
equations implementation.
"""
import os
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[name] = "1"

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad
from scipy.linalg import solve_triangular
from scipy.stats import norm
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
SN = ROOT / "sources/repos/CobayaSampler__sn_data/PantheonPlus"
BAO = ROOT / "sources/repos/CobayaSampler__bao_data/desi_bao_dr2"


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def curve(z, family, a, b, ng=128):
    nodes, weights = leggauss(ng)
    zz = z[:, None] * (nodes[None, :] + 1) / 2
    if family == "w":
        e = np.sqrt(a * (1 + zz) ** 3 + (1 - a) * (1 + zz) ** (3 * (1 + b)))
    else:
        e = (1 + zz) ** (1 + a + b) * np.exp(-b * zz / (1 + zz))
    integral = z / 2 * np.sum(weights[None, :] / e, axis=1)
    return 5 * np.log10(integral / z)


def quad_curve(z, family, a, b):
    def reciprocal(u):
        if family == "w":
            e = np.sqrt(a * (1 + u) ** 3 + (1 - a) * (1 + u) ** (3 * (1 + b)))
        else:
            e = (1 + u) ** (1 + a + b) * np.exp(-b * u / (1 + u))
        return 1 / e
    return 5 * np.log10(quad(reciprocal, 0, float(z), epsabs=1e-12, epsrel=1e-12)[0] / z)


def main():
    protocol = json.loads((OUT / "protocol.json").read_text())
    amendment = json.loads((OUT / "preexecution-amendment.json").read_text())
    result = json.loads((OUT / "result.json").read_text())
    assert sha(OUT / "protocol.json") == amendment["protocol_sha256"] == result["protocol_sha256"]
    assert sha(OUT / "preexecution-amendment.json") == result["preexecution_amendment_sha256"]
    assert sha(OUT / "frozen_source.py") == protocol["input_sha256"]["scripts/research_2026_09_26/sn_bao_design.py"] == amendment["old_source_sha256"]
    assert sha(OUT / "executed_source.py") == sha(ROOT / "scripts/research_2026_09_26/sn_bao_design.py") == amendment["new_source_sha256"]
    for rel, digest in protocol["input_sha256"].items():
        if not rel.endswith("sn_bao_design.py"):
            assert sha(ROOT / rel) == digest, rel
    for filename, digest in result["artifacts_sha256"].items():
        assert sha(OUT / filename) == digest, filename
    source = (OUT / "executed_source.py").read_text()
    assert "usecols=['CID','zHD','zHEL']" in source
    assert "m_b_corr" not in source[source.index("def design():"):source.index("if __name__")]

    # The metadata-only row selection keeps each observation, including repeats.
    metadata = pd.read_csv(SN / "Pantheon+SH0ES.dat", sep=r"\s+", usecols=["CID", "zHD", "zHEL"])
    raw_cov = np.loadtxt(SN / "Pantheon+SH0ES_STAT+SYS.cov")
    n = int(raw_cov[0]); assert len(metadata) == n and len(raw_cov) == 1 + n * n
    full = raw_cov[1:].reshape(n, n)
    asym = float(np.max(abs(full - full.T)))
    assert abs(asym - result["input_covariance_asymmetry_mag2"]) < 1e-15
    full = (full + full.T) / 2
    take = np.flatnonzero(metadata.zHD.to_numpy() > .01)
    sn = metadata.iloc[take].reset_index(drop=True)
    C = full[np.ix_(take, take)]
    z = sn.zHD.to_numpy()
    table = pd.read_csv(OUT / "design-table.csv")
    with np.load(OUT / "weights.npz", allow_pickle=False) as saved:
        archive = {k: saved[k] for k in saved.files}
    assert len(z) == result["sn_rows"] == len(archive["original_sn_rows"])
    assert sn.CID.nunique() == result["sn_distinct_CID"]
    assert np.array_equal(take, archive["original_sn_rows"])
    assert np.array_equal(z, archive["sn_z"])
    repeats = len(z) - sn.CID.nunique()

    bao = pd.read_csv(BAO / "desi_gaussian_bao_ALL_GCcomb_mean.txt", sep=r"\s+", comment="#", names=["z", "value", "kind"])
    raw_bao_cov = np.loadtxt(BAO / "desi_gaussian_bao_ALL_GCcomb_cov.txt")
    select = np.flatnonzero((bao.kind == "DM_over_rs") & (bao.z >= z.min()) & (bao.z <= z.max()))
    zb = bao.z.to_numpy()[select]; d = bao.value.to_numpy()[select]
    marginal = raw_bao_cov[np.ix_(select, select)]
    J = np.diag(5 / np.log(10) / d)
    cbmag = J @ marginal @ J
    assert np.array_equal(zb, archive["bao_z"])
    assert np.max(abs(cbmag - archive["bao_magnitude_cov"])) < 1e-15
    families = [("w", om, w) for om in (.15, .3, .5) for w in (-1.5, -1., -.5)]
    families += [("q", q0, q1) for q0 in (-1., 0., .5) for q1 in (0., 1., 2.)]
    zz = np.r_[z, zb]
    curves = np.array([curve(zz, *f) for f in families])
    sampled_z = np.r_[z[np.linspace(0, len(z) - 1, 20, dtype=int)], zb]
    max_quad_reference_difference = max(abs(curve(sampled_z, *f) - np.array([quad_curve(x, *f) for x in sampled_z])).max() for f in families)

    max_weight_difference = 0.0
    max_curve_bias_difference = 0.0
    max_sn_cov_difference = 0.0
    max_total_cov_difference = 0.0
    gate = np.ones(len(zb), dtype=bool)
    independent = {}
    for label, degree, width in (("primary", 2, .10), ("narrow", 2, .075), ("wide", 2, .125), ("cubic", 3, .10)):
        W = np.zeros((len(zb), len(z)))
        for j, node in enumerate(zb):
            local = np.flatnonzero(abs(np.log1p(z) - np.log1p(node)) <= width)
            x = (np.log1p(z[local]) - np.log1p(node)) / width
            X = np.column_stack([x ** k for k in range(degree + 1)])
            Cp = C[np.ix_(local, local)]
            L = np.linalg.cholesky(Cp)
            A = solve_triangular(L, X, lower=True)
            Q, R = np.linalg.qr(A, mode="reduced")
            direction = solve_triangular(R.T, np.eye(degree + 1)[0], lower=True)
            w = solve_triangular(L.T, Q @ direction, lower=False)
            W[j, local] = w
            ix = (table.design == label) & (abs(table.z_node - node) < 1e-12)
            assert ix.sum() == 1
            row = table.loc[ix].iloc[0]
            counts = [sn.CID.iloc[local].nunique(), sn.CID.iloc[local[x < 0]].nunique(), sn.CID.iloc[local[x > 0]].nunique()]
            bias = curves[:, :len(z)] @ W[j] - curves[:, len(z) + j]
            max_curve_bias_difference = max(max_curve_bias_difference, abs(abs(bias).max() - row.max_finite_family_bias_mag))
            cond = np.linalg.cond(A) ** 2
            passed = counts[0] >= 20 and min(counts[1:]) >= 5 and cond < 1e8 and abs(bias).max() <= .005
            assert bool(row.admitted) == bool(passed)
            assert counts == [int(row.distinct_CID), int(row.distinct_below), int(row.distinct_above)]
            assert abs(cond - row.condition) / row.condition < 1e-9
            assert abs(np.dot(w, Cp @ w) ** .5 - row.sn_sigma_mag) < 1e-11
            assert np.max(abs(w @ X - np.r_[1., np.zeros(degree)])) < 1e-10
            gate[j] &= passed
        max_weight_difference = max(max_weight_difference, float(np.max(abs(W - archive[label]))))
        csn = W @ C @ W.T
        V = csn + cbmag
        max_sn_cov_difference = max(max_sn_cov_difference, float(np.max(abs(csn - archive[label + "_sn_cov"]))))
        max_total_cov_difference = max(max_total_cov_difference, float(np.max(abs(V - archive[label + "_total_cov"]))))
        admitted = np.flatnonzero(archive["admitted"])
        if len(admitted) >= 2:
            lo, hi = admitted[0], admitted[-1]
            # With exactly two nodes, the free-intercept endpoint estimator is
            # the difference and its variance is the Helmert contrast of V.
            assert len(admitted) == 2
            contrast = np.zeros(len(zb)); contrast[lo] = -1; contrast[hi] = 1
            se = float(np.sqrt(contrast @ V @ contrast))
            assert abs(se - result["forecast"][label]["endpoint_drift_se_mag"]) < 1e-11
            for amplitude in (.02, .05, .10):
                critical = norm.ppf(.975)
                power = float(norm.cdf(-critical-amplitude/se) + norm.sf(critical-amplitude/se))
                assert abs(power - result["forecast"][label]["two_sided_5pct_power"][str(amplitude)]) < 1e-12
            independent[label] = {"free_intercept_contrast_se": se, "first_three_primary_weight_max_abs_error": float(np.max(abs(W[:3] - archive[label][:3]))) if label == "primary" else None}
    assert np.array_equal(gate, archive["admitted"])
    assert np.array_equal(zb[gate], np.asarray(result["admitted_nodes"]))
    assert max_weight_difference < 1e-9
    assert max_curve_bias_difference < 1e-9
    assert max_sn_cov_difference < 1e-10
    assert max_total_cov_difference < 1e-10
    assert max_quad_reference_difference < 1e-10
    verification = {
        "gate": "pass",
        "scope": "Outcome-free independent design verification; no observed SN magnitudes loaded or drift scored.",
        "input_and_output_hashes_pass": True,
        "preexecution_source_amendment_pass": True,
        "sn_rows": len(z),
        "distinct_cid": int(sn.CID.nunique()),
        "repeated_observation_rows": repeats,
        "marginal_bao_dm_nodes": zb.tolist(),
        "admitted_intersection": zb[gate].tolist(),
        "primary_only_admitted": table[(table.design == "primary") & table.admitted].z_node.tolist(),
        "max_qr_gls_weight_difference": max_weight_difference,
        "max_recomputed_curve_bias_difference_mag": max_curve_bias_difference,
        "max_quadrature_vs_scipy_quad_difference_mag": float(max_quad_reference_difference),
        "max_recomputed_sn_covariance_difference_mag2": max_sn_cov_difference,
        "max_recomputed_total_covariance_difference_mag2": max_total_cov_difference,
        "free_intercept_forecast": independent,
        "source_sha256": sha(Path(__file__)),
    }
    (OUT / "independent-check.json").write_text(json.dumps(verification, indent=2, allow_nan=False) + "\n")
    print(json.dumps(verification, indent=2))


if __name__ == "__main__":
    main()
